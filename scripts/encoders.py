###copyright Alan Moses 2021-2025
##permission to use and copy is granted for research purposes only
##software provided as is with no warrantee of any kind

import numpy as np
import tensorflow as tf


# ── One-hot encoding helpers ───────────────────────────────────────────────────

def standardize_dna(seq, target_len, repeat_pad=True, random_crop=True):
    """Pad or crop a DNA string to exactly target_len characters."""
    if len(seq) == target_len:
        return seq
    if len(seq) > target_len:
        if not random_crop:
            return seq[0:target_len]
        b = np.random.randint(0, len(seq) - target_len)
        return seq[b:b + target_len]
    # shorter: repeat-pad then random-crop
    if not repeat_pad:
        return seq + '0' * (target_len - len(seq))
    return standardize_dna(seq + seq, target_len)


def dna_one_hot2(seq, standardize=None):
    return _one_hot(seq, standardize, bases="ACGT")


def _one_hot(seq, standardize=None, bases="ACGT"):
    if standardize is not None:
        seq = standardize_dna(seq, standardize, repeat_pad=False)
    seq_length = len(seq)
    one_hot_seq = np.zeros((seq_length, len(bases)))
    index = {bases[b]: b for b in range(len(bases))}
    for j in range(seq_length):
        if seq[j] in index:
            one_hot_seq[j, index[seq[j]]] = 1
    return one_hot_seq


# ── Custom Keras layers ────────────────────────────────────────────────────────

class Reverse_Complement(tf.keras.layers.Layer):
    """Reverse-complement a batch of one-hot encoded DNA sequences."""
    def __init__(self, complement=None, **kwargs):
        super().__init__(**kwargs)
        self.com = complement if complement is not None else [3, 2, 1, 0]

    def call(self, seqs):
        rev = tf.reverse(seqs, axis=[1])
        return tf.gather(rev, indices=self.com, axis=-1)

    def get_config(self):
        cfg = super().get_config()
        cfg.update({"complement": self.com})
        return cfg


class TrainablePooling(tf.keras.layers.Layer):
    """
    Learns continuously between max and average pooling
    (from the intrinsically interpretable encoder paper, Balci et al. 2023).
    """
    def __init__(self, seqL, filterN, initial_alpha=1.0, trainable=True, **kwargs):
        super().__init__(**kwargs)
        self.fn = filterN
        self.sl = seqL
        self.ia = initial_alpha
        self.tr = trainable

    def build(self, input_shape):
        self.alpha = self.add_weight(
            "alpha", shape=(1, self.fn),
            initializer=tf.keras.initializers.Constant(self.ia),
            constraint=tf.keras.constraints.NonNeg(),
            trainable=self.tr)

    def call(self, inputs):
        A = tf.repeat(self.alpha, self.sl, axis=0)
        W = tf.nn.softmax(A * inputs, axis=1)
        return tf.reduce_sum(W * inputs, axis=1)

    def get_config(self):
        cfg = super().get_config()
        cfg.update({"fn": self.fn, "sl": self.sl,
                    "ia": self.ia, "tr": self.tr})
        return cfg


class TrainableScaling(tf.keras.layers.Layer):
    """
    Learnable per-motif scale and threshold
    (from the intrinsically interpretable encoder paper, Balci et al. 2023).
    """
    def __init__(self, seqL, filterN, trainable=True,
                 initial_center=0.0, initial_scale=1.0, **kwargs):
        super().__init__(**kwargs)
        self.fn = filterN
        self.sl = seqL
        self.ic = initial_center
        self.isc = initial_scale
        self.tr = trainable

    def build(self, input_shape):
        self.center = self.add_weight(
            "center", shape=(1, self.fn),
            initializer=tf.keras.initializers.Constant(self.ic),
            constraint=tf.keras.constraints.NonNeg(),
            trainable=self.tr)
        self.scale = self.add_weight(
            "scale", shape=(1, self.fn),
            initializer=tf.keras.initializers.Constant(self.isc),
            constraint=tf.keras.constraints.NonNeg(),
            trainable=self.tr)

    def call(self, inputs):
        C = tf.repeat(self.center, self.sl, axis=0)
        S = tf.repeat(self.scale,  self.sl, axis=0)
        return S * inputs - C

    def get_config(self):
        cfg = super().get_config()
        cfg.update({"fn": self.fn, "sl": self.sl,
                    "ic": self.ic, "isc": self.isc, "tr": self.tr})
        return cfg


class PWMConstraint(tf.keras.constraints.Constraint):
    """
    Projects conv1d weights onto the space of valid log-odds PWMs
    so that probabilities at each position sum to 1.
    See Moses et al. 2026 Methods for derivation.
    """
    def __init__(self, bg=None):
        self.bg = bg if bg is not None else tf.constant([0.25, 0.25, 0.25, 0.25])
        self._bg = tf.reshape(self.bg, (1, 1, -1))   # shape (1,1,4) for broadcasting

    def __call__(self, w):
        # w shape: (filter_len, 4, filter_num)
        bg_part    = tf.math.log(self._bg)                              # (1, 1, 4) broadcast
        motif_part = tf.math.log(
            tf.reduce_sum(tf.math.exp(w), axis=1, keepdims=True))      # (L, 1, N)
        return w - bg_part - motif_part

    def get_config(self):
        return {'bg': self.bg.numpy().tolist()}


# ── Main encoder class ─────────────────────────────────────────────────────────

class DNAencoder:
    def __init__(self, filterN=256, filterL=15, alphabet=4,
                 bases="ACGT", trainable=True, print_summary=True):
        self.filterN       = filterN
        self.filterL       = filterL   # default matches paper (15), not 12
        self.alphabet      = alphabet
        self.bases         = bases
        self.trainable     = trainable
        self.print_summary = print_summary
        self.com = [3, 2, 1, 0]
        if bases == "ATGC":
            self.com = [1, 0, 3, 2]
        self.RC = Reverse_Complement(complement=self.com)
        if self.print_summary:
            print(f"Building DNA encoder... bases={bases}, com={self.com}")

    def interpretable_encoder(self, seq_len, W=10, activation_type='sigmoid'):
        """
        Interpretable motif-based encoder (Balci et al. 2023 / Moses et al. 2026).

        Architecture:
          Input
            -> PWM conv scan (both strands, max-pool over strands)
            -> MaxPool1D (window W to avoid counting overlapping matches)
            -> TrainableScaling  (learnable threshold + scale per motif)
            -> Sigmoid activation
            -> TrainablePooling  (learnable between max and avg pooling)
            -> BatchNormalization  ← CORRECT POSITION (before interactions)
            -> Dense (motif interaction matrix)
            -> Sigmoid
            -> Element-wise multiply with pooled matches
          Output (256-dim representation)

        FIX vs original: BatchNorm is placed BEFORE the interaction Dense layer,
        not after it. Placing it after the sigmoid-weighted interactions caused
        motif collapse (many motifs converging to the same representation).
        """
        seq_in = tf.keras.layers.Input(
            shape=(seq_len, self.alphabet), name="seqs")

        # ── Strand-aware PWM scanning ──────────────────────────
        conv = tf.keras.layers.Conv1D(
            self.filterN, self.filterL,
            padding="valid",
            use_bias=False,
            kernel_constraint=PWMConstraint(),
            name="matrix_scan")

        fwd = conv(seq_in)
        rev = conv(self.RC(seq_in))
        rev = tf.keras.layers.Lambda(
            lambda t: tf.reverse(t, axis=[1]))(rev)

        fwd = tf.expand_dims(fwd, axis=1)
        rev = tf.expand_dims(rev, axis=1)
        both = tf.keras.layers.Concatenate(axis=1)([fwd, rev])

        # Max over the two strands at each position
        matches = tf.keras.layers.MaxPooling2D(
            pool_size=(2, 1), strides=None,
            padding="valid", data_format="channels_last",
            name="rc_max_pool")(both)
        matches = tf.keras.layers.Reshape(
            (matches.shape[2], self.filterN))(matches)

        # ── Non-overlapping match counting ─────────────────────
        matches = tf.keras.layers.MaxPooling1D(W)(matches)

        # ── Per-motif scaling and thresholding ─────────────────
        matches = TrainableScaling(
            matches.shape[1], self.filterN,
            initial_center=0.0, initial_scale=1.0,
            name="motif_scaling", trainable=True)(matches)
        matches = tf.keras.layers.Activation(activation_type)(matches)

        # ── Pooling over sequence positions ────────────────────
        matches = TrainablePooling(
            matches.shape[1], self.filterN,
            initial_alpha=2.0,
            name="match_pooling", trainable=True)(matches)

        # BatchNorm BEFORE interactions -
        # This stabilises training by normalising the per-motif scores
        # before the interaction matrix mixes them. Placing it AFTER
        # the interactions (as in the original buggy code) caused the
        # interaction matrix to dominate and drove motif collapse.
        matches = tf.keras.layers.BatchNormalization(
            center=True, scale=True)(matches)

        # ── Motif–motif interaction matrix ─────────────────────
        weights = tf.keras.layers.Dense(
            self.filterN,
            use_bias=False,
            activation=None,
            kernel_initializer=tf.keras.initializers.Identity(gain=1.0),
            trainable=True,
            name="motif_interactions")(matches)
        weights = tf.keras.layers.Activation('sigmoid')(weights)

        # ── Final representation ────────────────────────────────
        out = weights * matches

        model = tf.keras.models.Model(inputs=seq_in, outputs=out)
        if self.print_summary:
            model.summary()
        return model
