###copyright Alan Moses 2021-2025
##permission to use and copy is granted for research purposes only
##software provided as is with no warrantee of any kind

###trains an interpretable encoder using reverse homology (contrastive learning).
###outputs: model weights (.h5), matrix layer text file, training info csv.
###reads a directory of per-OG fasta files of orthologous promoters.

import sys
import os
import numpy as np
import tensorflow as tf
import mol_evol as me
import encoders as my_encoders

# ── Hyperparameters ────────────────────────────────────────────────────────────
# These match 00_config.sh. Keep them in sync if you change config.
learning_rate = 0.001
batch_size    = 256
num_epochs    = 100
encoder_type  = 'interpretable'
filter_num    = 256    # number of PWMs
filter_len    = 15     # PWM width 
fam_size      = 4      # homologs averaged per contrastive step

# SEQ_L=800 during training: the model is trained on 800bp sequences.
# compute_promoter_rep.py uses SEQ_L=500 for the final representations.
# This split is intentional and matches the paper (Methods section).
SEQ_L = 800

N          = fam_size + 4   # minimum sequences required per OG file
FILE_EXT   = ".fa"
REVCOM_AUG = True           # reverse-complement augmentation with p=0.5

# v_str: prefix of OG names held out for validation / early stopping.
# Set to None to train for the full num_epochs (used for Aspergillus
# where no canonical validation split exists).
v_str    = None
test_str = None

SAMPLE_BY_CLASS = False
learn_motifs    = True
alignment_QC    = True

SAVE_INIT       = False
SAVE_FINAL      = True
SAVE_MOTIF_DB   = True

# ── Input ──────────────────────────────────────────────────────────────────────
if len(sys.argv) < 2:
    print("usage: promoter_reverse_homology.py <path_to_og_up800_dir>")
    sys.exit(1)

input_dir = sys.argv[1]
filename  = (f"promoter_rh_{encoder_type}"
             f"N{filter_num}L{filter_len}_f{fam_size}_b{batch_size}_e{num_epochs}")

# ── Load data ──────────────────────────────────────────────────────────────────
print("Reading sequence data...")
DB = me.AlignmentDB()
DB.read_from_fa_dir(input_dir, ext=FILE_EXT)
print(f"{len(DB.names)} OG files loaded. Preparing training/validation splits...")


def good_DNA_seqs(ALN, max_gap_frac=0.9):
    """Return names of sequences that pass quality filters."""
    good = []
    for s in ALN.names:
        seq  = ALN.seq[s]
        gaps = seq.count('-')
        if len(seq) == gaps:
            continue
        ACGTs = sum(seq.count(b) for b in "ACGT")
        if ACGTs + gaps != len(seq):
            continue                          # no ambiguous bases
        if gaps > 0 and float(gaps) / float(gaps + ACGTs) > max_gap_frac:
            continue
        good.append(s)
    return good


seqs      = []
labels    = []
testseqs  = []
testlabels = []
i = 0

for a in DB.names:
    if test_str is not None and test_str in a:
        continue
    ALN = DB.data[a]
    if len(ALN.names) < N:
        continue
    for s in ALN.names:
        ALN.seq[s] = ALN.seq[s].upper().replace("-", "")
    sequences_to_use = good_DNA_seqs(ALN) if alignment_QC else ALN.names
    if len(sequences_to_use) < N:
        continue
    if v_str is not None and v_str in a:
        for j in range(len(sequences_to_use)):
            testseqs.append([a, sequences_to_use[j]])
            testlabels.append(int(i))
    else:
        for j in range(len(sequences_to_use)):
            seqs.append([a, sequences_to_use[j]])
            labels.append(int(i))
    i += 1

seqs       = np.array(seqs)
labels     = np.array(labels)
testseqs   = np.array(testseqs)
testlabels = np.array(testlabels)

print(f"Training sequences: {seqs.shape}, Validation sequences: {testseqs.shape}")

steps   = len(seqs)     // batch_size
v_steps = len(testseqs) // batch_size


# ── Batch generator ────────────────────────────────────────────────────────────

def RH_batch_pos_only(DB, seqs, labels, idx,
                      fam_size=4, batch_size=256,
                      sample_by_class=False, stand=None):
    refs      = []
    positives = []
    pos_classes = []

    for i in range(batch_size):
        if sample_by_class:
            pos_class    = idx[i]
            all_pos_i    = np.squeeze(np.argwhere(labels == pos_class))
            pos_sample   = np.random.choice(all_pos_i, size=fam_size + 1, replace=False)
            pos_i        = pos_sample[-1]
        else:
            pos_i        = idx[i]
            pos_class    = labels[idx[i]]
            all_pos_i    = np.squeeze(np.argwhere(labels == pos_class))
            all_pos_i    = np.delete(all_pos_i, np.where(all_pos_i == pos_i))
            pos_sample   = np.random.choice(all_pos_i, size=fam_size, replace=False)

        pos_classes.append(pos_class)
        fam_i = pos_sample[0:fam_size]
        positives.append(pos_i)

        these_refs = []
        for j in fam_i:
            dna = DB.data[seqs[j][0]].seq[seqs[j][1]]
            oh  = my_encoders.dna_one_hot2(dna, standardize=stand)
            if REVCOM_AUG and np.random.choice([True, False]):
                oh = oh[::-1, ::-1]   # reverse complement in one-hot space
            these_refs.append(oh)
        refs.append(np.array(these_refs))

    pos_classes = np.array(pos_classes)
    outputs     = np.equal.outer(pos_classes, pos_classes).astype(float)
    outputs    /= np.sum(outputs, axis=-1, keepdims=True)

    reorder   = np.random.permutation(batch_size)
    outputs   = np.take(outputs, reorder, axis=1)
    positives = [positives[j] for j in reorder]

    pseqs = []
    for j in positives:
        dna = DB.data[seqs[j][0]].seq[seqs[j][1]]
        oh  = my_encoders.dna_one_hot2(dna, standardize=stand)
        if REVCOM_AUG and np.random.choice([True, False]):
            oh = oh[::-1, ::-1]
        pseqs.append(oh)

    return [np.array(refs), np.array(pseqs)], outputs


def RH_generator(DB, seqs, labels, fam_size=4, batch_size=256,
                 sample_by_class=False, stand=None):
    i_list = range(len(seqs)) if not sample_by_class else np.unique(labels)
    if len(i_list) <= batch_size:
        print(f"Only {len(i_list)} indices for batch_size={batch_size}")
        sys.exit(1)
    while True:
        rand_ord = np.random.permutation(len(i_list))
        i = 0
        while i < len(i_list) - batch_size:
            yield RH_batch_pos_only(DB, seqs, labels,
                                    rand_ord[i:i + batch_size],
                                    fam_size=fam_size,
                                    batch_size=batch_size,
                                    sample_by_class=sample_by_class,
                                    stand=stand)
            i += batch_size


# ── Build model ────────────────────────────────────────────────────────────────
if SAMPLE_BY_CLASS:
    sample_i = np.unique(labels)[range(batch_size)]
else:
    sample_i = range(batch_size)

bseqs, _ = RH_batch_pos_only(DB, seqs, labels, sample_i,
                              fam_size=fam_size, batch_size=batch_size,
                              sample_by_class=SAMPLE_BY_CLASS, stand=SEQ_L)
seq_len = bseqs[1].shape[1]

enc = my_encoders.DNAencoder(
    filterN=filter_num, filterL=filter_len,
    trainable=learn_motifs
).interpretable_encoder(seq_len=seq_len)

if SAVE_INIT:
    enc.save_weights(filename + "_init_weights.h5", save_format='h5')


def build_rh_model(encoder, seq_len, fam_size=4, bsize=256, temperature=0.1):
    fam_in = tf.keras.layers.Input(
        shape=(fam_size, seq_len, 4), batch_size=bsize, name="fam")
    pos_in = tf.keras.layers.Input(
        shape=(seq_len, 4), batch_size=bsize, name="pos")

    fam4enc = tf.reshape(fam_in, (bsize * fam_size, seq_len, 4))
    fam_z   = encoder(fam4enc)
    fam_z   = tf.reshape(fam_z, (bsize, fam_size, -1))
    ref_z   = tf.math.l2_normalize(
        tf.reduce_mean(fam_z, axis=1), axis=-1)

    pos_z   = tf.math.l2_normalize(encoder(pos_in), axis=-1)
    pos_z   = tf.repeat(tf.reshape(pos_z, (1, bsize, -1)), bsize, axis=0)

    preds   = tf.keras.layers.Dot(axes=-1, name='pos_dot')([ref_z, pos_z])
    preds   = preds / temperature
    preds   = tf.keras.layers.Softmax(name='logits', axis=-1)(preds)
    return tf.keras.models.Model(inputs=[fam_in, pos_in], outputs=preds)


cmod = build_rh_model(enc, seq_len=seq_len,
                      fam_size=fam_size, bsize=batch_size)
cmod.compile(
    optimizer=tf.keras.optimizers.Adam(learning_rate),
    loss=tf.keras.losses.categorical_crossentropy,
    metrics=[tf.keras.metrics.categorical_accuracy])
cmod.summary()

# ── Train ──────────────────────────────────────────────────────────────────────
train_gen = RH_generator(DB, seqs, labels,
                         fam_size=fam_size, batch_size=batch_size,
                         sample_by_class=SAMPLE_BY_CLASS, stand=SEQ_L)

if v_steps > 0:
    val_gen = RH_generator(DB, testseqs, testlabels,
                           fam_size=fam_size, batch_size=batch_size,
                           sample_by_class=SAMPLE_BY_CLASS, stand=SEQ_L)
    init_acc = cmod.evaluate(val_gen, steps=v_steps)[1]
    print(f"Initial val accuracy: {init_acc * 100:.2f}%")
    callbacks  = [tf.keras.callbacks.EarlyStopping(
        patience=10, monitor="val_loss", restore_best_weights=True)]
    val_kwargs = dict(validation_data=val_gen, validation_steps=v_steps)
else:
    print("No validation set (v_str=None). Training for full num_epochs.")
    init_acc   = float('nan')
    callbacks  = []
    val_kwargs = {}

history = cmod.fit(
    train_gen,
    steps_per_epoch=steps,
    epochs=num_epochs,
    callbacks=callbacks,
    **val_kwargs)

# ── Save motif database ────────────────────────────────────────────────────────
if SAVE_MOTIF_DB:
    bases   = "ACGT"
    outname = filename + "_matrix_layer.txt"
    print(f"Saving matrix layer to {outname}")
    with open(outname, 'w') as out:
        M   = np.squeeze(enc.get_layer('matrix_scan').get_weights())
        int_w = enc.get_layer('motif_interactions').get_weights()[0]
        pool  = np.squeeze(enc.get_layer('match_pooling').get_weights())
        scale = np.squeeze(enc.get_layer('motif_scaling').get_weights())
        thr   = scale[0, :]
        sca   = scale[1, :]
        for i in range(filter_num):
            f = M[:, :, i]
            out.write(f">matrix{i} length={filter_len}"
                      f"\tthr={thr[i]}\tsca={sca[i]}\tpool={pool[i]}"
                      f"\tinteractions:\t"
                      + "\t".join(str(int_w[i, j]) for j in range(filter_num))
                      + "\n")
            for b in range(4):
                out.write(bases[b] + "\t"
                          + "\t".join(str(f[p, b]) for p in range(filter_len))
                          + "\n")

if SAVE_FINAL:
    print(f"Saving encoder weights to {filename}_final_weights.h5")
    enc.save_weights(filename + "_final_weights.h5", save_format='h5')

# ── Save training log ──────────────────────────────────────────────────────────
log_name = filename + "_training_info.csv"
with open(log_name, 'w') as out:
    out.write(f"{encoder_type} {filter_num} filters x {filter_len}nt, "
              f"fam_size={fam_size}, batch={batch_size}\n")
    out.write(f"Initial val accuracy: {init_acc * 100:.2f}%\n")
    for k, v in history.history.items():
        out.write(k + "," + ",".join(str(x) for x in v) + "\n")
print(f"Training log saved to {log_name}")
