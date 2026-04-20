import numpy as np
import tensorflow as tf


def standardize_dna(seq,target_len,repeat_pad=True,random_crop=True): ##this does repeat pad + random crop augmentation
	if len(seq)==target_len: return seq
	if len(seq) > target_len:
		if not random_crop: return seq[0:target_len] #truncate the end
		b=np.random.randint(0,len(seq)-target_len)
		return seq[b:b+target_len]
	elif len(seq) < target_len:
		if not repeat_pad:
			tail= '0' * (target_len-len(seq)) #0-pad the end
			return seq + tail
		return standardize_dna(seq + seq,target_len) # once the sequence gets long enough
													#it will get randomly cropped
	else:
		print("problem with standardizing ",seq)
		exit()


def dna_one_hot(seq,standardize=None,bases = "ACGT"):
	if standardize is not None: seq=standardize_dna(seq,standardize,repeat_pad=False)
	#byte_seq=tf.cast(tf.io.decode_raw(seq,tf.uint8),tf.int32) ##turns out this iis 5x slower!
	#return tf.nn.embedding_lookup(table,byte_seq) ##turns out this iis 5x slower!
	seq_length = len(seq)
	one_hot_seq = np.zeros((seq_length, len(bases)))
	for b in range(0, len(bases)):
		index = [j for j in range(0, seq_length) if seq[j] == bases[b]]
		one_hot_seq[index, b] = 1
	return one_hot_seq

def dna_one_hot2(seq,standardize=None):
	return aa_one_hot2(seq,standardize,bases="ACGT") ##should now do the faster version... ~3x faster for DNA

def aa_one_hot(seq,standardize=None): ##adapted from DNA one hot. probably there is a faster way to do it...
	return dna_one_hot(seq,standardize,bases = "ACDEFGHIKLMNPQRSTVWY")

def aa_one_hot2(seq,standardize=None,bases = "ACDEFGHIKLMNPQRSTVWY"): ##much faster yay!
	if standardize is not None: seq = standardize_dna(seq, standardize, repeat_pad=False)
	seq_length = len(seq)
	one_hot_seq = np.zeros((seq_length, len(bases)))
	#index = { b:[] for b in bases }
	index = { bases[b]:b for b in range(0,len(bases)) } ## simple lookup table
	for j in range(0, seq_length): ##only loop over the sequence length once.
		#if seq[j] in bases: index[seq[j]].append(j) ##list of indicies
		if seq[j] in bases: one_hot_seq[j,index[seq[j]]] = 1 ##just assigning them one by one is a bit faster
	#for j in range(0, len(bases)): one_hot_seq[index[bases[j]],j] = 1 ##vectorized indexing
	return one_hot_seq

#print(dna_one_hot('AAACCCGGGTTT',dna_table))
#print(np.shape(dna_one_hot('AAACCCGGGTTT',dna_table)))


class Reverse_Complement(tf.keras.layers.Layer):
	def __init__(self, complement=[3, 2, 1, 0], **kwargs):  ##default is ACGT order
		super(Reverse_Complement, self).__init__(**kwargs)
		self.com = complement

	# print("creating an RC layer with "+str(self.com))
	##expects a batch of one hot encoded sequences
	def call(self, seqs):
		rev = tf.reverse(seqs, axis=[1])
		revcom = tf.gather(rev, indices=self.com, axis=-1)
		return revcom

	def get_config(self):  ###not sure if this worked...
		config = super(Reverse_Complement, self).get_config()
		config.update({"complement": self.com})


class TrainablePooling(	tf.keras.layers.Layer):  # learns continuiously between max and average pooling, from intrinsically interpretable
	def __init__(self, seqL, filterN, initial_alpha=1.0, trainable=True, **kwargs):
		super(TrainablePooling, self).__init__(**kwargs)
		self.fn = filterN
		self.sl = seqL
		self.ia = initial_alpha
		self.tr = trainable

	def build(self, input_shape):
		self.alpha = self.add_weight("alpha", shape=(1, self.fn),
									 initializer=tf.keras.initializers.Constant(self.ia),
									 constraint=tf.keras.constraints.NonNeg(),
									 trainable=self.tr)

	def call(self, inputs):
		A = tf.repeat(self.alpha, self.sl, axis=0)
		# print ("A:",np.shape(A))
		W = tf.nn.softmax(A * inputs, axis=1)
		# print("W:",np.shape(W)," inputs",np.shape(inputs)," out:",np.shape(tf.reduce_sum(W*inputs,axis=1)) )
		return tf.reduce_sum(W * inputs, axis=1)

	def get_config(self):  ###doesn't seem to need this to save to h5...
		config = super(TrainablePooling, self).get_config()
		config.update({"fn": self.fn})
		config.update({"sl": self.sl})
		config.update({"ia": self.ia})
		config.update({"tr": self.tr})
		config.update({"alpha": self.alpha})


class TrainableScaling(tf.keras.layers.Layer):  # scales and subtracts from intrinsically interpretable
	def __init__(self, seqL, filterN, trainable=True, initial_center=0.0, initial_scale=1.0, min_scale=None, **kwargs):
		super(TrainableScaling, self).__init__(**kwargs)
		self.fn = filterN
		self.sl = seqL
		self.ic = initial_center
		self.isc = initial_scale
		self.tr = trainable
		self.minsc = min_scale

	def build(self, input_shape):
		self.center = self.add_weight("center", shape=(1, self.fn), initializer=tf.keras.initializers.Constant(self.ic),
									  constraint = tf.keras.constraints.NonNeg(),
									  trainable=self.tr)
		self.scale = self.add_weight("scale", shape=(1, self.fn), initializer=tf.keras.initializers.Constant(self.isc),
									 constraint=tf.keras.constraints.NonNeg(),
									 trainable=self.tr)
		if self.minsc is not None: self.min_sc = self.add_weight("min_scale", shape=(1, self.fn), initializer=tf.keras.initializers.Constant(self.minsc),
									  trainable=False) ##try to prevent the scale from going to 0.

	def call(self, inputs):
		C = tf.repeat(self.center, self.sl, axis=0)
		S = tf.repeat(self.scale, self.sl, axis=0)
		if self.minsc is not None: minS = tf.repeat(self.min_sc, self.sl, axis=0)
		# print ("A:",np.shape(A))
		if self.minsc is not None: S = tf.math.maximum(S,minS) ##try to prevent the scale from going to 0, because then no gradients for that filter...
		out = S * inputs
		out = out - C
		# out = S*out
		# print("W:",np.shape(W)," inputs",np.shape(inputs)," out:",np.shape(tf.reduce_sum(W*inputs,axis=1)) )
		return out


class PWMConstraint(tf.keras.constraints.Constraint):
	def __init__(self, bg=tf.constant([0.25, 0.25, 0.25, 0.25])):
		self.bg = tf.expand_dims(bg, axis=0)
		self.bg = tf.expand_dims(self.bg, axis=-1)

	def __call__(self, w):
		# print(w.numpy())
		# print(w.get_shape())
		bg_part = tf.repeat(tf.math.log(self.bg), repeats=w.get_shape()[0], axis=0)
		bg_part = tf.repeat(bg_part, repeats=w.get_shape()[-1], axis=-1)
		# print("bg:",bg_part.numpy())
		# print("bg:",bg_part.get_shape())
		motif_part = tf.math.reduce_sum(tf.math.exp(w), axis=-2, keepdims=True)
		motif_part = tf.repeat(tf.math.log(motif_part), repeats=w.get_shape()[-2], axis=-2)
		# print(motif_part.numpy())
		# print(motif_part.get_shape())
		return w - bg_part - motif_part

	def get_config(self):
		return {'bg': self.bg}


class DNAencoder():
	def __init__(self, filterN=256, filterL=12, alphabet=4, bases="ACGT", trainable=True, print_summary=True):
		super(DNAencoder, self).__init__()
		self.print_summary = print_summary
		self.filterN = filterN
		self.filterL = filterL
		self.alphabet = alphabet
		self.bases = bases
		self.trainable = trainable
		self.com = [3, 2, 1, 0]  ##needs to match bases
		if (self.bases == "ATGC"): self.com = [1, 0, 3, 2]
		self.RC = Reverse_Complement(complement=self.com)
		if self.print_summary: print('Bulding DNA encoder...' + self.bases + str(self.com))


	def interpretable_encoder(self, seq_len, W=10, activation_type='sigmoid'):  ##after the tiSFM intrinsically interpretable paper
		##seems like we need sigmoids to make the trainable pooling interpretable.
		seq_in = tf.keras.layers.Input(shape=(seq_len, self.alphabet), name="seqs")
		if self.trainable:
			l1 = tf.keras.layers.Conv1D(self.filterN, self.filterL, padding="valid", use_bias=False,
										kernel_constraint=PWMConstraint(),  # tf.keras.constraints.NonNeg(),
										name="matrix_scan")
			forward_scan = l1(seq_in)  # -1.0*l1(seq_in) #for non-negative constraint
			revcom_scan = l1(self.RC(seq_in))  # -1.0*l1(self.RC(seq_in))
			revcom_scan = tf.keras.layers.Lambda(lambda t: tf.reverse(t, axis=[1]))(revcom_scan)
			forward_scan = tf.expand_dims(forward_scan, axis=1)  ##if
			revcom_scan = tf.expand_dims(revcom_scan, axis=1)
			both_scan = tf.keras.layers.Concatenate(axis=1)([forward_scan, revcom_scan])
			layer2 = tf.keras.layers.MaxPooling2D(pool_size=(2, 1),  ##max pool over both the rc
												  strides=None,
												  padding="valid",
												  data_format="channels_last",
												  name="rc_max_pool")
			matches = layer2(both_scan)
			matches = tf.keras.layers.Reshape((matches.get_shape()[2], self.filterN))(matches)
			# matches = tf.keras.layers.BatchNormalization(center=True,scale=True)(matches) ##still interpretable?
			matches = tf.keras.layers.MaxPooling1D(W)(
				matches)  ##max pool over W to ensure overlapping motifs are not counted
			matches = TrainableScaling(matches.get_shape()[1], self.filterN, initial_center=0.0, initial_scale=1.0, name="motif_scaling", trainable=True)(matches)
			matches = tf.keras.layers.Activation(activation_type)(matches)
			matches = TrainablePooling(matches.get_shape()[1], self.filterN, initial_alpha=2.0, name="match_pooling", trainable=True)(matches)
			weights = tf.keras.layers.Dense(self.filterN,  ##create a trainable weight matrix of TF x TF interactions
											use_bias=False,
											activation=None,
											kernel_initializer=tf.keras.initializers.Identity(gain=1.0), #initially, ignore presence of others
											#kernel_regularizer=keras.regularizers.L1(l1=0.01) #could this help reliability of interactions?
											trainable=True,  ##gets to loss 3.95 with b64 with this set to identity.
											name="motif_interactions")(matches)
		# print(weights)
		else:
			print ("not available")
			exit()
		weights = tf.keras.layers.Activation('sigmoid')(weights)
		out = weights * matches
		out = tf.keras.layers.BatchNormalization()(
			out)  ##seems to help learning with GELUs and sigmoid and rh=1. sigmoid still not great.
		model = tf.keras.models.Model(inputs=seq_in, outputs=out)
		if self.print_summary:
			print(model.summary())
		return model
