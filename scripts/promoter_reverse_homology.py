###copyright Alan Moses 2021-2025
##permission to use and copy is granted for research purposes only
##software provided as is with no warrantee of any kind

###this is supposed to be 1.0 released with paper about inferring fungal regulatory networks

### trains an interpretable encoder
###outputs: model hd5, motifs and interactions in a txt, training info in a csv

###reads a directory of fasta files of orthologous promoters.
### gaps will be removed, but otherwise all non-ACGT will just be ignored


import sys
import numpy as np
import tensorflow as tf
import mol_evol as me ## molecular evolution functions
import tf_tf ## transcription factors for tensorflow
import encoders as my_encoders


learning_rate = 0.001  #also works for batch size 64
batch_size = 256 ##this is also the comp_size+1 if reusing negatives.
num_epochs = 100 ##100 used #2 was what was initially used
encoder_type= 'interpretable'
filter_num=256 ##this is the number of PWMs that will be learned
filter_len=15  ##this is the width of the PWM
fam_size=4 ##number of homologs to average for reverse homology. 1=simCLR, 8=Lu et al. PLoS Comp. Bio. 2022

#SEQ_L=500 #max sequence length. Longer sequences will be truncated, shorter will be 0- padded #Uncomment this part for using original script
SEQ_L=800 #This is something I am changing to train my data on 800bp first 
N=fam_size+4 ##at least N examples in the class/homolog set
FILE_EXT=".fa" ##this will be used to find the input files in the directory

#choose how to define the validation split used for early stopping only
#v_str="NCU00" #v_str="YB" ##chromosome 2 has 379 alignments  #v_str="YAR" ##only 15 alignments #v_str="NCU00" ##889 homolog sets for neurospora #v_str="C4_" ### for c. albicans
v_str=None ##just go through all the num_epochs
#test_str="None"##"chr2_" ## skip these files entirely (if you want a held out test set.) Not usually needed because the method is unsupervised
test_str=None #modified for Aspergillus

SAMPLE_BY_CLASS=False ##true is the way Alex did it in PLoS CB 2022. Avoids clashes #False is what Alan was doing
### True Doesn't see all the sequences in each epoch, so much faster.
##if True, number of classes in the validation set must be bigger than the batch size/comp_size.
learn_motifs=True
alignment_QC=True

SAVE_INIT=False ##save the random initial encoder weights to a file
SAVE_FINAL=True ##save the final encoder weights to a file
SAVE_MOTIF_DB=True ## save the PWMs and associated parameters to a human readible text file

if (len(sys.argv)<2):
	print ("usage: path to alignment directory")
	quit()	


filename="promoter_rh_"+encoder_type+"N"+str(filter_num)+"L"+str(filter_len)+"_f"+str(fam_size)+"_b"+str(batch_size)
filename=filename + "_e" + str(num_epochs)

print("reading sequence data...")
DB=me.AlignmentDB() # intialize the data structure to hold alignments
DB.read_from_fa_dir(sys.argv[1],ext=FILE_EXT) ##read in all the fasta files
print(str(len(DB.names)) + " alignments loaded. preparing training/test splits...")


def good_DNA_seqs(ALN,max_gap_frac=0.9,standardize=None,remove_sp=None): ##quality control filters on sequences, returns the names
	goodseqs=[]
	for s in ALN.names:
		if remove_sp is not None:
			if remove_sp in s: continue
		gaps=ALN.seq[s].count('-')
		#print (s+" "+str(gaps))
		if (len(ALN.seq[s])==gaps): continue
		ACGTs=0
		for b in list("ACGT"):
			ACGTs+=ALN.seq[s].count(b)
		if standardize is not None:
			if ACGTs < standardize: continue
			#stand_seq = ALN.seq[s][:standardize]  ##Take only the first standardize bp. # Used when upstream sequences are taken for the homolog sets, not actual promoters.
			#print (len(stand_seq))
			#ALN.seq[s]=stand_seq[:]
		if (ACGTs + gaps != len(ALN.seq[s])): continue ##no exotic letters allows
		if (float(gaps)/float(gaps+ACGTs) > max_gap_frac): continue
		goodseqs.append(s)
	return goodseqs


seqs =[]
labels=[]
testseqs=[]
testlabels=[]
i=0
for a in DB.names:
	if (test_str is not None) and (test_str in a): continue
	ALN=DB.data[a]
	#ALN.print_mfa()
	if len(ALN.names)<N: continue

	for s in ALN.names: ##1-hot encoding and evolutionary simulations assume upper case
		ALN.seq[s]=ALN.seq[s].upper().replace("-","")
		#ALN.seq[s]=ALN.seq[s][len(ALN.seq[s])-SEQ_L:] ##taking the last SEQ_L bp. Could not get this to work in a subroutine! pointers?!?!?
	if alignment_QC: sequences_to_use=good_DNA_seqs(ALN)
	else: sequences_to_use=ALN.names ##use allt he sequences
	#print(a)
	if len(sequences_to_use) < N: continue
	if v_str is not None and v_str in a:
		for j in range(0,len(sequences_to_use)):
			testseqs.append([a,sequences_to_use[j]])
			testlabels.append(int(i))		
	else:	
		for j in range(0,len(sequences_to_use)):
			seqs.append([a,sequences_to_use[j]])
			labels.append(int(i))	
	i=i+1
##convert everything to np	
seqs=np.array(seqs)
labels=np.array(labels)
testseqs=np.array(testseqs)
testlabels=np.array(testlabels)

print("done. training data...")
print(np.shape(seqs))
print(np.shape(labels))
print("valedation data...")
print(np.shape(testseqs))
print(np.shape(testlabels))
#exit()
if not SAMPLE_BY_CLASS:
##all sequences in each epoch
	steps = len(seqs) // batch_size
	v_steps = len(testseqs) // batch_size
else:
##all classes in each epoch, N ( = fam_size+4 ) x "coverage" to try to see more data in each epoch 
	steps = N * len(np.unique(labels)) // batch_size
	v_steps = N * len(np.unique(testlabels)) // batch_size


# this will generate a batch for reverse homology
def RH_batch_pos_only(DB,seqs,labels,idx,fam_size=8,batch_size=16,sample_by_class=True,stand=None,REVCOM_AUG=True): ##no comp size. use batch_size
	refs=[]
	positives=[]
	outputs=[]
	positive_classes=[]
	for i in range(0,batch_size):
		#pos_class=labels[idx[i]] ##take the next sequence as the positive class
		if sample_by_class:
			pos_class=idx[i] ##class indices are given directly, clashes are impossible
			all_pos_i = np.squeeze(np.argwhere(labels==pos_class)) ##get the indices for this class
			pos_sample = np.random.choice(all_pos_i,size=fam_size+1,replace=False)
			pos_i= pos_sample[-1] #this will be the one mixed in with the negatives if not using the example given
		else:
			pos_i=idx[i] ##actually use the sequence given as the positive example
			pos_class=labels[idx[i]]
			all_pos_i = np.squeeze(np.argwhere(labels==pos_class)) ##get the indices for this class 
			all_pos_i = np.delete(all_pos_i, np.where(all_pos_i == pos_i))
			pos_sample = np.random.choice(all_pos_i,size=fam_size,replace=False)	
		#if (pos_class in positive_classes): print("clash: "+str(pos_class)+ " is already in"+str(positive_classes))
		positive_classes.append(pos_class) ##with small datasets/validation datasets, we are getting the same class. This could hinder learning.
		#print(str(pos_i)+" is the positive index "+str(labels[idx[i]])+" is the positive class")
		#pos_sample = all_pos_i[0:fam_size+1] #take the first ones
		fam_i=pos_sample[0:fam_size] #these are the reference family
		#pos_i=pos_sample[0] ##positive control: should always be classified correctly
		positives.append(pos_i) ##store the positive 
		#print(pos_i, " was the positive index",fam_i, " were the references from ",all_pos_i)
		these_refs = []
		for j in fam_i:
			dna=DB.data[seqs[j][0]].seq[seqs[j][1]]
			#print(seqs[j][0],seqs[j][1],dna)
			one_hot_dna=my_encoders.dna_one_hot2(dna,standardize=stand)
			if REVCOM_AUG:
				revcom_this_seq = np.random.choice([True, False])
				if revcom_this_seq: one_hot_dna = tf_tf.single_seq_revcom_one_hot(one_hot_dna)
			#print(one_hot_dna)
			#print(np.shape(one_hot_dna))
			these_refs.append(one_hot_dna) ##get the sequences for the references
		#print(np.shape(these_refs))	
		refs.append(np.array(these_refs))
	#outputs = np.eye(batch_size)
	#print(np.shape(np.array(refs)))
	#loss_w = np.ones((batch_size,batch_size)) ##set clashes to weight =0 so that they are not included in the loss	
	
	positive_classes=np.array(positive_classes)
	outputs = np.equal.outer(positive_classes,positive_classes) ##allow multiple positive examples in each row. Let keras cce deal with it
	outputs = outputs.astype(float)
	outputs /= np.sum(outputs,axis=-1)
	
	reorder=np.random.choice(range(0,batch_size),size=batch_size,replace=False)
	outputs= np.take(outputs,reorder,axis=1) #reorder the labels for each row
	##sample_weights = np.take(loss_w,reorder,axis=1)

	#print(sample_weights)
	#u, indices = np.unique(positive_classes, return_index=True)
	#sample_weights[indices]=1.0 ##just set the ones with clashes to be 0
	#print(np.shape(sample_weights),np.shape(outputs))	
	positives = [positives[j] for j in reorder]
	pseqs=[]
	for j in positives:
		#print ("posisitves:")
		dna=DB.data[seqs[j][0]].seq[seqs[j][1]]
		#print(seqs[j][0],seqs[j][1],dna)
		one_hot_dna=my_encoders.dna_one_hot2(dna,standardize=stand)
		if REVCOM_AUG:
			revcom_this_seq = np.random.choice([True, False])
			if revcom_this_seq: one_hot_dna = tf_tf.single_seq_revcom_one_hot(one_hot_dna)
		#print(one_hot_dna)
		#print(np.shape(one_hot_dna))
		pseqs.append(one_hot_dna)
	inputs = [np.array(refs), np.array(pseqs)]
	return inputs, outputs
	#return inputs, outputs, sample_weights

#this is the generator that will generate endless batches
def RH_generator(DB,seqs,labels,fam_size=8,comp_size=2, batch_size=16, reuse_neg=False, sample_by_class=True,stand=False): ##default is triplet
	if comp_size is None: reuse_neg=True ##batch size is the comp size if reusing negatives.
	##each epoch should cycle through the "classes" or the sequences?? ##if we do classes, we won't use all of our data evenly
	if not sample_by_class:
		i_list = range(0,len(seqs)) ## if we do sequences, we see the same 'class' over and over in each epoch
	else: 
		i_list = np.unique(labels) ## if we do classes, we won't see all the sequences
	if len(i_list) <= batch_size: 
		print("only "+str(len(i_list))+" labels found for b="+str(batch_size))
		exit()
	#else:
	#	print("making batches from "+str(len(i_list))+" indices")	
	while True: ##never ending loop
		rand_ord = np.random.choice(i_list,size=len(i_list),replace=False) ##need to randomize if reusing negatives/positives
		i=0
		while i<(len(i_list)-batch_size):
			if reuse_neg:
				yield RH_batch_pos_only(DB,seqs,labels,rand_ord[i:i+batch_size],
					fam_size=fam_size,batch_size=batch_size, sample_by_class=sample_by_class,stand=stand)
			else:
				print("not implemented")
				exit()
				#yield RH_batch(seqs,labels,rand_ord[i:i+batch_size],fam_size=fam_size,comp_size=comp_size, batch_size=batch_size, sample_by_class=sample_by_class)
			i=i+batch_size

##create a batch to get the shape
if SAMPLE_BY_CLASS: sample_i=np.unique(labels)[range(0,batch_size)]
else: sample_i=range(0,batch_size)

bseqs,ans=RH_batch_pos_only(DB,seqs,labels,sample_i,fam_size=fam_size,batch_size=batch_size,
	sample_by_class=SAMPLE_BY_CLASS,stand=SEQ_L)
#print(np.shape(bseqs[0])[1],np.shape(bseqs[1])[1],np.shape(ans))
#exit()
##create the encoder


if encoder_type == 'interpretable':
    enc = my_encoders.DNAencoder(filterN=filter_num, filterL=filter_len, trainable=learn_motifs).interpretable_encoder(seq_len=np.shape(bseqs[1])[1])
else:
	print("only interpretable encoder is available")
	exit()

if SAVE_INIT:
	print("Saving initial encoder weights to "+filename+ "_init_weights.h5")
	enc.save_weights(filename+ "_init_weights.h5",save_format='h5')

##create the full reverse homology model.

def create_reverse_homology_pos_only(encoder, seq_len, fam_size=8, bsize=16, temper=0.1, alphabet=4):
	fam_in = tf.keras.layers.Input(shape=(fam_size, seq_len, alphabet), batch_size=bsize, name="fam")
	pos_in = tf.keras.layers.Input(shape=(seq_len, alphabet), batch_size=bsize, name="pos")
	fam4enc = tf.keras.layers.Lambda(lambda t:
									 tf.keras.backend.reshape(t, (bsize * fam_size, seq_len, alphabet))
									 )(fam_in)

	fam = encoder(fam4enc)
	fam = tf.keras.layers.Lambda(lambda t:
								 tf.keras.backend.reshape(t, (bsize, fam_size, -1))
								 )(fam)

	ref_z = tf.keras.layers.GlobalAveragePooling1D(data_format='channels_last')(
		fam)  ##average features over the input family
	ref_z = tf.math.l2_normalize(ref_z, axis=-1)
	pos_z = encoder(pos_in)
	pos_z = tf.math.l2_normalize(pos_z, axis=-1)

	pos_z = tf.keras.layers.Lambda(lambda t: tf.keras.backend.reshape(t, (1, bsize, -1)))(pos_z)
	pos_z = tf.keras.layers.Lambda(lambda t: tf.repeat(t, repeats=bsize, axis=0))(pos_z)

	preds = tf.keras.layers.Dot(name='pos_dot', axes=-1)([ref_z, pos_z])
	preds = tf.divide(preds, temper)  ##seems models can learn without this, but learn way faster and more with it...
	# print(np.shape(preds)) ##should be batch_size x batch_size
	preds = tf.keras.layers.Softmax(name='logits', axis=-1)(preds)  ##better to send logits to the cce if possible
	return tf.keras.models.Model(inputs=[fam_in, pos_in], outputs=preds)


cmod = create_reverse_homology_pos_only(enc, seq_len=np.shape(bseqs[1])[1], fam_size=fam_size,bsize=batch_size)

cmod.compile(
	optimizer=tf.keras.optimizers.Adam(learning_rate),
	loss=tf.keras.losses.categorical_crossentropy, ##seems to not require only 1s and 0s as the targets, as long as they can be normalized...
	metrics=[tf.keras.metrics.categorical_accuracy], ##uses argmax of the predictions, need to check behavior for multiple positives
)
print(cmod.summary())
print(f"evaluating  validation set with {v_steps} steps")
if v_steps > 0:
    val_gen = RH_generator(DB, testseqs, testlabels, fam_size=fam_size, comp_size=None,
                           batch_size=batch_size, sample_by_class=SAMPLE_BY_CLASS, stand=SEQ_L)
    accuracy = cmod.evaluate(val_gen, steps=v_steps)[1]
    print(f"Initial val accuracy: {round(accuracy * 100, 2)}%")
    callbacks = [tf.keras.callbacks.EarlyStopping(patience=10, monitor="val_loss", restore_best_weights=True)]
    val_kwargs = dict(validation_data=val_gen, validation_steps=v_steps)
else:
    print("No validation set (v_str=None). Skipping evaluation/early stopping.")
    accuracy = float('nan')
    callbacks = []
    val_kwargs = {}

history = cmod.fit(
    RH_generator(DB, seqs, labels, fam_size=fam_size, comp_size=None,
                 batch_size=batch_size, sample_by_class=SAMPLE_BY_CLASS, stand=SEQ_L),
    callbacks=callbacks,
    steps_per_epoch=steps,
    batch_size=batch_size,
    epochs=num_epochs,
    **val_kwargs
)

if SAVE_MOTIF_DB:
	bases="ACGT"
	print("saving matrix layer to "+filename+"_matrix_layer.txt")
	outfile = open(filename + "_matrix_layer.txt", 'w')
	M = enc.get_layer('matrix_scan').get_weights()
	M = np.squeeze(M)
	print("motifs:",np.shape(M))
	if encoder_type=='interpretable':
		int_w = enc.get_layer('motif_interactions').get_weights()[0]
		print("interaction weights:", np.shape(int_w))  # , "self-importance:",np.shape(int_b))
		#int_b = enc.get_layer('motif_interactions').get_weights()[1]
	if encoder_type == 'interpretable' or encoder_type == 'motif':
		pooling_param = enc.get_layer('match_pooling').get_weights()
		pooling_param = np.squeeze(pooling_param)
		print("motif pooling:", np.shape(pooling_param))
		motif_scale = enc.get_layer('motif_scaling').get_weights()
		motif_scale = np.squeeze(motif_scale)
		print("motif scaling and thresholds:", np.shape(motif_scale))
		thr = motif_scale[0, :]
		sca = motif_scale[1, :]

	if learn_motifs:
		motif_N=filter_num
	else:
		motif_N=len(motif_names)
	for i in range(motif_N):
		f = M[:, :, i]
		if learn_motifs:
			m = "matrix"+str(i)
		else:
			m=motif_names[i]
		outfile.write(">"+str(m)+" length="+str(len(f)))
		if encoder_type=='interpretable' or encoder_type == 'motif':
			#if not learn_motifs:
			outfile.write("\tthr="+str(thr[i])+"\tsca="+str(sca[i])+"\tpool="+str(pooling_param[i]))#+"\tself-imp="+str(int_b[i]))
		if encoder_type == 'interpretable':
			outfile.write("\tinteractions:\t")
			outfile.write("\t".join([str(int_w[i,j]) for j in range(motif_N)]))
		outfile.write("\n")
		if learn_motifs:
			for b in range(len(bases)):
				outfile.write(bases[b]+"\t"+"\t".join([str(f[p,b]) for p in range(0,len(f))])+"\n")
	outfile.close()

if SAVE_FINAL:
	print("Saving encoder weights to "+filename+ "_final_weights.h5")
	enc.save_weights(filename+ "_final_weights.h5",save_format='h5')

training_file_name=filename+ "_training_info.csv"
print("Saving training info to "+training_file_name)
outfile=open(training_file_name,'w')
outfile.write(encoder_type+ " with "+str(filter_num)+" conv. filters by " +str(filter_len)+"nts.\n "+str(fam_size)+" homologs by batch size "+str(batch_size)+" as negatives contrastive learning\n")
if SAVE_FINAL: 
	outfile.write("final encoder weights written to "+filename+ "_final_weights.h5"+"\n")
outfile.write(f"Initial val accuracy: {round(accuracy * 100, 2)}%"+"\n")	
for h in history.history.keys():
	outfile.write(h+","+",".join([str(val) for val in history.history[h]])+"\n")
	#print(history.history[h])
outfile.close()
	

	

