###copyright Alan Moses 2021-2025
##permission to use and copy is granted for research purposes only
##software provided as is with no warrantee of any kind

###this is supposed to be 1.0 released with paper about inferring fungal regulatory networks
###this code converts sequences into numerical representations.

###inputs: model hd5 and a directory of fasta files of orthologous promoters.
### gaps will be removed, but otherwise all non-ACGT will just be ignored

### outputs: a csv file where rows are fasta files and columsn are signals for each motif

import sys
import numpy as np
import mol_evol as me ## my molecular evolution functions
import tf_tf ## transcription factors in tensorflow
import encoders as my_encoders

filename="rep.csv" #name for the output file

##need to say how we are dealing with the input files
alignment_QC=True
alignment_dir=True
FILE_EXT="fa" ##if reading in files from a directory

SEQ_L=500 ##does not need to match what the interpretable encoder was trained with.
#SEQ_L=800

###encoder hyperparameters need to match the ones used to train the model in DNA_RH. ###
learn_motifs=True ##this needs to be set.
trained_weights=True #set to False to get results for a randomly initialized model
allowed_encoders=['interpretable']

##these are defaults, but we will try to parse them from the name of the h5 file
batch_size = 64 ##this is also the comp_size+1 if reusing negatives.
encoder_type='interpretable'##'bestmatch'##'transformer' ##'cnn' ##
filter_num=256
filter_len=15 ##deepstarr is 7



if (len(sys.argv)<3):
	print ("usage: weights_file.h5, sequence path")
	quit()	

for et in allowed_encoders: ##get the encoder from the filename
	if et in sys.argv[1]:
		encoder_type = et
		print("recognized encoder as:", encoder_type)
if "N" in sys.argv[1] and "L" in sys.argv[1]:
	tempstr = sys.argv[1].replace("N","_")
	tempstr = tempstr.replace("L", "_")
	vals = tempstr.split("_")
	filter_num=int(vals[3])
	filter_len=int(vals[4])
	batch_size=int(vals[6].replace("b",""))
	print("found hyperparameters",filter_num,filter_len,batch_size)


if encoder_type == 'interpretable':
	enc = my_encoders.DNAencoder(filterN=filter_num, filterL=filter_len, trainable=learn_motifs).interpretable_encoder(seq_len=SEQ_L)
else: 
	print("not available")

if trained_weights:
	print("initializing "+encoder_type+" from "+sys.argv[1])
	enc.load_weights(sys.argv[1])
REP = {}
if alignment_dir:
	print("reading alignment data...")
	DB=me.AlignmentDB() # intialize the data structure to hold alignments
	DB.read_from_fa_dir(sys.argv[2],ext=FILE_EXT) ##read in all the fasta files
	print(str(len(DB.names)) + " alignments loaded. encoding ...")
	for a in DB.names:
		if len(DB.data[a].names)>1:
			for s in DB.data[a].names: ##"promoter" sequences are taken "upstream" so that the start codon follows the end.
				DB.data[a].seq[s] = DB.data[a].seq[s].replace("-", "") ##need to make sure the sequences are unaligned
				if len(DB.data[a].seq[s])>SEQ_L: DB.data[a].seq[s] = DB.data[a].seq[s][-SEQ_L:] #avoid taking the first SEQ_L bases
			REP[a] = np.mean(enc.predict_on_batch(tf_tf.fasta2one_hot(DB.data[a],standardize=SEQ_L)),axis=0)
	dim=len(REP[DB.names[0]])	
	print("encoded "+str(len(REP))+ " alignments")
else:
	ALN=me.Alignment()
	ALN.read_mfa(sys.argv[2])
	print("read file of "+str(len(ALN.names))+" seqs")
	for s in ALN.names:
		ALN.seq[s]= ALN.seq[s].replace("-", "")
		if len(ALN.seq[s]) > SEQ_L: ALN.seq[s] = ALN.seq[s][-SEQ_L:]
	seqs=tf_tf.fasta2one_hot(ALN,standardize=SEQ_L)
	print("converted to one-hot")
	#quit()
	print(np.shape(seqs))
	for i in range(0,len(seqs)-batch_size,batch_size):
		z=enc.predict_on_batch(seqs[i:i+batch_size])
		#print (np.shape(z))
		#print([i,ALN.names[i],np.shape(z[i])])
		for j in range(0,batch_size):
			#print(i,j,ALN.names[i+j])
			REP[ALN.names[i+j]]=z[j].numpy()
	lastz = enc(seqs[i:len(seqs)])
	for j in range(0,len(lastz)):
		REP[ALN.names[i+j]]=lastz[j].numpy()
	print("encoded "+str(len(REP))+ " seqs")		
	dim=len(REP[ALN.names[0]])	
print (str(dim)+ " dimensional representation created. writing to "+filename)

outfile=open(filename,'w')
outfile.write("FILE,"+",".join([str(i) for i in range(0,dim)])+"\n")
for a in REP:
	#outfile.write(a+" n="+str(len(DB.data[a].names))+","+",".join([str(REP[a][i]) for i in range (0,len(REP[a]))]))
	#print(a,np.shape(REP[a]))
	outfile.write(a+","+",".join([str(REP[a][i]) for i in range (0,len(REP[a]))]))
	outfile.write("\n")
outfile.close()		
