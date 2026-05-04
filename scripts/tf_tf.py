
import numpy as np
##really should be called tfs in np

def fasta2one_hot(mfa, bases="ACGT", standardize=None):
	###after peter koo's munge.py function, returns a numpy array###
	###expects an mol_evol.Alignment object, sequence order is the order in names###
	###all the characters should be in the list bases, order determines the one hot order?###
	###standardize controls what do with sequences that are not all the same length. default is "None" which does nothing
	###if a number is given, sequences are truncated or 0-padded to the length given. 
	one_hot_data = [] #start with regular python list
	for s in mfa.names:
		seq = mfa.seq[s].upper()
		seq = seq.replace("-","")
		seq_length = len(seq)
		if (standardize is not None): 
			seq_length=int(standardize)
		one_hot_seq = np.zeros((seq_length,len(bases)))
		seq_length=min(seq_length,len(seq))
		for b in range(0,len(bases)):
			index = [j for j in range(0,seq_length) if seq[j] == bases[b]]
			one_hot_seq[index,b] = 1
		one_hot_data.append(one_hot_seq)
	one_hot_data = np.array(one_hot_data)	### convert to numpy array
	return one_hot_data

def seq2one_hot(seq,bases="ACGT",standardize=None):
	seq=seq.upper()
	seq=seq.replace("-","")
	seq_length=len(seq)
	if (standardize is not None): 
		seq_length=int(standardize)
	one_hot_seq = np.zeros((seq_length,len(bases)))
	seq_length=min(seq_length,len(seq))
	for b in range(0,len(bases)):
		index = [j for j in range(0,seq_length) if seq[j] == bases[b]]
		one_hot_seq[index,b] = 1
	return one_hot_seq		

def one_hot2fasta(seqs, bases="ACGT"): ##expects numpy of shape (n,seq_len,len(bases))
	##returns a list of strings
	mfa=[]
	N=seqs.shape[0]
	L=seqs.shape[1]
	A=len(bases)
	for s in range(0,N):
		this_seq=""
		for i in range(0,L):
			for b in range(0,A):
				if seqs[s,i,b]==1: 
					this_seq+=bases[b]
		mfa.append(this_seq)
	return mfa

def revcom_one_hot(seqs,complement=[3,2,1,0]): ### returns the reverse complement of an np matrix. 
							###assumes A,C,G,T is the one hot encoding order by default
							##expects numpy of shape (n,seq_len,4)
	
	N=seqs.shape[0]
	A=seqs.shape[2]
	revcom=np.zeros(seqs.shape)
	rev = np.flip(seqs,axis=1)
	for s in range(0,N):
		for b in range(0,A):
			revcom[s,:,b]=rev[s,:,complement[b]]
	return revcom		

def single_seq_revcom_one_hot(seq, complement=[3, 2, 1, 0]):  ### returns the reverse complement of an np matrix.
###assumes A,C,G,T is the one hot encoding order by default
##expects numpy of shape (seq_len,4)
	A = seq.shape[1]
	revcom = np.zeros(seq.shape)
	rev = np.flip(seq, axis=0)
	for b in range(0, A):
		revcom[:, b] = rev[:, complement[b]]
	return revcom


def read_pwms(path):  ##reading in a flat file of pwms and return a big dict of dicts
	f = open(path, 'r') #open the file
	Jid="None"
	PWMs = {}
	for line in f:
		line=line.rstrip("\n\r")
		if ('>' in line):
			vals = line.split()
			Jid=vals[0].replace(">","").strip()
			PWMs[Jid]={}
		else:
			vals = line.split()
			b = vals[0].strip()
			pos=0
			for c in vals[1:len(vals)]: ##corrected bug that Andrew found
				if not pos in PWMs[Jid]:
					PWMs[Jid][pos]= {}
				PWMs[Jid][pos][b]=float(c)
				pos+=1
	f.close()
	return PWMs			
            
		
def pwm2filter(pwm, k, bases="ACGT", middle=True): ### convert a pwm dict to a numpy conv filter matrix of kernel size k
	mat=np.zeros((k,len(bases)))
	##put the pwm in the middle? we can use padding=same? but positions of matches will be offset.
	offset=0
	if middle: offset = (k - len(pwm)) // 2
	if offset<0: print ("problem with "+str(offset))
	
	for pos in range(0,len(pwm)):
		for b in range(0,len(bases)):
			mat[pos+offset,b]=pwm[pos][bases[b]]
	return mat		

def sample_from_filter(f,N=10,bases="ACGT",alphabet=4): #sample sequences from the conv filter matrix
	motlen=len(f)
	seqs=[]
	for i in range(0,N):
		s=""
		for pos in range(0,motlen):
			pr=f[pos,:]
			pr=pr/np.sum(pr)
			b=np.random.choice(pr,1,p=pr)
			s+=bases[b]
		seqs.append(s)
	return seqs		
	
def pwm_rev_com(pwm,bases="ACGT"): ###returns the reverse complement of a pwm in dict format
	pwmrc={}
	complement={"A":"T","C":"G","G":"C","T":"A"}
	motlen=len(pwm)
	for pos in range(0,motlen):
		pwmrc[pos]={}
		for b in (list(bases)):
			pwmrc[pos][b]=pwm[motlen-pos-1][complement[b]]
	return pwmrc

def matrix_rev_com(mot, complement=[3,2,1,0]): ### returns the reverse complement of an np matrix. ###assumes A,C,G,T is the one hot encoding order by default
	motifrc=np.zeros(np.shape(mot))
	motlen=len(motifrc)
	for pos in range(0,motlen):
		for b in (range(0,4)):
			motifrc[pos,b]=mot[motlen-pos-1,complement[b]]
	return motifrc

def mean_pwm_score_given_bg(pwm,bases="ACGT",bg={"A":0.25,"C":0.25,"G":0.25,"T":0.25}):
	scr = float(0)
	for p in range(0, len(pwm)):
		for b in list(bases):
			scr += bg[b] * pwm[p][b]
	return scr

def mean_pwm_score_given_m(pwm,bases="ACGT",bg={"A":0.25,"C":0.25,"G":0.25,"T":0.25}): ##assumes pwm is natural log and properly normalized
	scr = float(0)
	for p in range(0, len(pwm)):
		for b in list(bases):
			scr += bg[b] * np.exp(pwm[p][b]) * pwm[p][b]
	return scr


def var_pwm_score_given_bg(pwm,bases="ACGT",bg={"A":0.25,"C":0.25,"G":0.25,"T":0.25}):
	var = float(0)
	for p in range(0, len(pwm)):
		posvar = float(0)
		posmean = float(0)
		for a in list(bases):
			posmean += bg[a] * pwm[p][a]
		for a in list(bases):
			posvar += bg[a] * (pwm[p][a] - posmean) ** 2
		var += posvar  ###assumes each position is independent
	return var
