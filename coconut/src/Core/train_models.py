
################################################################################
#   Copyright (C) 2016-2024 Neelesh Soni <neeleshsoni03@gmail.com>,
#   <neelesh.soni@alumni.iiserpune.ac.in>
#
#   This library is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with this library.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

from collections import OrderedDict
from numpy import array, newaxis, exp, sort, median, mean
from scipy import stats
import numpy as np
from sklearn.neighbors import KernelDensity
import sys
import joblib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import chi2

matplotlib.rcParams['font.sans-serif']='Arial'
matplotlib.rcParams['font.size']=14


branching_dict={'A':1,'C':3.5,'D':3,'E':2,'F':3,'G':0,'H':3,'I':3.25,'K':2,'L':2.5,'M':2.5,'N':3,'P':4,'Q':2,'R':2,'S':2,'T':3,'V':3,'W':3,'Y':3};

def get_BS(value):
	"""Summary
	
	Args:
	    value (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	return branching_dict[value]

def check_heptad_repeat(hep):
	"""Summary
	
	Args:
	    hep (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	flag=True
	ph=hep[0];
	for h in hep[1:]:
		ph_n = ((ord(ph)-97)+1)%7
		h_n = ord(h)-97
		flag = flag * (ph_n==h_n)
		ph=h
	return flag

def read_cc_file(fname, redundancy_cutoff, alignment, skip_redundancy_check=False):
	"""Summary
	
	Args:
	    fname (TYPE): Description
	    redundancy_cutoff (TYPE): Description
	    alignment (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	dimers=OrderedDict()
	inf=open(fname,'r')
	lines=inf.readlines();
	inf.close();

	for l in lines:
		toks=l.strip().split('\t')
		key=toks[0].split('_')
		dimerkey='_'.join(key[0:2])
		chain=toks[1].split(':')[0]

		if skip_redundancy_check==False:
			if redundancy_cutoff >= float(toks[4][:-1]):

				if dimerkey not in dimers.keys():
					dimers[dimerkey]=[alignment]
					dimers[dimerkey].append([toks[2],toks[3]])

				else:
					#Skip for database cases where multiple keys are present for a dimer
					if len(dimers[dimerkey])==3:
						continue
					dimers[dimerkey].append([toks[2],toks[3]])
		else:
			
			if dimerkey not in dimers.keys():
				dimers[dimerkey]=[alignment]
				dimers[dimerkey].append([toks[2],toks[3]])

			else:
				#Skip for database cases where multiple keys are present for a dimer
				if len(dimers[dimerkey])==3:
					continue
				dimers[dimerkey].append([toks[2],toks[3]])

	#anti-parallel and parallel dimers that dont have the length of both the sequence identical were skipped
	#dimer_dict_new = filter_dimer(dimers)
	return dimers

def filter_dimer(dimers, min_length_cutoff,max_length_cutoff):
	"""Summary
	
	Args:
	    dimers (TYPE): Description
	    min_length_cutoff (TYPE): Description
	    max_length_cutoff (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	dimer_dict_new=OrderedDict()

	for k,v in dimers.items():
			
		#minimum and maximum length cutoff
		if ( (len(v[1][0])<min_length_cutoff) or (len(v[1][0])>max_length_cutoff) ):
			continue
		elif ( (len(v[2][0])<min_length_cutoff) or (len(v[2][0])>max_length_cutoff) ):
			continue

		#skip if the two dimer sequence is not same. This could lead to wrong heptad repeat
		if len(v[1][0])!=len(v[2][0]):
			#print(k,"sequence length dont match. Skipping...")
			continue

		#Check if "X" is present in seq
		Hep1 = v[1][1];
		Hep2 = v[2][1];
		adeg_indexes_H1 = [n for n,x in enumerate(Hep1) if (x in ['a','d','e','g'])];
		adeg_indexes_H2 = [n for n,x in enumerate(Hep2) if (x in ['a','d','e','g'])];
		seq1 = array(list(v[1][0]))[adeg_indexes_H1];
		seq2 = array(list(v[2][0]))[adeg_indexes_H2];
		
		if ( ('X' in seq1) or ('X' in seq2) ):
			continue

		#Both the heptad repeats shouldn't be discontinuos
		if check_heptad_repeat(Hep1)*check_heptad_repeat(Hep2)==0:
			continue

		#Added to make a parallel dimer sequence compatible for anti-parallel dimer sequence or vice versa
		#If dimer is parallel; start with either a/d and end with d/a
		#If dimer is anti-parallel; start with either a/d and end with d/a

		if v[0]==1:

			if (v[1][1][0] == 'a'):
				if (v[1][1][-1] == 'd'):
					pass
				else:
					#seq1 starts/ends a-a
					#chop last four heptad (e, f, g, a)
					v[1][1] = v[1][1][:-4];
					#chop last four residues 
					v[1][0] = v[1][0][:-4];

					#seq2 starts/ends a-a
					#chop last four heptad (e, f, g, a)
					v[2][1] = v[2][1][:-4];
					#chop last four residues 
					v[2][0] = v[2][0][:-4];

			elif (v[1][1][0] == 'd'):
				if (v[1][1][-1] == 'a'):
					pass
				else:
					#seq1 starts/ends d-d
					#chop last three heptad (b,c,d)
					v[1][1] = v[1][1][:-3];
					#chop last three residues 
					v[1][0] = v[1][0][:-3];

					#seq2 starts/ends d-d
					#chop last three heptad (b,c,d)
					v[2][1] = v[2][1][:-3];
					#chop last three residues 
					v[2][0] = v[2][0][:-3];

			
		else:
			#Since the input sequences for anti-parallel dimers in the database files are already reversed. 
			#No need to change for parallel or anti-parllel. same condition holds true
			#Except the cropping site

			if (v[1][1][0] == 'a'):
				if (v[1][1][-1] == 'd'):
					pass
				else:
					#seq1 starts/ends a-a
					#chop last four heptad (e, f, g, a)
					v[1][1] = v[1][1][:-4];
					#chop last four residues 
					v[1][0] = v[1][0][:-4];

					#seq2 starts/ends d-d
					#chop first four heptad (d, e, f, g)
					v[2][1] = v[2][1][4:];
					#chop last four residues 
					v[2][0] = v[2][0][4:];

			elif (v[1][1][0] == 'd'):
				if (v[1][1][-1] == 'a'):
					pass
				else:
					#seq1 starts/ends d-d
					#chop last three heptad (b,c,d)
					v[1][1] = v[1][1][:-3];
					#chop last three residues 
					v[1][0] = v[1][0][:-3];

					#seq2 starts/ends a-a
					#chop last three heptad (a, b,c)
					v[2][1] = v[2][1][3:];
					#chop last three residues 
					v[2][0] = v[2][0][3:];

		dimer_dict_new[k]=v;

	return dimer_dict_new

def load_training_data(CCDDIR):
	"""Summary
	
	Args:
	    CCDDIR (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	redundancy_cutoff = 35
	min_length_cutoff = 11
	max_length_cutoff = 100
	
	#####
	# PARALLEL DIMERS
	#####
	alignment = 1;	#1 for parallel, -1 for anti-parallel

	P_HOMO11SC= read_cc_file(CCDDIR+'/parallel_dimers_redundant_homo_SC.txt',redundancy_cutoff,alignment)
	P_HOMO11DC= read_cc_file(CCDDIR+'/parallel_dimers_redundant_homo_DC.txt',redundancy_cutoff,alignment)
	P_HETE11SC= read_cc_file(CCDDIR+'/parallel_dimers_redundant_hete_SC.txt',redundancy_cutoff,alignment)
	P_HETE11DC= read_cc_file(CCDDIR+'/parallel_dimers_redundant_hete_DC.txt',redundancy_cutoff,alignment)

	#####
	# ANTIPARALLEL DIMERS
	#####
	alignment = -1;	#1 for parallel, -1 for anti-parallel

	AP_HOMO11SC= read_cc_file(CCDDIR+'/antiparallel_dimers_redundant_homo_SC.txt',redundancy_cutoff,alignment)
	AP_HOMO11DC= read_cc_file(CCDDIR+'/antiparallel_dimers_redundant_homo_DC.txt',redundancy_cutoff,alignment)
	AP_HETE11SC= read_cc_file(CCDDIR+'/antiparallel_dimers_redundant_hete_SC.txt',redundancy_cutoff,alignment)
	AP_HETE11DC= read_cc_file(CCDDIR+'/antiparallel_dimers_redundant_hete_DC.txt',redundancy_cutoff,alignment)

	P_DIMERS = {**P_HOMO11SC,**P_HOMO11DC, **P_HETE11SC, **P_HETE11DC}
	AP_DIMERS = {**AP_HOMO11SC,**AP_HOMO11DC, **AP_HETE11SC, **AP_HETE11DC}

	print("Parallel dimers:",len(P_DIMERS),"Anti-Parallel dimers:",len(AP_DIMERS))

	print("Raw Dimer Stats:")
	print(len(P_DIMERS.items())+len(AP_DIMERS.items()),len(P_DIMERS.items()),len(AP_DIMERS.items()))

	P_DIMERS = filter_dimer(P_DIMERS, min_length_cutoff, max_length_cutoff);
	AP_DIMERS = filter_dimer(AP_DIMERS, min_length_cutoff, max_length_cutoff);

	DIMERS = {**P_DIMERS, **AP_DIMERS}
	print("Filtered Dimer Stats:")
	print(len(DIMERS.items()),len(P_DIMERS.items()),len(AP_DIMERS.items()))

	HOMO_DIMERS = {**P_HOMO11SC,**P_HOMO11DC,**AP_HOMO11SC,**AP_HOMO11DC}
	HETERO_DIMERS = {**P_HETE11SC, **P_HETE11DC, **AP_HETE11SC, **AP_HETE11DC}
	
	print("Homodimers:",len(HOMO_DIMERS))
	print("Heterodimers:",len(HETERO_DIMERS))

	HOMO_DIMERS = filter_dimer(HOMO_DIMERS, min_length_cutoff, max_length_cutoff);
	HETERO_DIMERS = filter_dimer(HETERO_DIMERS, min_length_cutoff, max_length_cutoff);

	print("Homodimers Filtered:",len(HOMO_DIMERS))
	print("Heterodimers Filtered:",len(HETERO_DIMERS))

	return DIMERS


def load_training_data_singlefile(CCDDIR, skip_filtering=False):
	"""Summary
	
	Args:
	    CCDDIR (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	redundancy_cutoff = 35
	min_length_cutoff = 11
	max_length_cutoff = 100
	
	#####
	# PARALLEL DIMERS
	#####
	alignment = 1;	#1 for parallel, -1 for anti-parallel

	skip_redundancy_check=False
	skip_length_filtering=False

	if skip_filtering==True:
		skip_redundancy_check=True

	P_DIMERS = read_cc_file(CCDDIR+'/parallel_dimers_all.out',redundancy_cutoff,alignment, skip_redundancy_check)

	#####
	# ANTIPARALLEL DIMERS
	#####
	alignment = -1;	#1 for parallel, -1 for anti-parallel

	AP_DIMERS = read_cc_file(CCDDIR+'/antiparallel_dimers_all.out',redundancy_cutoff,alignment, skip_redundancy_check)

	#-------------------------------------

	print("Parallel dimers:",len(P_DIMERS),"Anti-Parallel dimers:",len(AP_DIMERS))

	print("Raw Dimer Stats:")
	print(len(P_DIMERS.items())+len(AP_DIMERS.items()),len(P_DIMERS.items()),len(AP_DIMERS.items()))

	P_DIMERS = filter_dimer(P_DIMERS, min_length_cutoff, max_length_cutoff);
	AP_DIMERS = filter_dimer(AP_DIMERS, min_length_cutoff, max_length_cutoff);

	DIMERS = {**P_DIMERS, **AP_DIMERS}
	print("Filtered Dimer Stats:")
	print(len(DIMERS.items()),len(P_DIMERS.items()),len(AP_DIMERS.items()))

	return DIMERS

def generate_dodec_pairs(MODDIR,dimers):
	"""Summary
	
	Args:
	    MODDIR (TYPE): Description
	    dimers (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	Pada_wlines,Pdad_wlines,APada_wlines,APdad_wlines="","","",""

	for k,dimer in dimers.items():
		
		#Parallel dimer
		for hx in dimer[3]:

			C1=hx[3];C2=hx[4];E1=hx[7];E2=hx[8];
			sC1=hx[1];sC2=hx[2];sE1=hx[5];sE2=hx[5];
				
			#ada and geg; gadega-gadega
			if hx[3]=='ada':
				dodec1 = E1[0]+C1[0]+C1[1]+E1[1]+E1[2]+C1[2]
				dodec2 = E2[0]+C2[0]+C2[1]+E2[1]+E2[2]+C2[2]
				sdodec1 = sE1[0]+sC1[0]+sC1[1]+sE1[1]+sE1[2]+sC1[2]
				sdodec2 = sE2[0]+sC2[0]+sC2[1]+sE2[1]+sE2[2]+sC2[2]

				Pada_wlines+=k+'\t'+sdodec1+'\t'+sdodec2+'\t'+dodec1+'\t'+dodec2+'\n'
			
			#dad and ege; degade-degade
			else:
				dodec1 = C1[0]+E1[0]+E1[1]+C1[1]+C1[2]+E1[2]
				dodec2 = C2[0]+E2[0]+E2[1]+C2[1]+C2[2]+E2[2]
				sdodec1 = sC1[0]+sE1[0]+sE1[1]+sC1[1]+sC1[2]+sE1[2]
				sdodec2 = sC2[0]+sE2[0]+sE2[1]+sC2[1]+sC2[2]+sE2[2]

				Pdad_wlines+=k+'\t'+sdodec1+'\t'+sdodec2+'\t'+dodec1+'\t'+dodec2+'\n'
		

		for hx in dimer[5]:

			C1=hx[3];C2=hx[4];E1=hx[7];E2=hx[8];
			sC1=hx[1];sC2=hx[2];sE1=hx[5];sE2=hx[5];
				
			#ada-dad and geg-ege; gadega-edaged
			if hx[3]=='ada':
				dodec1 = E1[0]+C1[0]+C1[1]+E1[1]+E1[2]+C1[2]
				dodec2 = E2[0]+C2[0]+C2[1]+E2[1]+E2[2]+C2[2]
				sdodec1 = sE1[0]+sC1[0]+sC1[1]+sE1[1]+sE1[2]+sC1[2]
				sdodec2 = sE2[0]+sC2[0]+sC2[1]+sE2[1]+sE2[2]+sC2[2]

				APada_wlines+=k+'\t'+sdodec1+'\t'+sdodec2+'\t'+dodec1+'\t'+dodec2+'\n'
			
			#dad-ada and ege-geg; degade-agedag
			else:
				dodec1 = C1[0]+E1[0]+E1[1]+C1[1]+C1[2]+E1[2]
				dodec2 = C2[0]+E2[0]+E2[1]+C2[1]+C2[2]+E2[2]
				sdodec1 = sC1[0]+sE1[0]+sE1[1]+sC1[1]+sC1[2]+sE1[2]
				sdodec2 = sC2[0]+sE2[0]+sE2[1]+sC2[1]+sC2[2]+sE2[2]

				APdad_wlines+=k+'\t'+sdodec1+'\t'+sdodec2+'\t'+dodec1+'\t'+dodec2+'\n'

	inf = open(MODDIR+'/P_ada.out','w')
	inf.writelines(Pada_wlines)
	inf.close()

	inf = open(MODDIR+'/P_dad.out','w')
	inf.writelines(Pdad_wlines)
	inf.close()

	inf = open(MODDIR+'/AP_ada.out','w')
	inf.writelines(APada_wlines)
	inf.close()

	inf = open(MODDIR+'/AP_dad.out','w')
	inf.writelines(APdad_wlines)
	inf.close()

	return

def dump_frequency_models(DISTDIR, Res_stats):
	"""Summary
	
	Args:
	    DISTDIR (TYPE): Description
	    Res_stats (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	write_dict(Res_stats[0],DISTDIR+'/pairwise_aa_residues.out')
	write_dict(Res_stats[1],DISTDIR+'/pairwise_dd_residues.out')
	write_dict(Res_stats[2],DISTDIR+'/pairwise_ad_residues.out')

	write_dict(Res_stats[3],DISTDIR+'/pairwise_ag_residues.out')
	write_dict(Res_stats[4],DISTDIR+'/pairwise_de_residues.out')
	write_dict(Res_stats[5],DISTDIR+'/pairwise_ae_residues.out')
	write_dict(Res_stats[6],DISTDIR+'/pairwise_dg_residues.out')

	write_dict(Res_stats[7],DISTDIR+'/pairwise_ee_residues.out')
	write_dict(Res_stats[8],DISTDIR+'/pairwise_gg_residues.out')
	write_dict(Res_stats[9],DISTDIR+'/pairwise_eg_residues.out')

	return

def write_dict(Res_pairs,fname):
	"""Summary
	
	Args:
	    Res_pairs (TYPE): Description
	    fname (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	alines="";
	for key,value in Res_pairs.items():
		alines+=str(key[0])+" "+str(key[1])+" "+str(value)+"\n";

	outf=open(fname,'w')
	outf.writelines(alines);
	outf.close()

	return

def read_dodec_pairs(fname):
	"""Summary
	
	Args:
	    fname (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	inf=open(fname,'r');
	lines=inf.readlines();
	inf.close()

	wlines=""
	dodecp={}

	for l in lines:
		if l=='#':
			continue
		toks = l.strip().split('\t')

		S1=toks[1]; S2=toks[2]; H1=toks[3]; H2=toks[4];
		
		#Storing only unique keys with respect to amino acids
		if ((H1=='gadega') and (H2=='gadega')):
			AA1 = [ S1[1],S1[2],S1[5] ]
			AA2 = [ S2[1],S2[2],S2[5] ]
			#addS = array(list(map(get_BS,AA1)))+array(list(map(get_BS,AA2)))

			AA3 = [ S1[0],S1[3],S1[4] ]
			AA4 = [ S2[0],S2[3],S2[4] ]
			#sbscore=find_saltbridge_score(AA3,AA4)
			
			key1=''.join(S1)+'-'+''.join(S2); key2=''.join(S2)+'-'+''.join(S1);
			if ((key1 not in dodecp.keys()) and (key2 not in dodecp.keys())):
				#density = addS[1] - (addS[0]+addS[2])/2.0
				dodecp[key1]=[H1, H2, AA1, AA2]

		if ((H1=='degade') and (H2=='degade')):
			AA1 = [ S1[0],S1[3],S1[4] ]
			AA2 = [ S2[0],S2[3],S2[4] ]
			#addS = array(list(map(get_BS,AA1)))+array(list(map(get_BS,AA2)))

			AA3 = [ S1[1],S1[2],S1[5] ]
			AA4 = [ S2[1],S2[2],S2[5] ]
			#sbscore=find_saltbridge_score(AA3,AA4)
			
			key1=''.join(S1)+'-'+''.join(S2); key2=''.join(S2)+'-'+''.join(S1);
			if ((key1 not in dodecp.keys()) and (key2 not in dodecp.keys())):
				#density = addS[1] - (addS[0]+addS[2])/2.0
				dodecp[key1]=[H1, H2, AA1, AA2]

		if ((H1=='gadega') and (H2=='edaged')):
			AA1 = [ S1[1],S1[2],S1[5] ]
			AA2 = [ S2[1],S2[2],S2[5] ]
			#addS = array(list(map(get_BS,AA1)))+array(list(map(get_BS,AA2)))
			
			AA3 = [ S1[0],S1[3],S1[4] ]
			AA4 = [ S2[0],S2[3],S2[4] ]
			#sbscore=find_saltbridge_score(AA3,AA4)

			key1=''.join(S1)+'-'+''.join(S2); key2=''.join(S2)+'-'+''.join(S1);
			if ((key1 not in dodecp.keys()) and (key2 not in dodecp.keys())):
				#density = addS[1] - (addS[0]+addS[2])/2.0
				dodecp[key1]=[H1, H2, AA1, AA2]

		if ((H1=='degade') and (H2=='agedag')):
			AA1 = [ S1[0],S1[3],S1[4] ]
			AA2 = [ S2[0],S2[3],S2[4] ]
			#addS = array(list(map(get_BS,AA1)))+array(list(map(get_BS,AA2)))

			AA3 = [ S1[1],S1[2],S1[5] ]
			AA4 = [ S2[1],S2[2],S2[5] ]
			#sbscore=find_saltbridge_score(AA3,AA4)
			
			key1=''.join(S1)+'-'+''.join(S2); key2=''.join(S2)+'-'+''.join(S1);
			if ((key1 not in dodecp.keys()) and (key2 not in dodecp.keys())):
				#density = addS[1] - (addS[0]+addS[2])/2.0
				dodecp[key1]=[H1, H2, AA1, AA2]

	return dodecp

def calculate_feature_values(dodecp):
	"""Summary
	
	Args:
	    dodecp (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	dodecF = {}
	for key, value in dodecp.items():

		#print(key,value)
		
		H1, H2, AA1, AA2 = value[0],value[1],value[2],value[3]

		addS = array(list(map(get_BS,AA1)))+array(list(map(get_BS,AA2)))
		density = addS[1] - (addS[0]+addS[2])/2.0;
		#sbscore=find_saltbridge_score(AA3,AA4)
		dodecF[key]=[density];

	return dodecF

def bin_values(values, num_bins):
    """Summary
    
    Args:
        values (TYPE): Description
        num_bins (TYPE): Description
    
    Returns:
        TYPE: Description
    """
    # Sort the values
    sorted_values = sorted(values)
    
    # Determine bin size
    bin_size = len(sorted_values) // num_bins

    # Initialize empty bins
    bins = []

    # Iterate over the sorted values and distribute them into bins
    for i in range(0, len(sorted_values), bin_size):
        bin_contents = sorted_values[i:i+bin_size]
        bins.append(bin_contents)

    return bins

def check_distribution_pval(pval, CI):
	"""Summary
	
	Args:
	    pval (TYPE): Description
	    CI (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	if pval.pvalue<CI:
		print("Reject Null Hypothesis. Different distribution")
	else:
		print("Can't Reject Null Hypothesis. Identical distribution")
	return

def get_kerneldensity(atomdensity, fcolor=None, PlotData=True):
	"""Summary
	
	Args:
	    atomdensity (TYPE): Description
	    fcolor (None, optional): Description
	    PlotData (bool, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	atomdensity = sort(atomdensity)
	atomdensity = atomdensity[:,newaxis]
	
	# Gaussian KDE
	#kde = KernelDensity(kernel='gaussian', bandwidth=0.18).fit(atomdensity)
	kde = KernelDensity(kernel='gaussian', bandwidth='scott').fit(atomdensity)
	
	log_dens = kde.score_samples(atomdensity)
	
	if PlotData==True:

		plt.fill(atomdensity[:, 0], exp(log_dens), fc=fcolor, alpha=0.4)
		plt.plot(atomdensity[:, 0], exp(log_dens), color='gray', linewidth=0.8)

	return kde, exp(log_dens), atomdensity[:,0]

def train_kde_models(DODEC_DIR, ESTM_DIR, PLOT_DIR, PlotData=True):
	"""Summary
	
	Args:
	    DODEC_DIR (TYPE): Description
	    ESTM_DIR (TYPE): Description
	    PLOT_DIR (TYPE): Description
	    PlotData (bool, optional): Description
	
	Returns:
	    TYPE: Description
	"""
	dodecp1 = read_dodec_pairs(DODEC_DIR+'/AP_ada.out')
	dodecp2 = read_dodec_pairs(DODEC_DIR+'/AP_dad.out')
	
	fig, ax1 = plt.subplots(1,1,figsize=(5,4))

	#FOR plotting purpose
	dodecF1 = calculate_feature_values(dodecp1)
	atomdensity = array(list(dodecF1.values()))[:,0]
	kde, log_dens_ap_ada, atomdens_ap_ada  = get_kerneldensity(atomdensity, fcolor='blue', PlotData=True)

	dodecF1 = calculate_feature_values(dodecp2)
	atomdensity = array(list(dodecF1.values()))[:,0]
	kde, log_dens_ap_dad, atomdens_ap_dad  = get_kerneldensity(atomdensity, fcolor='turquoise', PlotData=True)
	

	dodec_ap = {**dodecp1,**dodecp2}
	dodecF1 = calculate_feature_values(dodec_ap)
	atomdensity = array(list(dodecF1.values()))[:,0]
	kde, log_dens, atomdens = get_kerneldensity(atomdensity, PlotData=False)
	joblib.dump(kde, ESTM_DIR+'/KDE_atomicdensity_AP_ada.joblib')
	joblib.dump(kde, ESTM_DIR+'/KDE_atomicdensity_AP_dad.joblib')

	dodecp3 = read_dodec_pairs(DODEC_DIR+'/P_ada.out')
	dodecF3 = calculate_feature_values(dodecp3)
	atomdensity = array(list(dodecF3.values()))[:,0]
	kde, log_dens_p_ada , atomdens_p_ada = get_kerneldensity(atomdensity, fcolor='olivedrab', PlotData=True)
	joblib.dump(kde, ESTM_DIR+'/KDE_atomicdensity_P_ada.joblib')
	
	dodecp4 = read_dodec_pairs(DODEC_DIR+'/P_dad.out')
	dodecF4 = calculate_feature_values(dodecp4)
	atomdensity = array(list(dodecF4.values()))[:,0]
	kde, log_dens_p_dad, atomdens_p_dad = get_kerneldensity(atomdensity, fcolor='lightsalmon', PlotData=True)
	joblib.dump(kde, ESTM_DIR+'/KDE_atomicdensity_P_dad.joblib')

	if PlotData==True:
		blue_patch = mpatches.Patch(color='cornflowerblue',alpha=0.9, label='Anti-P ada')
		orange_patch = mpatches.Patch(color='turquoise',alpha=0.4, label='Anti-P dad')
		green_patch = mpatches.Patch(color='olivedrab',alpha=0.4, label='P ada')
		red_patch = mpatches.Patch(color='lightsalmon',alpha=0.4, label='P dad')

		plt.legend(handles=[blue_patch,orange_patch,green_patch,red_patch])
		#plt.legend(handles=[blue_patch,green_patch,red_patch])
		plt.xlabel('ada/dad Atomic Density')
		plt.ylabel('Density')
		plt.xlim(-4,4)
		plt.title('Coiled-coil Atomic Density Scores')
		plt.savefig(PLOT_DIR+'/AP_P_ada_dad_atomic_density.pdf')
		plt.show()
		plt.close()
	
	pval1 = stats.ttest_ind(log_dens_p_dad, log_dens_p_ada)
	print("\nKS TEST P:ada and P:dad:",pval1.pvalue)
	check_distribution_pval(pval1, 0.05)
	pval2 = stats.ttest_ind(log_dens_ap_dad, log_dens_ap_ada)
	print("\nKS TEST AP:ada and AP:dad:",pval2.pvalue)
	check_distribution_pval(pval2, 0.05)
	pval3 = stats.ttest_ind(log_dens_p_ada, log_dens_ap_ada)
	print("\nKS TEST P:ada and AP:ada:",pval3.pvalue)
	check_distribution_pval(pval3, 0.05)
	pval4 = stats.ttest_ind(log_dens_p_dad, log_dens_ap_dad)
	print("\nKS TEST P:dad and AP:dad:",pval4.pvalue)
	check_distribution_pval(pval4, 0.05)
	pval5 = stats.ttest_ind(log_dens_p_ada, log_dens_ap_dad)
	print("\nKS TEST P:ada and AP:dad:",pval3.pvalue)
	check_distribution_pval(pval5, 0.05)
	pval6 = stats.ttest_ind(log_dens_p_dad, log_dens_ap_ada)
	print("\nKS TEST P:dad and AP:ada:",pval4.pvalue)
	check_distribution_pval(pval6, 0.05)
	
	return



