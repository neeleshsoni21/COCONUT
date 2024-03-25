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


from read_CC_properties import generate_Hex_pairs_dimers
from numpy import array, pad, vstack, median, mean, dstack, prod, exp
from collections import OrderedDict
import os

class ScoreCC:

	"""Summary
	
	Attributes:
	    ALLPROBS (TYPE): Description
	    ALLSCORES (list): Description
	    AP_Prior_Prob (float): Description
	    AP_Prior_Prob_dimer (TYPE): Description
	    branching_dict (TYPE): Description
	    D_PREDS (list): Description
	    D_PROBS (list): Description
	    D_PROBS_mean (TYPE): Description
	    Distributions (TYPE): Description
	    P_Prior_Prob (float): Description
	    P_Prior_Prob_dimer (TYPE): Description
	"""
	
	def __init__(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		self.ALLSCORES=[]

		self.ALLPROBS=OrderedDict()
		
		self.branching_dict={'X':1,'A':1,'C':3.5,'D':3,'E':2,'F':3,'G':0,'H':3,'I':3.25,'K':2,'L':2.5,'M':2.5,'N':3,'P':4,'Q':2,'R':2,'S':2,'T':3,'V':3,'W':3,'Y':3};

		self.D_PREDS=[]

		self.D_PROBS=[]

		self.P_Prior_Prob, self.AP_Prior_Prob = None, None
		self.P_Prior_Prob_dimer, self.AP_Prior_Prob_dimer = None, None

		self._SCOR_ROOT = os.path.abspath(  os.path.join(os.path.join( os.path.dirname(__file__), os.pardir), os.pardir))

		self._COCO_CCS = self._SCOR_ROOT+'/src/Models/DatasetScores'

		return

	def get_BS(self, value):
		"""Summary
		
		Args:
		    value (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		return self.branching_dict[value]

	def get_Distributions(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.Distributions

	def Score_Coiled_coil(self, input_obj, model_obj):
		"""Summary
		
		Args:
		    input_obj (TYPE): Description
		    model_obj (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		from Core.test_decoy2 import predict_partner

		#Predictions
		predict_partner(input_obj.interface_dict, input_obj.interface_true_false, model_obj)

		return

	def calculate_prior_probabilities(self, DIMERS, outputdir):
		"""Summary
		
		Args:
		    DIMERS (TYPE): Description
		    outputdir (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		P_Dist = []
		AP_Dist = []

		P_HexPair_Dist = []
		AP_HexPair_Dist = []

		
		for k,d in DIMERS.items():

			if d[0]==-1:
				#Dimer length in residues
				AP_Dist.append(len(d[1][0]))
				#Dimer length in hexpairs
				AP_HexPair_Dist.append(len(d[3]))
				
			else:
				P_Dist.append(len(d[1][0]))
				P_HexPair_Dist.append(len(d[3]))
				

		Total_Hex_Pairs = sum(P_HexPair_Dist)+sum(AP_HexPair_Dist)
		
		#self.P_Prior_Prob = sum(P_HexPair_Dist)/Total_Hex_Pairs
		#self.AP_Prior_Prob = sum(AP_HexPair_Dist)/Total_Hex_Pairs
		self.P_Prior_Prob = 0.5
		self.AP_Prior_Prob = 0.5

		Total_Dimers = sum(P_Dist)+sum(AP_Dist)
		self.P_Prior_Prob_dimer = sum(P_Dist)/Total_Dimers
		self.AP_Prior_Prob_dimer = sum(AP_Dist)/Total_Dimers

		self.Distributions = [ P_Dist, AP_Dist, P_HexPair_Dist, AP_HexPair_Dist ]

		print("Total number of Parallel Hexpairs:",sum(P_HexPair_Dist))
		print("Total number of Anti-parallel Hexpairs:",sum(AP_HexPair_Dist))

		wlines="#Parallel_Prior\tAnti-parallel_Prior\n"
		wlines+=str(self.P_Prior_Prob)+"\t"+str(self.AP_Prior_Prob)+"\n"
		
		outf = open(outputdir+'/Prior_Probability_Dataset.out','w')
		outf.writelines(wlines)
		outf.close()

		return

	def load_prior_probabilities(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		inf = open(self._COCO_CCS+'/Prior_Probability_Dataset.out','r')
		lines=inf.readlines()
		inf.close()

		for l in lines:
			if ((l[0]=='#') or (len(l)==1)):
				continue
			toks = l.strip().split('\t')

			if len(toks)==2:
				self.P_Prior_Prob=float(toks[0])
				self.AP_Prior_Prob=float(toks[1])

		return

	def Generate_CC_Bayesian_Scores(self, DIMERS, model_obj, CCDSDIR):
		"""Summary
		
		Args:
		    DIMERS (TYPE): Description
		    model_obj (TYPE): Description
		    CCDSDIR (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		scoring_dicts = model_obj.get_scoring_dicts()
		Hexpairs = model_obj.get_KDE_Models()

		for key,dimer in DIMERS.items():

			#dscore = self.estimate_dimer_probabilities(scoring_dicts,dimer)
			dscore = self.score_dimer_bayesian(scoring_dicts,Hexpairs,dimer)

			dscore.append(key);
			dscore.append(dimer[0]);

			self.ALLPROBS[key]=dscore;

		wlines=""
		for key,dscore in self.ALLPROBS.items():
			wlines += ' '.join(map(str,dscore))+'\n'

		outf = open(CCDSDIR+'/CC_dataset_dimers_Bayes_probabilities.txt','w')
		outf.writelines(wlines)
		outf.close()

		return

	def score_dimer_bayesian(self,scoring_dicts,Hexpairs, dimer):
		"""Summary
		
		Args:
		    scoring_dicts (TYPE): Description
		    Hexpairs (TYPE): Description
		    dimer (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		AAscores = [scoring_dicts['aa'],scoring_dicts['dd'],scoring_dicts['ad'],scoring_dicts['ae'],scoring_dicts['ag'],scoring_dicts['de'],scoring_dicts['dg']]
		aa_scores, dd_scores, ad_scores,ae_scores, ag_scores, de_scores, dg_scores = AAscores[0], AAscores[1], AAscores[2],AAscores[3], AAscores[4], AAscores[5], AAscores[6]
		[kde_AP_ada,kde_AP_dad,kde_P_ada,kde_P_dad] = Hexpairs

		Aa_scores_a = []; Aa_scores_d = []
		Aa_ada_scores = []; Aa_dad_scores = [];
		ada_score = []; dad_score=[]
		
		#Append scores with padding zeros
		DSCORES_Prob_E_P, DSCORES_Prob_E_AP = 0,0
		DSCORES2 = [0,0];

		#Iterate over all AA pairs
		for aa_pairs in dimer[4]:

			AA1 = aa_pairs[0]; AA2 = aa_pairs[1];

			if ((aa_pairs[2]=='a') and (aa_pairs[3]=='a')):
				Aa_scores_a+= [aa_scores[tuple([AA1,AA2])]]
			elif ((aa_pairs[2]=='d') and (aa_pairs[3]=='d')):
				Aa_scores_d+= [dd_scores[tuple([AA1,AA2])]]

		#Iterate over all hex pairs
		for hexp in dimer[3]:

			AA1 = list(hexp[1]); AA2 = list(hexp[2]); AA3 = list(hexp[5]); AA4 = list(hexp[6]);

			addS1 = array(list(map(self.get_BS,AA1)))
			addS2 = array(list(map(self.get_BS,AA2)))
			
			density1 = addS1[1] - (addS1[0]+addS1[2])/2.0
			density2 = addS2[1] - (addS2[0]+addS2[2])/2.0

			dd = array(density1+density2).reshape(-1,1);

			#Parallel
			if ( (hexp[3]=='ada') and (hexp[4]=='ada') ):
				Aa_ada_scores+=self.find_score_ada_ada(AA1,AA2,AA3,AA4,AAscores);
				ada_score+=[exp(kde_P_ada.score(dd))]
				
			elif ( (hexp[3]=='dad') and (hexp[4]=='dad') ):
				Aa_dad_scores+=self.find_score_dad_dad(AA1,AA2,AA3,AA4,AAscores);
				dad_score+=[exp(kde_P_dad.score(dd))]


		DSCORES3 = []
		DSCORES3+=Aa_scores_a; DSCORES3+=Aa_scores_d;
		DSCORES3+=Aa_ada_scores;DSCORES3+=Aa_dad_scores;
		DSCORES3+=ada_score;DSCORES3+=dad_score
		DSCORES_Prob_E_P = prod(DSCORES3)

		Aa_scores_a = []; Aa_scores_d = []
		Aa_ada_scores = []; Aa_dad_scores = [];
		ada_score = []; dad_score=[]

		#Iterate over all AA pairs
		for aa_pairs in dimer[6]:

			AA1 = aa_pairs[0]; AA2 = aa_pairs[1];

			if ((aa_pairs[2]=='a') and (aa_pairs[3]=='d')):
				Aa_scores_a+= [ad_scores[tuple([AA1,AA2])]]
			elif ((aa_pairs[2]=='d') and (aa_pairs[3]=='a')):
				Aa_scores_d+= [ad_scores[tuple([AA2,AA1])]]

		#Iterate over all hex pairs
		for hexp in dimer[5]:

			AA1 = list(hexp[1]); AA2 = list(hexp[2]); AA3 = list(hexp[5]); AA4 = list(hexp[6]);

			addS1 = array(list(map(self.get_BS,AA1)))
			addS2 = array(list(map(self.get_BS,AA2)))
			
			density1 = addS1[1] - (addS1[0]+addS1[2])/2.0
			density2 = addS2[1] - (addS2[0]+addS2[2])/2.0

			dd = array(density1+density2).reshape(-1,1);
			
			#anti-Parallel
			if ( (hexp[3]=='ada') and (hexp[4]=='dad') ):
				Aa_ada_scores+=self.find_score_ada_dad(AA1,AA2,AA3,AA4,AAscores);
				ada_score+=[exp(kde_AP_ada.score(dd))]
				
			elif ( (hexp[3]=='dad') and (hexp[4]=='ada') ):
				Aa_dad_scores+=self.find_score_dad_ada(AA1,AA2,AA3,AA4,AAscores);
				dad_score+=[exp(kde_AP_dad.score(dd))]
		
		
		DSCORES3 = []
		DSCORES3+=Aa_scores_a;DSCORES3+=Aa_scores_d;
		DSCORES3+=Aa_ada_scores;DSCORES3+=Aa_dad_scores;
		DSCORES3+=ada_score;DSCORES3+=dad_score
		DSCORES_Prob_E_AP = prod(DSCORES3)

		TotalProb = self.P_Prior_Prob*DSCORES_Prob_E_P + self.AP_Prior_Prob*DSCORES_Prob_E_AP
		DSCORES2[0] = self.AP_Prior_Prob*DSCORES_Prob_E_AP/TotalProb
		DSCORES2[1] = self.P_Prior_Prob*DSCORES_Prob_E_P/TotalProb

		return DSCORES2

	def Generate_CC_RF_Scores(self, DIMERS, model_obj, CCDSDIR):
		"""Summary
		
		Args:
		    DIMERS (TYPE): Description
		    model_obj (TYPE): Description
		    CCDSDIR (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		#get the maximum length of the hex pairs. for padding
		max_hexpairs_length=0;
		for key, value in DIMERS.items():
			if max_hexpairs_length < len(value[3]):
				max_hexpairs_length=len(value[3])

		outf = open(CCDSDIR+'/Padding_Length.out','w')
		outf.writelines("#Padding Length in the RandomForest Model\n"+str(max_hexpairs_length))
		outf.close()

		#make default dict for b/c/e/f/g positions that is also a default dict with score 0
		scoring_dicts = model_obj.get_scoring_dicts()
		Hexpairs = model_obj.get_KDE_Models()

		print("max_hexpairs_length_test",max_hexpairs_length)

		for key,dimer in DIMERS.items():
			dscore = self.score_dimer_RF(scoring_dicts,Hexpairs, dimer,max_hexpairs_length)

			dscore.append(key);
			dscore.append(dimer[0])

			self.ALLSCORES.append(dscore);
		
		wlines=""
		for dscore in self.ALLSCORES:
			wlines += ' '.join(map(str,dscore))+'\n'

		outf = open(CCDSDIR+'/CC_dataset_dimers_RF_scores.txt','w')
		outf.writelines(wlines)
		outf.close()

		return

	def score_dimer_RF(self,scoring_dicts,Hexpairs, dimer, max_hexpairs_length):
		"""Summary
		
		Args:
		    scoring_dicts (TYPE): Description
		    Hexpairs (TYPE): Description
		    dimer (TYPE): Description
		    max_hexpairs_length (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		AAscores = [scoring_dicts['aa'],scoring_dicts['dd'],scoring_dicts['ad'],scoring_dicts['ae'],scoring_dicts['ag'],scoring_dicts['de'],scoring_dicts['dg']]
		aa_scores, dd_scores, ad_scores,ae_scores, ag_scores, de_scores, dg_scores = AAscores[0], AAscores[1], AAscores[2],AAscores[3], AAscores[4], AAscores[5], AAscores[6]
		[kde_AP_ada,kde_AP_dad,kde_P_ada,kde_P_dad] = Hexpairs

		Aa_scores_a = []; Aa_scores_d = []
		Aa_ada_scores = []; Aa_dad_scores = [];
		ada_score = []; dad_score=[]
		
		#Append scores with padding zeros
		DSCORES=[];
		DSCORES_Prob_E_P, DSCORES_Prob_E_AP = 0,0
		DSCORES2 = [0,0];

		#Iterate over all AA pairs
		for aa_pairs in dimer[4]:

			AA1 = aa_pairs[0]; AA2 = aa_pairs[1];

			if ((aa_pairs[2]=='a') and (aa_pairs[3]=='a')):
				Aa_scores_a+= [aa_scores[tuple([AA1,AA2])]]
			elif ((aa_pairs[2]=='d') and (aa_pairs[3]=='d')):
				Aa_scores_d+= [dd_scores[tuple([AA1,AA2])]]

		#Iterate over all hex pairs
		for hexp in dimer[3]:

			AA1 = list(hexp[1]); AA2 = list(hexp[2]); AA3 = list(hexp[5]); AA4 = list(hexp[6]);

			addS1 = array(list(map(self.get_BS,AA1)))
			addS2 = array(list(map(self.get_BS,AA2)))
			
			density1 = addS1[1] - (addS1[0]+addS1[2])/2.0
			density2 = addS2[1] - (addS2[0]+addS2[2])/2.0

			dd = array(density1+density2).reshape(-1,1);

			#Parallel
			if ( (hexp[3]=='ada') and (hexp[4]=='ada') ):
				Aa_ada_scores+=self.find_score_ada_ada(AA1,AA2,AA3,AA4,AAscores);
				ada_score+=[exp(kde_P_ada.score(dd))]
				
			elif ( (hexp[3]=='dad') and (hexp[4]=='dad') ):
				Aa_dad_scores+=self.find_score_dad_dad(AA1,AA2,AA3,AA4,AAscores);
				dad_score+=[exp(kde_P_dad.score(dd))]


		paddingvalue=0;

		DSCORES+=list(pad(Aa_scores_a,(0,max_hexpairs_length+2-len(Aa_scores_a)),mode='constant',constant_values=(0, paddingvalue)))
		DSCORES+=list(pad(Aa_scores_d,(0,max_hexpairs_length+2-len(Aa_scores_d)),mode='constant',constant_values=(0, paddingvalue)))
		
		DSCORES3 = []
		DSCORES3+=Aa_scores_a;DSCORES3+=Aa_scores_d;
		DSCORES3+=Aa_ada_scores;DSCORES3+=Aa_dad_scores;
		DSCORES3+=ada_score;DSCORES3+=dad_score
		#print("P",DSCORES3)
		DSCORES_Prob_E_P = prod(DSCORES3)

		Aa_scores_a = []; Aa_scores_d = []
		Aa_ada_scores = []; Aa_dad_scores = [];
		ada_score = []; dad_score=[]

		#Iterate over all AA pairs
		for aa_pairs in dimer[6]:

			AA1 = aa_pairs[0]; AA2 = aa_pairs[1];

			if ((aa_pairs[2]=='a') and (aa_pairs[3]=='d')):
				Aa_scores_a+= [ad_scores[tuple([AA1,AA2])]]
			elif ((aa_pairs[2]=='d') and (aa_pairs[3]=='a')):
				Aa_scores_d+= [ad_scores[tuple([AA2,AA1])]]

		#Iterate over all hex pairs
		for hexp in dimer[5]:

			AA1 = list(hexp[1]); AA2 = list(hexp[2]); AA3 = list(hexp[5]); AA4 = list(hexp[6]);

			addS1 = array(list(map(self.get_BS,AA1)))
			addS2 = array(list(map(self.get_BS,AA2)))
			
			density1 = addS1[1] - (addS1[0]+addS1[2])/2.0
			density2 = addS2[1] - (addS2[0]+addS2[2])/2.0

			dd = array(density1+density2).reshape(-1,1);

			#anti-Parallel
			if ( (hexp[3]=='ada') and (hexp[4]=='dad') ):
				Aa_ada_scores+=self.find_score_ada_dad(AA1,AA2,AA3,AA4,AAscores);
				ada_score+=[exp(kde_AP_ada.score(dd))]
				
			elif ( (hexp[3]=='dad') and (hexp[4]=='ada') ):
				Aa_dad_scores+=self.find_score_dad_ada(AA1,AA2,AA3,AA4,AAscores);
				dad_score+=[exp(kde_AP_dad.score(dd))]
		
		paddingvalue=0;

		DSCORES+=list(pad(Aa_scores_a,(0,max_hexpairs_length+2-len(Aa_scores_a)),mode='constant',constant_values=(0, paddingvalue)))
		DSCORES+=list(pad(Aa_scores_d,(0,max_hexpairs_length+2-len(Aa_scores_d)),mode='constant',constant_values=(0, paddingvalue)))
		
		DSCORES3 = []
		DSCORES3+=Aa_scores_a;DSCORES3+=Aa_scores_d;
		DSCORES3+=Aa_ada_scores;DSCORES3+=Aa_dad_scores;
		DSCORES3+=ada_score;DSCORES3+=dad_score
		DSCORES_Prob_E_AP = prod(DSCORES3)

		TotalProb = self.P_Prior_Prob*DSCORES_Prob_E_P + self.AP_Prior_Prob*DSCORES_Prob_E_AP
		
		DSCORES2[0] = self.AP_Prior_Prob*DSCORES_Prob_E_AP/TotalProb
		DSCORES2[1] = self.P_Prior_Prob*DSCORES_Prob_E_P/TotalProb

		DSCORES_NEW = []
		DSCORES_NEW+=DSCORES2; DSCORES_NEW+=DSCORES

		return DSCORES_NEW

	def Generate_CC_Scores_Bayes(self, input_obj, model_obj):
		"""Summary
		
		Args:
		    input_obj (TYPE): Description
		    model_obj (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		max_hexpairs_length = model_obj.get_padding_size()

		Hexpairs = model_obj.get_KDE_Models()

		#make default dict for b/c/e/f/g positions that is also a default dict with score 0
		#scoring_dicts = defaultdict(lambda: defaultdict(lambda: 0))
		scoring_dicts = model_obj.get_scoring_dicts()

		for CCid, CC_obj in input_obj.interface_dict.items():
		
			dimer = CC_obj.get_CC()
				
			dimer_hex = generate_Hex_pairs_dimers({CCid:dimer})
			dimer = dimer_hex[CCid]
			
			#dimer_hex = generate_Hex_pairs_reverse_dimers({key:dimer})
			#dimer = dimer_hex[key]
				
			dscore = self.score_dimer_bayesian(scoring_dicts, Hexpairs, dimer)

			self.ALLSCORES.append(dscore);
		
		return

	def Generate_CC_Scores_RF(self, input_obj, model_obj):
		"""Summary
		
		Args:
		    input_obj (TYPE): Description
		    model_obj (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		max_hexpairs_length = model_obj.get_padding_size()

		Hexpairs = model_obj.get_KDE_Models()

		#make default dict for b/c/e/f/g positions that is also a default dict with score 0
		#scoring_dicts = defaultdict(lambda: defaultdict(lambda: 0))
		scoring_dicts = model_obj.get_scoring_dicts()

		#These are sliding windows of dimers
		for CCid, CC_obj in input_obj.interface_dict.items():
		
			dimer = CC_obj.get_CC()
				
			dimer_hex = generate_Hex_pairs_dimers({CCid:dimer})
			dimer = dimer_hex[CCid]
			
			#dimer_hex = generate_Hex_pairs_reverse_dimers({key:dimer})
			#dimer = dimer_hex[key]
				
			dscore = self.score_dimer_RF(scoring_dicts, Hexpairs, dimer,max_hexpairs_length)

			self.ALLSCORES.append(dscore);
		
		return

	def Write_Raw_Scores(self,input_obj, fname):
		"""Summary
		
		Args:
		    input_obj (TYPE): Description
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		wlines = "#Id Seq1 Hep1 Seq2 Hep2 Prob_AP Prob_P RawScores\n"
		

		for i,dscores in enumerate(self.ALLSCORES):

			sumdscores = sum(dscores);
			dimervals = list(input_obj.interface_dict.values())[i].get_CC()
			wlines += str(i)+" "
			wlines += dimervals[1][0] +" "+ dimervals[1][1] +" "+ dimervals[2][0] +" "+ dimervals[2][1] + " "
			wlines += str(round(self.D_PROBS_mean[i][0],6)) +" "+str(round(self.D_PROBS_mean[i][1],6))+" "
			wlines+= str(round(sumdscores,6))+" "
			wlines+= dimervals[0]+"\n"
			
		outf = open(fname,'w')
		outf.writelines(wlines)
		outf.close()

		return

	def Write_Bayesian_Probabilities(self,input_obj, fname):
		"""Summary
		
		Args:
		    input_obj (TYPE): Description
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		wlines = "#Id Seq1 Hep1 Seq2 Hep2 Prob_AP Prob_P\n"
		
		for i,dscores in enumerate(self.ALLSCORES):

			sumdscores = sum(dscores);
			dimervals = list(input_obj.interface_dict.values())[i].get_CC()
			
			wlines += str(i)+" "
			wlines += dimervals[1][0] +" "+ dimervals[1][1] +" "+ dimervals[2][0] +" "+ dimervals[2][1] + " "
			wlines+= str(round(dscores[0],6))+" "+str(round(dscores[1],6))+" "
			wlines+= str(round(sumdscores,6))+" "
			wlines+= dimervals[0]+"\n"
			
		outf = open(fname,'w')
		outf.writelines(wlines)
		outf.close()

		return

	def Predict_CC_Stability(self, model_obj):
		"""Summary
		
		Args:
		    model_obj (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		#Get maximal predicted output label
		clf = model_obj.get_RF_Models()[0];
		self.D_PREDS = clf.predict(self.ALLSCORES);
		self.D_PROBS = clf.predict_proba(self.ALLSCORES);

		for i in range(1,len(model_obj.get_RF_Models())):
			clf = model_obj.get_RF_Models()[i];
			d_preds = clf.predict(self.ALLSCORES);
			d_probs = clf.predict_proba(self.ALLSCORES);
			
			self.D_PREDS = vstack((self.D_PREDS,d_preds))
			self.D_PROBS = dstack((self.D_PROBS,d_probs))

		self.D_PREDS = self.D_PREDS.T
		self.D_PROBS = self.D_PROBS

		#self.D_PROBS_median = median(self.D_PROBS, axis=0)
		self.D_PROBS_mean = mean(self.D_PROBS, axis=2)

		return

	def Get_Predictions(self):
		"""Summary
		"""
		print("\n-------------------------------------------")
		for i, dimer_prob in enumerate(self.D_PROBS_mean):
			
			if dimer_prob[0]>=0.5:
				#return -1.0	#antiparallel
				#return "Anti-Parallel"
				print("dimer",str(i),"Anti-Par ","Prob: ",dimer_prob)
			else:
				#return +1.0	#parallel
				#return "Parallel"
				print("dimer",str(i),"Parallel ","Prob: ",dimer_prob)

	def Get_Prediction_Probs(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.D_PROBS_mean

	def find_score_ada_ada(self, AA1,AA2,AA3,AA4,AAscores):
		"""Summary
		
		Args:
		    AA1 (TYPE): Description
		    AA2 (TYPE): Description
		    AA3 (TYPE): Description
		    AA4 (TYPE): Description
		    AAscores (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		aa_scores, dd_scores, ad_scores,ae_scores, ag_scores, de_scores, dg_scores = AAscores[0], AAscores[1], AAscores[2],AAscores[3], AAscores[4], AAscores[5], AAscores[6]

		#default frequnecy is 1.0, For 'X' amino acids
		#s1 = aa_scores[tuple([AA1[0],AA2[0]])]
		#s2 = dd_scores[tuple([AA1[1],AA2[1]])]
		#s3 = aa_scores[tuple([AA1[2],AA2[2]])]

		#add e/g residue to complete 4 amino acid hole and knob packing
		s4 = de_scores[tuple([AA1[1],AA4[1]])]
		s5 = de_scores[tuple([AA2[1],AA3[1]])]
		
		return [s4,s5]
		
	def find_score_dad_dad(self, AA1,AA2,AA3,AA4,AAscores):
		"""Summary
		
		Args:
		    AA1 (TYPE): Description
		    AA2 (TYPE): Description
		    AA3 (TYPE): Description
		    AA4 (TYPE): Description
		    AAscores (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		aa_scores, dd_scores, ad_scores,ae_scores, ag_scores, de_scores, dg_scores = AAscores[0], AAscores[1], AAscores[2],AAscores[3], AAscores[4], AAscores[5], AAscores[6]

		#default frequnecy is 1.0, For 'X' amino acids
		#s1 = dd_scores[tuple([AA1[0],AA2[0]])]
		#s2 = aa_scores[tuple([AA1[1],AA2[1]])]
		#s3 = dd_scores[tuple([AA1[2],AA2[2]])]

		#add e/g residue to complete 4 amino acid hole and knob packing
		s4 = ag_scores[tuple([AA1[1],AA4[1]])]
		s5 = ag_scores[tuple([AA2[1],AA3[1]])]

		return [s4,s5]

	def find_score_ada_dad(self, AA1,AA2,AA3,AA4,AAscores):
		"""Summary
		
		Args:
		    AA1 (TYPE): Description
		    AA2 (TYPE): Description
		    AA3 (TYPE): Description
		    AA4 (TYPE): Description
		    AAscores (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		aa_scores, dd_scores, ad_scores,ae_scores, ag_scores, de_scores, dg_scores = AAscores[0], AAscores[1], AAscores[2],AAscores[3], AAscores[4], AAscores[5], AAscores[6]

		#default frequnecy is 1.0, For 'X' amino acids
		#s1 = ad_scores[tuple([AA1[0],AA2[0]])]
		#s2 = ad_scores[tuple([AA2[1],AA1[1]])]
		#s3 = ad_scores[tuple([AA1[2],AA2[2]])]

		#add e/g residue to complete 4 amino acid hole and knob packing
		s4 = dg_scores[tuple([AA1[1],AA4[1]])]
		s5 = ae_scores[tuple([AA2[1],AA3[1]])]

		return [s4,s5]

	def find_score_dad_ada(self, AA1,AA2,AA3,AA4,AAscores):
		"""Summary
		
		Args:
		    AA1 (TYPE): Description
		    AA2 (TYPE): Description
		    AA3 (TYPE): Description
		    AA4 (TYPE): Description
		    AAscores (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		aa_scores, dd_scores, ad_scores,ae_scores, ag_scores, de_scores, dg_scores = AAscores[0], AAscores[1], AAscores[2],AAscores[3], AAscores[4], AAscores[5], AAscores[6]

		#default frequnecy is 1.0, For 'X' amino acids
		#s1 = ad_scores[tuple([AA2[0],AA1[0]])]
		#s2 = ad_scores[tuple([AA1[1],AA2[1]])]
		#s3 = ad_scores[tuple([AA2[2],AA1[2]])]

		#add e/g residue to complete 4 amino acid hole and knob packing
		s4 = ae_scores[tuple([AA1[1],AA4[1]])]
		s5 = dg_scores[tuple([AA2[1],AA3[1]])]

		return [s4,s5]

