
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


from collections import defaultdict
import numpy as np
import os
from collections import OrderedDict, defaultdict
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import statistics

class DefaultDictNoAdd(defaultdict):

	"""Summary
	"""
	
	def __missing__(self, key):
		"""Summary
		
		Args:
		    key (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		return ['-','-',0]

class ExtendCoiledCoil:

	"""Summary
	
	Attributes:
	    DistDict (TYPE): Description
	    DistMat (TYPE): Description
	    p1_PC (TYPE): Description
	    p2_PC (TYPE): Description
	    pc_prot1_file (TYPE): Description
	    pc_prot2_file (TYPE): Description
	    phep1 (str): Description
	    phep2 (str): Description
	    ppath (TYPE): Description
	    PROB_CUTOFF (TYPE): Description
	    prot1 (TYPE): Description
	    prot1_max (TYPE): Description
	    prot1_min (TYPE): Description
	    prot1res (TYPE): Description
	    prot2 (TYPE): Description
	    prot2_max (TYPE): Description
	    prot2_min (TYPE): Description
	"""
	
	def __init__(self, ppath, prot1_prot2, pc_prot1_file, pc_prot2_file, PROB_CUTOFF):
		"""Summary
		
		Args:
		    ppath (TYPE): Description
		    prot1_prot2 (TYPE): Description
		    pc_prot1_file (TYPE): Description
		    pc_prot2_file (TYPE): Description
		    PROB_CUTOFF (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		self.PROB_CUTOFF = PROB_CUTOFF

		self.ppath = ppath

		self.pc_prot1_file = pc_prot1_file
		self.pc_prot2_file = pc_prot2_file

		self.phep1 = 'a'
		self.phep2 = 'd'

		plot_True=True;

		toks = prot1_prot2.split('_')
		
		self.prot1=toks[0]
		self.prot2=toks[1]

		#This is a anti-parallel flag. For parallel use False
		for APFlag in [True, False]:
		#for APFlag in [False]:
			
			Segments, self.DistDict, self.DistMat = self.Sample_Plot_Combine( APFlag)

			##APFlag = True ;#This is a anti-parallel flag. For parallel use False
			self.Generate_Alignments_Scores(Segments, APFlag, plot_True)
		
		return

	def def_value(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return None

	def load_paircoil_alignment(self, pc_file):
		"""Summary
		
		Args:
		    pc_file (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		inf = open(pc_file,'r')
		lines = inf.readlines()
		inf.close()

		prot_seq = DefaultDictNoAdd()

		for l in lines:
			if l[0]=='#':
				continue

			if len(l)<2:
				continue

			toks = l.strip().split()
		
			if float(toks[3]) >= self.PROB_CUTOFF:
				prot_seq[int(toks[0])]=[toks[1],toks[2],float(toks[3])]
			else:
				prot_seq[int(toks[0])]=[toks[1],'-',float(toks[3])]

		return prot_seq

	def Sample_Plot_Combine(self, APFlag):
		"""Summary
		
		Args:
		    APFlag (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		scoreidx=-1
		if APFlag==False:
			scoreidx = 6
		else:
			scoreidx = 5

		self.p1_PC = self.load_paircoil_alignment(self.pc_prot1_file)
		self.p2_PC = self.load_paircoil_alignment(self.pc_prot2_file)

		self.DistMat = np.zeros((len(self.p1_PC.keys()),len(self.p2_PC.keys())))

		
		inf = open(self.ppath+self.prot1+'_'+self.prot2+'_CoiledCoil_Scores_'+self.phep1+'.dat','r')
		lines=inf.readlines()
		inf.close()


		self.DistDict=defaultdict(self.def_value)

		G_SWDs1=[];G_SWDs2=[]

		for l in lines:
			if l[0]=='#':
				continue

			toks = l.strip().split()
			res = toks[8].split('-')

			sw1 = tuple([int(res[0]),int(res[1])])
			sw2 = tuple([int(res[2]),int(res[3])])

			#SWDs1.append(sw1);SWDs2.append(sw2);
			G_SWDs1.append(sw1);G_SWDs2.append(sw2);

			#Probability of parallel dimer
			#self.DistDict[tuple([sw1,sw2])] = -1*np.log(float(toks[scoreidx]))
			self.DistDict[tuple([sw1,sw2])] = float(toks[scoreidx])

		Segments=[]

		inf = open(self.ppath+self.prot1+'_'+self.prot2+'_CoiledCoil_Scores_'+self.phep2+'.dat','r')
		lines=inf.readlines()
		inf.close()

		SWDs1=[];SWDs2=[]

		for l in lines:
			if l[0]=='#':
				continue
			toks = l.strip().split()
			res = toks[8].split('-')
			sw1 = tuple([int(res[0]),int(res[1])])
			sw2 = tuple([int(res[2]),int(res[3])])

			#SWDs1.append(sw1);SWDs2.append(sw2);
			G_SWDs1.append(sw1);G_SWDs2.append(sw2);

			#Probability of parallel dimer
			#self.DistDict[tuple([sw1,sw2])] = -1*np.log(float(toks[scoreidx]))
			self.DistDict[tuple([sw1,sw2])] = float(toks[scoreidx])

		G_SWDs1 = list(set(G_SWDs1))
		G_SWDs2 = list(set(G_SWDs2))

		G_SWDs1.sort(key=lambda k: (k[0], k[1]))
		G_SWDs2.sort(key=lambda k: (k[0], k[1]))

		#Iterate over all sliding windows
		for sw1 in G_SWDs1:
			for sw2 in G_SWDs2:

				#Skip that are not present or gets scored
				swval = self.DistDict[tuple([sw1,sw2])]
				if swval==None:
					continue

				resrange1 = sw1
				resrange2 = sw2

				#Residue numbers starts with index 1 not zero
				reslist1 = list(range(resrange1[0],resrange1[1]+1))
				reslist2 = list(range(resrange2[0],resrange2[1]+1))

				swcount=1;
				Segments.append([reslist1,reslist2,swval,swcount,None])
		
		return Segments, self.DistDict, self.DistMat

	def plot_individual_cc_scores(self, Segments_Temp):
		"""Summary
		
		Args:
		    Segments_Temp (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		Homodimer_Scores = OrderedDict()

		#This loop extends the overlapping segments
		for i,segs1 in tqdm.tqdm(enumerate(Segments_Temp)):
			if segs1[0]==segs1[1]:

				for resid in segs1[0]:
					if resid in Homodimer_Scores.keys():
						Homodimer_Scores[resid].append(segs1[2])
					else:
						Homodimer_Scores[resid]=[]
						Homodimer_Scores[resid].append(segs1[2])

		Homo_Residue_Scores = []
		Homo_Residue_HighScores = []
		Homo_Residue_HighScores2 = []
		Homo_Residue_Scores_SD = []
		Homo_Residues = []
		Homodimer_CCProb = []

		for k,v in self.p1_PC.items():
			
			Homo_Residues.append(k)
			Homodimer_CCProb.append(self.p1_PC[k][2])

			if k in Homodimer_Scores.keys():
				
				mean_val = np.mean(Homodimer_Scores[k])
				sd_val = np.std(Homodimer_Scores[k])
				
				Homo_Residue_Scores.append(mean_val)
				Homo_Residue_Scores_SD.append(sd_val)

				if mean_val >= self.PROB_CUTOFF:
					Homo_Residue_HighScores2.append(-0.25) #Position where to put a dot
				else:
					Homo_Residue_HighScores2.append(None)

			else:
				Homo_Residue_Scores.append(None)
				Homo_Residue_Scores_SD.append(None)
				Homo_Residue_HighScores2.append(None)

		ccwlines='';ccwlines_dgt=''

		fig = plt.figure(figsize=(16,4))
		sns.set(font_scale=1)
		sc2 = plt.plot(Homo_Residues, Homodimer_CCProb,label='Pcoil Scores')
		sc2 = plt.plot(Homo_Residues, Homo_Residue_Scores,label='Coconut Scores')
		sc2 = plt.scatter(Homo_Residues, Homo_Residue_HighScores2, s=2, label='CoilCoiled')
		plt.legend()
		plt.savefig(self.ppath+self.prot1+'_'+self.prot2+'_Coconut_Scores_Raw.pdf')
		plt.xlabel("Amino acid sequence homodimer")
		plt.ylabel("Amino acid CC Scores")
		#plt.show()

		#Digitize the probabilities at 0.25
		digit_prob = 0.25

		for i, val in enumerate(Homodimer_CCProb):
			if val==None:
				continue
			if val >= digit_prob:
				Homodimer_CCProb[i] = 1.0
			else:
				Homodimer_CCProb[i] = 0.0

		for i, val in enumerate(Homo_Residue_Scores):
			if val==None:
				ccwlines+=str(Homo_Residues[i])+'\t'+'-'+'\n'
				continue
			if val >= digit_prob:
				ccwlines+=str(Homo_Residues[i])+'\t'+str(round(Homo_Residue_Scores[i],2))+'\n'
				Homo_Residue_Scores[i] = 0.9
			else:
				ccwlines+=str(Homo_Residues[i])+'\t'+str(round(Homo_Residue_Scores[i],2))+'\n'
				Homo_Residue_Scores[i] = 0.1
		

		with open(self.ppath+self.prot1+'_'+self.prot2+'_Coconut_Scores.txt','w') as outf:
			outf.writelines(ccwlines)

		fig = plt.figure(figsize=(16,4))
		sns.set(font_scale=1)
		sc2 = plt.plot(Homo_Residues, Homodimer_CCProb,label='Pcoil Scores')
		sc2 = plt.plot(Homo_Residues, Homo_Residue_Scores,label='Coconut Scores')
		sc2 = plt.scatter(Homo_Residues, Homo_Residue_HighScores2, s=2, label='CoilCoiled')
		plt.legend()
		plt.savefig(self.ppath+self.prot1+'_'+self.prot2+'_Coconut_Scores.pdf')
		plt.xlabel("Amino acid sequence homodimer")
		plt.ylabel("Amino acid CC Scores")
		#plt.show()

		return

	def extend_segment(self, Segments, APFlag):
		"""Summary
		
		Args:
		    Segments (TYPE): Description
		    APFlag (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		Segments_Temp=[]
		for i,segs in enumerate(Segments):
			reslist1 = list(sorted(set(segs[0])));
			reslist2 = list(sorted(set(segs[1]),reverse=APFlag))
			Segments_Temp.append([reslist1,reslist2,segs[2],segs[3],list(zip(reslist1,reslist2))])


		self.plot_individual_cc_scores(Segments_Temp)

		Segments_Extended=[];
		Segments_New3_Dict=OrderedDict()


		segs1 = Segments_Temp[0]
		s1_p1p2_rl = segs1[4]
		s1_res_s = s1_p1p2_rl[0];
		s1_res_e = s1_p1p2_rl[-1];
		s1_key = tuple([s1_res_s,s1_res_e])

		keys_to_delete=[]

		#This loop extends the overlapping segments
		for i,segs1 in tqdm.tqdm(enumerate(Segments_Temp)):
			
			s1_p1_rl = segs1[0]
			s1_p2_rl = segs1[1]
			s1_swval = segs1[2]
			s1_count   = segs1[3]
			s1_p1p2_rl = segs1[4]
			s1_res_s = s1_p1p2_rl[0];
			s1_res_e = s1_p1p2_rl[-1];
			s1_key = tuple([s1_res_s,s1_res_e])
			
			flag=0
			for k, v in Segments_New3_Dict.items():
				if s1_res_s in v[4]:
					k_key = tuple([k[0],s1_res_e])
					Segments_New3_Dict[k_key] =[[],[],0,0,[],[]]
					Segments_New3_Dict[k_key][0] = Segments_New3_Dict[k][0] + s1_p1_rl;
					Segments_New3_Dict[k_key][1] = Segments_New3_Dict[k][1] + s1_p2_rl;
					Segments_New3_Dict[k_key][2] = Segments_New3_Dict[k][2] + s1_swval;
					Segments_New3_Dict[k_key][3] = Segments_New3_Dict[k][3] + s1_count;
					Segments_New3_Dict[k_key][4] = Segments_New3_Dict[k][4] + s1_p1p2_rl;
					Segments_New3_Dict[k_key][4] = list(sorted(set(Segments_New3_Dict[k_key][4]),reverse=APFlag));
					Segments_New3_Dict[k_key][5]+= Segments_New3_Dict[k][5][:]
					Segments_New3_Dict[k_key][5]+= [s1_swval]

					flag=1
					del Segments_New3_Dict[k]
					break
			
			if flag==0:
				Segments_New3_Dict[s1_key]=[[],[],0,0,[],[]]
				Segments_New3_Dict[s1_key][0] = s1_p1_rl;
				Segments_New3_Dict[s1_key][1] = s1_p2_rl;
				Segments_New3_Dict[s1_key][2] = s1_swval;
				Segments_New3_Dict[s1_key][3] = s1_count;
				Segments_New3_Dict[s1_key][4] = s1_p1p2_rl;
				Segments_New3_Dict[s1_key][4] = list(sorted(set(Segments_New3_Dict[s1_key][4]),reverse=APFlag));
				Segments_New3_Dict[s1_key][5].append(s1_swval)
		

		Segments_Extended=list(set(Segments_Extended))
		Segments_New3=[]
		
		Segments_New3 = list(Segments_New3_Dict.values())
			
		return Segments_New3

	def Generate_Alignments_Scores(self, Segments, APFlag, plot_True):
		"""Summary
		
		Args:
		    Segments (TYPE): Description
		    APFlag (TYPE): Description
		    plot_True (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		#APFlag = True #This is a anti-parallel flag. For parallel use False
		Segments_New = self.extend_segment(Segments, APFlag)


		self.prot1res=[];self.prot2res=[]
		for segs in Segments_New:
			self.prot1res+=segs[0]
			self.prot2res+=segs[1]
		self.prot1_min=np.min(self.prot1res)
		self.prot2_min=np.min(self.prot2res)
		self.prot1_max=np.max(self.prot1res)
		self.prot2_max=np.max(self.prot2res)
		
		diag_cutoff1 = 0.3*(self.prot1_max-self.prot1_min)
		diag_cutoff2 = 0.3*(self.prot2_max-self.prot2_min)

		print("cutoff",diag_cutoff1)
		segslines=""

		Diagonal_Score = 0
		Diagonal_Count = 0

		for segs in Segments_New:
			
			reslist1 = segs[0];
			reslist2 = segs[1];
			reslist12 = segs[4]
			segs[5] = list(set(segs[5]))

			segm1 = str(reslist1[0])+"-"+str(reslist1[-1])
			segm2 = str(reslist2[0])+'-'+str(reslist2[-1])
			
			segslines+=segm1+'\t'+segm2+'\t'+str(segs[2]/segs[3])+"\n"
			
			for i, [resnum1, resnum2] in enumerate(reslist12):

				self.DistMat[resnum1-1,resnum2-1]=segs[2]/segs[3]

		fileAPname=""
		if APFlag==False:
			fileAPname='parallel'
			#self.get_parallel_diagonal_alignment(Segments_New, fileAPname)
		else:
			fileAPname='antiparallel'

		outf = open(self.ppath + self.prot1+'_'+self.prot2+"_Coiled-Coil_Segments_"+fileAPname+".dat",'w')
		outf.writelines(segslines)
		outf.close()

		if plot_True:
			figsize1=(4,3)
			self.plot_alignment( APFlag, figsize1)

		return

	def get_parallel_diagonal_alignment(self, Segments_New, fileAPname):
		"""Summary
		
		Args:
		    Segments_New (TYPE): Description
		    fileAPname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		segslines=""
		Diagonal_Score = 0
		Diagonal_Count = 0

		for segs in Segments_New:
			reslist1 = segs[0];
			reslist2 = segs[1];
			reslist12 = segs[4]

			segm1 = str(reslist1[0])+"-"+str(reslist1[-1])
			segm2 = str(reslist2[0])+'-'+str(reslist2[-1])

			if segm1 == segm2:
				Diagonal_Score+=float(segs[2]/segs[3])
				Diagonal_Count+=1
				segslines+=segm1+'\t'+segm2+'\t'+str(segs[2]/segs[3])+"\n"

		Avg_Diagonal_Score = Diagonal_Score / Diagonal_Count
		if Avg_Diagonal_Score>=0.5:
			print("Parallel Coiled-coil is predicted. Total Score:",Avg_Diagonal_Score)
		else:
			print("Anti-Parallel Coiled-coil is predicted. Total Score:",Avg_Diagonal_Score)

		outf = open(self.ppath + self.prot1+'_'+self.prot2+"_Coiled-Coil_Diagonal_Segments_"+fileAPname+".dat",'w')
		outf.writelines(segslines)
		outf.close()

		return

	def plot_alignment(self, APFlag, figsize1=(6,3), fscale=1):
		"""Summary
		
		Args:
		    APFlag (TYPE): Description
		    figsize1 (tuple, optional): Description
		    fscale (int, optional): Description
		
		Returns:
		    TYPE: Description
		"""
		p1_res = sorted(list(self.p1_PC.keys()))
		p2_res = sorted(list(self.p2_PC.keys()))

		self.DistMat = self.DistMat.transpose()

		score_nonzero_x=[]
		score_nonzero_y=[]
		score_color=[]
		for i,row in enumerate(self.DistMat):
			for j,val in enumerate(row):
				#Skiiping with all zero score values, they will remian blank in the figure.
				if val==0:
					continue
				score_nonzero_y.append(i)
				score_nonzero_x.append(j)
				score_color.append(val)


		fig = plt.figure(figsize=figsize1) #For 2ztaA-2ztaB
		sns.set(font_scale=fscale)
		sc2 = plt.scatter(x=score_nonzero_x,y=score_nonzero_y,s=1,c=score_color,cmap="cool",alpha=1.0)
		plt.colorbar(sc2)
		plt.tight_layout()
		plt.xlabel("Amino acid sequence 1")
		plt.ylabel("Amino acid sequence 2")
		plt.xlim(p1_res[0]-1,p1_res[-1]+1)
		plt.ylim(p2_res[0]-1,p2_res[-1]+1)
		
		
		fileAPname=""
		if APFlag==False:
			fileAPname='parallel'
		else:
			fileAPname='antiparallel'
		plt.title(fileAPname+' alignment')
		plt.savefig(self.ppath+ self.prot1+'_'+self.prot2+'_Aligned_score_'+fileAPname+'.pdf')
		plt.close()
		
		return








