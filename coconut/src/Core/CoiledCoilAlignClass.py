"""Summary
"""
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


from collections import defaultdict, OrderedDict
from sys import setrecursionlimit
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import networkx as nx
from sklearn.cluster import DBSCAN
import pickle
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib import collections  as mc
from networkx.drawing.nx_pydot import graphviz_layout
import sys
import os
from itertools import islice
import copy

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

class CoiledCoilAlign:

	"""Summary
	
	Attributes:
	    BaseAlignmentFile (TYPE): Description
	    constant_score (float): Description
	    fsseq1 (TYPE): Description
	    fsseq2 (TYPE): Description
	    pc_prot1_file (TYPE): Description
	    pc_prot2_file (TYPE): Description
	    ppath (TYPE): Description
	    PROB_CUTOFF (TYPE): Description
	    prot1 (TYPE): Description
	    prot2 (TYPE): Description
	    topk (int): Description
	

	"""
	
	#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7122542/
	#RESIDUES_PER_TURN=3.5

	#https://pubmed.ncbi.nlm.nih.gov/14222887/
	#AXIAL_RISE_PER_RESIDUE=1.485

	def __init__(self, ppath, prot1_prot2, pc_prot1_file, pc_prot2_file, PROB_CUTOFF, BaseAlignmentFile=None):
		"""Summary
		
		Args:
		    ppath (TYPE): Description
		    prot1_prot2 (TYPE): Description
		    pc_prot1_file (TYPE): Description
		    pc_prot2_file (TYPE): Description
		    PROB_CUTOFF (TYPE): Description
		    BaseAlignmentFile (None, optional): Description
		
		Returns:
		    TYPE: Description
		"""
		#global __TOP_N_SEGMENTS
		self.__TOP_N_SEGMENTS=100

		prots = prot1_prot2.split('_')

		self.prot1 = prots[0]
		self.prot2 = prots[1]
		self.ppath = ppath
		self.PROB_CUTOFF = PROB_CUTOFF

		self.topk = 1 #Number of alignments to write 
		
		self.constant_score=0.1  #Base score added to the alignment based on alignment from cliustal omega. Infuture this will be replced by blosum matrix
		self.pc_prot1_file = pc_prot1_file
		self.pc_prot2_file = pc_prot2_file

		self.fsseq1 = self.read_fasta(self.ppath+prots[0]+'.fasta')
		self.fsseq2 = self.read_fasta(self.ppath+prots[1]+'.fasta')

		self.BaseAlignmentFile=BaseAlignmentFile

		alignment = self.coconut_align(self.fsseq1, self.fsseq2, self.ppath, prot1_prot2, plotchart=True)
		
		#CHECK THESE FUNCTIONS NOT WORKING PROPERLY. NEED to assing source node properly
		##alignment = self.read_coco_alignment_graph(self.fsseq1, self.fsseq2, self.ppath, prot1_prot2)
		#alignment = self.read_coco_alignment_graph_new(self.fsseq1, self.fsseq2, self.ppath, prot1_prot2, self.topk)

		return

	def create_base_alignment(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		lines = ''
		lines+='>'+self.prot1+'\n'
		lines+=self.fsseq1+'\n'
		lines+='>'+self.prot2+'\n'
		lines+=self.fsseq2+'\n'
		with open(self.BaseAlignmentFile,'w') as outf:
			outf.writelines(lines)

		return

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

	def read_fasta(self, fname):
		"""Summary
		
		Args:
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		inf = open(fname, 'r')
		lines = inf.readlines()
		inf.close()
		Qryseq = ""
		for l in lines:
			if l[0] == '>':
				continue
			Qryseq += l.strip()
		return Qryseq

	def coconut_align(self, seq1, seq2, ppath, prot1_prot2, plotchart = False):
		"""Summary
		
		Args:
		    seq1 (TYPE): Description
		    seq2 (TYPE): Description
		    ppath (TYPE): Description
		    prot1_prot2 (TYPE): Description
		    plotchart (bool, optional): Description
		
		Returns:
		    TYPE: Description
		"""
		#get the length of the sequences in two variables, n and m
		n = len(seq1)
		m = len(seq2)
		
		print("Length",n,m)

		#set the recursion limit according to the longest sequence + offset (normally 1000)
		setrecursionlimit(max(n,m)+1000)
		
		print("Loading Coiled-coil segments")
		Segment_Dict_full = self.load_Coiled_Coil_Segments_full(ppath+prot1_prot2+'_Coiled-Coil_Segments_parallel.dat')
		#Segment_Dict_full = self.load_Coiled_Coil_Segments_full(ppath+prot1_prot2+'_Coiled-Coil_Segments_antiparallel.dat')

		#Add default alignment in Segment_dict for base scores of BLOSUM.
		alignment_file = self.BaseAlignmentFile
		Segment_Dict_full = self.add_alignment_scores_segmentdict(Segment_Dict_full, seq1, seq2, alignment_file, self.constant_score )
		
		
		figscale_const = int(float(np.max([n,m])/800)+1)
		fig,ax = plt.subplots(figsize=((4*figscale_const,3*figscale_const)))

		score_cutoff=0.0
		
		print("Filtering Coiled-coil segments")
		prot1,prot2 = prot1_prot2.split('_')
		
		#Get the diagonal segments with bandwidth 20%
		bandwidth1=1;bandwidth2=1
		
		Segment_Dict_full, Segment_Data, max_pos_key, min_pos_key = self.Filter_OffDiagonal_Segments(Segment_Dict_full,n,m,bandwidth1,bandwidth2,score_cutoff)
		
		min_pos_key_copy,max_pos_key_copy = copy.copy(min_pos_key),copy.copy(max_pos_key)

		#Add Source and Sink nodes segments
		Segment_Dict_full[tuple([min_pos_key[1],min_pos_key[3]])]=[1,min_pos_key[0],min_pos_key[2],0.0000001];
		Segment_Dict_full[tuple([max_pos_key[1],max_pos_key[3]])]=[1,max_pos_key[0],max_pos_key[2],0.0000001];
		
		Segment_Dict_full = dict(sorted(Segment_Dict_full.items(), key=lambda x: (x[0][0],x[0][1]), reverse=False ) )

		#Add a default stop value for alignments
		Segment_Data.append([min_pos_key[1],min_pos_key[3]])
		Segment_Data.append([max_pos_key[1],max_pos_key[3]])

		assert len(Segment_Data)==len(Segment_Dict_full.keys())

		#Create KDTree of all the segment filtered from the off diagonal elements using their end points
		KDtree = KDTree(Segment_Data)

		#plot_segments(Segment_Dict_full,ax,0.5,None,score_cutoff)
		#plt.savefig(ppath+prot1_prot2+'_Aligned_score.pdf')
		#plt.show()

		
		Total_Segments=len(Segment_Dict_full.keys())
		print("Total Num of Segment:",Total_Segments)

		if Total_Segments < self.__TOP_N_SEGMENTS:
			print("Setting Maximum segment search!")
			self.__TOP_N_SEGMENTS = Total_Segments

		print("Starting Graph creation and traceback...")
		alignments = nx.DiGraph()
		

		###---------------------------------------------------------
		#
		###---------------------------------------------------------
		self.__NODE_INDEX=-1
		#Create Graph nodes; Since dictionary is sorted, Source Node (Maximum Seq + XX residue) INDEX will be 1.
		for k, v in Segment_Dict_full.items():
			
			self.__NODE_INDEX+=1

			seq1seq2 = tuple([v[1],k[0],v[2],k[1]])
			
			node_score = None
			if v[0]==0:
				node_score = 0.0000001
			else:
				node_score = v[0]

			#FOR TESTING
			#USE BLOSUM SCORES
			#alignments.add_node(self.__NODE_INDEX, seq1seq2=seq1seq2, score=-1*np.log(node_score)-1*np.log(v[3]) , alignflag=0)
			#DO NOT USE BLOSUM SCORES. 
			#LOWER SCORE the BETTER, RANGE: 0 + inf
			##alignments.add_node(self.__NODE_INDEX, seq1seq2=seq1seq2, score=-1*np.log(node_score), alignflag=0)
			#LOWER SCORE the BETTER, node_score range: 0 to 1 (probability), ALL SCORE RANGE: 0 to 1, 0 best
			
			alignments.add_node(self.__NODE_INDEX, seq1seq2=seq1seq2, score=1-node_score, alignflag=0)
			Segment_Dict_full[k].append(self.__NODE_INDEX)

		###
		###---------------------------------------------------------
		#
		###---------------------------------------------------------
		GLOB_SOURCE=self.__NODE_INDEX
		GLOB_TARGET=0

		Prev_TopNLeafs=[]; TopNLeafs=[self.__NODE_INDEX];Depth_Level=0; max_depth=1

		iteration=0
		while TopNLeafs!=Prev_TopNLeafs:
			print("Running Depth Level:",Depth_Level, len(TopNLeafs))
			
			templeafs=[]
			for nodeid in TopNLeafs:

				currentid=nodeid;sourceid=nodeid;
				max_pos_key = alignments.nodes()[nodeid]['seq1seq2']
				
				alignments = self.segment_traceback( max_pos_key, alignments, Segment_Dict_full, 
					KDtree, Segment_Data, max_depth,sourceid,currentid)

				sourceid=nodeid; #cluster and Merge similar scoring segments
				
				
			Prev_TopNLeafs=TopNLeafs

			#topn=1000;sourceid=0;
			topn=1;sourceid=0;
			TopNLeafs = [x for x in alignments.nodes() if alignments.out_degree(x)==0 and alignments.in_degree(x)!=0]
			
			Depth_Level+=max_depth

		
		#write_dot(alignments,ppath+prot1_prot2+'_alignment_graph.dot')
		with open(ppath+prot1_prot2+'_alignment_graph.gpickle', 'wb') as f:
			pickle.dump(alignments, f, pickle.HIGHEST_PROTOCOL)

		#self.plot_segments(Segment_Dict_full,ax,1.0,2.0,None,score_cutoff)
		#self.plot_segments(Segment_Dict_full,ax,2.0,5.0,None,score_cutoff)
		self.plot_segments(Segment_Dict_full,ax,2.0,5.0,None,score_cutoff)

		##plt.scatter(ALL_XLINKS_PAIRS[:,0],ALL_XLINKS_PAIRS[:,1],s=1,zorder=3)
		
		print("Ploting alignment")

		TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE = self.get_tracepath_alignments(alignments, GLOB_SOURCE,GLOB_TARGET, seq1,seq2, self.topk)
		
		print("Total Alignments:",len(TRACE))
		for i,tracepath in enumerate(TRACE):
			print("Alignment:",i,ALISCORE[i])
			#self.plot_alignment_segments(TRACE_SEGMENT[i],ax,0.5,"gray")
			self.plot_alignment_segments(TRACE_SEGMENT[i],ax,0.7,"gray")

		#ax.set_xlim(min_pos_key_copy[0]-10,max_pos_key_copy[0]+10)
		#ax.set_ylim(min_pos_key_copy[2]-10,max_pos_key_copy[2]+10)
		ax.set_xlim(min_pos_key_copy[0]-1,max_pos_key_copy[0]+1)
		ax.set_ylim(min_pos_key_copy[2]-1,max_pos_key_copy[2]+1)


		print(ppath+prot1_prot2+'_Aligned_score.pdf')
		plt.title(prot1+' '+prot2+' all alignments')
		plt.xlabel(prot1+' amino acids')
		plt.ylabel(prot2+' amino acids')
		plt.savefig(ppath+prot1_prot2+'_Aligned_score.pdf')
		#plt.show()
		

		fig,ax = plt.subplots(figsize=(25,12))
		#self.plot_tree(alignments,ppath+prot1_prot2+'_alignment_graph.pdf')
		#self.plot_tree_simple(alignments,ppath+prot1_prot2+'_alignment_graph.pdf')

		
		oaname=self.ppath+self.prot1+'_'+self.prot2+'_alignment.cali'
		topn=1
		self.Write_Topn_Alignments(oaname, topn, self.prot1, self.prot2, seq1, seq2, TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE)

		TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE = self.get_tracepath_alignments_fullseq(alignments, GLOB_SOURCE, GLOB_TARGET, seq1, seq2, self.topk)
		oaname=self.ppath+self.prot1+'_'+self.prot2+'_alignment.xali'
		topn=1
		self.Write_Topn_Alignments(oaname, topn, self.prot1, self.prot2, seq1, seq2, TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE)
		
		return

	def load_Coiled_Coil_Segments_full(self, fname):
		"""Summary
		
		Args:
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		inf = open(fname,'r')
		lines = inf.readlines()
		inf.close()

		#Segment_Dict=defaultdict(defvalseg2)
		Segment_Dict=OrderedDict()
		
		for l in lines:
			toks = l.strip().split('\t')
			prot_res1=toks[0].split('-'); 
			prot_res2=toks[1].split('-')
			#Segment_Dict[tuple([int(prot_res1[0]),int(prot_res1[1]),int(prot_res2[0]),int(prot_res2[1])])]=float(toks[2])
			Segment_Dict[tuple([int(prot_res1[1]),int(prot_res2[1])])]=[float(toks[2]),int(prot_res1[0]),int(prot_res2[0])]

		return Segment_Dict

	def add_alignment_scores_segmentdict(self, Segment_Dict_full, seq1, seq2, alignment_file, constant_score ):
		"""Summary
		
		Args:
		    Segment_Dict_full (TYPE): Description
		    seq1 (TYPE): Description
		    seq2 (TYPE): Description
		    alignment_file (TYPE): Description
		    constant_score (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		from BaseScores import BaseScoresAlign

		_EXP1_ROOT = os.path.abspath(os.path.dirname(__file__))
		matrix_file = _EXP1_ROOT+'/../../DATASET/BLOSUM62.txt'

		BSA =  BaseScoresAlign(matrix_file)
		BaseScoreMat = BSA.get_scoremat()

		if alignment_file!=None:

			cwseq1, cwseq2  = self.read_aligned_fasta(alignment_file)

			MLP_CC = {'P1_CW':cwseq1,'P2_CW':cwseq2}

			Seq_dict = {**MLP_CC}
			
			#Gettin ght e maximum length of the coverage
			MaxLengthAlignment = max(len(Seq_dict['P1_CW']),len(Seq_dict['P2_CW']))

			ResNumMap = OrderedDict()
			ResNumMap = {'P1_MAP':[],'P2_MAP':[]}

			mlp1_itr=0; mlp2_itr=0;

			for i in range(0,MaxLengthAlignment):
				
				if Seq_dict['P1_CW'][i]!='-':
					mlp1_itr+=1
					ResNumMap['P1_MAP'].append(mlp1_itr)
				else:
					ResNumMap['P1_MAP'].append(None)

				if Seq_dict['P2_CW'][i]!='-':
					mlp2_itr+=1
					ResNumMap['P2_MAP'].append(mlp2_itr)
				else:
					ResNumMap['P2_MAP'].append(None)

			Pairing_score=constant_score
			PairedResidues=OrderedDict()
			for i in range(0,MaxLengthAlignment):
				if ((ResNumMap['P1_MAP'][i]!=None) and (ResNumMap['P2_MAP'][i]!=None)):
					PairedResidues[tuple([int(ResNumMap['P1_MAP'][i]),int(ResNumMap['P2_MAP'][i])])]=Pairing_score

			for k,v in Segment_Dict_full.items():
				if k in PairedResidues.keys():
					Segment_Dict_full[k].append(PairedResidues[k])
					
				else:
					Segment_Dict_full[k].append(0.0000001)

		else:

			for k,v in Segment_Dict_full.items():

				#Segment_Dict_full[k].append(constant_score)
				segment_residues = zip( list(range(v[1],k[0]+1)), list(range(v[2],k[1]+1)) )

				TotalBlosumScore=[]
				for respair in segment_residues:
					respair_key = tuple([seq1[respair[0]],seq2[respair[1]]])
					TotalBlosumScore.append(BaseScoreMat[respair_key])

				TotalBlosumScore = np.mean(TotalBlosumScore)

				Segment_Dict_full[k].append(TotalBlosumScore)
				

		return Segment_Dict_full

	def read_aligned_fasta(self, fname):
		"""Summary
		
		Args:
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		inf = open(fname, 'r')
		lines = inf.readlines()
		inf.close()

		ALL_SEQS=[]
		Qryseq="";
		for l in lines:
			if l[0] == '>':
				ALL_SEQS.append(Qryseq)
				Qryseq=""
				continue
			else:
				Qryseq += l.strip()
		ALL_SEQS.append(Qryseq)

		return ALL_SEQS[1:]

	def Filter_OffDiagonal_Segments(self, Segment_Dict_full, n, m, bandwidth1, bandwidth2, score_cutoff):
		"""Summary
		
		Args:
		    Segment_Dict_full (TYPE): Description
		    n (TYPE): Description
		    m (TYPE): Description
		    bandwidth1 (TYPE): Description
		    bandwidth2 (TYPE): Description
		    score_cutoff (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		Segment_Dict_OffDiag=OrderedDict()

		#COnsidering largest bandwith for each diagonal side to avoid missing good models
		B1 = int(bandwidth1*n)
		B2 = int(bandwidth2*n)

		Segment_Data = []
		max_pos_key = None;
		min_pos_key = None
		seq1_max_2 = -1;
		seq2_max_2 = -1;
		seq1_min_2 = 100000;
		seq2_min_2 = 100000;


		for k,v in Segment_Dict_full.items():

			#If segments are in the band matrix
			if not ((B1 < k[0]*(m/n) - k[1]) or ( B2 < k[1] - k[0]*(m/n) )):
				
				Segment_Dict_OffDiag[k]=v
				Segment_Data.append([k[0],k[1]])

				seq1_res = k[0]
				seq2_res = k[1]

				seq1_res_s = v[1]
				seq2_res_s = v[2]

				#FInd the maximum residue number in each sequence
				if seq1_res>seq1_max_2:
					seq1_max_2=seq1_res

				if seq2_res>seq2_max_2:
					seq2_max_2=seq2_res

				#FInd the minimum residue number in each sequence
				if seq1_res_s<seq1_min_2:
					seq1_min_2=seq1_res_s

				if seq2_res_s<seq2_min_2:
					seq2_min_2=seq2_res_s

		#Make key for pseudo node which is greater than every point. start position of traceback
		max_pos_key=tuple([seq1_max_2+1,seq1_max_2+1,seq2_max_2+1,seq2_max_2+1])
		min_pos_key=tuple([seq1_min_2-1,seq1_min_2-1,seq2_min_2-1,seq2_min_2-1])

		return Segment_Dict_OffDiag, Segment_Data, max_pos_key, min_pos_key

	def segment_traceback(self, tracekey, alignments, Segment_Dict_full,KDtree, Segment_Data, max_depth,sourceid,currentid):
		"""Summary
		
		Args:
		    tracekey (TYPE): Description
		    alignments (TYPE): Description
		    Segment_Dict_full (TYPE): Description
		    KDtree (TYPE): Description
		    Segment_Data (TYPE): Description
		    max_depth (TYPE): Description
		    sourceid (TYPE): Description
		    currentid (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		#Tracekey: [seq1_min, seq1_max, seq2_min, seq2_max]
		
		if len(nx.shortest_path(alignments,source=sourceid,target=currentid))>max_depth:
			return alignments

		seglist = self.get_next_segment_lists(tracekey,Segment_Dict_full,KDtree, Segment_Data)

		current_node_index=currentid
			
		for i, pos_key in enumerate(seglist):

			succ_seg_node_idx = Segment_Dict_full[tuple([pos_key[1],pos_key[3]])][4]

			succ_seg_node_score = alignments.nodes()[succ_seg_node_idx]['score']

			residue_gap_score = pos_key[5]
			
			#add distance weight
			alignments.add_edge(current_node_index, succ_seg_node_idx, weight = (succ_seg_node_score+residue_gap_score)/2.0, directed = True)
			currentid=succ_seg_node_idx

			alignments = self.segment_traceback(pos_key, alignments, Segment_Dict_full, KDtree, Segment_Data, max_depth,sourceid,currentid)

		return alignments

	def get_next_segment_lists(self, tracekey,Segment_Dict_full,KDtree, Segment_Data):
		"""Summary
		
		Args:
		    tracekey (TYPE): Description
		    Segment_Dict_full (TYPE): Description
		    KDtree (TYPE): Description
		    Segment_Data (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		query = KDtree.query([tracekey[0],tracekey[2]], k=self.__TOP_N_SEGMENTS)

		ValidKeyList=[]
		
		#For all the indexes of the nearest segments
		for i, q in enumerate(query[1]):
			
			#seg_key contains upper limit of segment, seg_value is the lower limit of segment
			seg_key = tuple(Segment_Data[q])
			seg_value = Segment_Dict_full[seg_key]

			#All segments that are within lower limit of current segment
			if seg_key[0]<tracekey[0] and seg_key[1]<tracekey[2]:

				#Add absolute value of gaps for a given succeesive segment
				gap = abs((tracekey[0] - seg_key[0]) - (tracekey[2] - seg_key[1]))/1.0   
				
				#Add distance at the end
				#seq1_min, seq1_max, seq2_min, seq2_max
				ValidKeyList.append([seg_value[1], seg_key[0], seg_value[2], seg_key[1], query[0][i], gap ])

		
		#This is for removing the succesive segments from the list of all segments to search for
		#Following code prevents the skipping of the segments
		ValidKeyList2=[]
		for seg1 in ValidKeyList:
			flag=1
			for seg2 in ValidKeyList:
				if seg1[1]<seg2[0] and seg1[3]<seg2[2]:
					flag=0
					break
			if flag==1:
				ValidKeyList2.append(seg1)
		
		return ValidKeyList2
	
	def plot_segments(self, Segment_Dict_full,ax,lwd,ssize,seg_color, score_cutoff):
		"""Summary
		
		Args:
		    Segment_Dict_full (TYPE): Description
		    ax (TYPE): Description
		    lwd (TYPE): Description
		    ssize (TYPE): Description
		    seg_color (TYPE): Description
		    score_cutoff (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		Data_Max=[]
		Data_Min=[]
		Data_Line=[];Data_Color=[]
		for k,v in Segment_Dict_full.items():

			##x,y of segment high
			Data_Max.append([k[0],k[1]])
			##x,y of segment low
			Data_Min.append([v[1],v[2]])

			Data_Line.append([tuple([v[1],v[2]]),tuple([k[0],k[1]])])
			Data_Color.append((v[0]+v[3])/2.0)

		Data_Max = np.array(Data_Max)
		Data_Min = np.array(Data_Min)
		Data_Line = np.array(Data_Line)
		Data_Color = np.array(Data_Color)
		
		plt.scatter(Data_Max[:,0],Data_Max[:,1],s=ssize)
		plt.scatter(Data_Min[:,0],Data_Min[:,1],s=ssize)
		
		lines=None
		if seg_color==None:

			norm = colors.Normalize(Data_Color.min(), Data_Color.max())
			coolmap = cm.get_cmap("cool")
			#seg_color2 = coolmap(norm(Data_Color))
			seg_color2 = coolmap(Data_Color)
			lc = mc.LineCollection(Data_Line,color=seg_color2, linewidths=lwd)
			lines = ax.add_collection(lc)
			sm = plt.cm.ScalarMappable(cmap=coolmap, norm=plt.Normalize(vmin=0, vmax=1))
			sm.set_array([])
			plt.colorbar(sm, ax=ax)
				
		else:
			#lc = mc.LineCollection(Data_Line,color=seg_color, linewidths=lwd,linestyle='dotted')
			lc = mc.LineCollection(Data_Line,color=seg_color, linewidths=lwd,linestyle='dashed')
			lines = ax.add_collection(lc)

		return

	def k_shortest_paths(self, G, source, target, k, weight=None):
		"""Summary
		
		Args:
		    G (TYPE): Description
		    source (TYPE): Description
		    target (TYPE): Description
		    k (TYPE): Description
		    weight (None, optional): Description
		
		Returns:
		    TYPE: Description
		"""
		return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))

	def get_tracepath_alignments(self, alignments, sourceid, targetid, seq1,seq2, topk):
		"""Summary
		
		Args:
		    alignments (TYPE): Description
		    sourceid (TYPE): Description
		    targetid (TYPE): Description
		    seq1 (TYPE): Description
		    seq2 (TYPE): Description
		    topk (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		TRACE=[]
		TRACE_SEGMENT=[]
		ALI1=[]
		ALI2=[]
		ALISCORE=[]

		Alignment_Paths = self.k_shortest_paths(alignments,source=sourceid,target=targetid, k=topk, weight='weight')

		print("Total Alignments:",len(Alignment_Paths))

		for path in Alignment_Paths:

			Rpath = path[::-1]
			
				
			ali1="";ali2="";tracepath=[];Aliscore=0;ali1e=0;ali2e=0;tracesegment=[]
			#Iterate over all the segments that constitute an alignment
			for n in Rpath:

				segment=alignments.nodes[n]['seq1seq2']
				
				ali1+="-"*(len(seq1[ali1e:segment[0]])-1)
				ali2+="-"*(len(seq2[ali2e:segment[2]])-1)
					
				ali1+=seq1[segment[0]-1:segment[1]]
				ali2+=seq2[segment[2]-1:segment[3]]

				ali1e = segment[1]
				ali2e = segment[3]
					
				seq1resnos = range(segment[0],segment[1]+1)
				seq2resnos = range(segment[2],segment[3]+1)
					
				tracepath+=list(zip(seq1resnos,seq2resnos))
				tracesegment.append([tuple([segment[0],segment[2]]),tuple([segment[1],segment[3]])])

			#Add end of alignments
			ali1+="-"*len(seq1[ali1e:])
			ali2+="-"*len(seq2[ali2e:])

			Aliscore = 0
			for i in range(len(Rpath)-1):
				Aliscore += alignments.edges[Rpath[i+1], Rpath[i]]['weight']
		
			tracepath=np.array(tracepath)

			TRACE.append(tracepath)
			TRACE_SEGMENT.append(tracesegment)
			ALI1.append(ali1)
			ALI2.append(ali2)
			ALISCORE.append(Aliscore)

		return TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE

	def plot_alignment_segments(self, TRACE_SEGMENTi,ax,lwd,seg_color):
		"""Summary
		
		Args:
		    TRACE_SEGMENTi (TYPE): Description
		    ax (TYPE): Description
		    lwd (TYPE): Description
		    seg_color (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		TRACE_SEGMENT2=[]
		#Do not use the pseudo max point (TRACE_SEGMENTi[:-1]) as a starting point for the alignment in the plots
		for i,val in enumerate(TRACE_SEGMENTi[:-1]):
			TRACE_SEGMENT2.append([tuple([TRACE_SEGMENTi[i][1][0],TRACE_SEGMENTi[i][1][1]]),
				tuple([TRACE_SEGMENTi[i+1][0][0],TRACE_SEGMENTi[i+1][0][1]])])
		
		
		#lc = mc.LineCollection(Data_Line,color=seg_color, linewidths=lwd,linestyle='dotted')
		lc = mc.LineCollection(TRACE_SEGMENTi,color=seg_color, linewidths=lwd,linestyle='solid')
		lines = ax.add_collection(lc)
		lc = mc.LineCollection(TRACE_SEGMENT2,color=seg_color, linewidths=lwd,linestyle='dotted')
		lines = ax.add_collection(lc)

		return

	def plot_tree(self, alignments,fname):
		"""Summary
		
		Args:
		    alignments (TYPE): Description
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		Leafs = [x for x in alignments.nodes() if alignments.out_degree(x)==0 and alignments.in_degree(x)==1]

		labeldict={};values=[]
		for n, d in alignments.nodes(data=True):
			
			if n in Leafs:
				labeldict[n]=[round(d['score'],2)]
				#labeldict[n]=[d['seq1seq2'][1],d['seq1seq2'][3]]
				values.append(d['score'])
			else:
				labeldict[n]=[d['seq1seq2'][1],d['seq1seq2'][3]]
				values.append(d['score'])


		pos = graphviz_layout(alignments, prog="dot") #prog="twopi" #prog="circo"
		#nx.draw(alignments, pos, labels=labeldict,with_labels=True,font_size=3)
		nx.draw_networkx(alignments, pos, node_size=30, arrowsize=10, labels=labeldict,
			with_labels=True,font_size=3,cmap=plt.get_cmap('cool'),node_color=values)
		
		values = np.array(values)
		coolmap = cm.get_cmap("cool")
		sm = plt.cm.ScalarMappable(cmap=coolmap, norm=plt.Normalize(vmin=values.min(), vmax=values.max()))
		plt.colorbar(sm)
		#plt.tight_layout()
		plt.savefig(fname)
		#plt.show()

		return

	def plot_tree_simple(self, alignments,fname):
		"""Summary
		
		Args:
		    alignments (TYPE): Description
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		pos = graphviz_layout(alignments, prog="dot")
		#nx.draw(alignments, pos, labels=labeldict,with_labels=True,font_size=3)
		nx.draw_networkx(alignments, pos)
		plt.savefig(fname)
		#plt.show()
		return

	def read_coco_alignment_graph(self, seq1, seq2, ppath,prot1_prot2):
		"""Summary
		
		Args:
		    seq1 (TYPE): Description
		    seq2 (TYPE): Description
		    ppath (TYPE): Description
		    prot1_prot2 (TYPE): Description
		"""
		#get the length of the sequences in two variables, n and m
		n = len(seq1)
		m = len(seq2)
		
		print("Length",n,m)

		print("Loading Coiled-coil segments")
		Segment_Dict_full = self.load_Coiled_Coil_Segments_full(ppath+prot1_prot2+'_Coiled-Coil_Segments_parallel.dat')
		Segment_Dict_full[tuple([0,0])]=[1,0,0];
		
		with open(ppath+prot1_prot2+'_alignment_graph.gpickle', 'rb') as f:
			alignments = pickle.load(f)

		print("Ploting alignment")
		sourceid=0
		TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE = self.get_tracepath_alignments(alignments,sourceid, seq1,seq2)

		fig,ax = plt.subplots()
		score_cutoff=0.0
		#plot_segments(Segment_Dict_full,ax,0.5,None,score_cutoff)

		plt.scatter(0,0)
		plt.scatter(len(seq1),len(seq2))
		print("Total Alignments:",len(TRACE))
		for i,tracepath in enumerate(TRACE):
			self.plot_alignment_segments(TRACE_SEGMENT[i],ax,0.7,"gray")

		plt.savefig(ppath+prot1_prot2+'_Aligned_score.pdf')
		#plt.show()
		
		fig,ax = plt.subplots(figsize=(200,12))
		self.plot_tree(alignments,ppath+prot1_prot2+'_alignment_graph.pdf')

	def get_hamming_distances(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		HammDistances=[]
		for i,tracepath_i in enumerate(TRACE):
			"""<FRESHLY_INSERTED>"""
			print("SCORE:",ALISCORE[i], "CCDom:",tracepath_i[:,0].shape[0])
			print("S1: ",ALI1[i])
			print("S2: ",ALI2[i])
			print("")
			
			temp_hamm=[]
			for j,tracepath_j in enumerate(TRACE):

				hamm_dist = hamming([ALI1[i],ALI1[j]],[ALI2[i],ALI2[j]])

				temp_hamm.append(hamm_dist)

				lines1,lines2="",""
				for s1,s2 in zip(ALI1[i],ALI2[i]):
					if s1==s2:
						lines1+='-'
					else:
						lines1+='*'
				for s1,s2 in zip(ALI1[j],ALI2[j]):
					if s1==s2:
						lines2+='-'
					else:
						lines2+='*'
			HammDistances.append(temp_hamm)
		
		plt.imshow(HammDistances,cmap="gist_earth")
		plt.savefig(ppath+prot1_prot2+'_alignment_hamming.pdf')
		#plt.show()
		

		return

	def write_alignment(self, alignment,prot1, prot2, outfilename):
		"""Summary
		
		Args:
		    alignment (TYPE): Description
		    prot1 (TYPE): Description
		    prot2 (TYPE): Description
		    outfilename (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		outf = open(outfilename,'w')
		outf.write('\n>'+prot1+'\n')
		outf.write(alignment[0])
		outf.write('\n')
		outf.write('\n>'+prot2+'\n')
		outf.write(alignment[1])
		outf.close()
		return

	def write_alignment_with_heptad(self, alignment, prot1, prot2, p1_PC, p2_PC, outfilename):
		"""Summary
		
		Args:
		    alignment (TYPE): Description
		    prot1 (TYPE): Description
		    prot2 (TYPE): Description
		    p1_PC (TYPE): Description
		    p2_PC (TYPE): Description
		    outfilename (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		resnum=0;alilines1,heplines1="",""
		
		for ch in alignment[0]:
			if ch=='-':
				alilines1+=ch
				heplines1+='-'
			else:

				resnum+=1;
				assert p1_PC[resnum][0]==ch
				alilines1+=ch
				heplines1+=p1_PC[resnum][1]
		
		resnum=0;alilines2,heplines2="",""
		for ch in alignment[1]:
			if ch=='-':
				alilines2+=ch
				heplines2+='-'
			else:
				resnum+=1
				assert p2_PC[resnum][0]==ch
				alilines2+=ch
				heplines2+=p2_PC[resnum][1]


		outf = open(outfilename,'w')
		outf.write('\n>'+prot1+'\n')
		outf.writelines(alilines1)
		outf.write('\n')
		outf.write('\n>'+prot1+'_hep\n')
		outf.writelines(heplines1)
		outf.write('\n')
		outf.write('\n>'+prot2+'\n')
		outf.writelines(alilines2)
		outf.write('\n')
		outf.write('\n>'+prot2+'_hep\n')
		outf.writelines(heplines2)
		outf.write('\n')
		
		outf.close()
		return

	def read_coco_alignment_graph_new(self, seq1, seq2, ppath,prot1_prot2, topk):
		"""Summary
		
		Args:
		    seq1 (TYPE): Description
		    seq2 (TYPE): Description
		    ppath (TYPE): Description
		    prot1_prot2 (TYPE): Description
		    topk (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		#get the length of the sequences in two variables, n and m
		n = len(seq1)
		m = len(seq2)
		
		print("Length",n,m)

		print("Loading Coiled-coil segments")
		Segment_Dict_full = self.load_Coiled_Coil_Segments_full(ppath+prot1_prot2+'_Coiled-Coil_Segments_parallel.dat')
		Segment_Dict_full[tuple([0,0])]=[1,0,0];
		
		with open(ppath+prot1_prot2+'_alignment_graph.gpickle', 'rb') as f:
			alignments = pickle.load(f)

		
		ALLnodes = alignments.nodes()
		sourceid=max(ALLnodes)
		targetid=0

		print("Ploting alignment")
		sourceid=0
		TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE = self.get_tracepath_alignments(alignments,sourceid, targetid, seq1,seq2, topk)

		
		fig,ax = plt.subplots()
		score_cutoff=0.0

		plt.scatter(0,0)
		plt.scatter(len(seq1),len(seq2))
		print("Total Alignments:",len(TRACE))
		for i,tracepath in enumerate(TRACE):
			self.plot_alignment_segments(TRACE_SEGMENT[i],ax,0.7,"gray")

		plt.savefig(ppath+prot1_prot2+'_Aligned_score.pdf')
		#plt.show()
		
		fig,ax = plt.subplots(figsize=(200,12))
		#self.plot_tree(alignments,ppath+prot1_prot2+'_alignment_graph.pdf')
		
		oaname=self.ppath+self.prot1+'_'+self.prot2+'_alignment.cali'
		topn=1
		self.Write_Topn_Alignments(oaname, topn, self.prot1, self.prot2, seq1, seq2, TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE)

		TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE = self.get_tracepath_alignments_fullseq(alignments,sourceid, targetid, seq1,seq2, topk)
		oaname=self.ppath+self.prot1+'_'+self.prot2+'_alignment.xali'
		topn=1
		self.Write_Topn_Alignments(oaname, topn, self.prot1, self.prot2, seq1, seq2, TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE)

		return

	def Write_Topn_Alignments(self, oaname, topn, prot1, prot2, seq1, seq2, TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE):
		"""Summary
		
		Args:
		    oaname (TYPE): Description
		    topn (TYPE): Description
		    prot1 (TYPE): Description
		    prot2 (TYPE): Description
		    seq1 (TYPE): Description
		    seq2 (TYPE): Description
		    TRACE (TYPE): Description
		    TRACE_SEGMENT (TYPE): Description
		    ALI1 (TYPE): Description
		    ALI2 (TYPE): Description
		    ALISCORE (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		owlines=""

		#Obtain pseudo-crosslinks possible for every alignment
		for i,tracepath in enumerate(TRACE[0:topn]):

			SeqSegments = TRACE_SEGMENT[i]
			Alimentscore = ALISCORE[i]

			owlines+="#"+prot1+"\n"
			owlines+=ALI1[i]+"\n"
			owlines+="#"+prot2+"\n"
			owlines+=ALI2[i]+"\n\n"
			owlines+="#Segments as rigid bodies \n"

			zcoord = 0.0
			for seg in SeqSegments:

				seg1 = [seg[0][0] , seg[1][0]]
				seg2 = [seg[0][1] , seg[1][1]]

				#If it is a pseudo start/end segments (added for alignment)
				if seg1[0]==seg1[1] and seg2[0]==seg2[1]:
					continue

				
				pline1 = prot1+"\t"+str(seg1[0])+"-"+str(seg1[1])
				pline2 = prot2+"\t"+str(seg2[0])+"-"+str(seg2[1])
				owlines+= pline1 + "\t" + pline2 +"\t"
				owlines+= seq1[seg1[0]-1:seg1[1]] +"\t"+seq2[seg2[0]-1:seg2[1]]
				owlines+= "\n"

		outf = open(oaname,'w')
		outf.writelines(owlines)
		outf.close()

		return

	def add_sequence_at_cterm(self, ali1,ali2,seq1,seq2,ali1e,ali2e,segment):
		"""Summary
		
		Args:
		    ali1 (TYPE): Description
		    ali2 (TYPE): Description
		    seq1 (TYPE): Description
		    seq2 (TYPE): Description
		    ali1e (TYPE): Description
		    ali2e (TYPE): Description
		    segment (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		lenccseg1 = len(seq1[ali1e:])
		lenccseg2 = len(seq2[ali2e:])

		if lenccseg1<lenccseg2:
			ali1+=seq1[ali1e:]
			ali2+=seq2[ali2e:]
			ali1+="-"*(lenccseg2-lenccseg1)
			

		elif lenccseg1>lenccseg2:
			ali1+=seq1[ali1e:]
			ali2+=seq2[ali2e:]
			ali2+="-"*(lenccseg1-lenccseg2)
			

		else:
			ali1+=seq1[ali1e:]
			ali2+=seq2[ali2e:]

		return ali1,ali2

	def add_sequence_with_gaps(self, ali1,ali2,seq1,seq2,ali1e,ali2e,segment):
		"""Summary
		
		Args:
		    ali1 (TYPE): Description
		    ali2 (TYPE): Description
		    seq1 (TYPE): Description
		    seq2 (TYPE): Description
		    ali1e (TYPE): Description
		    ali2e (TYPE): Description
		    segment (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		lenccseg1 = len(seq1[ali1e:segment[0]])
		lenccseg2 = len(seq2[ali2e:segment[2]])

		if lenccseg1<lenccseg2:
			ali1+="-"*(lenccseg2-lenccseg1)
			ali1+=seq1[ali1e:segment[0]-1]
			ali2+=seq2[ali2e:segment[2]-1]

		elif lenccseg1>lenccseg2:
			ali2+="-"*(lenccseg1-lenccseg2)
			ali1+=seq1[ali1e:segment[0]-1]
			ali2+=seq2[ali2e:segment[2]-1]

		else:

			if ((ali1e<segment[0]-1) and (ali2e<segment[2]-1)):

				ali1+=seq1[ali1e:segment[0]-1]
				ali2+=seq2[ali2e:segment[2]-1]

		return ali1,ali2

	def get_tracepath_alignments_fullseq(self, alignments, sourceid, targetid, seq1,seq2, topk):
		"""Summary
		
		Args:
		    alignments (TYPE): Description
		    sourceid (TYPE): Description
		    targetid (TYPE): Description
		    seq1 (TYPE): Description
		    seq2 (TYPE): Description
		    topk (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		TRACE=[]
		TRACE_SEGMENT=[]
		ALI1=[]
		ALI2=[]
		ALISCORE=[]

		Alignment_Paths = self.k_shortest_paths(alignments,source=sourceid,target=targetid, k=topk, weight='weight')

		print("Total Alignments:",len(Alignment_Paths))

		for path in Alignment_Paths:

			Rpath = path[::-1]
			
				
			ali1="";ali2="";tracepath=[];Aliscore=0;ali1e=0;ali2e=0;tracesegment=[]
			#Iterate over all the segments that constitute an alignment
			for n in Rpath:

				segment=alignments.nodes[n]['seq1seq2']
				
				#ali1+="-"*(len(seq1[ali1e:segment[0]])-1)
				#ali2+="-"*(len(seq2[ali2e:segment[2]])-1)

				ali1,ali2 = self.add_sequence_with_gaps(ali1,ali2,seq1,seq2,ali1e,ali2e,segment)

				#ali1+=seq1[ali1e:segment[0]]
				#ali2+=seq2[ali2e:segment[2]]
				

				ali1+=seq1[segment[0]-1:segment[1]]
				ali2+=seq2[segment[2]-1:segment[3]]

				#ali1+=str(segment[0])+"-"+str(segment[1])+"\t"
				#ali2+=str(segment[2])+"-"+str(segment[3])+"\t"

				ali1e = segment[1]
				ali2e = segment[3]
					
				seq1resnos = range(segment[0],segment[1]+1)
				seq2resnos = range(segment[2],segment[3]+1)
					
				tracepath+=list(zip(seq1resnos,seq2resnos))
				tracesegment.append([tuple([segment[0],segment[2]]),tuple([segment[1],segment[3]])])

			#Add end of alignments
			ali1,ali2 = self.add_sequence_at_cterm(ali1,ali2,seq1,seq2,ali1e,ali2e,segment)

			Aliscore = 0
			for i in range(len(Rpath)-1):
				Aliscore += alignments.edges[Rpath[i+1], Rpath[i]]['weight']

			#Leaf has the full score of the alignment
			#Aliscore=alignments.nodes[Rpath[0]]['score']
					
			tracepath=np.array(tracepath)

			TRACE.append(tracepath)
			TRACE_SEGMENT.append(tracesegment)
			ALI1.append(ali1)
			ALI2.append(ali2)
			ALISCORE.append(Aliscore)

		return TRACE, TRACE_SEGMENT, ALI1, ALI2, ALISCORE


