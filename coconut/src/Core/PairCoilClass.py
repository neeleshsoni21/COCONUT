
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


from collections import OrderedDict, defaultdict

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

class PairCoilParse:

	"""Summary
	
	Attributes:
	    header (str): Description
	    p1_PC (TYPE): Description
	    pc_prot1_file (TYPE): Description
	    phep1 (str): Description
	    phep2 (str): Description
	    PROB_CUTOFF (TYPE): Description
	"""
	
	def __init__(self, pc_prot1_file, PROB_CUTOFF, probability_type = 'pcoil'):
		"""Summary
		
		Args:
		    pc_prot1_file (TYPE): Description
		    PROB_CUTOFF (TYPE): Description
		    probability_type (str, optional): Description
		
		Returns:
		    TYPE: Description
		"""
		self.PROB_CUTOFF = PROB_CUTOFF

		self.pc_prot1_file = pc_prot1_file

		self.phep1 = 'a'
		self.phep2 = 'd'

		self.header =''

		if probability_type == 'paircoil':

			self.p1_PC = self.load_paircoil_alignment(self.pc_prot1_file)

		elif probability_type == 'pcoil':
			self.p1_PC = self.load_Pcoil_alignment(self.pc_prot1_file)

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

		prot_seq = defaultdict(self.defaultaa)

		for l in lines:
			if l[0]=='#':
				self.header+=l
				continue

			if len(l)<2:
				continue

			toks = l.strip().split()
		
			if float(toks[3]) >= self.PROB_CUTOFF:
				prot_seq[int(toks[0])]=[toks[1],toks[2],float(toks[3])]
			else:
				prot_seq[int(toks[0])]=[toks[1],'-',float(toks[3])]

		return prot_seq

	def load_Pcoil_alignment(self, pc_file):
		"""Summary
		
		Args:
		    pc_file (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		prob_col1=1;prob_col2=2;prob_col3=3;
		
		inf = open(pc_file,'r')
		lines = inf.readlines()
		inf.close()

		prot_seq = DefaultDictNoAdd()

		for l in lines:

			if l[0]=='>':
				continue
			toks = l.strip().split()

			if len(toks)<8:
				self.header+='#'+l
				continue
			
			probidx1 = 2*prob_col1+1;heptadidx1 = 2*prob_col1
			probidx2 = 2*prob_col2+1;heptadidx2 = 2*prob_col2
			probidx3 = 2*prob_col3+1;heptadidx3 = 2*prob_col3

			if float(toks[probidx3]) >= self.PROB_CUTOFF:
				prot_seq[int(toks[0])]=[toks[1],toks[heptadidx3],float(toks[probidx3])]

			elif float(toks[probidx2]) >= self.PROB_CUTOFF:
				prot_seq[int(toks[0])]=[toks[1],toks[heptadidx2],float(toks[probidx2])]

			elif float(toks[probidx1]) >= self.PROB_CUTOFF:
				prot_seq[int(toks[0])]=[toks[1],toks[heptadidx1],float(toks[probidx1])]
			
			else:
				prot_seq[int(toks[0])]=[toks[1],'-',float(toks[probidx1])]

		return prot_seq

	def load_Pcoil_alignment_backup(self, pc_file, prob_col=3):
		"""Summary
		
		Args:
		    pc_file (TYPE): Description
		    prob_col (int, optional): Description
		
		Returns:
		    TYPE: Description
		"""
		prob_col1=1;prob_col2=1;prob_col3=1;

		inf = open(pc_file,'r')
		lines = inf.readlines()
		inf.close()

		prot_seq = DefaultDictNoAdd()

		for l in lines:

			if l[0]=='>':
				continue

			toks = l.strip().split()

			if len(toks)<8:
				self.header+='#'+l
				continue
			
			probidx = 2*prob_col+1
			heptadidx = 2*prob_col
			if float(toks[probidx]) >= self.PROB_CUTOFF:
				prot_seq[int(toks[0])]=[toks[1],toks[heptadidx],float(toks[probidx])]
			else:
				prot_seq[int(toks[0])]=[toks[1],'-',float(toks[probidx])]

		return prot_seq

	def get_pc_seq(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.p1_PC

	def hamming_distance(self, s1, s2, resrange):
		"""Summary
		
		Args:
		    s1 (TYPE): Description
		    s2 (TYPE): Description
		    resrange (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		HamDist = 0
		DiffPos = [];
		DiffHep = [];
		hepnames = ['a','b','c','d','e','f','g']
		
		for i, (c1, c2) in enumerate(zip(s1, s2)):
			if c1!=c2:
				HamDist+=1
				DiffPos.append(resrange[i])
				DiffHep.append(hepnames[i])

		#Hamming Distance is len(DiffPos), which can be calculated at later stages from Diff Pos
		return DiffPos, DiffHep

	def get_heptad_repeats_pairs(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		heptad_seqs = []

		#Get all pairs started by heptad 'a'
		cc_flag=0
		for k,v in self.p1_PC.items():
			register = v[1]
			v.append('-')
			#v.append(None)
			
			if register == '-':
				cc_flag=0
				heptad_repeat=[]
				continue
			else:
				cc_flag=1
	
			if cc_flag==1:
				if register=='a':
					heptemp = []
					for i in range(k,k+7):
						heptemp.append(self.p1_PC[i][1])
					heptad_seqs.append([tuple([k,k+6]), ''.join(heptemp)])
		
		return heptad_seqs

	def correct_incorrectly_Indentified_heptad(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		Heptad_pairs = OrderedDict()

		#heptad_repeat=[]
		hep_marker = 0

		#Get all heptad repeat sequences 'abcdefg' starting with a and length 7
		heptad_seqs = self.get_heptad_repeats_pairs()

		#------------------------------------------------------
		# #Get all heptad repeats that are not contigous
		# Change those with hamming distance == 1. All heptad with more changes skipped
		#------------------------------------------------------
		for i,v in enumerate(heptad_seqs):
			
			#Skip all seqs that are complete heptad repeats or gaps identified by '-'
			if ((v[1]=='abcdefg') or ('-' in v[1])):
				continue

			Residue_range = list(range(v[0][0],v[0][1]+1))
			#Get hamming distance
			DiffPos, DiffHep = self.hamming_distance(v[1],'abcdefg',Residue_range)

			#Ski heptad seqs with hamming distNCE >1
			if len(DiffPos)>1:
				continue

			#Modify the heptad position such that it completes the heptad repeat
			for i,changepos in enumerate(DiffPos):
				self.p1_PC[changepos][1] = DiffHep[i]

		#------------------------------------------------------
		# Find and Add markers for Handecad repeat
		#------------------------------------------------------
		heptad_seqs_new = self.get_heptad_repeats_pairs()
		Complete_HeptadSeqs = []
		for i,v in enumerate(heptad_seqs_new):
			if (v[1]=='abcdefg'):
				Complete_HeptadSeqs.append(v)

		for i in range(0,len(Complete_HeptadSeqs)-1):
			vi = Complete_HeptadSeqs[i]
			vi1 = Complete_HeptadSeqs[i+1]

			resdiff = vi1[0][0] - vi[0][0]
			if resdiff==11:
				
				for resnum in range(vi[0][0],vi[0][1]+1 + 4):

					self.p1_PC[resnum][3] = 'Handecad'
			else:
				for resnum in range(vi[0][0],vi[0][1]+1):
					self.p1_PC[resnum][3] = 'Heptad'
		return

	def write_Updated_Paircoil_Probs(self, fname):
		"""Summary
		
		Args:
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		wlines = '#File Written by the COCONUT code\n'
		wlines += self.header

		for k, v in self.p1_PC.items():

			try:
				marker = ' '+str(v[3]).rjust(8)
				wlines+=str(k).rjust(8)+' '+str(v[0]).rjust(8)+' '+str(v[1]).rjust(14)+' '+"{:.5f}".format(v[2]).rjust(14)+marker+'\n'
			except:
				wlines+=str(k).rjust(8)+' '+str(v[0]).rjust(8)+' '+str(v[1]).rjust(14)+' '+"{:.5f}".format(v[2]).rjust(14)+'\n'

		outf = open(fname,'w')
		outf.writelines(wlines)
		outf.close()

		return








