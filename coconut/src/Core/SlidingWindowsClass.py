
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

class SlidingWindows:

	"""Summary
	
	Attributes:
	    pc_prot1_file (TYPE): Description
	    pc_prot2_file (TYPE): Description
	    PROB_CUTOFF (TYPE): Description
	    prot1_seq (TYPE): Description
	    prot2_seq (TYPE): Description
	    sw_file_a (TYPE): Description
	    sw_file_d (TYPE): Description
	    swlen (TYPE): Description
	    swlen_d (TYPE): Description
	"""
	
	def __init__(self, pc_prot1_file, pc_prot2_file, sw_file_a, sw_file_d, PROB_CUTOFF=0.2, swlen=18):
		"""Summary
		
		Args:
		    pc_prot1_file (TYPE): Description
		    pc_prot2_file (TYPE): Description
		    sw_file_a (TYPE): Description
		    sw_file_d (TYPE): Description
		    PROB_CUTOFF (float, optional): Description
		    swlen (int, optional): Description
		
		Returns:
		    TYPE: Description
		"""
		self.pc_prot1_file = pc_prot1_file
		self.pc_prot2_file = pc_prot2_file

		self.sw_file_a = sw_file_a
		self.sw_file_d = sw_file_d

		self.swlen = swlen

		self.PROB_CUTOFF = PROB_CUTOFF

		self.prot1_seq = self.load_paircoil_alignment(self.pc_prot1_file)
		self.prot2_seq = self.load_paircoil_alignment(self.pc_prot2_file)

		self.GLOBAL_SeqHep1_a, self.GLOBAL_SeqHep1_d = self.Get_Sliding_Windows(self.prot1_seq)
		self.GLOBAL_SeqHep2_a, self.GLOBAL_SeqHep2_d = self.Get_Sliding_Windows(self.prot2_seq)

		self.Write_CC_InputFile(self.GLOBAL_SeqHep1_a, self.GLOBAL_SeqHep2_a, self.sw_file_a)
		self.Write_CC_InputFile(self.GLOBAL_SeqHep1_d, self.GLOBAL_SeqHep2_d,self.sw_file_d)

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

	def check_heptad_repeat(self, hep):
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
			#print(ph,h, ph_n, h_n, ph_n==h_n,flag )
			ph=h
		return flag

	def Get_Sliding_Windows(self, prot_seq):
		"""Summary
		
		Args:
		    prot_seq (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		MLPSeq=[]
		MLPHep=[]
		MLPRes=[]

		for k, v in prot_seq.items():
			MLPRes.append(k)
			MLPSeq.append(v[0])
			MLPHep.append(v[1])
		
		Hepas = [];Hepds = [];
		for i, h in enumerate(MLPHep):
			if h=='a':
				Hepas.append(i)
			if h=='d':
				Hepds.append(i)

		GLOBAL_SeqHep_a=[]
		GLOBAL_SeqHep_d=[]
		
		for idx in Hepas:
			
			seq1= ''.join(MLPSeq[idx:idx+self.swlen]);
			hep1= ''.join(MLPHep[idx:idx+self.swlen]);
			
			if idx+self.swlen>=len(MLPRes):
				continue
			start_res = MLPRes[idx]
			end_res = MLPRes[idx+self.swlen-1]

			if self.check_heptad_repeat(hep1)==True and len(hep1)==self.swlen:

				#print([ start_res, end_res ])
				GLOBAL_SeqHep_a.append([ seq1, hep1, start_res, end_res ])
				"""
				#Add e/f/g in the extended length, NOT in the sequence
				#if idx+self.swlen-1+3<len(MLPHep):
				if idx+self.swlen<len(MLPHep):

					#if MLPHep[idx+self.swlen-1+3]!='-':
					if MLPHep[idx+self.swlen]!='-':

						#end_res_ext = idx+self.swlen-1+3;
						end_res_ext = idx+self.swlen;
						GLOBAL_SeqHep_a.append([ seq1, hep1, start_res, end_res_ext ])
				else:
					GLOBAL_SeqHep_a.append([ seq1, hep1, start_res, end_res ])
				"""

		self.swlen_d = self.swlen + 1
		for idx in Hepds:
			
			seq1= ''.join(MLPSeq[idx:idx+self.swlen_d]);
			hep1= ''.join(MLPHep[idx:idx+self.swlen_d]);

			if idx+self.swlen_d>=len(MLPRes):
				continue
			start_res = MLPRes[idx]
			end_res = MLPRes[idx+self.swlen_d-1]

			if self.check_heptad_repeat(hep1)==True and len(hep1)==self.swlen_d:

				#print([ start_res, end_res ])
				GLOBAL_SeqHep_d.append([ seq1, hep1, start_res, end_res ])
				"""
				
				#Add b/c in the extended length, NOT in the sequence
				#if idx+self.swlen_d-1+2<len(MLPHep):
				if idx+self.swlen_d<len(MLPHep):

					#if MLPHep[idx+self.swlen_d-1+2]!='-':
					if MLPHep[idx+self.swlen_d]!='-':

						#end_res_ext = idx+self.swlen_d-1+2
						end_res_ext = idx+self.swlen_d
						GLOBAL_SeqHep_d.append([ seq1, hep1, start_res, end_res_ext ])
				else:
					GLOBAL_SeqHep_d.append([ seq1, hep1, start_res, end_res ])
				"""

		return GLOBAL_SeqHep_a, GLOBAL_SeqHep_d

	def Write_CC_InputFile(self, SeqHep1, SeqHep2, fname):
		"""Summary
		
		Args:
		    SeqHep1 (TYPE): Description
		    SeqHep2 (TYPE): Description
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		wlines=""
		for i, S1 in enumerate(SeqHep1):
			for j, S2 in enumerate(SeqHep2):
				
				label = str(S1[2])+'-'+str(S1[3])+'-'+str(S2[2])+'-'+str(S2[3])
				wlines+=S1[0]+'\t'+S1[1]+'\t'+S2[0]+'\t'+S2[1]+'\t'+label+'\n'

		outf = open(fname,'w')
		outf.writelines(wlines)
		outf.close()

		return




