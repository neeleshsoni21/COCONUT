
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

from CoiledCoilClass import CoiledCoil
from collections import OrderedDict
import sys

class InputCC:

	"""Summary
	
	Attributes:
	    input_filename (TYPE): Description
	    interface_dict (TYPE): Description
	    num_monomers (int): Description
	    valid_aas (TYPE): Description
	    valid_hep (list): Description
	"""
	
	def __init__(self, input_filename):
		"""Summary
		
		Args:
		    input_filename (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		#TODO: Add helper and sampling functions in the class for a given input

		self.input_filename = input_filename

		self.valid_aas=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];

		self.valid_hep=['a','b','c','d','e','f','g']

		#hard coded for dimer for now
		self.num_monomers = 2
		
		return

	def raise_error(self,flag):
		"""Summary
		
		Args:
		    flag (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		print("Error",flag)
		sys.exit(1)
		return

	def check_inp_len(self, aa_seq, aa_hep):
		"""Summary
		
		Args:
		    aa_seq (TYPE): Description
		    aa_hep (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		if len(aa_seq)<11:
			print(aa_seq)
			print("WARNING: Sequence should be greater than 11 residues")
			self.raise_error(1)

		if len(aa_seq)!=len(aa_hep):
			print("WARNING: Heptads should be equal")
			self.raise_error(2)			

		return

	def check_aa_seq(self, aa_seq):
		"""Summary
		
		Args:
		    aa_seq (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		for aa in aa_seq:
			if aa not in self.valid_aas:
				
				print("WARNING: Heptad should be from a/b/c/d/e/f/g")
				self.raise_error(3)

		return

	def check_aa_hep(self, aa_hep):
		"""Summary
		
		Args:
		    aa_hep (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		for hp in aa_hep:
			if hp not in self.valid_hep:
				self.raise_error(4)

		return

	def get_cc_orientation(self,value):
		"""Summary
		
		Args:
		    value (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		if value == 'parallel':
			return 1
		elif value == 'antiparallel':
			return -1
		else:
			print("WARNING: Orientation should be either 'parallel' or 'antiparallel'")
			self.raise_error(5)

	def correct_P_AP_compatibility(self, dseq, dhep):
		"""Summary
		
		Args:
		    dseq (TYPE): Description
		    dhep (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		if dhep[0]=='a' and dhep[-1]=='d':
			pass
		elif dhep[0]=='d' and dhep[-1]=='a':
			pass
		else:
			print("WARNING: Heptad repeat should start with 'a' and ends with 'd' or vice-versa")
			self.raise_error(6)

		return

	def print_input_info(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		print("Input Sequences:")
		#Iterator for the input variables
		for CCid, CC_obj in self.interface_dict.items():
			print("CC id:",CCid,"#Monomers:",CC_obj.get_num_monomers());
			print("Details:",CC_obj.get_CC());
			print()
			#CC_obj.display_CC()

		return

	def load_CCData(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		self.interface_dict=OrderedDict()

		inf = open(self.input_filename,'r')
		lines = inf.readlines()
		inf.close()

		P_counter=0;AP_counter=0;

		iteration=0;

		#Read the input file line by line
		for l in lines:
			
			#skip all lines starting with '#'
			if ((l[0]=='#') or (len(l)==1)):
				continue

			#split the line with <space> or tabs
			toks = l.strip().split()

			iteration+=1

			CCinp = CoiledCoil(iteration,self.num_monomers)

			CCinp.set_CC_orientation(toks[-1])

			for monomer_i in range(self.num_monomers):

				idx_seq = monomer_i*2
				idx_hep = monomer_i*2 + 1

				#Check the length of monomers
				self.check_inp_len(toks[idx_seq],toks[idx_hep]);
				
				#Validate the amino acid sequece
				self.check_aa_seq(toks[idx_seq]);

				#Validate the amino acid heptads
				self.check_aa_hep(toks[idx_hep]);

				self.correct_P_AP_compatibility(toks[idx_seq],toks[idx_hep])
				
				#Set the values of coiled-coil seq and heptad
				CCinp.set_CC_values(monomer_i, toks[idx_seq],toks[idx_hep])
				
			self.interface_dict[iteration]=CCinp
			
		return 

