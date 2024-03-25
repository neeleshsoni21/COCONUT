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


class CoiledCoil:

	"""Summary
	
	Attributes:
	    CC (list): Description
	"""
	
	def __init__(self, CCid, num_monomers=2):
		"""Summary
		
		Args:
		    CCid (TYPE): Description
		    num_monomers (int, optional): Description
		
		Returns:
		    TYPE: Description
		"""
		self.CC = []

		self.__CC_num_monomers = num_monomers;

		#Add default 1, 1 for parallel dimer, -1 for anti-paralell dimers
		self.CC.append(1)

		for monomer in range(0,self.__CC_num_monomers):

			monomer = []
			#1st element for sequence
			monomer.append("")
			#2nd element for hep
			monomer.append("")

			self.CC.append(monomer)
		
		#Add default label 'true' to the coiled-coil. Useful for testing/validating datasets
		#self.CC.append("true")

		return

	def set_CC_orientation(self, orientation):
		"""Summary
		
		Args:
		    orientation (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		self.CC[0]=orientation

		return

	def set_CC_values(self,monomer_idx, monomer_seq,monomer_hep):
		"""Summary
		
		Args:
		    monomer_idx (TYPE): Description
		    monomer_seq (TYPE): Description
		    monomer_hep (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		#Increment the index by one to accomodate the orientation int he list
		self.CC[monomer_idx+1][0]=monomer_seq
		
		self.CC[monomer_idx+1][1]=monomer_hep

		return

	def get_CC(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.CC

	def get_num_monomers(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.__CC_num_monomers

	def display_CC(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		print("Orientation:",self.CC[0])
		for monomer in range(1,self.__CC_num_monomers+1):
			print(self.CC[monomer])
		return