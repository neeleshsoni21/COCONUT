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


class BaseScoresAlign:

	"""Summary
	
	Attributes:
	    scoremat (dict): Description
	"""
	
	def __init__(self, fname):
		"""Summary
		
		Args:
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		#Score matrix variable for storing substitution scores
		self.scoremat={}

		self.read_matrix_file(fname)

		return

	def read_matrix_file(self, fname):
		"""Summary
		
		Args:
		    fname (TYPE): Description
		
		Returns:
		    TYPE: Description
		"""
		
		#Open the file in reading mode
		inf = open(fname,'r')
		lines = inf.readlines();
		inf.close();

		#temperory list for storing
		tempmat=[];

		#iterate over all lines
		for l in lines:
			#skip if first character in the file is '#'
			if l[0]=='#':
				continue;

			#first remove the trailing '\n' and then split the lines 
			#for removing multiple spaces as delimeter
			toks = l.strip().split();

			#append the lins converted to a list
			tempmat.append(toks)

		
		#Get all rows except first that contains matrix headers
		submat = tempmat[1:]

		#Iterate over all rows
		for i in range(0,len(submat)):
			#Iterate over columns
			for j in range(0,len(tempmat[0])):

				#get column res
				colid = tempmat[0][j]
				#get row res
				rowid = submat[i][0]
				#get the value of res-res pair from substitution matrix
				value = submat[i][j+1]
				
				#Subtract -11.0 to make the scores all negative -15 to 0. divide by 15 (maximum number)
				#to make best scores 1 and worst score 0
				self.scoremat[(rowid,colid)]= 1 + 1*(float(value)-11)/15.0;

		return

	def get_scoremat(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.scoremat