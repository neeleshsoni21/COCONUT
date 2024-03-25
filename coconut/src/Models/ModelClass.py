
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

from glob import glob
from joblib import load
from collections import OrderedDict
from read_frequency_scores import load_AA_Frequency_scores

import os

_MDL_ROOT = os.path.abspath(os.path.dirname(__file__))

class Models:

	"""Summary
	
	Attributes:
	    Hexscores (TYPE): Description
	    RF_Models (list): Description
	    RF_PaddingSize (TYPE): Description
	    scoring_dicts (TYPE): Description
	"""
	
	def __init__(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		self.RF_Models= []

		self.scoring_dicts = OrderedDict()

		self.Hexscores = None

		self.RF_PaddingSize=None;

		return

	def load_RF_models(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		#Load and sort the models using name. NOT NUMERICALLY
		mdl_files = sorted(glob(_MDL_ROOT+'/trained_network/*.pkl'))

		#Load the models
		for mdl in mdl_files:
			self.RF_Models.append( load(mdl) )

		return

	def load_AA_distributions(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		#make default dict for b/c/e/f/g positions that is also a default dict with score 0
		#scoring_dicts = defaultdict(lambda: defaultdict(lambda: 0))
		

		self.scoring_dicts['aa'] = load_AA_Frequency_scores(_MDL_ROOT+'/Distributions/pairwise_aa_residues.out')
		self.scoring_dicts['dd'] = load_AA_Frequency_scores(_MDL_ROOT+'/Distributions/pairwise_dd_residues.out')
		self.scoring_dicts['ad'] = load_AA_Frequency_scores(_MDL_ROOT+'/Distributions/pairwise_ad_residues.out')

		self.scoring_dicts['ag'] = load_AA_Frequency_scores(_MDL_ROOT+'/Distributions/pairwise_ag_residues.out')
		self.scoring_dicts['de'] = load_AA_Frequency_scores(_MDL_ROOT+'/Distributions/pairwise_de_residues.out')
		self.scoring_dicts['ae'] = load_AA_Frequency_scores(_MDL_ROOT+'/Distributions/pairwise_ae_residues.out')
		self.scoring_dicts['dg'] = load_AA_Frequency_scores(_MDL_ROOT+'/Distributions/pairwise_dg_residues.out')

		return

	def load_kde_models(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		#Load the model
		kde_AP_ada = load(_MDL_ROOT+'/estimators/KDE_atomicdensity_AP_ada.joblib')
		kde_AP_dad = load(_MDL_ROOT+'/estimators/KDE_atomicdensity_AP_dad.joblib')
		kde_P_ada = load(_MDL_ROOT+'/estimators/KDE_atomicdensity_P_ada.joblib')
		kde_P_dad = load(_MDL_ROOT+'/estimators/KDE_atomicdensity_P_dad.joblib')

		self.Hexscores = [kde_AP_ada,kde_AP_dad,kde_P_ada,kde_P_dad]

		return

	def get_padding_size(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		inf=open(_MDL_ROOT+'/DatasetScores/Padding_Length.out','r')
		lines = inf.readlines()
		inf.close()
		for l in lines:
			if ((l[0]=='#') or (len(l)==1)):
				continue
			self.RF_PaddingSize = int(l.strip())
		
		return self.RF_PaddingSize

	def get_KDE_Models(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.Hexscores

	def get_RF_Models(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.RF_Models

	def get_scoring_dicts(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.scoring_dicts

	def plot_KDE_Models(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return