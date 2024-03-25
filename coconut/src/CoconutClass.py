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


from sys import exit
import os
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import zeros, flip
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
import matplotlib.colors as colors
import numpy as np
import seaborn as sns


font = matplotlib.rcParams['font.sans-serif']='Arial'
font = matplotlib.rcParams['font.size']=14

# Class definition
class COCONUT:

	"""Summary
	"""
	
	def __init__(self):
		"""Summary
		
		Returns:
			TYPE: Description
		"""
		from coconut.version import __version__
		print("\nCOCONUT Version:", __version__,"\n")

		self._COCO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

		self._COCO_CCD = os.path.join(self._COCO_ROOT,'DATASET','CC_Dataset')

		self._COCO_DIST = os.path.join(self._COCO_ROOT,'src','Models','Distributions')

		self._COCO_CCS = os.path.join(self._COCO_ROOT,'src','Models','DatasetScores')

		self._COCO_NET = os.path.join(self._COCO_ROOT,'src','Models','trained_network')

		self._COCO_PLOT = os.path.join(self._COCO_ROOT,'src','Models','Plots')

		self._COCO_DODEC = os.path.join(self._COCO_ROOT,'src','Models','dodec_pairs')

		self._COCO_ESTM = os.path.join(self._COCO_ROOT,'src','Models','estimators')

		self._DIMERS = None

		return

	def train_coconut_models(self):
		"""Summary
		
		Returns:
			TYPE: Description
		"""
		from Core.train_models import load_training_data
		self._DIMERS =  load_training_data(self._COCO_CCD)

		from Core.read_CC_properties import generate_Hex_pairs_dimers
		self._DIMERS = generate_Hex_pairs_dimers(self._DIMERS)

		from Core.train_models import generate_dodec_pairs
		generate_dodec_pairs(self._COCO_DODEC, self._DIMERS)

		from Core.train_models import train_kde_models
		train_kde_models(self._COCO_DODEC, self._COCO_ESTM, self._COCO_PLOT, PlotData=True)

		from Core.read_CC_properties import generate_AA_pairs_frequency
		Res_stats = generate_AA_pairs_frequency(self._DIMERS)	

		self.plot_pairwise_distributions(Res_stats)
		
		from Core.train_models import dump_frequency_models
		dump_frequency_models(self._COCO_DIST, Res_stats)

		return

	def test_coconut_models(self):
		"""Summary
		
		Returns:
			TYPE: Description
		"""
		from ModelClass import Models
		#Create Models class object to load all the pre-Modeled Data
		model_obj = Models()

		#Load Amino acid Distributions in the model object
		model_obj.load_AA_distributions()

		#Load the kde models of hex pairs in the model object
		model_obj.load_kde_models()

		from ScoringClass import ScoreCC
		
		#Create the ScoreCC obj to score Given input sequences
		score_obj = ScoreCC()

		#Get the hexpair distributions and calculate the prior probabilities
		score_obj.calculate_prior_probabilities(self._DIMERS, self._COCO_CCS)

		print("Prior Probability Parallel Dimers:",round(score_obj.P_Prior_Prob,2))
		print("Prior Probability Anti-Parallel Dimers:",round(score_obj.AP_Prior_Prob,2))

		#Plot the dimer distribution
		self.plot_dimer_length_distribution(score_obj.get_Distributions())

		from ml_prediction import train_RF_models, train_RF_models_plot, train_bayesian_models
		
		print("Training RandomForest Models...")
		score_obj.Generate_CC_RF_Scores(self._DIMERS, model_obj, self._COCO_CCS)
		train_RF_models(self._COCO_CCS, self._COCO_NET, self._COCO_PLOT, self._DIMERS)
		
		#FOR PLOTTING USE BELOW
		train_RF_models_plot(self._COCO_CCS, self._COCO_NET, self._COCO_PLOT, self._DIMERS, model_obj)
		
		#print("Getting Bayesian-Inference Models...")
		#score_obj.Generate_CC_Bayesian_Scores(self._DIMERS, model_obj, self._COCO_CCS)
		#train_bayesian_models(self._COCO_CCS, self._COCO_NET, self._COCO_PLOT, self._DIMERS)

		return

	def predict_coiled_coil(self):
		"""Summary
		
		Returns:
			TYPE: Description
		"""
		#from Core.predict_cc import predict_coiled_coil
		#predict_coiled_coil()
		pass

		return

	def prediction_example(self):
		"""Summary
		
		Returns:
			TYPE: Description
		"""
		from pipeline1.load_pipeline import example1

		data_dir = '2zta'
		protein1, protein1_pcoils = '2zta','2zta_pcoils.out'
		protein2, protein2_pcoils = '2zta','2zta_pcoils.out'
		example1(data_dir, protein1, protein1_pcoils, protein2, protein2_pcoils)


		data_dir = 'k5k14'
		protein1, protein1_pcoils = 'k5','pcoils_k5.out'
		protein2, protein2_pcoils = 'k14','pcoils_k14.out'
		example1(data_dir, protein1, protein1_pcoils, protein2, protein2_pcoils)
		
		return

	def plot_pairwise_distributions(self, Res_stats):
		"""Summary
		
		Args:
			Res_stats (TYPE): Description
		
		Returns:
			TYPE: Description
		"""
		midnorm = MidpointNormalize(vmin=0., vcenter=0.04, vmax=1.0)
		colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 256))
		colors_land = plt.cm.terrain(np.linspace(0.25, 1, 256))
		all_colors = np.vstack((colors_undersea, colors_land))
		terrain_map = colors.LinearSegmentedColormap.from_list('terrain_map', all_colors)

		AA =['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
		heppairs = ['aa','dd','ad','ag','de','ae','dg','ee','gg','eg']

		DistMat = zeros((len(AA),len(AA)))
		
		for hidx, htype in enumerate(heppairs):

			for i,aa1 in enumerate(AA[::-1]):

				for j,aa2 in enumerate(AA):

					DistMat[i,j]=Res_stats[hidx][tuple([aa1,aa2])]

			fig = plt.figure()
			xlabs=AA;ylabs=AA[::-1];
			ax = sns.heatmap(DistMat, cmap="YlGnBu",
			 xticklabels=xlabs,yticklabels=ylabs,
			 square=True,cbar=True,
			 norm=SymLogNorm(linthresh=0.0001),
			 cbar_kws={'label': 'Probability'})

			ax.set_xlabel("Residue at Heptad "+htype[0])
			ax.set_ylabel("Residue at Heptad "+htype[1])
			
			plt.savefig(self._COCO_PLOT+'/Residue_pair_'+htype+'.pdf')
			#plt.show()
		
		
		return

	def plot_dimer_length_distribution(self, Distributions):
		"""Summary
		
		Args:
			Distributions (TYPE): Description
		
		Returns:
			TYPE: Description
		"""
		[ P_Dist, AP_Dist, P_HexPair_Dist, AP_HexPair_Dist ] = Distributions

		fig, ax1 = plt.subplots(1,2,figsize=(8,4))

		sns.histplot(AP_Dist,ax=ax1[0], binwidth=2,color="salmon")
		sns.histplot(P_Dist,ax=ax1[0],binwidth=2, color="teal")
		ax1[0].set_xlabel("Number of Residues in Dimers")
		ax1[0].set_xlim(0,40)
		ax1[0].legend(labels=["Anti-Parallel","Parallel"])
		ax1[0].set_title("Residue Distributions in CoiledCoils")

		sns.histplot(AP_HexPair_Dist,ax=ax1[1], binwidth=1, color="salmon")
		sns.histplot(P_HexPair_Dist,ax=ax1[1], binwidth=1,  color="teal")
		ax1[1].set_xlabel("Number of Hex-pairs in Dimers")
		ax1[1].set_xlim(0,20)
		ax1[1].legend(labels=["Anti-Parallel","Parallel"])
		ax1[1].set_title("Hexpair Distributions in CoiledCoils")
			
		plt.savefig(self._COCO_PLOT+'/Dimer-Hexpair_Length_Distribution.pdf')
		#plt.show()


		return

#TODO CHANGE THIS
class MidpointNormalize(colors.Normalize):

	"""Summary
	
	Attributes:
		vcenter (TYPE): Description
	"""
	
	def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
		"""Summary
		
		Args:
			vmin (None, optional): Description
			vmax (None, optional): Description
			vcenter (None, optional): Description
			clip (bool, optional): Description
		"""
		self.vcenter = vcenter
		super().__init__(vmin, vmax, clip)

	def __call__(self, value, clip=None):
		"""Summary
		
		Args:
			value (TYPE): Description
			clip (None, optional): Description
		
		Returns:
			TYPE: Description
		"""
		x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1.]
		return np.ma.masked_array(np.interp(value, x, y,
											left=-np.inf, right=np.inf))

	def inverse(self, value):
		"""Summary
		
		Args:
			value (TYPE): Description
		
		Returns:
			TYPE: Description
		"""
		y, x = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
		return np.interp(value, x, y, left=-np.inf, right=np.inf)




