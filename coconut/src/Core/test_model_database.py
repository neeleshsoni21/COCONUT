
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

from collections import OrderedDict
from numpy import array

def test_model():
	"""Summary
	
	Returns:
	    TYPE: Description
	"""
	from InputClass import InputCC
	from ModelClass import Models
	from ScoringClass import ScoreCC
	from read_CC_properties import generate_Hex_pairs_dimers
	from train_models import load_training_data, load_training_data_singlefile


	#DIMERS = load_training_data('./coconut/DATASET/CC_Dataset')
	DIMERS = load_training_data_singlefile('./coconut/DATASET/CC_Dataset', skip_filtering=False)
	
	DIMERS = generate_Hex_pairs_dimers(DIMERS)

	#Create Models class object to load all the pre-Modeled Data
	model_obj1 = Models()
	
	#Load the models
	model_obj1.load_RF_models()
	
	#Load Distributions
	model_obj1.load_AA_distributions()

	#Load the kde models of hex pairs in the model object
	model_obj1.load_kde_models()

	Hexpairs = model_obj1.get_KDE_Models()
	
	scoring_dicts1 = model_obj1.get_scoring_dicts()
	
	nof1 = model_obj1.get_RF_Models()[0].n_features_in_#n_features_
	
	#This is obtained by subtracting the paddings. See ScoringClass for details
	# 6 (prod(P_probs), prod(AP_Probs),2 P_Aa_scores, 2 AP_Aa_scores 
	# 8 (4 parallel and 4 antiparall) additional zeros were added. 
	max_hexpairs_length1 = int((nof1-6)/6)+8; 

	#Create the ScoreCC obj to score Given input sequences
	score_obj1 = ScoreCC()
	score_obj1.load_prior_probabilities()

	for key,dimer in DIMERS.items():
		dscore = score_obj1.score_dimer_RF(scoring_dicts1, Hexpairs, dimer, max_hexpairs_length1)
		score_obj1.ALLSCORES.append(dscore);
		dimer.append(sum(dscore))

	#Predict the input data using pre-tranined models
	score_obj1.Predict_CC_Stability(model_obj1)
	#score_obj1.Get_Predictions()

	Right_AP,Right_P,Wrong_AP,Wrong_P=0,0,0,0
	print("\n-------------------------------------------")
	for i,k in enumerate(DIMERS.keys()):

		dimer1 = DIMERS[k]
		dimer_prob1 =score_obj1.D_PROBS_mean[i]

		#If Dimer label is antiparallel
		if DIMERS[k][0]==-1:
			#if antiparallel probability is greater than parallel prob
			if dimer_prob1[0]>=dimer_prob1[1]:
				Right_AP+=1
			else:
				Wrong_AP+=1
		
		#If Dimer label is parallel
		elif DIMERS[k][0]==1:
			#if antiparallel probability is less than parallel prob
			if dimer_prob1[0]<=dimer_prob1[1]:
				Right_P+=1
			else:
				Wrong_P+=1

	print("Total AP Dimers:",(Right_AP+Wrong_AP))
	print("Total P Dimers:",(Right_P+Wrong_P))
	
	print("\nCorrect AP Predictions:", Right_AP)
	print("InCorrect AP Predictions:", Wrong_AP)
	print("Correct AP%",int(Right_AP*100/(Right_AP+Wrong_AP)))

	print("\nCorrect P Predictions:", Right_P)
	print("InCorrect P Predictions:", Wrong_P)
	print("Correct P%",int(Right_P*100/(Right_P+Wrong_P)))
	
	return

import sys
import copy
sys.path.insert(1, './coconut/src/Models/')
sys.path.insert(1, './coconut/src')
sys.path.insert(1, './coconut/src/Core')

if __name__ == '__main__':
	test_model()