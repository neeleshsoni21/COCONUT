
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


import os

_EXP1_ROOT = os.path.abspath(os.path.dirname(__file__))

from InputClass import InputCC
from ModelClass import Models
from ScoringClass import ScoreCC
from SlidingWindowsClass import SlidingWindows
from ExtendCCClass import ExtendCoiledCoil
from CoiledCoilAlignClass import CoiledCoilAlign
from PairCoilClass import PairCoilParse

def example1(data_dir, protein1, protein1_pcoils, protein2, protein2_pcoils):
	"""Summary
	
	Args:
	    data_dir (TYPE): Description
	    protein1 (TYPE): Description
	    protein1_pcoils (TYPE): Description
	    protein2 (TYPE): Description
	    protein2_pcoils (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	merlen = 11;
	PROB_CUTOFF=0.25

	Sample_List = []
	
	Sample_List.append([data_dir, protein1, protein1_pcoils, protein2, protein2_pcoils])
	
	#------------------------------------------------------
	# Correct Paircoil file by doing minimal changes to the
	# missing/incorrect heptad repeats
	#------------------------------------------------------
	for ppath, prot1, prot1_paircoil, prot2, prot2_paircoil in Sample_List:

		prot1_prot2 = prot1+'_'+prot2
		prot1_paircoil_new = prot1_paircoil[:-4]+'_new.out'
		prot2_paircoil_new = prot2_paircoil[:-4]+'_new.out'

		
		PC_Prot1_File = os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_paircoil
		#PC1_Obj = PairCoilParse(PC_Prot1_File, PROB_CUTOFF)
		PC1_Obj = PairCoilParse(PC_Prot1_File, PROB_CUTOFF, probability_type = 'pcoil')
		p1_PC = PC1_Obj.get_pc_seq()
		PC1_Obj.correct_incorrectly_Indentified_heptad()

		PC_Prot1_File_New = os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_paircoil_new
		PC1_Obj.write_Updated_Paircoil_Probs(PC_Prot1_File_New)

		
		PC_Prot2_File = os.path.join(_EXP1_ROOT, ppath)+'/'+prot2_paircoil
		#PC2_Obj = PairCoilParse(PC_Prot2_File, PROB_CUTOFF)
		PC2_Obj = PairCoilParse(PC_Prot2_File, PROB_CUTOFF, probability_type = 'pcoil')
		p2_PC = PC2_Obj.get_pc_seq()
		PC2_Obj.correct_incorrectly_Indentified_heptad()

		PC_Prot2_File_New = os.path.join(_EXP1_ROOT, ppath)+'/'+prot2_paircoil_new
		PC2_Obj.write_Updated_Paircoil_Probs(PC_Prot2_File_New)
	
	
	#------------------------------------------------------
	# Generate Sliding Windows
	#------------------------------------------------------
	for ppath, prot1, prot1_paircoil, prot2, prot2_paircoil in Sample_List:

		print("Generating SLiding Windows for ",prot1_prot2)

		prot1_prot2 = prot1+'_'+prot2

		prot1_paircoil_new = prot1_paircoil[:-4]+'_new.out'
		prot2_paircoil_new = prot2_paircoil[:-4]+'_new.out'

		PC_Prot1_File = os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_paircoil_new
		PC_Prot2_File = os.path.join(_EXP1_ROOT, ppath)+'/'+prot2_paircoil_new

		phep = 'a'
		SlidingWindowFile_a = os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_prot2+'_CoiledCoil_a.out'
		phep = 'd'
		SlidingWindowFile_d = os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_prot2+'_CoiledCoil_d.out'

		slidew_obj = SlidingWindows(PC_Prot1_File, PC_Prot2_File, SlidingWindowFile_a, SlidingWindowFile_d, PROB_CUTOFF, merlen)
	
	

	#------------------------------------------------------
	# Generate Coiled-Coil segment scores from COCONUT Base scores
	#------------------------------------------------------	
	for ppath, prot1, prot1_paircoil, prot2, prot2_paircoil in Sample_List:

		prot1_prot2 = prot1+'_'+prot2

		for phep in ['a','d']:

			print("Generating Coiled-Coil Segment Scores for ",prot1_prot2," for Heptad ",phep)
		
			#input_file = _EXP1_ROOT+'/MLP2_MLP2_CoiledCoil.out'
			input_file = os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_prot2+'_CoiledCoil_'+phep+'.out'
			input_obj = InputCC(input_file)

			input_obj.load_CCData()
			
			#Display input info
			#input_obj.print_input_info()
			
			#----------------------------------------------#
			# Model load can go outside the loop
			#----------------------------------------------#
			#Create Models class object to load all the pre-Modeled Data
			model_obj = Models()

			#Load the models
			model_obj.load_RF_models()

			#Load Distributions
			model_obj.load_AA_distributions()

			#Load the kde models of hex pairs in the model object
			model_obj.load_kde_models()

			#Create the ScoreCC obj to score Given input sequences
			score_obj = ScoreCC()

			score_obj.load_prior_probabilities()

			#PREDICTION using RF Models
			score_obj.Generate_CC_Scores_RF(input_obj, model_obj)
			score_obj.Predict_CC_Stability(model_obj)
			score_obj.Write_Raw_Scores(input_obj, os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_prot2+'_CoiledCoil_Scores_'+phep+'.dat')
			##score_obj.Get_Predictions()

			#PREDICTION using Bayesian Models
			#score_obj.Generate_CC_Scores_Bayes(input_obj, model_obj)
			#score_obj.Write_Bayesian_Probabilities(input_obj, os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_prot2+'_CoiledCoil_Scores_'+phep+'.dat')

			#check_accuracy()
	
	#------------------------------------------------------
	# Extend Coiled-Coil segments that are overlapping
	#------------------------------------------------------
	for ppath, prot1, prot1_paircoil, prot2, prot2_paircoil in Sample_List:

		prot1_prot2 = prot1+'_'+prot2

		prot1_paircoil_new = prot1_paircoil[:-4]+'_new.out'
		prot2_paircoil_new = prot2_paircoil[:-4]+'_new.out'

		PC_Prot1_File = os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_paircoil_new
		PC_Prot2_File = os.path.join(_EXP1_ROOT, ppath)+'/'+prot2_paircoil_new

		ExtendCC = ExtendCoiledCoil(os.path.join(_EXP1_ROOT, ppath)+'/', prot1_prot2, PC_Prot1_File, PC_Prot2_File, PROB_CUTOFF)
	
	
	#------------------------------------------------------
	# Align two sequences using Coconut scores of extended segments, only parallel alignment is supported at this moment
	#------------------------------------------------------
	for ppath, prot1, prot1_paircoil, prot2, prot2_paircoil in Sample_List:

		prot1_prot2 = prot1+'_'+prot2

		prot1_paircoil_new = prot1_paircoil[:-4]+'_new.out'
		prot2_paircoil_new = prot2_paircoil[:-4]+'_new.out'

		PC_Prot1_File = os.path.join(_EXP1_ROOT, ppath)+'/'+prot1_paircoil_new
		PC_Prot2_File = os.path.join(_EXP1_ROOT, ppath)+'/'+prot2_paircoil_new

		
		#BaseAlignmentFile can be clustal Omega
		#k5_k14_basealignment=None;#k5_k14_basealignment=_EXP1_ROOT+'/k5k14/k5_k14_aligned.fasta'
		#AlignCC = CoiledCoilAlign(_EXP1_ROOT+'/'+ppath, prot1_prot2, PC_Prot1_File, PC_Prot2_File, PROB_CUTOFF, BaseAlignmentFile=k5_k14_basealignment)
		AlignCC = CoiledCoilAlign( os.path.join(_EXP1_ROOT,ppath)+'/', prot1_prot2, PC_Prot1_File, PC_Prot2_File, PROB_CUTOFF, BaseAlignmentFile=None)
	

	return

