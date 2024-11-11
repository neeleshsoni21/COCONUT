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


#-------------------------------#
#Code for re-building COCONUT Models
#-------------------------------#

import coconut as cc

cc_t_obj = cc.COCONUT()

'''
cc_t_obj.train_coconut_models()


cc_t_obj.test_coconut_models()


#-------------------------------#
#Code for running input example
#-------------------------------#

cc_obj = cc.COCONUT()
cc_obj.prediction_example()


#-------------------------------#
#Alternative way for running input example
#-------------------------------#

from pipeline1.load_pipeline import example1

data_dir = '2zta'
protein1, protein1_pcoils = '2zta','2zta_pcoils.out'
protein2, protein2_pcoils = '2zta','2zta_pcoils.out'
example1(data_dir, protein1, protein1_pcoils, protein2, protein2_pcoils)


data_dir = 'k5k14'
protein1, protein1_pcoils = 'k5','pcoils_k5.out'
protein2, protein2_pcoils = 'k14','pcoils_k14.out'
example1(data_dir, protein1, protein1_pcoils, protein2, protein2_pcoils)


#-------------------------------#
#Code for running on entire dataset
#-------------------------------#

from test_model_database import test_model

test_model()
'''



#-------------------------------#
#Code for generating PDB files for the coiled-coils
#using alignment files. Standard coiled-coils can also be constructed
#Some examples
#-------------------------------#

from cc_writer import *

ExampleRoot = os.path.join(cc_t_obj._COCO_ROOT,'./Externals/pdb_cc_writer/')
#Modeling of 2zta pdb
run_example1(ExampleRoot)

#Modeling of 2zta long
run_example2(ExampleRoot)

#Modeling for cc builder comparison
run_example3(ExampleRoot)

#Modeling of standard trimer
run_example4(ExampleRoot)

#This is for Nup60 modeling
cc_segments=[ [[27,47],[27,47],'g'],  [[91,104],[91,104],'g'],[[106,119],[106,119],'g'],  [[121,140],[121,140],'g'],[[142,162],[142,162],'g']]
run_example7(ExampleRoot, cc_segments)

#This is for Nup1 modeling
cc_segments=[ [[1,32],[1,32],'g'],  [[85,104],[85,104],'g'],[[106,123],[106,123],'g']]
run_example8(ExampleRoot, cc_segments)

#This is for generic or unified MLPs model
run_example9(ExampleRoot)

#This is for Mouse TPR Modeling
run_example12(ExampleRoot)

#This is for Nup153 modeling mouse
cc_segments=[ [[36,57],[36,57],'g'] ]
run_example13(ExampleRoot, cc_segments)

#This is for Keratin dimer modeling
run_example15(ExampleRoot)


