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

