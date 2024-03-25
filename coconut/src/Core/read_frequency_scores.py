
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
from numpy import log, sqrt

def load_AA_Frequency_scores(fname):
	"""Summary
	
	Args:
	    fname (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	#default -log(frequency) is 100.0
	Res_pair_frequency = defaultdict(lambda: 0.0)

	inf = open(fname,'r')
	lines = inf.readlines()
	inf.close()

	max_value=-1;

	for l in lines:
		toks = l.strip().split();
		Res_pair_frequency[tuple([toks[0],toks[1]])]=float(toks[2]);

	return Res_pair_frequency

def load_AA_Frequency_scores_logs(fname):
	"""Summary
	
	Args:
	    fname (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	#default -log(frequency) is 100.0
	Res_pair_frequency = defaultdict(lambda: 100.0)

	inf = open(fname,'r')
	lines = inf.readlines()
	inf.close()

	max_value=-1;

	for l in lines:
		toks = l.strip().split();
		if float(toks[2])!=0.0:
			Res_pair_frequency[tuple([toks[0],toks[1]])]=-1*log(float(toks[2]));

			#getting the max scores among the distribution
			if -1*log(float(toks[2])) > max_value:
				max_value = -1*log(float(toks[2]));
				
		else:
			#Since max_value for aa/dd/ad dicts are 6.7, 7.3, 6.8. Thus to prevent too much penalizing, log(0) was set as 100.
			Res_pair_frequency[tuple([toks[0],toks[1]])]=100;

	return Res_pair_frequency


