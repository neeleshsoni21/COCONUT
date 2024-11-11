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

def read_fasta(fname):

	inf = open(fname,'r')
	lines = inf.readlines()
	inf.close()

	Qryseq=""
	SeqId=""
	for l in lines:
		if len(l)==1 or l[0]=='':
			continue
		if l[0]=='>':
			SeqId=l[1:].strip()
			continue
		Qryseq += l.strip()

	return Qryseq, SeqId