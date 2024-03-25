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


from copy import copy

def max_normalize_dicts(Res_stats):
	"""Summary
	
	Args:
	    Res_stats (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	for i in range(0,len(Res_stats)):
		Res_pairs = Res_stats[i]

		#Max normalize all dictionaries
		maxval=0
		for key,value in Res_pairs.items():
			if maxval<value:
				maxval=value;
		for key,value in Res_pairs.items():
			Res_pairs[key]=Res_pairs[key]/float(maxval);

		Res_stats[i] = Res_pairs;

	return Res_stats

def total_normalize_dicts(Res_stats):
	"""Summary
	
	Args:
	    Res_stats (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	for i in range(0,len(Res_stats)):
		Res_pairs = Res_stats[i]

		#total normalize all dictionaries
		totval=0
		for key,value in Res_pairs.items():
			totval+=value;

		for key,value in Res_pairs.items():
			Res_pairs[key]=Res_pairs[key]/float(totval);

		Res_stats[i] = Res_pairs;

	return Res_stats

def write_unique_hexpairs(Realdimers):
	"""Summary
	
	Args:
	    Realdimers (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	hp=[];

	for key,value in Realdimers.items():
		for hexp in value[3]:
			hp.append(tuple([hexp[1],hexp[2]]))
	Uhp=set(hp)
	wlines=""
	for hexp in Uhp:
		wlines+=hexp[0]+" "+hexp[1]+"\n";

	outf=open("Unique_hex_pairs.out",'w')
	outf.writelines(wlines);
	outf.close()

	return

def shifting_align_dimer(value):
	"""Summary
	
	Args:
	    value (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	seq1,hep1,seq2,hep2 = value[1][0],value[1][1],value[2][0],value[2][1]
	
	#If anti-parallel
	if value[0]==-1:

		#if seq1 starts with a then remove 'abc' from seq 1, and last 3 from seq2
		if (hep1[0]==hep1[-1]):

			#Case 1
			#if hep1 starts with 'a', remove first 3 letters from hep1, and last 3 from 'd'
			if hep1[0]=='a':
				hep1 = hep1[3:];seq1=seq1[3:]
				hep2 = hep2[:-3];seq2=seq2[:-3]

			#Case 2
			#if hep1 starts with 'd', means hep2 starts with 'a', remove first 3 letters from hep2
			elif hep1[0]=='d':
				hep1 = hep1[:-3];seq1=seq1[:-3]
				hep2 = hep2[3:];seq2=seq2[3:]
				
		else:
			#Case 3, Case 4
			#No change required for converting anti-parallel to parallel
			pass
	
	#If parallel
	else:

		#if seq1 starts with a then remove 'abc' from seq 1, and last 3 from seq2
		if (hep1[0]==hep1[-1]):

			#Case 5
			#if hep1 starts with 'a', remove first 3 letters from hep1, and first 3 from hep2
			if hep1[0]=='a':
				hep1 = hep1[3:];seq1=seq1[3:]
				hep2 = hep2[3:];seq2=seq2[3:]

			#Case 6
			#if hep1 starts with 'd', means hep2 starts with 'a', remove first 3 letters from hep2
			elif hep1[0]=='d':
				hep1 = hep1[:-3];seq1=seq1[:-3]
				hep2 = hep2[:-3];seq2=seq2[:-3]
				
		else:
			#Case 7, Case 8
			#No change required for converting parallel to anti-parallel
			pass
	valueR = [value[0],[seq1,hep1],[seq2,hep2]];

	return valueR

def generate_Hex_pairs_dimers(dimers):
	"""Summary
	
	Args:
	    dimers (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	#dimers is a dict
	#key is the dimer id obtained from the database
	#value is a list of 2 elements. 1st element indicates whether dimer is paralle or not.
	#2nd element provides seq and corresponding heptads.
	#3rd element adds (in this function) corresponding ada-dad pairs in the dimers dict in a list data type
	funcdict={1:generate_Hex_pairs_parallel_dimers,-1:generate_Hex_pairs_antiparallel_dimers}

	#for key,value in list(dimers.items())[3:4]:
	for key,value in dimers.items():
		
		seq1,hep1,seq2,hep2 = value[1][0],value[1][1],value[2][0],value[2][1]

		#Add parallel dodec pairs
		#Get dodec pairs considering parallel dimers
		GDODECPAIRS1 , AA_pairs1= funcdict[1](key,value);
		#add the diemrs ada-ada and dad-dad sequence pairs in the dict
		dimers[key].append(GDODECPAIRS1);
		dimers[key].append(AA_pairs1);

		#Add parallel dodec pairs
		#Get dodec pairs considering anti-parallel dimers
		GDODECPAIRS2, AA_pairs2 = funcdict[-1](key,value);
		#add the diemrs ada-dad sequence pairs in the dict
		dimers[key].append(GDODECPAIRS2);
		dimers[key].append(AA_pairs2);

	return dimers

def generate_Hex_pairs_parallel_dimers(key,v):
	"""Summary
	
	Args:
	    key (TYPE): Description
	    v (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	GDODECPAIRS=[]; AA_pairs=[];

	Seq1,Seq2,Hep1,Hep2=v[1][0],v[2][0],v[1][1],v[2][1]
		
	DodecPairIdx=[]
	DodecPairSeq=[]
	DodecPairHep=[]

	HexPairHep=[]
	HexPairIdx=[]

	for i in range(0,len(Seq1)):
		
		if Hep1[i] in ['a','d']:
			#No need to check if the two heptad repeat are same or not. already checked
			HexPairHep.append(Hep1[i])
			HexPairIdx.append(i)
			#print(HexPairHep)
			if len(HexPairHep)==3:
				if ''.join(HexPairHep)=="dad":
					
					#For parallel dimer indexes will be same for both helices
					#For a,d,e,g positons
					#helix1seq=[HexPairIdx[0],HexPairIdx[0]+1,HexPairIdx[1]-1,HexPairIdx[1],HexPairIdx[2],HexPairIdx[2]+1]
					#helix2seq=[HexPairIdx[0],HexPairIdx[0]+1,HexPairIdx[1]-1,HexPairIdx[1],HexPairIdx[2],HexPairIdx[2]+1]

					#For a and d only
					helix1seq=[HexPairIdx[0],HexPairIdx[1],HexPairIdx[2]]
					helix2seq=[HexPairIdx[0],HexPairIdx[1],HexPairIdx[2]]

					#For e,g positons
					helix3seq=[HexPairIdx[0]+1,HexPairIdx[1]-1,HexPairIdx[2]+1]
					helix4seq=[HexPairIdx[0]+1,HexPairIdx[1]-1,HexPairIdx[2]+1]

					DodecPairIdx.append([helix1seq,helix2seq,helix3seq,helix4seq])
						
				elif ''.join(HexPairHep)=="ada":
					
					#For parallel dimer indexes will be same for both helices
					#For a,d,e,g positons
					#helix1seq=[HexPairIdx[0]-1,HexPairIdx[0],HexPairIdx[1],HexPairIdx[1]+1,HexPairIdx[2]-1,HexPairIdx[2]]
					#helix2seq=[HexPairIdx[0]-1,HexPairIdx[0],HexPairIdx[1],HexPairIdx[1]+1,HexPairIdx[2]-1,HexPairIdx[2]]
						
					#For a and d only
					helix1seq=[HexPairIdx[0],HexPairIdx[1],HexPairIdx[2]]
					helix2seq=[HexPairIdx[0],HexPairIdx[1],HexPairIdx[2]]

					#For e,g positons
					helix3seq=[HexPairIdx[0]-1,HexPairIdx[1]+1,HexPairIdx[2]-1]
					helix4seq=[HexPairIdx[0]-1,HexPairIdx[1]+1,HexPairIdx[2]-1]
											

					DodecPairIdx.append([helix1seq,helix2seq,helix3seq,helix4seq])
						
				HexPairHep.pop(0)
				HexPairIdx.pop(0)
			else:
				continue

	for dp in DodecPairIdx:
		temp1=[]
		htemp1=[]
		
		#Code for helix 1 sequence (a/d positions)
		#This code is useless since ada or dad always exists. '-' can never happen. IMP for e and g positions
		for i in dp[0]:
			if (i<0) | (i>=len(Seq1)):
				temp1.append('-')
				htemp1.append('-')
			else:
				temp1.append(Seq1[i])
				htemp1.append(Hep1[i])

		temp2=[]
		htemp2=[]

		#Code for helix 2 sequence (a/d positions)
		#This code is useless since ada or dad always exists. '-' can never happen. IMP for e and g positions
		for i in dp[1]:
			if (i<0) | (i>=len(Seq2)):
				temp2.append('-')
				htemp2.append('-')
			else:
				temp2.append(Seq2[i])
				htemp2.append(Hep2[i])

		temp3=[]
		htemp3=[]
		#Code for helix 1 sequence (e/g positions)
		#IMP for e and g positions
		for i in dp[2]:
			if (i<0) | (i>=len(Seq1)):
				temp3.append('-')
				htemp3.append('-')
			else:
				temp3.append(Seq1[i])
				htemp3.append(Hep1[i])

		temp4=[]
		htemp4=[]
		#Code for helix 2 sequence (e/g positions)
		#IMP for e and g positions
		for i in dp[3]:
			if (i<0) | (i>=len(Seq2)):
				temp4.append('-')
				htemp4.append('-')
			else:
				temp4.append(Seq2[i])
				htemp4.append(Hep2[i])

		DodecPairSeq.append([temp1,temp2,temp3,temp4])
		DodecPairHep.append([htemp1,htemp2,htemp3,htemp4])

	for i in range(0,len(DodecPairIdx)):
		
		templ=[]; templ.append(str(key)+'_'+str(i));
		templ.append(''.join(DodecPairSeq[i][0])); templ.append(''.join(DodecPairSeq[i][1]))
		templ.append(''.join(DodecPairHep[i][0])); templ.append(''.join(DodecPairHep[i][1]))
		
		templ.append(''.join(DodecPairSeq[i][2])); templ.append(''.join(DodecPairSeq[i][3]))
		templ.append(''.join(DodecPairHep[i][2])); templ.append(''.join(DodecPairHep[i][3]))

		GDODECPAIRS.append(templ);

	for i in range(0,len(Seq1)):
		
		if ((Hep1[i]=='a') and (Hep2[i]=='a')):
			AA_pairs.append([Seq1[i],Seq2[i],Hep1[i],Hep2[i]])

		if ((Hep1[i]=='d') and (Hep2[i]=='d')):
			AA_pairs.append([Seq1[i],Seq2[i],Hep1[i],Hep2[i]])

	return GDODECPAIRS, AA_pairs

def generate_Hex_pairs_antiparallel_dimers(key,v):
	"""Summary
	
	Args:
	    key (TYPE): Description
	    v (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	GDODECPAIRS=[]; AA_pairs=[]
	
	Seq1,Seq2,Hep1,Hep2=v[1][0],v[2][0],v[1][1],v[2][1]

	#get indices to trim the sequences
	s1flag=0
	for i in range(0,len(Seq1)):
		if Hep1[i] in ['a','d']:
			s1flag=i;
			break
	e1flag=0
	for i in range(0,len(Seq1)):
		j=len(Seq1)-i-1;
		if Hep1[j] in ['a','d']:
			e1flag=j+1;
			break
	s2flag=0
	for i in range(0,len(Seq2)):
		if Hep2[i] in ['a','d']:
			s2flag=i;
			break
	e2flag=0
	#Reverse indexes for antiparallel dimers
	for i in range(0,len(Seq2)):
		j=len(Seq2)-i-1;
		if Hep2[j] in ['a','d']:
			e2flag=j+1;
			break;
		
	#Readjust the both helix Seqs so that they start with a/d.
	Seq1,Seq2,Hep1,Hep2 = Seq1[s1flag:e1flag],Seq2[s2flag:e2flag],Hep1[s1flag:e1flag],Hep2[s2flag:e2flag]
		
	DodecPairIdx=[]
	DodecPairSeq=[]
	DodecPairHep=[]

	HexPairHep=[]
	HexPairHep2=[]
	HexPairIdx=[]
	HexPairIdx2=[]
	

	for i in range(0,len(Seq1)):
		#Reverse indexes for antiparallel dimers
		j=len(Seq2)-i-1
		#print (Seq1[i],Seq2[i],Hep1[i],Hep2[i])
		if Hep1[i] in ['a','d']:
			#No need to check if the two heptad repeat are same or not. already checked
			HexPairHep.append(Hep1[i])
			HexPairIdx.append(i)
			HexPairIdx2.append(j)
			HexPairHep2.append(Hep2[j])
			#print(HexPairHep)
			if len(HexPairHep)==3:
				#Helix 1 corresponds to degade, partner helix2 corresponds to agedag
				if ''.join(HexPairHep)=="dad":
					
					#For parallel dimer indexes will be same for both helices
					#For a,d,e,g positons
					#helix1seq=[HexPairIdx[0],HexPairIdx[0]+1,HexPairIdx[1]-1,HexPairIdx[1],HexPairIdx[2],HexPairIdx[2]+1]
					#helix2seq=[HexPairIdx2[0],HexPairIdx2[0]-1,HexPairIdx2[1]+1,HexPairIdx2[1],HexPairIdx2[2],HexPairIdx2[2]-1]
						
					#For a and d only
					helix1seq=[HexPairIdx[0],HexPairIdx[1],HexPairIdx[2]]
					helix2seq=[HexPairIdx2[0],HexPairIdx2[1],HexPairIdx2[2]]

					#For e,g positons
					helix3seq=[HexPairIdx[0]+1,HexPairIdx[1]-1,HexPairIdx[2]+1]
					helix4seq=[HexPairIdx2[0]-1,HexPairIdx2[1]+1,HexPairIdx2[2]-1]
					
					DodecPairIdx.append([helix1seq,helix2seq,helix3seq,helix4seq])
					
				#Helix 1 corresponds to gadega, partner helix2 corresponds to edaged
				elif ''.join(HexPairHep)=="ada":
					
					#For parallel dimer indexes will be same for both helices
					#For a,d,e,g positons
					#helix1seq=[HexPairIdx[0]-1,HexPairIdx[0],HexPairIdx[1],HexPairIdx[1]+1,HexPairIdx[2]-1,HexPairIdx[2]]
					#helix2seq=[HexPairIdx2[0]+1,HexPairIdx2[0],HexPairIdx2[1],HexPairIdx2[1]-1,HexPairIdx2[2]+1,HexPairIdx2[2]]

					#For a and d only
					helix1seq=[HexPairIdx[0],HexPairIdx[1],HexPairIdx[2]]
					helix2seq=[HexPairIdx2[0],HexPairIdx2[1],HexPairIdx2[2]]					

					#For e,g positons
					helix3seq=[HexPairIdx[0]-1,HexPairIdx[1]+1,HexPairIdx[2]-1]
					helix4seq=[HexPairIdx2[0]+1,HexPairIdx2[1]-1,HexPairIdx2[2]+1]

					DodecPairIdx.append([helix1seq,helix2seq,helix3seq,helix4seq])
						
				HexPairHep.pop(0)
				HexPairIdx.pop(0)
				HexPairIdx2.pop(0)
				HexPairHep2.pop(0)
			else:
				continue

	for dp in DodecPairIdx:
		#Code for helix 1 sequence (a/d positions)
		#Useless for a/d positions
		temp1=[]
		htemp1=[]
			
		for k in dp[0]:
			if (k<0) | (k>=len(Seq1)):
				temp1.append('-')
				htemp1.append('-')
			else:
				temp1.append(Seq1[k])
				htemp1.append(Hep1[k])

		#Code for helix 2 sequence (a/d positions)
		#Useless for a/d positions
		temp2=[]
		htemp2=[]
		for k in dp[1]:
			if (k<0) | (k>=len(Seq2)):
				temp2.append('-')
				htemp2.append('-')
			else:
				temp2.append(Seq2[k])
				htemp2.append(Hep2[k])

		temp3=[]
		htemp3=[]
		#Code for helix 1 sequence (e/g positions)
		#IMP for e and g positions
		for i in dp[2]:
			if (i<0) | (i>=len(Seq1)):
				temp3.append('-')
				htemp3.append('-')
			else:
				temp3.append(Seq1[i])
				htemp3.append(Hep1[i])

		temp4=[]
		htemp4=[]
		#Code for helix 2 sequence (e/g positions)
		#IMP for e and g positions
		for i in dp[3]:
			if (i<0) | (i>=len(Seq2)):
				temp4.append('-')
				htemp4.append('-')
			else:
				temp4.append(Seq2[i])
				htemp4.append(Hep2[i])

		DodecPairSeq.append([temp1,temp2,temp3,temp4])
		DodecPairHep.append([htemp1,htemp2,htemp3,htemp4])

	for i in range(0,len(DodecPairIdx)):
		
		templ=[]; templ.append(str(key)+'_'+str(i));
		templ.append(''.join(DodecPairSeq[i][0])); templ.append(''.join(DodecPairSeq[i][1]))
		templ.append(''.join(DodecPairHep[i][0])); templ.append(''.join(DodecPairHep[i][1]))
		
		templ.append(''.join(DodecPairSeq[i][2])); templ.append(''.join(DodecPairSeq[i][3]))
		templ.append(''.join(DodecPairHep[i][2])); templ.append(''.join(DodecPairHep[i][3]))

		GDODECPAIRS.append(templ);

	for i in range(0,len(Seq1)):
		#Reverse indexes for antiparallel dimers
		j=len(Seq2)-i-1

		if ((Hep1[i]=='a') and (Hep2[j]=='d')):
			AA_pairs.append([Seq1[i],Seq2[j],Hep1[i],Hep2[j]])

		if ((Hep1[i]=='d') and (Hep2[j]=='a')):
			#AA_pairs.append([Seq2[j],Seq1[i],Hep2[j],Hep1[i]])
			AA_pairs.append([Seq1[i],Seq2[j],Hep1[i],Hep2[j]])

	return GDODECPAIRS, AA_pairs

def generate_AA_pairs_frequency(dimers):
	"""Summary
	
	Args:
	    dimers (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	#dimers is a dict
	#key is the dimer id obtained from the database
	#value is a list of 2 elements. 1st element indicates whether dimer is paralle or not.
	#2nd element provides seq and corresponding heptads.
	#3rd element adds (in this function) corresponding ada-dad pairs in the dimers dict in a list data type
	funcdict={1:generate_AA_pairs_parallel_dimers,-1:generate_AA_pairs_antiparallel_dimers}

	#Since the two succeesive dodec pairs are overlapped. The frequency counts will be three times except for the terminal ones
	Res_pairs_aa,Res_pairs_dd,Res_pairs_ad,Res_pairs_ag,Res_pairs_de,Res_pairs_ae,Res_pairs_dg,Res_pairs_ee,Res_pairs_gg,Res_pairs_eg={},{},{},{},{},{},{},{},{},{}
	aa_list=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];

	Defaultval=0.00000000000000000000000000001
	for i in range(0,len(aa_list)):
		for j in range(0,len(aa_list)):
			key1=aa_list[i];
			key2=aa_list[j];
			Res_pairs_aa[tuple([key1,key2])]=Defaultval
			Res_pairs_dd[tuple([key1,key2])]=Defaultval
			Res_pairs_ad[tuple([key1,key2])]=Defaultval

			Res_pairs_ag[tuple([key1,key2])]=Defaultval
			Res_pairs_de[tuple([key1,key2])]=Defaultval
			Res_pairs_ae[tuple([key1,key2])]=Defaultval
			Res_pairs_dg[tuple([key1,key2])]=Defaultval

			Res_pairs_ee[tuple([key1,key2])]=Defaultval
			Res_pairs_gg[tuple([key1,key2])]=Defaultval
			Res_pairs_eg[tuple([key1,key2])]=Defaultval

	Res_stats = [Res_pairs_aa,Res_pairs_dd,Res_pairs_ad,Res_pairs_ag,Res_pairs_de,Res_pairs_ae,Res_pairs_dg,Res_pairs_ee,Res_pairs_gg,Res_pairs_eg]

	for key,value in dimers.items():
		
		seq1,hep1,seq2,hep2 = value[1][0],value[1][1],value[2][0],value[2][1]
		
		#Both the heptad repeats shouldn't be discontinuos
		if check_heptad_repeat(hep1)*check_heptad_repeat(hep2)==0:
			continue

		#In dimer[key], initially parallel dodec pairs and then anti-parallel dodec pairs
		if value[0]==1:
			#Add parallel dodec pairs
			#Get dodec pairs considering parallel dimers
			Res_stats = funcdict[1](key,value,Res_stats);
			
		else:

			#Get dodec pairs considering anti-parallel dimers
			Res_stats  = funcdict[-1](key,value,Res_stats);


	#Sort according to key  then value
	#Res_pairs_aa_sorted=sorted(Res_pairs_aa.items(),key=itemgetter(0,0,1), reverse=False)
	#Res_pairs_dd_sorted=sorted(Res_pairs_dd.items(),key=itemgetter(0,0,1), reverse=False)
	#Res_pairs_ad_sorted=sorted(Res_pairs_ad.items(),key=itemgetter(0,0,1), reverse=False)

	#First normalize according to the frequency
	Res_stats = total_normalize_dicts(Res_stats)

	#Then scale with the max value for scaling
	#Res_stats = max_normalize_dicts(Res_stats)

	#Make the Res_pairs_aa and Res_pairs_dd symmetric, since aa/ and dd is symmetric.
	#Make ee and gg and eg pairs also symetric
	#Symmetry is done by adding non-diagonal elements
	Res_pairs_aa = Res_stats[0]
	Res_pairs_aa1 = copy(Res_stats[0])
	for key,value in Res_pairs_aa.items():
		key1=key[0]; key2=key[1];
		k1 = tuple([key1,key2]); k2 = tuple([key2,key1]);
		if key1==key2:
			continue
		temp = Res_pairs_aa[k1] + Res_pairs_aa[k2]
		Res_pairs_aa1[k1] = temp;

	Res_stats[0] = Res_pairs_aa1
	

	Res_pairs_dd = Res_stats[1]
	Res_pairs_dd1 = copy(Res_stats[1])
	for key,value in Res_pairs_dd.items():
		key1=key[0]; key2=key[1];		
		if key1==key2:
			continue
		k1 = tuple([key1,key2]); k2 = tuple([key2,key1]);
		temp = Res_pairs_dd[k1] + Res_pairs_dd[k2]
		Res_pairs_dd1[k1] = temp;

	Res_stats[1] = Res_pairs_dd1
	
	
	Res_pairs_ee = Res_stats[7]
	Res_pairs_ee1 = copy(Res_stats[7])
	for key,value in Res_pairs_ee.items():
		key1=key[0]; key2=key[1];
		if key1==key2:
			continue
		k1 = tuple([key1,key2]); k2 = tuple([key2,key1]);
		temp = Res_pairs_ee[k1] + Res_pairs_ee[k2]
		Res_pairs_ee1[k1] = temp; 
		
	Res_stats[7] = Res_pairs_ee1
	
	
	Res_pairs_gg = Res_stats[8]
	Res_pairs_gg1 = copy(Res_stats[8])
	for key,value in Res_pairs_gg.items():
		key1=key[0]; key2=key[1];
		if key1==key2:
			continue
		k1 = tuple([key1,key2]); k2 = tuple([key2,key1]);
		temp = Res_pairs_gg[k1] + Res_pairs_gg[k2]
		Res_pairs_gg1[k1] = temp;
	
	Res_stats[8] = Res_pairs_gg1


	Res_pairs_eg = Res_stats[9]
	Res_pairs_eg1 = copy(Res_stats[9])
	for key,value in Res_pairs_eg.items():
		key1=key[0]; key2=key[1];
		if key1==key2:
			continue
		k1 = tuple([key1,key2]); k2 = tuple([key2,key1]);
		temp = Res_pairs_eg[k1] + Res_pairs_eg[k2]
		Res_pairs_eg1[k1] = temp;
	
	Res_stats[9] = Res_pairs_eg1
	
	return Res_stats

def generate_AA_pairs_parallel_dimers(key,v,Res_stats):
	"""Summary
	
	Args:
	    key (TYPE): Description
	    v (TYPE): Description
	    Res_stats (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	Res_pairs_aa,Res_pairs_dd,Res_pairs_ad,Res_pairs_ag,Res_pairs_de,Res_pairs_ae,Res_pairs_dg,Res_pairs_ee,Res_pairs_gg,Res_pairs_eg = Res_stats[0],Res_stats[1],Res_stats[2], Res_stats[3], Res_stats[4], Res_stats[5], Res_stats[6], Res_stats[7], Res_stats[8], Res_stats[9]

	Seq1,Seq2,Hep1,Hep2=v[1][0],v[2][0],v[1][1],v[2][1]

	#Both the heptad repeats shouldn't be discontinuos
	if check_heptad_repeat(Hep1)*check_heptad_repeat(Hep2)==0:

		return Res_pairs_aa,Res_pairs_dd,Res_pairs_ad,Res_pairs_ag,Res_pairs_de,Res_pairs_ae,Res_pairs_dg,Res_pairs_ee,Res_pairs_gg,Res_pairs_eg
	
	for i in range(0,len(Seq1)):
		
		if ((Seq1[i]=='X') or (Seq2[i]=='X')):
			continue

		if ((Hep1[i]=='a') and (Hep2[i]=='a')):
			Res_pairs_aa[tuple([Seq1[i],Seq2[i]])]+=1

			if (i-1) >= 0:
				#Seq1[i] ==> a, Seq2[i-1] ==> g
				Res_pairs_ag[tuple([Seq1[i],Seq2[i-1]])]+=1
				Res_pairs_ag[tuple([Seq2[i],Seq1[i-1]])]+=1 ; #print(Hep1[i],Hep2[i-1],Hep2[i],Hep1[i-1])

			if (i-3) >= 0:
				#Seq1[i-1] ==> g, Seq2[i-3] ==> e
				Res_pairs_eg[tuple([Seq2[i-3],Seq1[i-1]])]+=1
				#Res_pairs_eg[tuple([Seq1[i-3],Seq2[i-1]])]+=1

			if ((i-1) >= 0) and ((i+4) < len(Seq1)):
				Res_pairs_eg[tuple([Seq2[i+4],Seq1[i-1]])]+=1 ; #print(Hep1);print(Hep2);print(Hep1[i+4],Hep1[i-1],Hep1[i+4],Hep2[i-1])
				#Res_pairs_eg[tuple([Seq1[i+4],Seq2[i-1]])]+=1 ;

		if ((Hep1[i]=='d') and (Hep2[i]=='d')):

			Res_pairs_dd[tuple([Seq1[i],Seq2[i]])]+=1

			if (i+1) < len(Seq1):
				#Seq1[i] ==> d, Seq2[i+1] ==> e
				Res_pairs_de[tuple([Seq1[i],Seq2[i+1]])]+=1
				Res_pairs_de[tuple([Seq2[i],Seq1[i+1]])]+=1; #print(Hep1[i],Hep2[i+1],Hep2[i],Hep1[i+1])

			if (i+3) < len(Seq1):
				#Seq1[i+1] ==> e, Seq2[i+3] ==> g
				Res_pairs_eg[tuple([Seq1[i+1], Seq2[i+3]])]+=1
				#Res_pairs_eg[tuple([Seq2[i+1], Seq1[i+3]])]+=1; #print(Hep2[i+3],Hep1[i+1],Hep1[i+3],Hep2[i+1])

			if ((i-4) >= 0) and ((i+1) < len(Seq1)):
				Res_pairs_eg[tuple([Seq1[i+1], Seq2[i-4]])]+=1 ; #print(Hep2[i+3],Hep1[i+1],Hep1[i+3],Hep2[i+1])
				#Res_pairs_eg[tuple([Seq2[i+1], Seq1[i-4]])]+=1 ;

	Res_stats = [Res_pairs_aa,Res_pairs_dd,Res_pairs_ad,Res_pairs_ag,Res_pairs_de,Res_pairs_ae,Res_pairs_dg,Res_pairs_ee,Res_pairs_gg,Res_pairs_eg]

	return Res_stats

def check_heptad_repeat(hep):
	"""Summary
	
	Args:
	    hep (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	flag=True
	ph=hep[0];
	for h in hep[1:]:
		ph_n = ((ord(ph)-97)+1)%7
		h_n = ord(h)-97
		flag = flag * (ph_n==h_n)
		#print(ph,h, ph_n, h_n, ph_n==h_n,flag )
		ph=h
	return flag

def generate_AA_pairs_antiparallel_dimers(key,v,Res_stats):
	"""Summary
	
	Args:
	    key (TYPE): Description
	    v (TYPE): Description
	    Res_stats (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	Res_pairs_aa,Res_pairs_dd,Res_pairs_ad,Res_pairs_ag,Res_pairs_de,Res_pairs_ae,Res_pairs_dg,Res_pairs_ee,Res_pairs_gg,Res_pairs_eg = Res_stats[0],Res_stats[1],Res_stats[2], Res_stats[3], Res_stats[4], Res_stats[5], Res_stats[6], Res_stats[7], Res_stats[8], Res_stats[9]
	
	Seq1,Seq2,Hep1,Hep2=v[1][0],v[2][0],v[1][1],v[2][1]

	#Both the heptad repeats shouldn't be discontinuos
	if check_heptad_repeat(Hep1)*check_heptad_repeat(Hep2)==0:
		return Res_pairs_aa,Res_pairs_dd,Res_pairs_ad,Res_pairs_ag,Res_pairs_de,Res_pairs_ae,Res_pairs_dg,Res_pairs_ee,Res_pairs_gg,Res_pairs_eg

	#get indices to trim the sequences
	s1flag=0
	for i in range(0,len(Seq1)):
		if Hep1[i] in ['a','d']:
			s1flag=i;
			break
	e1flag=0
	for i in range(0,len(Seq1)):
		j=len(Seq1)-i-1;
		if Hep1[j] in ['a','d']:
			e1flag=j+1;
			break
	s2flag=0
	for i in range(0,len(Seq2)):
		if Hep2[i] in ['a','d']:
			s2flag=i;
			break
	e2flag=0
	#Reverse indexes for antiparallel dimers
	for i in range(0,len(Seq2)):
		j=len(Seq2)-i-1;
		if Hep2[j] in ['a','d']:
			e2flag=j+1;
			break;
		
	#Readjust the both helix Seqs so that they start with a/d.
	Seq1,Seq2,Hep1,Hep2 = Seq1[s1flag:e1flag],Seq2[s2flag:e2flag],Hep1[s1flag:e1flag],Hep2[s2flag:e2flag]

	for i in range(0,len(Seq1)):
		#Reverse indexes for antiparallel dimers
		j=len(Seq2)-i-1

		if ((Seq1[i]=='X') or (Seq2[j]=='X')):
			continue

		if ((Hep1[i]=='a') and (Hep2[j]=='d')):

			Res_pairs_ad[tuple([Seq1[i],Seq2[j]])]+=1
			
			if (i-1) >= 0:
				#Seq2[j] ==> d, Seq1[i-1] ==> g
				Res_pairs_dg[tuple([Seq2[j],Seq1[i-1]])]+=1 ; #print(Hep2[j],Hep1[i-1])

			if (j+1) < len(Seq1):
				#Seq1[i] ==> a, Seq2[j+1] ==> e
				Res_pairs_ae[tuple([Seq1[i],Seq2[j+1]])]+=1; #print(Hep1[i],Hep2[j+1])

			if ((i-1) >= 0) and ((j+3) < len(Seq1)):
				#Seq1[i-1] ==> g, Seq2[j+3] ==> g
				Res_pairs_gg[tuple([Seq1[i-1],Seq2[j+3]])]+=1; #print(Hep1[i-1],Hep2[j+3])

			if ((i-1) >= 0) and ((j-4) >= 0):
				#Seq1[i-1] ==> g, Seq2[j-4] ==> g
				Res_pairs_gg[tuple([Seq1[i-1],Seq2[j-4]])]+=1; #print(Hep1[i-1],Hep2[j-4])

			"""
			#ee pairs are added by below d-a pairings
			if ((i-3) >= 0) and ((j+1) < len(Seq1)):
				#Seq1[i-3] ==> e, Seq2[j+1] ==> e
				#Res_pairs_ee[tuple([Seq2[j+1],Seq1[i-3]])]+=1; #print(Hep2[j+1],Hep1[i-3])

			if ((i+4) < len(Seq1)) and ((j+1) < len(Seq1)):
				#Seq1[i-1] ==> e, Seq2[j-4] ==> e
				#Res_pairs_ee[tuple([Seq2[j+1],Seq1[i+4]])]+=1; #print(Hep2[j+1],Hep1[i+4])
			"""

		if ((Hep1[i]=='d') and (Hep2[j]=='a')):

			Res_pairs_ad[tuple([Seq2[j],Seq1[i]])]+=1
			
			if (i+1) < len(Seq1):
				#Seq2[j] ==> a, Seq1[i+1] ==> e
				Res_pairs_ae[tuple([Seq2[j],Seq1[i+1]])]+=1; #print(Hep2[j],Hep1[i+1])

			if (j-1) >= 0:
				#Seq1[i] ==> d, Seq2[j-1] ==> g
				Res_pairs_dg[tuple([Seq1[i],Seq2[j-1]])]+=1; #print(Hep1[i],Hep2[j-1])

			"""
			#gg pairs are added by above a-d pairings
			if ((j-1) >= 0) and ((i+3) < len(Seq1)):
				#Seq1[i-1] ==> g, Seq2[j+3] ==> g
				#Res_pairs_gg[tuple([Seq2[j-1],Seq1[i+3]])]+=1; #print(Hep2[j-1],Hep1[i+3])
				pass

			if ((j-1) >= 0) and ((i-4) >= 0):
				#Seq1[i-1] ==> g, Seq2[j-4] ==> g
				#Res_pairs_gg[tuple([Seq2[j-1],Seq1[i-4]])]+=1; #print(Hep2[j-1],Hep1[i-4])
				pass
			"""
			if ((j-3) >= 0) and ((i+1) < len(Seq1)):
				#Seq1[i-3] ==> e, Seq2[j+1] ==> e
				Res_pairs_ee[tuple([Seq1[i+1],Seq2[j-3]])]+=1; #print(Hep1[i+1],Hep2[j-3])
				

			if ((i+1) < len(Seq1)) and ((j+4) < len(Seq1)):
				#Seq1[i-1] ==> e, Seq2[j-4] ==> e
				Res_pairs_ee[tuple([Seq1[i+1],Seq2[j+4]])]+=1; #print(Hep1[i+1],Hep2[j+4])
				

	Res_stats = [Res_pairs_aa,Res_pairs_dd,Res_pairs_ad,Res_pairs_ag,Res_pairs_de,Res_pairs_ae,Res_pairs_dg,Res_pairs_ee,Res_pairs_gg,Res_pairs_eg]
	return Res_stats



