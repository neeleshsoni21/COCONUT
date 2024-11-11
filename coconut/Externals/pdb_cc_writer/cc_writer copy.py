
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

# PDB Coiled-coil writer
# Author: Neelesh Soni

#import IMP
#import IMP.core

import random
import os
import sys
import numpy as np
from collections import OrderedDict, defaultdict

from src.generate_helix import generate_helix_locus
from src.generate_helix import generate_nstranded_cc_locus

from src.pdbclass import PDB
from src.read_seq import read_fasta
from src.generate_pdb import generate_helix_coords

class DefaultDictNoAdd(defaultdict):
    def __missing__(self, key):
        """Returns a default value when a key is missing.

        Args:
            key: The key that is missing.

        Returns:
            list: A default list value ['-', '-', 0].
        """
        return ['-', '-', 0]

def read_coiled_coil_parameter_file(fname):
    """Reads the coiled-coil parameter file and returns the parameters as a dictionary.

    Args:
        fname (str): The name of the file containing the coiled-coil parameters.

    Returns:
        dict: A dictionary containing the coiled-coil parameters.
    """
    inf = open(fname, 'r')
    lines = inf.readlines()
    inf.close()

    params = {}
    for l in lines:
        if l[0] == '#' or len(l) == 1:
            continue
        toks = l.strip().split()
        assert len(toks) == 3

        params[toks[0]] = toks[2]

    try:
        # Check all keys
        params['NumStrand'] = int(params['NumStrand'])
        params['Radius_Major_Helix'] = float(params['Radius_Major_Helix'])
        params['Radius_Minor_Helix'] = float(params['Radius_Minor_Helix'])
        params['h'] = float(params['h'])
        params['Pitch'] = float(params['Pitch'])
        params['Periodicity'] = float(params['Periodicity'])
        params['Interface_angle'] = float(params['Interface_angle'])
        params['Register'] = str(params['Register'])
    except:
        print("Incorrect Spelling of parameters")
        sys.exit()

    return params

def get_coiled_coil_PDB(fcname, Helices, QrySeqs, ccparams, cc_segments, LINEAR_TRANS=5, LINEAR_SEGMENTS=True, singlestrandonly=False):
    """Generates the PDB representation of the coiled-coil structure.

    Args:
        fcname (str): File name to save the PDB.
        Helices (list): Helices data.
        QrySeqs (list): List of query sequences.
        ccparams (dict): Coiled-coil parameters.
        cc_segments (list): Coiled-coil segments.
        LINEAR_TRANS (int, optional): Linear translation distance. Defaults to 5.
        LINEAR_SEGMENTS (bool, optional): Whether to use linear segments. Defaults to True.
        singlestrandonly (bool, optional): If True, generates a single strand. Defaults to False.

    Returns:
        list: Helices information.
    """
    MDL_LISTS = []
    mdl_pairs = []

    ccnumstrands = ccparams['NumStrand']
    # This is useful when you want to write just one helix instead of the whole coiled-coil
    if singlestrandonly:
        ccnumstrands = 1

    for i in range(0, ccnumstrands):
        TEMP_MDL_LISTS = []
        init_zcoord = 0

        for j, ccseg in enumerate(cc_segments):
            x1, y1, z1 = Helices[i][j]

            ##############################
            mdl_new = PDB("")
            chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
            Coords = [x1, y1, z1]
            chain = chains[i]

            QrySegment = QrySeqs[i][int(ccseg[i][0]) - 1:ccseg[i][1]]
            initresnum = int(ccseg[i][0])

            mdl_new = generate_helix_coords(mdl_new, QrySegment, initresnum, init_zcoord, Coords, chain)

            if LINEAR_SEGMENTS:
                # Translate next segment by LINEAR_TRANS distance
                init_zcoord += z1[-1] + LINEAR_TRANS

            TEMP_MDL_LISTS.append(mdl_new)

        MDL_LISTS.append(TEMP_MDL_LISTS)

    if not LINEAR_SEGMENTS:
        TRandInitCoords = []
        for j in range(len(MDL_LISTS[0])):
            new_x = random.randint(0, 300)
            new_y = random.randint(0, 300)
            new_z = random.randint(0, 300)
            icoord = [new_x, new_y, new_z]

            TRandInitCoords.append(icoord)

        for j in range(len(MDL_LISTS[0])):
            for i in range(len(MDL_LISTS)):
                mdl_new = MDL_LISTS[i][j]
                mdl_new.translate_model(TRandInitCoords[j])

    nfcname = fcname + '.pdb'

    if os.path.exists(nfcname):
        os.remove(nfcname)

    for i in range(len(MDL_LISTS)):
        for j in range(len(MDL_LISTS[0])):
            mdl_new = MDL_LISTS[i][j]
            mdl_new.write_all_chains(nfcname)

    return Helices

def get_coiled_coil_Locus(Helices, cc_segments, ccparams):
    """Generates the locus of the coiled-coil structure.

    Args:
        Helices (list): List of helices.
        cc_segments (list): List of coiled-coil segments.
        ccparams (dict): Coiled-coil parameters.

    Returns:
        list: List representing the locus of the coiled-coil structure.
    """
    Locus = []
    start_hep = ccparams['Register']
    heplist = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    hepid = heplist.index(start_hep)

    AtomHeptad_List = []
    AtomRadius_List = []

    ccnumstrands = ccparams['NumStrand']

    for i in range(0, ccnumstrands):
        for j, ccseg in enumerate(cc_segments):
            x1, y1, z1 = Helices[i][j]

        AtomHeptad_List.append(heplist[1])
        AtomRadius_List.append(1)
        Locus.append(list(zip(x1[:], y1[:], z1[:], AtomRadius_List, AtomHeptad_List)))

    return Locus

def write_pymol_visualization(PyW, frameno, ParticlesList):
    """Writes visualization data for PyMol.

    Args:
        PyW: PyMol Writer object.
        frameno (int): Frame number for the visualization.
        ParticlesList (list): List of particles to visualize.
    """
    PyW.set_frame(frameno)
    for p in ParticlesList:
        # Display core::XYZR particle as a ball
        g = IMP.core.XYZRGeometry(p)
        PyW.add_geometry(g)
    return

def change_radius(m, ParticlesList):
    """Changes the radius of particles.

    Args:
        m: Model object.
        ParticlesList (list): List of particles.
    """
    for p in ParticlesList:
        # Since particle is already setup, just use the decorator with model
        d = IMP.core.XYZR(m, p)
        # Set new random radius
        d.set_radius(d.get_radius() + random.randint(-1, 1))
    return

def move_particles(m, ParticlesList):
    """Moves particles in random directions.

    Args:
        m: Model object.
        ParticlesList (list): List of particles.
    """
    for p in ParticlesList:
        # Since particle is already setup, just use the decorator with model
        d = IMP.core.XYZR(m, p)

        randtrans_x = random.randint(-3, 3)
        randtrans_y = random.randint(-3, 3)
        randtrans_z = random.randint(-3, 3)

        # Get coordinates
        x, y, z = d.get_coordinates()
        # Set the coordinates again
        d.set_coordinates((x + randtrans_x, y + randtrans_y, z + randtrans_z))
    return

def model_coiled_coil(fname, Helices):
    """Creates a coiled-coil model and writes visualization data for PyMol.

    Args:
        fname (str): File name to save the PyMol visualization.
        Helices (list): List of helices data.
    """
    # Create a model
    m = IMP.Model()

    ParticlesList = []

    for ResidueList in Helices:

        # Add particles for each residue
        for i, res in enumerate(ResidueList):

            # Create a generalized IMP particle
            p = IMP.Particle(m)
            p.set_name("P" + str(i))

            # Get appropriate decorators for a particle
            # Decorate with 3D vector
            pcoord = IMP.algebra.Vector3D(res[0], res[1], res[2])
            # Decorate with 3D sphere
            psphere = IMP.algebra.Sphere3D(pcoord, res[3])

            # Decorate IMP particle with 3D vector and its radius
            d = IMP.core.XYZR.setup_particle(p, psphere)

            # Don't let the local coordinate of the particle change by the optimizer
            d.set_coordinates_are_optimized(False)

            # Select color according to heptad
            hepcolor = [1, 2, 3, 4, 5, 6, 7]
            selcolor = hepcolor[i % 7]

            # Add color decorator for particle
            IMP.display.Colored.setup_particle(p, IMP.display.get_display_color(selcolor))

            ParticlesList.append(p)

    # Write to PyMol
    PyW = IMP.display.PymolWriter(fname)
    write_pymol_visualization(PyW, 0, ParticlesList)

    """
    for i in range(0,50):
        
        frameno = i
        write_pymol_visualization(PyW, frameno, ParticlesList)

        change_radius(m, ParticlesList)

        move_particles(m, ParticlesList)
    """
    del PyW

    return

def load_paircoil_alignment(fname):
    """Loads paircoil alignment data from a file, filtering based on a probability cutoff.

    Args:
        fname (str): The name of the file containing paircoil alignment data.

    Returns:
        OrderedDict: An ordered dictionary where keys are residue positions and values are lists containing
        residue information, heptad designation, and probability.
    """
    inf = open(fname,'r')
    lines = inf.readlines()
    inf.close()


    prot_seq = OrderedDict()

    for l in lines:
        if l[0]=='#':
            continue

        if len(l)<2:
            continue

        toks = l.strip().split()
    
        if float(toks[3]) >=PROB_CUTOFF:
            prot_seq[int(toks[0])]=[toks[1],toks[2],float(toks[3])]
        else:
            prot_seq[int(toks[0])]=[toks[1],'-',float(toks[3])]

    return prot_seq

def load_paircoil_alignment_nofilter(fname):
    """Loads paircoil alignment data from a file without applying a probability cutoff filter.

    Args:
        fname (str): The name of the file containing paircoil alignment data.

    Returns:
        OrderedDict: An ordered dictionary where keys are residue positions and values are lists containing
        residue information, heptad designation, and probability.
    """
    inf = open(fname,'r')
    lines = inf.readlines()
    inf.close()


    prot_seq = OrderedDict()

    for l in lines:
        if l[0]=='#':
            continue

        if len(l)<2:
            continue

        toks = l.strip().split()
    
        
        prot_seq[int(toks[0])]=[toks[1],toks[2],float(toks[3])]

    return prot_seq

def load_Pcoil_alignment(pc_file):
    """Loads Pcoil alignment data from a file, filtering based on a probability cutoff for multiple columns.

    Args:
        pc_file (str): The name of the file containing Pcoil alignment data.

    Returns:
        DefaultDictNoAdd: A default dictionary where keys are residue positions and values are lists containing
        residue information, heptad designation, and probability, with values determined by multiple probability columns.
    """
    prob_col1=1;prob_col2=2;prob_col3=3;
    
    inf = open(pc_file,'r')
    lines = inf.readlines()
    inf.close()

    prot_seq = DefaultDictNoAdd()

    header = ""
    for l in lines:

        if l[0]=='>':
            continue
        #toks = l.strip().split('\t')
        toks = l.strip().split()

        if len(toks)<8:
            header+='#'+l
            continue
        
        probidx1 = 2*prob_col1+1;heptadidx1 = 2*prob_col1
        probidx2 = 2*prob_col2+1;heptadidx2 = 2*prob_col2
        probidx3 = 2*prob_col3+1;heptadidx3 = 2*prob_col3

        if float(toks[probidx3]) >= PROB_CUTOFF:
            prot_seq[int(toks[0])]=[toks[1],toks[heptadidx3],float(toks[probidx3])]

        elif float(toks[probidx2]) >= PROB_CUTOFF:
            prot_seq[int(toks[0])]=[toks[1],toks[heptadidx2],float(toks[probidx2])]

        elif float(toks[probidx1]) >= PROB_CUTOFF:
            prot_seq[int(toks[0])]=[toks[1],toks[heptadidx1],float(toks[probidx1])]
        
        else:
            prot_seq[int(toks[0])]=[toks[1],'-',float(toks[probidx1])]

    return prot_seq

def read_alignment_file(fname, prot1, MLP1_PC, prot2, MLP2_PC):
    """Reads an alignment file and extracts coiled-coil segments based on provided protein data.

    Args:
        fname (str): The name of the alignment file.
        prot1 (str): The identifier for the first protein.
        MLP1_PC (dict): Paircoil data for the first protein.
        prot2 (str): The identifier for the second protein.
        MLP2_PC (dict): Paircoil data for the second protein.

    Returns:
        list: A list of coiled-coil segments, each represented by residue pairs and heptad designation.
    """
    inf = open(fname,'r')
    lines = inf.readlines()
    inf.close()

    cc_segments=[]

    for l in lines:
        
        if ((l[0]=='-') or (l[0]=='#') or len(l)==1):
            continue

        toks = l.strip().split('\t')
        
        respair1 = list(map(int,toks[1].split('-')))
        respair2 = list(map(int,toks[3].split('-')))

        assert toks[0]==prot1
        assert toks[2]==prot2

        #Assertion not working for mlp2 seq for one segment due to incoorect mlp2 paircoil assignment
        #assert MLP1_PC[respair1[0]][1]==MLP2_PC[respair2[0]][1]
        
        starthep = MLP1_PC[respair1[0]][1]

        cc_segments.append([ respair1, respair2, starthep ])

    return cc_segments

def run_example1(RootDir='./', IMP_Model=False):
    """Runs a coiled-coil modeling example based on a specified input sequence and writes the output to a PDB file.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.

    Returns:
        None
    """

    fname = os.path.join(RootDir,'./test_example/2zta/2ztaorg.fasta')
    oname = os.path.join(RootDir,'./test_example/2zta/2ztaorg_pseudo_chain')

    print("Input File:",fname)
    print("Output File:",oname)

    QrySeq, Seqid = read_fasta(fname)
    no_of_residues = len(QrySeq)
    print("Number of Residues/helix:",no_of_residues)

    ccparams = read_coiled_coil_parameter_file(os.path.join(RootDir,'./input/standard_cc_dimer.par'))

    cc_segments = [[[1,len(QrySeq)],[1,len(QrySeq)], ccparams['Register']]]
    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[ccparams['NumStrand']]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[ccparams['NumStrand']])

        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)
        
        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])

    get_coiled_coil_PDB(oname, Helices, [QrySeq,QrySeq], ccparams, cc_segments, LINEAR_SEGMENTS=True)

    if IMP_Model == True:
        Locus = get_coiled_coil_Locus(Helices, ccparams)
        model_coiled_coil('./test_example/2zta/Chain-particles.pym', Locus)

    return

def run_example2(RootDir='./', IMP_Model=False):
    """Runs a coiled-coil modeling example based on a specified input sequence and writes the output to a PDB file.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.

    Returns:
        None
    """
    fname = './test_example/2zta_long/2zta_long_17mer.fasta'
    oname = './test_example/2zta_long/2zta_long_pseudo'

    print("Input File:",fname)
    print("Output File:",oname)

    QrySeq, Seqid = read_fasta(fname)
    no_of_residues = len(QrySeq)
    print("Number of Residues/helix:",no_of_residues)

    ccparams = read_coiled_coil_parameter_file('./input/standard_cc_dimer_long.par')

    cc_segments = [[[1,len(QrySeq)],[1,len(QrySeq)], ccparams['Register']]]
    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[ccparams['NumStrand']]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[ccparams['NumStrand']])

        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)
        
        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])

    get_coiled_coil_PDB(oname, Helices, [QrySeq,QrySeq], ccparams, cc_segments, LINEAR_SEGMENTS=True)

    if IMP_Model == True:
        Locus = get_coiled_coil_Locus(Helices, ccparams)
        model_coiled_coil('./test_example/2zta_long/Chain-particles.pym', Locus)

    return

def run_example3(RootDir='./', IMP_Model=False):
    """Runs a coiled-coil modeling example based on a specified input sequence and writes the output to a PDB file.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.

    Returns:
        None
    """
    fname = './test_example/ccbuilder_dimer/ccbuilder.fasta'
    oname = './test_example/ccbuilder_dimer/ccbuilder_pseudo'

    print("Input File:",fname)
    print("Output File:",oname)

    QrySeq, Seqid = read_fasta(fname)
    no_of_residues = len(QrySeq)
    print("Number of Residues/helix:",no_of_residues)

    ccparams = read_coiled_coil_parameter_file('./input/standard_cc_dimer_ccbuilder.par')

    cc_segments = [[[1,len(QrySeq)],[1,len(QrySeq)], ccparams['Register']]]
    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[ccparams['NumStrand']]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[ccparams['NumStrand']])

        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)
        
        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])

    get_coiled_coil_PDB(oname, Helices, [QrySeq,QrySeq], ccparams, cc_segments, LINEAR_SEGMENTS=True)

    if IMP_Model == True:
        Locus = get_coiled_coil_Locus(Helices, ccparams)
        model_coiled_coil('./test_example/ccbuilder_dimer/Chain-particles.pym', Locus)

    return

def run_example4(RootDir='./', IMP_Model=False):
    """Runs a coiled-coil modeling example based on a specified input sequence and writes the output to a PDB file.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.

    Returns:
        None
    """
    fname = './test_example/trimer/1flk_ccdomain.fasta'
    oname = './test_example/trimer/Standard_trimer'

    print("Input File:",fname)
    print("Output File:",oname)

    QrySeq, Seqid = read_fasta(fname)
    no_of_residues = len(QrySeq)
    print("Number of Residues/helix:",no_of_residues)

    ccparams = read_coiled_coil_parameter_file('./input/standard_cc_trimer.par')

    cc_segments = [[[1,len(QrySeq)],[1,len(QrySeq)],[1,len(QrySeq)],ccparams['Register']]]
    
    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[ccparams['NumStrand']]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[ccparams['NumStrand']])

        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)

        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])

    get_coiled_coil_PDB(oname, Helices, [QrySeq,QrySeq,QrySeq], ccparams, cc_segments, LINEAR_SEGMENTS=True)

    if IMP_Model == True:
        Locus = get_coiled_coil_Locus(Helices, ccparams)
        model_coiled_coil('./test_example/trimer/Chain-particles.pym', Locus)

    return

def run_example7(RootDir='./', IMP_Model=False, cc_segments=[]):
    """Runs a coiled-coil modeling example for the Nup60 protein, generating a PDB file representation.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.
        cc_segments (list, optional): List of coiled-coil segments to model. Defaults to an empty list.

    Returns:
        None
    """
    fname1 = './test_example/nup60/Nup60.fasta'
    fname2 = './test_example/nup60/Nup60.fasta'

    oname = './test_example/nup60/Nup60'

    print("Input File:",fname1)
    print("Input File:",fname2)
    print("Output File:",oname)

    QrySeq1, Seqid1 = read_fasta(fname1)
    QrySeq2, Seqid2 = read_fasta(fname2)

    no_of_residues1 = len(QrySeq1)
    no_of_residues2 = len(QrySeq2)

    print("Number of Residues/helix:",no_of_residues1,no_of_residues2)

    ccparams = read_coiled_coil_parameter_file('./input/nup60_helix.par')

    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[ccparams['NumStrand']]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[ccparams['NumStrand']])

        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)
        
        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])
    
    get_coiled_coil_PDB(oname, Helices, [QrySeq1,QrySeq2], 
        ccparams, cc_segments, LINEAR_SEGMENTS=True, singlestrandonly=True)

    if IMP_Model == True:
        Locus = get_coiled_coil_Locus(Helices, ccparams)
        model_coiled_coil('./test_example/nup60/Chain-particles.pym', Locus)

    return

def run_example8(RootDir='./', IMP_Model=False, cc_segments=[]):
    """Runs a modeling example for the Nup1 protein, generating a PDB file representation.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.
        cc_segments (list, optional): List of coiled-coil segments to model. Defaults to an empty list.

    Returns:
        None
    """
    fname1 = './test_example/nup1/Nup1.fasta'
    fname2 = './test_example/nup1/Nup1.fasta'

    oname = './test_example/nup1/Nup1'

    print("Input File:",fname1)
    print("Input File:",fname2)
    print("Output File:",oname)

    QrySeq1, Seqid1 = read_fasta(fname1)
    QrySeq2, Seqid2 = read_fasta(fname2)

    no_of_residues1 = len(QrySeq1)
    no_of_residues2 = len(QrySeq2)

    print("Number of Residues/helix:",no_of_residues1,no_of_residues2)

    ccparams = read_coiled_coil_parameter_file('./input/nup1_helix.par')

    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[ccparams['NumStrand']]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[ccparams['NumStrand']])

        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)
        
        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])
        
    get_coiled_coil_PDB(oname, Helices, [QrySeq1,QrySeq2], 
        ccparams, cc_segments, LINEAR_SEGMENTS=True, singlestrandonly=True)

    if IMP_Model == True:
        Locus = get_coiled_coil_Locus(Helices, ccparams)
        model_coiled_coil('./test_example/nup1/Chain-particles.pym', Locus)

    return

def run_example9(RootDir='./', IMP_Model=False, cc_segments=[]):
    """Runs a coiled-coil modeling example for the generic MLP protein, generating a PDB file representation.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.
        cc_segments (list, optional): List of coiled-coil segments to model. Defaults to an empty list.

    Returns:
        None
    """
    fname1 = './test_example/mlp_cc_cw_unified/MLP.fasta'
    fname2 = './test_example/mlp_cc_cw_unified/MLP.fasta'
    oname = './test_example/mlp_cc_cw_unified/MLP_MLP_Random'

    prot1='MLP'
    prot2='MLP'

    alifile='./test_example/mlp_cc_cw_unified/'+prot1+'_'+prot2+'_CW_alignment.cali'
    paircoilfile1 = './test_example/mlp_cc_cw_unified/'+prot1+'_paircoil2.out'
    paircoilfile2 = './test_example/mlp_cc_cw_unified/'+prot2+'_paircoil2.out'
    
    #Sequence has both amino acid sequence and heptad repeat alignment form the pair coil
    MLP1_PC = load_paircoil_alignment(paircoilfile1)
    MLP2_PC = load_paircoil_alignment(paircoilfile2)
    cc_segments = read_alignment_file(alifile, prot1, MLP1_PC, prot2, MLP2_PC)
    
    
    print("Input File:",fname1)
    print("Input File:",fname2)
    print("Output File:",oname)

    QrySeq1, Seqid1 = read_fasta(fname1)
    QrySeq2, Seqid2 = read_fasta(fname2)

    no_of_residues1 = len(QrySeq1)
    no_of_residues2 = len(QrySeq2)

    print("Number of Residues/helix:",no_of_residues1,no_of_residues2)

    ccparams = read_coiled_coil_parameter_file('./input/mlp_ccdimer.par')

    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[2]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[2])
        
        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)
        
        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])

    get_coiled_coil_PDB(oname, Helices, [QrySeq1,QrySeq2], 
        ccparams, cc_segments, LINEAR_TRANS = 20, LINEAR_SEGMENTS=False)

    if IMP_Model == True:
        Locus = get_coiled_coil_Locus(Helices, ccparams)
        model_coiled_coil('./test_example/mlp_cc_cw_unified/Chain-particles.pym', Locus)

    return

def run_example12(RootDir='./', IMP_Model=False, cc_segments=[]):
    """Runs a coiled-coil modeling example for the Mouse TPR protein, generating a PDB file representation.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.
        cc_segments (list, optional): List of coiled-coil segments to model. Defaults to an empty list.

    Returns:
        None
    """
    fname1 = './test_example/MuTPR/MuTPR.fasta'
    fname2 = './test_example/MuTPR/MuTPR.fasta'
    oname = './test_example/MuTPR/MuTPR_MuTPR'

    prot1='MuTPR'
    prot2='MuTPR'

    alifile='./test_example/MuTPR/'+prot1+'_'+prot2+'_alignment.cali'
    paircoilfile1 = './test_example/MuTPR/'+prot1+'_paircoil2_new.out'
    paircoilfile2 = './test_example/MuTPR/'+prot2+'_paircoil2_new.out'
    
    #Sequence has both amino acid sequence and heptad repeat alignment form the pair coil
    MLP1_PC = load_paircoil_alignment_nofilter(paircoilfile1)
    MLP2_PC = load_paircoil_alignment_nofilter(paircoilfile2)
    cc_segments = read_alignment_file(alifile, prot1, MLP1_PC, prot2, MLP2_PC)
    
    print("Input File:",fname1)
    print("Input File:",fname2)
    print("Output File:",oname)

    QrySeq1, Seqid1 = read_fasta(fname1)
    QrySeq2, Seqid2 = read_fasta(fname2)

    no_of_residues1 = len(QrySeq1)
    no_of_residues2 = len(QrySeq2)

    print("Number of Residues/helix:",no_of_residues1,no_of_residues2)

    ccparams = read_coiled_coil_parameter_file('./input/mu_tpr_ccdimer.par')

    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[ccparams['NumStrand']]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[ccparams['NumStrand']])

        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)
        
        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])
        
    get_coiled_coil_PDB(oname, Helices, [QrySeq1,QrySeq2], 
        ccparams, cc_segments, LINEAR_TRANS = 20, LINEAR_SEGMENTS=True)

    if IMP_Model == True:
        Locus = get_coiled_coil_Locus(Helices, ccparams)
        model_coiled_coil('./test_example/numa_cc/Chain-particles.pym', Locus)

    return

def run_example13(RootDir='./', IMP_Model=False, cc_segments=[]):
    """Runs a coiled-coil modeling example for the Mouse Nup153 protein, generating a PDB file representation.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.
        cc_segments (list, optional): List of coiled-coil segments to model. Defaults to an empty list.

    Returns:
        None
    """
    fname1 = './test_example/munup153/Nup153.fasta'
    fname2 = './test_example/munup153/Nup153.fasta'

    oname = './test_example/munup153/Nup153'

    print("Input File:",fname1)
    print("Input File:",fname2)
    print("Output File:",oname)

    QrySeq1, Seqid1 = read_fasta(fname1)
    QrySeq2, Seqid2 = read_fasta(fname2)

    no_of_residues1 = len(QrySeq1)
    no_of_residues2 = len(QrySeq2)

    print("Number of Residues/helix:",no_of_residues1,no_of_residues2)

    ccparams = read_coiled_coil_parameter_file('./input/nup153_helix.par')

    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[ccparams['NumStrand']]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[ccparams['NumStrand']])

        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)
        
        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])
        
    get_coiled_coil_PDB(oname, Helices, [QrySeq1,QrySeq2], 
        ccparams, cc_segments, LINEAR_SEGMENTS=True, singlestrandonly=True)

    if IMP_Model == True:
        #Locus = get_coiled_coil_Locus(Helices, ccparams)
        #model_coiled_coil('./test_example/nup60/Chain-particles.pym', Locus)
        pass

    return

def run_example15(RootDir='./', IMP_Model=False, cc_segments=[]):
    """Runs a coiled-coil modeling example for the K5K14 dimer protein, generating a PDB file representation.

    Args:
        IMP_Model (bool): Flag indicating whether to use IMP modeling. If True, additional modeling steps are performed.
        cc_segments (list, optional): List of coiled-coil segments to model. Defaults to an empty list.

    Returns:
        None
    """
    fname1 = './test_example/k5k14/k5.fasta'
    fname2 = './test_example/k5k14/k14.fasta'
    oname = './test_example/k5k14/k5k14'

    prot1='k5'
    prot2='k14'
    #alifile='./test_example/mlp_cc/'+prot1+'_'+prot2+'_alignment.cali'

    alifile='./test_example/k5k14/'+prot1+'_'+prot2+'_alignment.cali'
    paircoilfile1 = './test_example/k5k14/'+prot1+'_pcoils.out'
    paircoilfile2 = './test_example/k5k14/'+prot2+'_pcoils.out'
    
    #Sequence has both amino acid sequence and heptad repeat alignment form the pair coil
    MLP1_PC = load_Pcoil_alignment(paircoilfile1)
    MLP2_PC = load_Pcoil_alignment(paircoilfile2)
    cc_segments = read_alignment_file(alifile, prot1, MLP1_PC, prot2, MLP2_PC)
    
    print("Input File:",fname1)
    print("Input File:",fname2)
    print("Output File:",oname)

    QrySeq1, Seqid1 = read_fasta(fname1)
    QrySeq2, Seqid2 = read_fasta(fname2)

    no_of_residues1 = len(QrySeq1)
    no_of_residues2 = len(QrySeq2)

    print("Number of Residues/helix:",no_of_residues1,no_of_residues2)

    ccparams = read_coiled_coil_parameter_file('./input/k5k14_ccdimer.par')

    Helices=[]
    for stid in range(0,ccparams['NumStrand']):
        Helices.append([])

    for ccseg in cc_segments:

        no_of_residues1 = (int(ccseg[0][1])-int(ccseg[0][0]))+1
        no_of_residues2 = (int(ccseg[1][1])-int(ccseg[1][0]))+1

        assert no_of_residues1==no_of_residues2

        assert( str(ccseg[2]) in ['a', 'b', 'c', 'd', 'e', 'f', 'g'] )

        ccparams['Register']=str(ccseg[2])

        Temp_Helices = generate_nstranded_cc_locus(no_of_residues1, ccparams)
        
        for stid in range(0,ccparams['NumStrand']):
            Helices[stid].append(Temp_Helices[stid])
        
    Helices = get_coiled_coil_PDB(oname, Helices, [QrySeq1,QrySeq2], 
        ccparams, cc_segments, LINEAR_TRANS = 20, LINEAR_SEGMENTS=True)

    if IMP_Model == True:
        
        Locus = get_coiled_coil_Locus(Helices, cc_segments, ccparams)
        model_coiled_coil('./test_example/k5k14/Chain-particles.pym', Locus)

    return


PROB_CUTOFF=0.20

if __name__ == '__main__':

    IMP_Model=False #NOT WORKING with TRUE

    #Modeling of 2zta pdb
    run_example1(IMP_Model)

    #Modeling of 2zta long
    run_example2(IMP_Model)

    #Modeling for cc builder comparison
    run_example3(IMP_Model)

    #Modeling of standard trimer
    run_example4(IMP_Model)

    #This is for Nup60 modeling
    cc_segments=[ [[27,47],[27,47],'g'],  [[91,104],[91,104],'g'],[[106,119],[106,119],'g'],  [[121,140],[121,140],'g'],[[142,162],[142,162],'g']]
    run_example7(IMP_Model, cc_segments)
    
    #This is for Nup1 modeling
    cc_segments=[ [[1,32],[1,32],'g'],  [[85,104],[85,104],'g'],[[106,123],[106,123],'g']]
    run_example8(IMP_Model, cc_segments)

    #This is for generic or unified MLPs model
    run_example9(IMP_Model)

    #This is for Mouse TPR Modeling
    run_example12(IMP_Model)

    #This is for Nup153 modeling mouse
    cc_segments=[ [[36,57],[36,57],'g'] ]
    run_example13(IMP_Model, cc_segments)

    #This is for Keratin dimer modeling
    run_example15(IMP_Model)


    

#################################################
#Benchmarking
#Results of comparisons of chain A superimpositions
    # ccbuilder 0.401
    # 2zta 0.566
    # 2zta_long_17mer 0.780
    # standard trimer 0.519
#################################################



