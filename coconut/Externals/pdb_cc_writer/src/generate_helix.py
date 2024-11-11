"""Summary

Attributes:
    deltaa (TYPE): Description
"""
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


import numpy as np

# w1 is calculated from the vector pointing outwards (np.pi), 
# a/d positions points towards the core (2*np.pi/7)/2
deltaa = np.pi - (2*np.pi/7)/2

__REGISTER_ADJUST = {
    'a': 0 * np.pi / 180 - deltaa,
    'b': 1*(720/7) * np.pi/ 180 - deltaa,
    'c': 2*(720/7) * np.pi/ 180 - deltaa,
    'd': 3*(720/7) * np.pi/ 180 - deltaa,
    'e': 4*(720/7) * np.pi/ 180 - deltaa,
    'f': 5*(720/7) * np.pi/ 180 - deltaa,
    'g': 6*(720/7) * np.pi/ 180 - deltaa

    #After shifting by 6pi/7, absolute number will be as follows
    #a=0, b=4pi/7, c=8pi/7, d=12pi/7, e=16pi/7, f=20pi/7, g=24pi/7
}

def rad2deg(v):
    """Convert radians to degrees.

    Args:
        v (float or np.array): Angle in radians.

    Returns:
        float or np.array: Angle converted to degrees.
    """
    return v*180/np.pi

def deg2rad(v):
    """Convert degrees to radians.

    Args:
        v (float or np.array): Angle in degrees.

    Returns:
        float or np.array: Angle converted to radians.
    """
    return v*np.pi/180.0

def generate_nstranded_cc_locus(no_of_residues, ccparams):
    """Generate coordinates for an n-stranded coiled-coil locus.

    Args:
        no_of_residues (int): The number of residues in the coiled-coil.
        ccparams (dict): Parameters for coiled-coil generation, including:
            - 'Interface_angle' (float): Angle between strands at the interface (degrees).
            - 'Radius_Major_Helix' (float): Radius of the major helix.
            - 'h' (float): Helical rise per turn.
            - 'Pitch' (float): Pitch of the superhelix.
            - 'Radius_Minor_Helix' (float): Radius of the minor helix.
            - 'Periodicity' (float): Number of residues per helical turn.
            - 'Register' (str): Starting register ('a', 'b', 'c', 'd', 'e', 'f', 'g').
            - 'NumStrand' (int): Number of strands in the coiled-coil.

    Returns:
        list: A list of helices, each represented by [x-coordinates, y-coordinates, z-coordinates].
    """
    t = np.arange(0,no_of_residues);

    #check this if true for trimers and others
    alpha = ccparams['Interface_angle']/2.0;
    alpha = alpha*np.pi/180 #Convert to radians

    r0 = ccparams['Radius_Major_Helix']
    w0 = 2*np.pi*ccparams['h']/ccparams['Pitch']    #Super helix is left handed (leucine zipper)., +ive p0 makes it right handed
    
    r1 = ccparams['Radius_Minor_Helix']
    w1 = 2*np.pi*(1/ccparams['Periodicity'])    #w1 +ive for right handed helix

    register = ccparams['Register']

    NumStrand = ccparams['NumStrand']

    p0 = ccparams['Pitch']

    Helices = []

    for i in range(1,NumStrand+1):

        xis, yis, zis = generate_istrand_cc_locus(t,r1,w1,r0,w0,p0,alpha, i, NumStrand, register )

        Helices.append([xis,yis,zis])

    return Helices

def generate_istrand_cc_locus(t,r1,w1,r0,w0,p0 ,alpha, istrand, NumStrand, register):
    
    """Generate coordinates for an individual strand of a coiled-coil.

    Args:
        t (np.array): Array representing residue positions.
        r1 (float): Radius of the minor helix.
        w1 (float): Angular frequency of the minor helix.
        r0 (float): Radius of the major helix.
        w0 (float): Angular frequency of the major helix.
        p0 (float): Pitch of the superhelix.
        alpha (float): Interface angle in radians.
        istrand (int): Index of the current strand.
        NumStrand (int): Total number of strands.
        register (str): Starting register for the strand ('a', 'b', 'c', 'd', 'e', 'f', 'g').

    Returns:
        tuple: x, y, z coordinates for the generated strand.
    """
    #start register adjustment
    srs = 1*__REGISTER_ADJUST[register]
    
    #ith strand adjustment
    psi = 1*2*np.pi*(istrand-1)/NumStrand; #0, or pi for 2 stranded
    #psi will be added to w0*t big helix rotation NOT in w1*t

    #Equations after adjustments
    x1 = r0*np.cos(w0*t+psi) + r1*np.cos(w0*t+psi)*np.cos(w1*t+srs) - r1*np.cos(alpha)*np.sin(w0*t+psi)*np.sin(w1*t+srs)
    y1 = r0*np.sin(w0*t+psi) + r1*np.sin(w0*t+psi)*np.cos(w1*t+srs) + r1*np.cos(alpha)*np.cos(w0*t+psi)*np.sin(w1*t+srs)
    z1 = (p0*( (w0*t)/(2*np.pi)) - r1*np.sin(alpha)*np.sin(w1*t))
    
    return x1,y1,z1


#Ideal helix locus (NOT coiled-coil)
#Good for strainght helix, in disctransformation
#No superhelical turn/pitch. Straight helix
#x1,y1,z1 = generate_helix_locus(t,r1,p1,w1, InitPos, theta1)

def generate_helix_locus(t,r1,p1,w1, InitPos, theta1):
    """Generate coordinates for a standard (non-coiled-coil) helix.

    Args:
        t (np.array): Array representing residue positions.
        r1 (float): Radius of the helix.
        p1 (float): Pitch of the helix.
        w1 (float): Angular frequency of the helix.
        InitPos (tuple): Initial position (xi, yi, zi) of the helix.
        theta1 (float): Initial angular shift for the helix.

    Returns:
        tuple: x, y, z coordinates for the generated helix.
    """
    xi,yi,zi = InitPos
    t = t + theta1
    x1 = r1 * np.cos(w1*t) +xi
    y1 = r1 * np.sin(w1*t) +yi
    z1 = 1*p1*(w1*t/(2*np.pi)) +zi
    
    return x1,y1,z1


