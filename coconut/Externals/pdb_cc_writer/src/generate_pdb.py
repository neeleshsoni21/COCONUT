"""Summary

Attributes:
    ONETHREE (dict): Description
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


ONETHREE={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU','S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS','G':'GLY','P':'PRO','A':'ALA','V':'VAL','I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR','W':'TRP'}


def generate_helix_coords2(mdl_new,QrySeq,initresnum,init_coord, Coords,chain):
    """Generate atomic coordinates for a helical region in a protein structure.

    Args:
        mdl_new (object): The model instance that stores the atomic details of the generated helix.
        QrySeq (str): The amino acid sequence of the helix represented as one-letter codes.
        initresnum (int): The starting residue number for the helix.
        init_coord (tuple): Initial (x, y, z) coordinates for positioning the helix in space.
        Coords (tuple): A tuple containing three arrays (x, y, z) that define the relative coordinates for the helix.
        chain (str): Chain ID to assign to the helix.

    Returns:
        object: Updated model instance containing the atomic details for the generated helix.

    Description:
    ------------
    This function generates the atomic coordinates for a given amino acid sequence (`QrySeq`) of a helix.
    It uses the provided initial coordinates (`init_coord`) to determine the starting point of the helix.
    The coordinates for each residue are calculated using the relative position data (`Coords`) and stored
    in `mdl_new`. The function appends information for each atom in the residue including the atom serial number,
    atom name, residue name, chain ID, and spatial coordinates.
    """
    x1,y1,z1 = Coords[0],Coords[1],Coords[2]
    init_xcoord, init_ycoord, init_zcoord = init_coord[0],init_coord[1],init_coord[2]
    iteri = 0
    for residx, resi in enumerate(QrySeq):
        resnum = residx+initresnum

        #for resnm in ['N','CA','C','O']:
        for resname in ['CA']:
            mdl_new.serial().append(iteri)
            mdl_new.name().append(resname)
            mdl_new.altLoc().append(" ")
            mdl_new.resName().append(ONETHREE[resi])
            mdl_new.chainID().append(chain)
            mdl_new.resSeq().append(resnum)
            #mdl_new.iCode().append(" ")
            mdl_new.x().append(x1[residx]+init_xcoord)
            mdl_new.y().append(y1[residx]+init_ycoord)
            mdl_new.z().append(z1[residx]+init_zcoord)
            mdl_new.T().append(1.0)

        iteri+=1
    #mdl_new.iter_length() = iteri
    mdl_new.set_iter_length(iteri)

    return mdl_new

def generate_helix_coords(mdl_new,QrySeq,initresnum,init_zcoord, Coords,chain):
    """Generate atomic coordinates for a helical region in a protein structure along the z-axis.

    Args:
        mdl_new (object): The model instance that stores the atomic details of the generated helix.
        QrySeq (str): The amino acid sequence of the helix represented as one-letter codes.
        initresnum (int): The starting residue number for the helix.
        init_zcoord (float): The initial z-coordinate for positioning the helix along the z-axis.
        Coords (tuple): A tuple containing three arrays (x, y, z) that define the relative coordinates for the helix.
        chain (str): Chain ID to assign to the helix.

    Returns:
        object: Updated model instance containing the atomic details for the generated helix.

    Description:
    ------------
    This function generates the atomic coordinates for a helical segment along the z-axis. It uses the given amino
    acid sequence (`QrySeq`) and shifts the helix along the z-axis by the provided `init_zcoord`. Coordinates are
    calculated from the input arrays (`Coords`) and stored in the model instance (`mdl_new`). The function appends
    information for each atom in the residue including serial number, residue name, chain ID, and spatial coordinates.
    """
    x1,y1,z1 = Coords[0],Coords[1],Coords[2]
    iteri = 0
    for residx, resi in enumerate(QrySeq):
        resnum = residx+initresnum

        #for resnm in ['N','CA','C','O']:
        for resname in ['CA']:
            mdl_new.serial().append(iteri)
            mdl_new.name().append(resname)
            mdl_new.altLoc().append(" ")
            mdl_new.resName().append(ONETHREE[resi])
            mdl_new.chainID().append(chain)
            mdl_new.resSeq().append(resnum)
            #mdl_new.iCode().append(" ")
            mdl_new.x().append(x1[residx])
            mdl_new.y().append(y1[residx])
            mdl_new.z().append(z1[residx]+init_zcoord)
            mdl_new.T().append(1.0)

        iteri+=1
    #mdl_new.iter_length() = iteri
    mdl_new.set_iter_length(iteri)

    return mdl_new