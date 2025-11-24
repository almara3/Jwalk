#!/usr/bin/env python

# ===============================================================================
#     This file is part of Jwalk.
#
#     Jwalk - A tool to calculate the solvent accessible surface distance (SASD)
#     between crosslinked residues.
#
#     Copyright 2016 Jwalk Inventor and Birkbeck College University of London.
#                          The Jwalk Inventor is: Josh Bullock
#
#
#     Jwalk is available under Public Licence.
#     This software is made available under GPL V3
#
#     Please cite your use of Jwalk in published work:
#
#     J.Bullock, J. Schwab, K. Thalassinos, M. Topf (2016)
#     The importance of non-accessible crosslinks and solvent accessible surface distance
#     in modelling proteins with restraints from crosslinking mass spectrometry.
#     Molecular and Cellular Proteomics (15) pp.2491-2500
#
# ===============================================================================

import os
import sys
import argparse
import multiprocessing as mp
from Jwalk import PDBTools, GridTools, SurfaceTools, SASDTools

# default parameters

max_dist = 60
vox = 1
surface = False
xl_path = ''
aa1 = "LYS"
aa2 = "LYS"
ncpus = mp.cpu_count()
pdb = ''

amino_acids = {
    "LYS": "lysines",
    "CYS": "cysteines",
    "ASP": "acidic residues",
    "GLU": "acidic residues",
    "VAL": "valines",
    "ILE": "isoleucines",
    "LEU": "leucines",
    "ARG": "arginines",
    "PRO": "prolines",
    "GLY": "glycines",
    "ALA": "alanines",
    "TRP": "tryptophans",
    "PHE": "phenylalanines",
    "SER": "serines",
    "GLN": "glutamines",
    "HIS": "histidines",
    "MET": "methionines",
    "THR": "threonines",
    "ASN": "asparagines",
    "TYR": "tyrosines"}

parser = argparse.ArgumentParser(description='JWALK: Calculate SASDs on your target PDB files')

parser.add_argument('-lys', action="store_true",
                    help='calculate lysine crosslinks (default)')
parser.add_argument('-xl_path', nargs=1,
                    help='calculate crosslinks from input list')
parser.add_argument('-i', nargs=1,
                    help='specify input pdb: -i <inputfile.pdb>')
parser.add_argument('-aa1', nargs=1,
                    help='specify start amino acid (three letter code e.g. LYS)')
parser.add_argument('-aa2', nargs=1,
                    help='specify end amino acid (three letter code e.g. LYS)')
parser.add_argument('-max_dist', nargs=1,
                    help='specify maximum crosslink distance in Angstroms')
parser.add_argument('-vox', nargs=1,
                    help='specify voxel size of grid')
parser.add_argument('-surface', action="store_true",
                    help='use higher accuracy method to calculate solvent accessibility - requires '
                    'Freesasa installation')
parser.add_argument('-ncpus', nargs=1,
                    help='specify number of cpus to use')

args = parser.parse_args()

if args.lys:
    aa1 = "LYS"
    aa2 = "LYS"

elif args.xl_path:
    xl_path = args.xl_path[0]

elif args.aa1:
    if args.aa2:
        aa1 = args.aa1[0].upper()
        aa2 = args.aa2[0].upper()
        # catch any dodgy typing
        if aa1 not in amino_acids or aa2 not in amino_acids:
            print("ERROR: Please type amino acid in three letter code format")
            sys.exit(2)
    else:
        print("Please specify both aa1 AND aa2 if you want to use this option")
        sys.exit(2)

if args.max_dist:
    max_dist = int(args.max_dist[0])

if args.i:
    pdb_path = args.i[0]

if args.vox:
    vox = int(args.vox[0])

if args.surface:
    surface = True

if args.ncpus:
    ncpus = int(args.ncpus[0])

#######################################


def runJwalk(pdb_path, max_dist, vox, surface, xl_path='', aa1='', aa2='', ncpus=mp.cpu_count()):
    """
        Execute Jwalk with processed command line options

        Args:
            pdb_path: pdb file to calculate crosslinks on
            max_dist: maximum distance Jwalk will search
            vox: angstoms per voxel in grid
            surface: if True use higher resolution surface method
            xl_path: list of specific crosslinks to calculate
            aa1: starting residue type
            aa2: ending residues type
            ncpus: number of cpus to use
    """
    print("calculating crosslinks on", pdb_path)
    # load pdb into Jwalk
    structure_instance = PDBTools.read_PDB_file(pdb_path)
    # generate grid of voxel size (vox) that encapsulates pdb
    grid = GridTools.makeGrid(structure_instance, vox)

    # mark C-alpha positions on grid
    if xl_path:  # if specific crosslinks need to be calculated
        crosslink_pairs, aa1_CA, aa2_CA = GridTools.mark_CAlphas_pairs(
            grid, structure_instance, xl_path)
        xl_path_basename = os.path.basename(xl_path)[:-4]
    elif aa1 and aa2:
        crosslink_pairs = []  # na if searching every combination between residue types
        aa1_CA, aa2_CA = GridTools.markCAlphas(grid, structure_instance, aa1, aa2)
        xl_path_basename = ""
    else:
        print("ERROR: please specify either amino acid types or crosslink list")
        sys.exit(2)

    if surface:
        # check more rigorously if residues are solvent accessible or not
        aa1_CA = SurfaceTools.check_solvent_accessibility_freesasa(
            pdb_path, aa1_CA, xl_list=xl_path)
        if aa1 != aa2 or xl_path:
            aa2_CA = SurfaceTools.check_solvent_accessibility_freesasa(
                pdb_path, aa2_CA, xl_list=xl_path)
        else:
            aa2_CA = aa1_CA.copy()

    dens_map = GridTools.generate_solvent_accessible_surface(
        grid, structure_instance, aa1_CA, aa2_CA)
    # identify which residues are on the surface
    aa1_voxels, remove_aa1 = GridTools.find_surface_voxels(aa1_CA, dens_map, surface, xl_path)
    aa2_voxels, remove_aa2 = GridTools.find_surface_voxels(aa2_CA, dens_map, surface, xl_path)

    crosslink_pairs = SurfaceTools.update_crosslink_pairs(
        crosslink_pairs, aa1_CA, aa2_CA, remove_aa1, remove_aa2)

    # calculate sasds
    sasds = SASDTools.parallel_BFS(
        aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, crosslink_pairs,
        max_dist, vox, ncpus, xl_path)

    # remove duplicates
    sasds = GridTools.remove_duplicates(sasds)
    sasds = SASDTools.get_euclidean_distances(sasds, pdb_path, aa1, aa2)

    # output sasds to .pdb file and .txt file
    PDBTools.write_sasd_to_txt(sasds, pdb_path, xl_path_basename)
    PDBTools.write_sasd_to_pdb(dens_map, sasds, pdb_path, xl_path_basename)
    print(len(sasds), "SASDs calculated")

    aas_on_path = PDBTools.get_aas_on_path(sasds, structure_instance, dens_map)

    return sasds, aas_on_path


if __name__ == "__main__":
    runJwalk(pdb_path, max_dist, vox, surface, xl_path, aa1, aa2, ncpus)
