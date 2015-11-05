#Still need to:
    #call stretch.py to pull surfactant (inputfile already made below as AOT_pull.csv)
    #put bond, angle, and dihedral definitions for AOT into forcefields file


#Ashley Kiemen
#This script appropriates code from PDBTools.py and pull.py by Kyle Huston
import MDAnalysis as MDA
import importlib
from subprocess import call
from subprocess import Popen, PIPE
import pandas as pd
import os
from PDBTools import pdbcat, RemoveBrokenSolvent, AppendToTopFromPDB, NDX
from decimal import *

#Characterize system
water_oil_system = False
water_air_system = True
surfactant = "AOT_53a6oxy.pdb"
surfactant_smiles = "[O-]S(=O)(=O)C(C(=O)OCC(CC)CCCC)CC(=O)OCC(CC)CCCC"
surfactant_itp_filename = "AOT_53a6oxy.itp"
surfactant_top_filename = "AOT_53a6oxy.top"
surfactant_resname = "AOT"

#oil
oil = "octane.pdb"
number_of_oil_molecules = 600   #INTEGER value from desired oil density in molecules/nm3 * oil slab volume
oil_inp_filename = "octane.inp"

#water
water_resname = "SOL"
water_sel_str = "resname SOL"
water_crop_dist = 5             #in angstroms (for removing water which is too close to the oil slab or surfactant lattice)
number_of_water_beads = 3       #3 or 1 depending on whether or not hydrogens are independent beads

#Box dimensions
oil_slab_length = 0
water_slab_length = 5
air_slab_length = 10
#box_length = sum of slab lengths
box_width = 5                   #in nm
box_height = 5                  #in nm
#surfactant lattice dimensions
number_of_molecules_wide = 5
number_of_molecules_tall = 5

min_crop_z = 0
max_crop_z = 60

#Define input packmol file lines
line1 = "tolerance 2.0\nfiletype pdb\nadd_box_sides 2.0\n\noutput "
line2_oil = "lat+oil.pdb\n"
line2_air = "lat.pdb\n"
line3 = "structure lattice.pdb\n    center\n    fixed "
line4_oil = str(5*box_width) + " " + str(5*box_height) + " " + str(5*water_slab_length - 2) + " 0 0 0\nend structure\n"
line4_air = str(5*box_width) + " " + str(5*box_height) + " " + str(5*air_slab_length - 2) + " 3.145926 0 0\nend structure\n"
line5 = "structure lattice.pdb\n    center\n    fixed "
line6_oil = str(5*box_width) + " " + str(5*box_height) + " " + str(5*water_slab_length + 10*oil_slab_length + 2) + " 3.145926 0 0\nend structure\n"
line6_air = str(5*box_width) + " " + str(5*box_height) + " " + str(5*air_slab_length + 10*water_slab_length + 2) + " 0 0 0\nend structure\n"
p = (5*water_slab_length) + (10*oil_slab_length)
line7 = "structure " + oil + "\n    number " +str(number_of_oil_molecules) + "\n     inside box 0.5 0.5 " + str(5*water_slab_length) + " " + str(10*box_width) + " " + str(10*box_height) + " " + str(p) + "\nend structure"

#Define template topology file lines
a = open("template.top", "w")
if water_air_system == True:
    a.write('; Include forcefield parameters \n#include ' + '"' + 'gromos53a6oxy+D_furan.ff/forcefield.itp' + '"' + '\n#include ' + '"' + 'ions.itp' + '"' '\n#include ' + '"' + surfactant_itp_filename + '"' + '\n')
if water_oil_system == True:
    a.write('; Include forcefield parameters\n#include ' + '"' + 'ffG53a6.itp' '"' + '\n#include ' + '"' + '53a6oxy_O.itp' + '"' + '\n#include ' + '"' + 'ions.itp' + '"' + '\n#include ' + '"' + surfactant_inp_filename + '"' + '\n#include ' + '"' + oil_inp_filename + '"' + '\n')
a.write('; Include water topology\n#ifdef POSRES\n#include ' + '"' + 'headgr_posre.itp' + '"' + '\n#endif\n#ifdef FLEXIBLE\n#include ' + '"' + 'flexspc.itp' + '"' + '\n#else\n#include ' '"' + 'spc.itp' + '"' + '\n#endif\n[ system ]\n; Name\nThree Isomers in 0% NaCl in water\n\n[ molecules ]\n; Compound        #mols')
a.close()

#stretch molecule and orient so that polar end points towards positive z
surfactant_pull = surfactant
#determine pull length from smiles
i = 0
smiles = ""
while i < len(surfactant_smiles):
    while  surfactant_smiles[i] == "(" or surfactant_smiles[i] == "[":
        while surfactant_smiles[i] != ")" and surfactant_smiles[i] != "]":
            i = i + 1
            if surfactant_smiles[i] == "(":
                while surfactant_smiles[i] != ")":
                    i = i + 1
                if surfactant_smiles[i] == ")":
                    i = i + 1
        i = i + 1
    if surfactant_smiles[i] != "=" and surfactant_smiles[i] != "+" and surfactant_smiles[i] != "-" and surfactant_smiles[i] != "." and surfactant_smiles[i] != ",":
        smiles = smiles + surfactant_smiles[i]
    i = i + 1
pull_length = Decimal(len(smiles)) + Decimal(1.05)

#determine head and tail pull atoms
a = open(surfactant, 'r')
atom_names = ""
for line in a:
    words = line.split()
    if words[0] == "ATOM":
        atom_names = atom_names + words[2] + " "
a.close()
words = atom_names.split()
head_1 = words[0]
head_2 = words[1]
head_3 = words[2]
tail_1 = words[-1]
tail_2 = words[-2]
tail_3 = words[-3]

#build input file
linein1 = "top                 = " + surfactant_top_filename
linein2 = "\npdb                 = " + surfactant
linein3 = "\nresname             = " + surfactant_resname
linein4 = "\ndt                  = 0.002"
linein5 = "\n\nsurf_sel_str        = resname " + surfactant_resname
linein6 = "\npull_head_sel_str   = name " + head_1 + " or name " + head_2 + " or name " + head_3
linein7 = "\npull_tail_sel_str   = name " + tail_1 + " or name " + tail_2 + " or name " + tail_3
linein8 = "\n\nstretch_k           = 6000\nstretch_rate        = 0.001\nstretch_target      = " + str(pull_length)
p = open(surfactant_resname + "_pull.csv", "w")
p.write(linein1+linein2+linein3+linein4+linein5+linein6+linein7+linein8)
p.close()

#pull molecule


#Build surfactant lattice from stretched molecule
mol_width = box_width/number_of_molecules_wide
mol_height = box_height/number_of_molecules_tall
call(["genbox", "-cp", surfactant_pull, "-box", str(mol_width), str(mol_height), "8", "-o", "surfactant_pull_box.pdb"])
call(["genconf", "-f", "surfactant_pull_box", "-nbox", str(number_of_molecules_wide), str(number_of_molecules_tall), "1", "-o", "lattice.pdb"])

#water-air system setup
if water_air_system == True:
    #Build box
    x = open("lat.inp", "w")
    x.write(line1+line2_air+line3+line4_air+line5+line6_air)
    x.close()
    with open('lat.inp') as inp:
        call(["packmol"], stdin=inp)
    call(["editconf", "-f", "lat.pdb", "-o", "lat_box.pdb", "-box", str(box_width), str(box_height), str(water_slab_length)])
    #Add water and resize box
    call(["genbox", "-cp", "lat_box.pdb", "-cs", "-o", "lat+water_1.pdb"])
    call(["editconf", "-f", "lat+water_1.pdb", "-o", "lat+water.pdb", "-box", str(box_width), str(box_height), str(water_slab_length+air_slab_length)])
    #Create topology and energy minimize
    AppendToTopFromPDB("template.top", "lat+water.top", "lat+water.pdb")
    call(["grompp", "-f", "em.mdp", "-c", "lat+water.pdb", "-p", "lat+water.top", "-o", "lat+water_em.tpr"])
    #mdrun(["mdrun", "-deffnm", "lat+water_em", "-c", "lat+water_em.pdb"])
    #Run simulation


#Water-oil system setup
if water_oil_system == True:
    #generate oil-slab using packmol
    x = open("lat+oil.inp", "w")
    x.write(line1+line2_oil+line3+line4_oil+line5+line6_oil+line7)
    x.close()
    call(["packmol", "<", "lat+oil.inp"])
    call(["editconf", "-f", "lat+oil.pdb", "-o", "lat+oil_box.pdb", "-box", str(box_width), str(box_height), str(oil_slab_length+water_slab_lenth)])    
    #Add water and remove water from center
    call(["genbox", "-cp", "lat+oil_box.pdb", "-cs", "-o", "lat+oil+water_1.pdb"])
    u = MDA.Universe('lat+oil+water_1.pdb')
    not_water = u.selectAtoms('not (%s)'%water_sel_str)
    if oil:
        water_save = u.selectAtoms('(%s) and not around %f not (%s)'%(water_sel_str,water_crop_dist,water_sel_str))
    elif crop:
        water_save = u.selectAtoms('(%s) and (prop z < %f or prop z > %f)'%(water_sel_str,min_crop_z,max_crop_z))
    W = MDA.Writer('lat+oil+water_broken.pdb')
    W.write(not_water+water_save)
    W.close()
    RemoveBrokenSolvent("lat+oil+water_broken.pdb", "lat+oil+water.pdb", water_resname, number_of_water_beads)
    #Create topology and energy minimize
    AppendToTopFromPDB("template.top", "lat+oil+water.top", "lat+oil+water.pdb")
    call(["grompp", "-f", "em.mdp", "-c", "lat+oil+water.pdb", "-p", "lat+oil+water.top", "-o", "lat+oil+water_em.tpr"])
    mdrun(["mdrun", "-deffnm", "lat+water_em", "-c", "lat+water_em.pdb"])
    #Z-squeeze to avoid vacuum space and frozen oil
    call(["grompp", "-f", "z-squeeze.mdp", "-c", "lat+oil+water_em.pdb", "-p", "lat+oil+water.top", "-o", "lat+oil+water_z.tpr"])
    mdrun(["mdrun", "-deffnm", "lat+oil+water_z", "-c", "lat+oil+water_z.pdb"])
    #Run simulation

#pull surfactant
#umbrella sampling
