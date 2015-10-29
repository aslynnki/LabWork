import subprocess
from subprocess import Popen, PIPE
import pandas as pd
import os
# made by Ashley Kiemen

#Nomenclature
smiles = "[O-]S(=O)(=O)C(C(=O)OCC(CC)CCCC)C(C(=O)OCC(CC)CCCC"
force_field = "53a6oxy"
output_filename = 'AOT_53a6oxy'
three_letter_abbreviation = "AOT"
#for complex molecules:
partial_charges_filename = "AOT_charges.csv"


#generate smiles and pdb files
pdb_filename = output_filename + ".pdb"
smiles_filename = output_filename + ".smi"
x = open(smiles_filename, 'w')
x.write(smiles + '\n')
x.close()
subprocess.call(['babel', '-ismi', smiles_filename, '-opdb', pdb_filename, '--gen3d'])

#if molecule begins or ends in an alcohol, append an H
H1 = 0
H2 = 0
if smiles.find("OC") == 0:
    smiles = "H" + smiles
    H1 = 1
if smiles[-2] == "C" and smiles[-1] == "O":
    smiles = smiles + "H"
    H2 = 1
if len(smiles) > 5:
    if smiles[-6] == "C" and smiles[-5] == "(" and smiles[-4] == "=" and smiles[-3] == "O" and smiles[-2] == ")" and smiles[-1] == "O":
        smiles = smiles + "H"
        H2 = 1
copy_smiles = smiles
Case = 0

#Reconfigure .pdb file and remove nonessential hydrogens
v = open(pdb_filename, 'r')
y = open("new.pdb", 'w')
z = open("bonds.csv", 'w')
oxygens = " "
hydrogens = " "
hydrogen_lines = " "
hydrogen_keep = " "
i = 0
for line in v:
    words = line.split()
    if words[0] == "HETATM":
        if words[2] == "O":
            oxygens = oxygens + words[1] + " "
        if words[2] == "H":
            hydrogens = hydrogens + words[1] + " "
            hydrogen_lines = hydrogen_lines + str(i+1) + " "
    i = i + 1
    if words[0] == "CONECT":
        if oxygens.find(words[1]) != -1:
            a = -1
            b = -1
            if len(words) > 2:
                a = hydrogens.find(" " + words[2] + " ")
            if len(words) > 3:
                b = hydrogens.find(" " + words[3] + " ")
            if a != -1 or b != -1:
                if a != -1:
                    hydrogen_keep = hydrogen_keep + words[2] + " "
                if b != -1:
                    hydrogen_keep = hydrogen_keep + words[3] + " "
v.close()
v = open(pdb_filename, 'r')
i = 0
j = 0
h_count = len(hydrogen_keep.split())
h_replace = ["" for x in range(h_count)]
h_replace_count = 0
h_replace_fix = 0
for line in v:
    k = 1
    i = i + 1
    line = line.replace(" LIG", three_letter_abbreviation)
    line = line.replace("HETATM", "ATOM  ")
    if hydrogen_lines.find(" " + str(i) + " ") == -1:
        j = j + 1
        if line.find("CONECT") == -1:
            if line.find("MASTER") != -1:
                line = "TER\n"
            y.write(line)
    if hydrogen_lines.find(" " + str(i) + " ") != -1:
        if hydrogen_keep.find(" " + str(i-2) + " ") != -1:
            line = line.replace(str(i-2), str(j-1), 1)
            h_replace[h_replace_count] = str(j-1)
            h_replace_count = h_replace_count + 1
            j = j + 1
            y.write(line)
    if line.find("CONECT") != -1:
        words = line.split()
        if hydrogens.find(" " + words[1] + " ") == -1:
            while k < len(words) - 1:
                k = k + 1
                if hydrogens.find(" " + words[k] + " ") != -1 and hydrogen_keep.find(" " + words[k] + " ") == -1:
                    line = line.replace(" " + words[k] + " ", " ")
                if hydrogen_keep.find(" " + words[k] + " ") != -1:
                    line = line.replace(" " + words[k] + " ", " " + h_replace[h_replace_fix] + " ")
                    h_replace_fix = h_replace_fix + 1
            z.write(line)
z.close()
v.close()
y.close()

#Assign atom types for forcefield 53a6oxy
atom_types = ["NULL" for x in range(len(smiles))]
while smiles.find("[O-]S(=O)(=O)") != -1:
    num = smiles.find("[O-]S(=O)(=O)")
    smiles = smiles.replace("[O-]S(=O)(=O)","GXGGXGGXGGGXG", 1)
    atom_types[num+1] = " OM"
    atom_types[num+4] = "  S"
    atom_types[num+7] = " OM"
    atom_types[num+11] = " OM"
    Case = 4
while smiles.find("C(=O)O") != -1:
    num = smiles.find("C(=O)O") 
    smiles = smiles.replace("C(=O)O", "XGGXGX", 1)
    atom_types[num] = "  C"
    atom_types[num+3] = "  O"
    atom_types[num+5] = " OM"
    if Case < 3:
        Case = 3
while smiles.find("COC") != -1:
    num = smiles.find("COC")
    smiles = smiles.replace("COC", "XXX", 1)
    atom_types[num] = "CH2"
    atom_types[num+1] = "OE2"
    atom_types[num+2] = "CH2"
    if Case < 2:
        Case = 2
while smiles.find("C(=O)") != -1:
    num = smiles.find("C(=O)")
    smiles = smiles.replace("C(=O)", "XGGXG", 1)
    atom_types[num] = "  C"
    atom_types[num+3] = "  O"
    if Case < 3:
        Case = 3
if smiles.find("H") == 0:
    smiles = smiles.replace("H", "X", 1)
    atom_types[0] = "  H"
    if Case < 2:
        Case = 2
    if smiles.find("O") == 1:
        smiles = smiles.replace("O", "X", 1)
        atom_types[1] = " OA"
if smiles.find("H") == (len(smiles) - 1):
    smiles = smiles.replace("H", "X", 1)
    atom_types[-1] = "  H"
    if Case < 2:
        Case = 2
    if smiles.find("O") == (len(smiles) - 2):
        smiles = smiles.replace("O", "X", 1)
        atom_types[-2] = " OA"
while smiles.find("C") != -1:
    num = smiles.find("C")
    smiles = smiles.replace("C", "X", 1)
    if num == 0 or num == (len(smiles) - 1):
        atom_types[num] = "CH3"
    else:
        atom_types[num] = "CH2"
while smiles.find("[Na+].") != -1:
    num = smiles.find("[Na+].")
    smiles = smiles.replace("[Na+].", "GXGGGG", 1)
    atom_types[num] = "NA+"
while smiles.find("[Cl-].") != -1:
    num = smiles.find("[Cl-].")
    smiles = smiles.replace("[Cl-].", "GXGGGG", 1)
    atom_types[num] = "CL-"
num1 = smiles.find("(")
while num1 != -1:
    smiles = smiles.replace("(", "G")
    smiles = smiles.replace(")", "G")
    num1 = smiles.find("(")
i = 0
atom_number = 0
while i < len(smiles):
    j = smiles.find("X")
    smiles = smiles.replace("X", "Y", 1)
    if j != -1:
        atom_number = atom_number + 1
    i = i + 1
atom_types_final = ["" for x in range(atom_number)]
i = 0
j = 0
while i < len(smiles):
    if atom_types[i] != "NULL":
        atom_types_final[j] = atom_types[i]
        j = j + 1
    i = i + 1

#Assign atomic mass
if force_field == "53a6oxy":
    atom_data = pd.read_csv('53a6oxy_masses.csv',keep_default_na=False,na_values='')
    atom_data = atom_data.drop_duplicates().set_index('type',drop=False)

#miscellaneous itp arrays
nr = ["0" for x in range(len(atom_types_final))]
cgnr = nr
resnr = ["1" for x in range(len(atom_types_final))]
atom = ["0" for x in range(len(atom_types_final))]
resid = [three_letter_abbreviation for x in range(len(atom_types_final))]
i = 0
while i < len(atom_types_final):
    nr[i] = str(i + 1)
    atom[i] = "A" + str(i + 1)
    i = i + 1

#rename atoms in pdb file
x = open("new.pdb", 'r')
y = open("new2.pdb", 'w')
i = -1
for line in x:
    words = line.split()
    if i > 0:
        words = line.split()
        if len(words) > 2:
            a = "A" + str(i)
            if len(words[2]) == 2:
                a = " " + a
            if i < 10:
                a = " " + a
            line = line.replace(words[2] + " ", a, 1)
    y.write(line)
    i = i + 1
x.close()
y.close()


#Assign partial charges
smiles = copy_smiles
#Case 1: only carbon
atomic_charges = ["0.00" for x in range(len(atom_types_final))]
bonds = ["0" for x in range(len(atom_types_final) - 1)]
#Case 2: only carbon and oxygen - no double bonds
if Case == 2:
    while smiles.find("COC") != -1:
        num = smiles.find("COC")
        smiles = smiles.replace("COC", "XXX", 1)
        atomic_charges[num] = " 0.29"
        atomic_charges[num+1] = "-0.58"
        atomic_charges[num+2] = " 0.29"
    while smiles.find("COH") != -1:
        num = smiles.find("COH")
        smiles = smiles.replace("COH", "XXX", 1)
        atomic_charges[num] = " 0.29"
        atomic_charges[num+1] = "-0.70"
        atomic_charges[num+2] = " 0.41"
    while smiles.find("HOC") != -1:
        num = smiles.find("HOC")
        smiles = smiles.replace("HOC", "XXX", 1)
        atomic_charges[num] = " 0.41"
        atomic_charges[num+1] = "-0.70"
        atomic_charges[num+2] = " 0.29"
    while smiles.find("C") != -1:
        num = smiles.find("C")
        smiles = smiles.replace("C", "X", 1)
        atomic_charges[num] = " 0.00"
#Case 3: carbon and oxygen with double bonds
if Case == 3:
    atomic_charges = ["x" for x in range(len(atom_types))]
    while smiles.find("C(=O)OH") != -1:
        num = smiles.find("C(=O)OH")
        smiles = smiles.replace("C(=O)OH", "XXXXXXX", 1)
        atomic_charges[num] = " 0.63"
        atomic_charges[num+3] = "-0.55"
        atomic_charges[num+5] = "-0.37"
        atomic_charges[num+6] = " 0.29"
    while smiles.find("HOC(=O)") != -1:
        num = smiles.find("HOC(=O)")
        smiles = smiles.replace("HOC(=O)", "XXXXXXX", 1)
        atomic_charges[num] = " 0.29"
        atomic_charges[num+1] = "-0.37"
        atomic_charges[num+2] = "-0.55"
        atomic_charges[num+5] = " 0.63"
    while smiles.find("C(=O)C") != -1:
        num = smiles.find("C(=O)C")
        smiles = smiles.replace("C(=O)C", "XXXXXX", 1)
        atomic_charges[num] = "0.29"
        atomic_charges[num+3] = "-0.58"
        atomic_charges[num+5] = "0.29"
    while smiles.find("COC") != -1:
        num = smiles.find("COC")
        smiles = smiles.replace("COC", "XXX", 1)
        atomic_charges[num] = " 0.29"
        atomic_charges[num+1] = "-0.58"
        atomic_charges[num+2] = " 0.29"
    while smiles.find("COH") != -1:
        num = smiles.find("COH")
        smiles = smiles.replace("COH", "XXX", 1)
        atomic_charges[num] = " 0.29"
        atomic_charges[num+1] = "-0.70"
        atomic_charges[num+2] = " 0.41"
    while smiles.find("HOC") != -1:
        num = smiles.find("HOC")
        smiles = smiles.replace("HOC", "XXX", 1)
        atomic_charges[num] = " 0.41"
        atomic_charges[num+1] = "-0.70"
        atomic_charges[num+2] = " 0.29"
    while smiles.find("C") != -1:
        num = smiles.find("C")
        smiles = smiles.replace("C", "X", 1)
        atomic_charges[num] = " 0.00"
    i = 0
    j = 0
    atomic_charges_final = ["x" for x in range(len(atom_types_final))]
    while i < len(smiles):
        if atomic_charges[i] != "x":
            atomic_charges_final[j] = atomic_charges[i]
            j = j + 1
        i = i + 1
    atomic_charges = atomic_charges_final
#Case4: complex molecule
if Case == 4:
    try:
        x = open(partial_charges_filename, 'r')
        i = 0
        while i < len(atom_types_final):
            atomic_charges[i] = x.readline()
            atomic_charges[i] = atomic_charges[i].replace("\n", "")
            i = i + 1
    except:
        print "\n\nFile containing partial charges not found.\n\nPlease make a .csv file containing atomic partial charges for the smiles string (including hydrogen charge for alcohols).\n\n Example file for CCO (partial charges listed in smiles order C, C, O, H):\n\n 0.00\n 0.29\n-0.70\n 0.41\n\n"
        raise SystemExit
#Create atomic mass array
atomic_masses = [atom_data['mass'][x.strip().upper()] for x in atom_types_final]

#Get bond sequences from .pdb file
smi = open("bonds.csv", 'r')
bond_list = ""
for line in smi:
    if line.find("CONECT") != -1:
        words = line.split()
        if len(words) >= 3:
            i = 0
            if len(words) == 3 and words[2] < words[1]:
                i = 4
            while i < len(words) - 1:
                i = i + 1
                if int(words[i]) >= int(words[1]):
                    if i == 2:
                        bond_list = bond_list + "A" + words[1] 
                        if int(words[1]) < 10:
                            bond_list = bond_list + " "
                        bond_list = bond_list + "      A" + words[2] + "\n"
                    if i > 2:
                        bond_list = bond_list + "A" + words[1]
                        if int(words[1]) < 10:
                            bond_list = bond_list + " "
                        bond_list = bond_list + "      A" + words[i] + "\n"
os.rename('new2.pdb', pdb_filename)
smi.close()

#Get improper dihedral sequences from .pdb file
smi = open("bonds.csv", 'r')
dihedral_list = ""
for line in smi:
    if line.find("CONECT") != -1:
        words = line.split()
        if len(words) == 5:
            i = 1
            while i < len(words):
                if int(words[i]) < 10:
                    words[i] = words[i] + " "
                i = i + 1
            dihedral_list = dihedral_list + "A" + words[1] + "      A" + words[2] + "      A" + words[3] + "      A" + words[4] + "\n"

#create .rtp file
rtp_filename = output_filename + ".rtp"
a = open(rtp_filename, 'w')
a.write("[ " + three_letter_abbreviation + " ]\n[ atoms ]\n")
a.write(";AtmName  AtmType  Charge  Count\n")
i = 0
while i < len(atom_types_final):
    if i < 9:
        a.write(" ")
    a.write(atom[i] + "      " + atom_types_final[i] + "      " + str(atomic_charges[i]) + "    ")
    if i < 9:
        a.write(" ")
    a.write(str(int(cgnr[i]) - 1) + "\n")
    i = i + 1
a.write("\n[ bonds ]\n;AtmName AtmName\n")
a.write(bond_list)
a.write("\n\n[ impropers ]\n;AtmName AtmName AtmName AtmName\n")
a.write(dihedral_list)
a.close()
smi.close()
p = subprocess.call(['mv', rtp_filename, 'gromos53a6oxy+D_furan.ff'])

#Create top and itp files
topology_filename_1 = output_filename + "1.top"
topology_filename = output_filename + ".top"
itp_filename = output_filename + ".itp"
p = Popen(["pdb2gmx", "-f", pdb_filename, "-p", topology_filename_1], stdin=PIPE)
p.communicate("1\n1\n")
a = open(topology_filename_1, 'r')
b = open(itp_filename, 'w')
for line in a:
    if line.find("Include Position restraint file") == -1:
        line = line.replace("UNNAMED", output_filename, 1)
        line = line.replace("Other", three_letter_abbreviation, 1)
        b.write(line)
    else:
        b.close()
        b = open(topology_filename, 'w')
        b.write(line)
        b.write("#include " + '"' + itp_filename + '"')
a.close()
b.close()
os.remove("bonds.csv")
os.remove("new.pdb")
os.remove(smiles_filename)
os.remove("conf.gro")
