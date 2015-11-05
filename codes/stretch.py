import sys
import math
import numpy as np
import MDAnalysis as MDA
from glob import glob
from PDBTools import RemoveBrokenSolvent,IncrementResnumAtPDBHeader,\
        MakeNDXFileFromPDB,AppendIndicesToNDXFile,pdbcat,NDX
from subprocess import call, PIPE, Popen

# USER INPUT PARAMETERS
# top file
# structure file
# template top file
# head selection
# tail selection
# pull_head selection
# pull_tail selection
# head-tail separation distance
# head-tail separation spring constant
# head-tail separation rate
# pre-baller energy threshold (excess energy over starting energy)
# baller spring constant
# baller spring constant shrink rate

# e^{-[k(x-x0)**2]/k_BT}
# k = (k_BT)/(x-x0)**2
# Shrink (x-x0) with time constant

k_BT = 2.494 # kJ/mol nm**2

class MicelleMaker():
    def __init__(self,inputfilepath):
        self.inputparams = LoadInputParams(inputfilepath)
        self.N = int(self.inputparams['n'])

    def LoadInputParams(self,inputfilepath):
        inputparams = {}
        with open(inputfilepath) as inputfile:
            for line in inputfile.readlines():
                try:
                    words = line.split('=')
                    key = words[0].strip().replace('-','').replace('_','').lower()
                    valuestr = words[1].strip()
                    if key in inputparams.keys():
                        sys.exit('Parameter '+key+' must only be specified once')
                    inputparams[key] = valuestr
                except:
                    pass
        return inputparams


def gro_grompp(*args):     call(['grompp']+list(args))
def gro_mdrun(*args):      call(['mdrun']+list(args))
def gro_editconf(*args):   call(['editconf']+list(args))
def gro_genbox(*args):     call(['genbox']+list(args))


def AppendMDASelectionStrToNDXFile(u,sel_str,ndx,label):
    sel = u.selectAtoms(sel_str)
    indices = list(np.array(sel.indices())+1)
    AppendIndicesToNDXFile(ndx,label,indices)

def AppendToTopFromPDB(templatetop,outtop,pdbin):
    recordname = ''
    lastresid = ''
    reslist = []
    resnumlist = []
    with open(pdbin,'r') as pdbfile:
        for line in pdbfile.readlines():
            recordname = line[0:6].strip()
            if recordname == 'ATOM':
                resname = line[17:20].strip()
                resid = line[22:26].strip()
                if lastresid == resid and lastresname == resname:
                    pass
                else:
                    lastresid = resid
                    if reslist and reslist[-1] == resname:
                        resnumlist[-1] = resnumlist[-1] + 1
                    else:
                        reslist.append(resname)
                        resnumlist.append(1)
                lastresname = resname

    linestowrite = []
    with open(templatetop,'r') as topfile:
        linestowrite = topfile.readlines()
    linestowrite.append('\n')

    for i in range(len(reslist)):
        linestowrite.append(reslist[i] + '    '+str(resnumlist[i])+'\n')

    with open(outtop,'w') as topfile:
        topfile.writelines(linestowrite)

def MDPFromTemplate(templatemdp,outmdp,replacements):
    with open(templatemdp,'r') as mdpfile:
        mdptext = mdpfile.read()
    for key in replacements:
        mdptext = mdptext.replace(key,replacements[key])
    with open(outmdp,'w') as mdpfile:
        mdpfile.write(mdptext)

class MicelleMakerGromacs(MicelleMaker):
    molecule_ndx = 'molecule.ndx'
    sphere_pdb = 'sphere.pdb'
    micelle_top = 'micelle.top'
    micelle_ndx = 'micelle.ndx'
    restr_relax_ndx_path = 'restr_relax.ndx'

    def __init__(self,inputfilepath):
        self.inputparams = self.LoadInputParams(inputfilepath)
        self.N = int(self.inputparams['n'])
        self.stretch_pdb = 'stretch.pdb'

    def StretchMoleculeGromacs(self):
        ##### Parameters #######################################
        dt = float(self.inputparams['dt'])
        stretch_k = float(self.inputparams['stretchk'])
        stretch_rate = float(self.inputparams['stretchrate'])
        stretch_target = float(self.inputparams['stretchtarget'])
        pull_head_sel_str = self.inputparams['pullheadselstr']
        pull_tail_sel_str = self.inputparams['pulltailselstr']
        pdb = self.inputparams['pdb']
        top = self.inputparams['top']
        stretch_pdb = self.stretch_pdb
        resname = self.inputparams['resname']
        ########################################################

        u = MDA.Universe(pdb,pdb)
        pull_head = u.selectAtoms(pull_head_sel_str)
        pull_tail = u.selectAtoms(pull_tail_sel_str)
        stretch_target_now = np.linalg.norm(pull_head.centroid()-pull_tail.centroid())

        moleculetop = 'molecule.top'
        AppendToTopFromPDB(top,moleculetop,pdb)

        self.molecule_ndx_path = 'molecule.ndx'
        MakeNDXFileFromPDB(pdb,self.molecule_ndx_path)
        AppendMDASelectionStrToNDXFile(u,pull_head_sel_str,self.molecule_ndx_path,'pull_head')
        AppendMDASelectionStrToNDXFile(u,pull_tail_sel_str,self.molecule_ndx_path,'pull_tail')

        stretch_steps = (stretch_target-stretch_target_now)/(10*stretch_rate*dt)
        if stretch_steps < 0:
            stretch_rate = -stretch_rate
            stretch_steps = -stretch_steps
        replacements = {}
        replacements['DT'] = str(dt)
        replacements['STRETCH_K'] = str(stretch_k)
        replacements['STRETCH_RATE'] = str(stretch_rate)
        replacements['STRETCH_STEPS'] = str(stretch_steps)
        MDPFromTemplate('stretch_template.mdp','stretch.mdp',replacements)

        gro_grompp('-f','stretch.mdp','-c',pdb,'-p',moleculetop,'-o','stretchrun.tpr','-n',self.molecule_ndx_path)
        gro_mdrun('-deffnm','stretchrun','-c',stretch_pdb)
        self.stretch_universe = MDA.Universe(stretch_pdb)
        #call(['rm',stretch_pdb])
            
    def Make(self):
        self.StretchMoleculeGromacs()

inputfilepath = 'inputfile2.txt'
m = MicelleMakerGromacs(inputfilepath)
m.StretchMoleculeGromacs()

#default_dt = 0.04 # 40 fs ... Dry martini
#default_stretch_k = 1247 # kJ/mol nm**2
#default_stretch_rate = 1 # nm/ns??
