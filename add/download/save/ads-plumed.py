# coding: utf-8
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from numpy import abs, tanh
from scipy.optimize import differential_evolution
from gmx_wrapper import Gromacs
from PDBTools import AppendToTopFromPDB, pdbcat
from glob import glob
import os
from shutil import copyfile
import pandas as pd

gmx = Gromacs(verbosity=0)
num_surfs_ads_rem = 0

def TanhProfile(z1,dz,w,rho1,x_edges):
    Lx = x_edges[-1]-x_edges[0]
    dx = x_edges[1]-x_edges[0]
    x = (x_edges[:-1]+x_edges[1:])/2.
    x_ext = np.linspace(x[0]-Lx,x[-1]+Lx,3*len(x))
    profile = rho1/2.*(tanh(2*(x_ext-z1)/w)-tanh(2*(x_ext-z1-abs(dz))/w))
    stack = np.vstack((profile[:len(x)],profile[len(x):-len(x)],profile[-len(x):]))
    return np.amax(stack,axis=0)

def TanhProfileResidual(z1_dz_w_rho1,x,rho_to_fit):
    z1   = z1_dz_w_rho1[0]
    dz   = z1_dz_w_rho1[1]
    w    = z1_dz_w_rho1[2]
    rho1 = z1_dz_w_rho1[3]
    return np.sum((TanhProfile(z1,dz,w,rho1,x)-rho_to_fit)**2.)

slab_sel_str = 'resname 50'
# Angstroms
dstart = 40
douter = 2*dstart

def identify_slab(top,traj,slab_sel_str):
    """
        Return slab edge positions, edge width, and bulk density assuming symmetric, intact slab (not pbc-broken).
        Performs this just for the first frame.
    """
    # Use first frame of trajectory to find interface
    u = mda.Universe(top, traj)
    slab = u.select_atoms(slab_sel_str)
    
    slab.pack_into_box()

    bins = 200
    hist, z_edges = np.histogram(slab.positions[:,2],weights=slab.masses/(u.dimensions[0]*u.dimensions[1]*u.dimensions[2]/bins),
                                 bins=bins,range=(0,u.dimensions[2]))
    z_centers = (z_edges[:-1]+z_edges[1:])/2.

    sol = differential_evolution(TanhProfileResidual,bounds=((0,u.dimensions[2]),(0,u.dimensions[2]),(0,u.dimensions[2]/2.),(0.1,4.9)),args=(z_edges,hist))
    z1f   = sol.values()[3][0]
    dzf   = sol.values()[3][1]
    wf    = sol.values()[3][2]
    rho1f = sol.values()[3][3]
    z2f = z1f + dzf
    return z1f,z2f,wf,rho1f

#plt.plot(z_centers,TanhProfile(z1f,dzf,wf,rho1f,z_centers))
#plt.plot(z_centers,hist)
#plt.show()

class Simulate_Surf(object):
    def __init__(self,surf_pdb,mdp,template_top,surf_sel_str=None):
        """ Generate random surfactant configurations"""
        self.surf_pdb = surf_pdb
        self.mdp = mdp
        self.template_top = template_top
        self.deffnm = 'surfactant'
        self.simulate()
    @property
    def tpr(self):
        return '%s.tpr'%self.deffnm
    @property
    def top(self):
        return '%s.top'%self.deffnm
    @property
    def traj(self):
        return '%s.xtc'%self.deffnm
    @property
    def gro(self):
        return '%s.gro'%self.deffnm

    def simulate(self):
        AppendToTopFromPDB(self.template_top,'%s.top'%(self.deffnm),self.surf_pdb)
        gmx.grompp('-f %s -c %s -p %s -o surfactant.tpr'%(self.mdp,self.surf_pdb,self.top))
        gmx.mdrun_mpi(1, '-deffnm surfactant')

class Surfer(object):
    def __init__(self,surf_pdb,mdp,template_top,surf_sel_str=None):
        """
            surf_pdb -> structure file of surfactant
            surf_sel_str -> (opt.) selection string for surfactant.
                            If not provided, it will be inferred from surf_pdb
        """
        self.deffnm_history = []
        self.interval = 40 # ps between checks
        self.surf_pdb = surf_pdb
        self.mdp = mdp
        self.template_top = template_top
        self.count = 0
        if surf_sel_str is None:
            residues = mda.Universe(surf_pdb).residues
            assert(len(residues) == 1)
            self.surf_sel_str = 'resname %s'%residues.resnames[0]
        else:
            self.surf_sel_str = surf_sel_str

    @property
    def deffnm(self):
        return '{:d}'.format(self.count)
    @property
    def tpr(self):
        return '{:s}.tpr'.format(self.deffnm)
    @property
    def top(self):
        return '{:s}.top'.format(self.deffnm)
    @property
    def traj(self):
        return '{:s}.xtc'.format(self.deffnm)
    @property
    def xtc(self):
        return self.traj
    @property
    def gro(self):
        return '{:s}.gro'.format(self.deffnm)

    def ready_move(self,out_frame=None,x=0.,y=0.,z=0.):
        #gmx.trjcat('-f {xtcs:s} -o combined_{i:d}.xtc'.format(
        #           xtcs=' '.join(['{:d}_{:d}.xtc'.format(self.count2, i) for i in range(1, self.count+1)]),
        #           i=self.count2))

        #copyfile(self.tpr, 'combined_{i:d}.xtc'.format(i=self.count2))
        # load and deal with pbc issues
        #gmx.editconf('-f %s -o %s'%(self.gro,'self.pdb'))
        u = mda.Universe(self.tpr, self.xtc)
        #u = mda.Universe(self.tpr, self.gro)
        #u.trajectory[out_frame]
        if out_frame is not None:
            u.trajectory[out_frame]
        else:
            out_frame = len(u.trajectory)-1
            u.trajectory[out_frame]
        new_molecule = u.select_atoms(self.surf_sel_str).residues[0]
        mda.lib.mdamath.make_whole(new_molecule)
        COM = new_molecule.center_of_mass()
        xyz = np.array((x,y,z))
        new_molecule.positions += (xyz-COM)
        # Ashley: maybe try moving surrouding solvent along with surfactant to avoid leaving gap immediately around solvent
        # Try to re-add intersecting solvent to random positions inside bounding box of where moved surfactant


        u.atoms.pack_into_box()
        COM = u.atoms.center_of_mass()
        Lx, Ly, Lz, _, _, _ = u.dimensions
        xyz = np.array((Lx/2., Ly/2., Lz/2.))
        u.atoms.positions += (xyz-COM)

        with mda.Writer('moved.pdb') as W:
            W.write(u.atoms)
        self.count += 1
        # pre-process simulation
        AppendToTopFromPDB(self.template_top,'%s.top'%self.deffnm,'moved.pdb')
        gmx.grompp('-f {mdp:s} -c {pdb:s} -p {top:s} -o {tpr:s}.tpr'.format(
            mdp=self.mdp,
            pdb='moved.pdb',
            top=self.top,
            tpr=self.deffnm))
        # set self state
        self.mdrun_args = ['-deffnm',self.deffnm, '-plumed', 'plumed.dat']
        self.state = 'start'

    def ready_new(self,x=0.,y=0.,z=0.):
        """Add a new molecule to x, y, z
        """
        # merge old box with new surfactant
        # Delete intersecting solvent
        # Try to re-add intersecting solvent to random positions inside solvent slab
        # write new structure
        u_box = mda.Universe(self.gro)
        u_surf = mda.Universe(self.surf_pdb)
        print('z: {:f}'.format(z))
        print('surf COM: {}'.format(u_surf.atoms.center_of_mass()))
        u_new = mda.Merge(u_surf.atoms, u_box.atoms)
        u_new.dimensions = u_box.dimensions
        print('u dimensions: {}'.format(u_new.dimensions))
        new_molecule = u_new.select_atoms(self.surf_sel_str).residues[0]
        COM = new_molecule.center_of_mass()
        xyz = np.array((x, y, z))
        new_molecule.positions += (xyz-COM)
        print('new_molecule z: {:f}'.format(new_molecule.center_of_mass()[2]))
        u_new.atoms.pack_into_box()
        COM = u_new.atoms.center_of_mass()
        Lx, Ly, Lz, _, _, _ = u_new.dimensions
        xyz = np.array((Lx/2., Ly/2., Lz/2.))
        u_new.atoms.positions += (xyz-COM)
        print('new_molecule z: {:f}'.format(new_molecule.center_of_mass()[2]))
        with mda.Writer('new.pdb') as W:
            W.write(u_new.atoms)

        self.count += 1
        # pre-process simulation
        AppendToTopFromPDB(self.template_top,'%s.top'%self.deffnm,'new.pdb')
        gmx.grompp('-f %s -c %s -p %s -o %s.tpr'%(self.mdp,'new.pdb',self.top,self.deffnm))
        # set self state
        self.mdrun_args = ['-deffnm',self.deffnm, '-plumed', 'plumed.dat']
        self.state = 'start'

    def run(self):
        if self.state:
            mdrun_args = self.mdrun_args
        elif self.state is None:
            raise Exception('No simulation is ready to run')
        else:
            raise Exception('not self.state but self.state is not None')
        gmx.mdrun_mpi(1, *mdrun_args)



#Inputs
#simulates surfactant for 10ns
surfactant_pdb = 'c12e8.pdb'
surfactant = Simulate_Surf(surfactant_pdb,'md1.mdp','template.top')
surfactant.deffnm = 'surfactant'
auto_corr_time = 100
print 'adding surfactant and starting trials'


surfer = Surfer('surfactant.gro','md_plumed.mdp','template.top')
N = 1

z1, z2, w, rho = identify_slab(surfer.gro, surfer.gro, slab_sel_str)

start = z2 + dstart

u = mda.Universe(surfer.gro)
Lz = u.dimensions[2]
if z2 + douter >= Lz:
    raise ValueError('Lz={:f} must be greater than z2+douter ({:f})'.format(Lz, z2+douter))

def write_plumed_dat(z2, douter, Lz):
    with open('plumed.dat', 'w') as f:
        f.write('c_surf: COM ATOMS=1-11\n'
                'p: POSITION ATOM=c_surf\n'
                'PRINT ARG=p.z STRIDE=10 FILE=COLVAR\n'
                'COMMITTOR ...\n'
                '  ARG=p.z\n'
                '  STRIDE=10\n'
                '  BASIN_A_LOWER=-{Lz:f}\n'
                '  BASIN_A_UPPER={lower:f}\n'
                '  BASIN_B_LOWER={upper:f}\n'
                '  BASIN_B_UPPER=0\n'
                '  FILE=COMMITTOR\n'
                '... COMMITTOR \n'.format(
                    lower=(z2-Lz)/10.,
                    upper=(z2 + douter - Lz)/10.,
                    Lz=Lz/10.))

# WARNING: IT IS ASSUMED THAT NOTHING WILL DESORB SPONTANEOUSLY
surfer.ready_new(z=start)
z1,z2,w,rho = identify_slab(surfer.tpr, 'new.pdb', slab_sel_str)
write_plumed_dat(z2, douter, Lz)
surfer.run()
print 'start = ' + str(start)
#print 'z2 = ' + str(z2)

while surfer.count < 100:
    start = z2 + dstart
    outer = z2 + douter
    surfactant_pull_time = auto_corr_time * N
    #gmx.trjconv('-s surfactant.tpr -f %s -o %s -b %s -e %s'%('surfactant.xtc',surfactant_pdb,surfactant_pull_time,surfactant_pull_time), stdin_text='System\n')
    
    state = None

    #Check COLVAR script for z positions
    committor = pd.read_csv('COMMITTOR', comment='#', sep=' ', skipinitialspace=True,
                         names=['t', 'Z'])
    colvar = pd.read_csv('COLVAR', comment='#', sep=' ', skipinitialspace=True,
                         names=['t', 'Z'])
    z_final = committor.Z.values[0]
    if colvar.Z[0] <= (z2-Lz)/10. or colvar.Z[0] >= (outer-Lz)/10.:
        raise Exception('ERROR: Simulation started outside of bounds')
    elif z_final <= z2/10. or z_final >= outer/10.:
        surfer.ready_move(z=start)
        z1,z2,w,rho = identify_slab(surfer.tpr, 'moved.pdb' ,slab_sel_str)
    else:
        raise Exception('ERROR: Simulation stopped but surfactant is not adsorbed or out')

    write_plumed_dat(z2, douter, Lz)
    surfer.run()
    for i in glob('*#*'):
        os.remove(i)
