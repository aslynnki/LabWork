# coding: utf-8
import pickle
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from glob import glob
#u.select_atoms('resname 44')
#u.select_atoms('resname 44').fragments
#u.select_atoms('resname 44').residues
#u.select_atoms('resname 44').residues.center_of_mass()

#create bins
#length_box = 243
length_bins = 8.1
bins_max = []
bin = 0
while bin < 400:
    bin += length_bins
    bins_max.append(bin)
num_bins = len(bins_max)
bin_count = np.zeros([num_bins,num_bins])
print bins_max

lag = 3

#checks center of mass every 10ps
xtcs = set(glob('*.xtc')) - set(['surfactant.xtc'])
tprs = [x.replace('xtc', 'tpr') for x in xtcs]
for tpr, xtc in zip(tprs, xtcs):
    u = mda.Universe(tpr, xtc)
    target = u.select_atoms('resname 44')

    x = np.array([target.center_of_mass()[2] for ts in u.trajectory])
    from_xs = x[:-lag]
    to_xs = x[lag:]

    for from_x, to_x in zip(from_xs, to_xs):
        from_bin = 0.1
        to_bin = 0.1
        i = 0
        print 'from = ' + str(from_x) 
        print 'to = ' + str(to_x)
        while from_bin == 0.1:
            if from_x < bins_max[i]:
                from_bin = i
            i += 1
        i = 0
        while to_bin == 0.1:    
            if to_x < bins_max[i]:
                to_bin = i
            i += 1

        bin_count[to_bin][from_bin] += 1
    
print bin_count
with open('bin.pck','w') as pck:
    pickle.dump(bin_count, pck)

with open('bin.pck') as pck:
    N = pickle.load(pck)
    rgb = np.dstack((1-np.clip(N/np.max(N),0.,1.),1-np.clip(N/np.max(N),0.,1.),1-np.clip(N/np.max(N),0.,1.)))
    plt.imshow(rgb)
    plt.savefig('example.png')

