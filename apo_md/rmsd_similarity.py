import sys
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms

ref_top = sys.argv[1]
ref_traj = sys.argv[2]
mobile_top = sys.argv[3]
mobile_traj = sys.argv[4]
name = sys.argv[5]
index_ref = sys.argv[6]
index_mob = sys.argv[7]
rmsd_val = sys.argv[8]

u = mda.Universe(ref_top,ref_traj)
v = mda.Universe(mobile_top,mobile_traj)

u.transfer_to_memory(step=10,verbose=True)
v.transfer_to_memory(step=10,verbose=True)

prmsd = np.zeros([len(u.trajectory), len(v.trajectory)])

for i in range(0,len(u.trajectory)):
    r = rms.RMSD(v, u, select='backbone', ref_frame=i).run()
    prmsd[i,:] = r.rmsd[:,2]

match = np.where((prmsd<float(rmsd_val)).any(axis=0))[0].shape[0]/len(v.trajectory) #Shouldn't axis=1  ?
print(name,index_ref,index_mob,rmsd_val,str(match)[:4],sep='\t')
