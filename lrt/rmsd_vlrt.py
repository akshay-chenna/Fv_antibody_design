import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align
import sys

def residue_displacement(pdb_mdenm,xtc_mdenm,pdb_ref,selection,v_lrt):
    u= mda.Universe(pdb_mdenm,xtc_mdenm)
    v = mda.Universe(pdb_ref,pdb_ref)

    v_lrt = np.load(v_lrt)
    v_lrt = v_lrt/np.linalg.norm(v_lrt)
    u_fvca = u.select_atoms(selection)
    v_fvca = v.select_atoms(selection)

    dp = []
    for _ in u.trajectory:
        _ = align.alignto(u_fvca,v_fvca)
        deformation = (u_fvca.positions - v_fvca.positions).flatten()
        dp.append((np.dot(deformation,v_lrt)**2/u_fvca.n_residues)**0.5)
    return np.array(dp)

z_v = [residue_displacement(str(i)+ '/mdframe.pdb',str(i)+ '/md_corrected.xtc','../15/20/mdframe_solute.pdb','(chainID A or chainID B) and name CA','inputs/v_lrt_15-20.npy') for i in range(5,35,5)]

plt.figure(figsize=(8,6))
plt.rcParams['font.size'] = '16'
for i in z_v:
    plt.plot(i, linewidth='2')
plt.xlabel('Time (ps)')
plt.ylabel('RMSD-v_lrt $(\AA)$')
plt.legend([str(i) for i in range(5,35,5)], ncol=3)
plt.savefig('rmsd_vlrt_15-20.png',dpi=300, facecolor='white',bbox_inches='tight')
