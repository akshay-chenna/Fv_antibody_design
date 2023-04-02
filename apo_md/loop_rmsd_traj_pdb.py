'''
Calculate the RMSD of the loops between a trajectory and a pdb. Supported loops are h/l1-3, heavy,light, all.
Use this when selections are not equivalent.
'''

import sys
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
import pyrosetta as py
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
from pyrosetta.rosetta.protocols import antibody
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
py.init('-out:level 0')

ref_top = sys.argv[1]
ref_traj = sys.argv[2]
mobile_top = sys.argv[3]
mobile_traj = sys.argv[4]
name = sys.argv[5]
index_ref = sys.argv[6]
index_mob = sys.argv[7]
rmsd_val = sys.argv[8]
renumbered_pdb = sys.argv[9]
cdr = sys.argv[10]

pose = py.pose_from_pdb(renumbered_pdb)
ab_info = antibody.AntibodyInfo(pose, antibody.Chothia_Scheme, antibody.North)

if cdr == 'all':
    vec = ab_info.get_all_cdrs()
elif cdr == 'heavy':
    vec = py.rosetta.utility.vector1_bool(6)
    vec[antibody.h1] = True
    vec[antibody.h2] = True
    vec[antibody.h3] = True

elif cdr == 'light':
    vec = py.rosetta.utility.vector1_bool(6)
    vec[antibody.l1] = True
    vec[antibody.l2] = True
    vec[antibody.l3] = True

else:
    vec = py.rosetta.utility.vector1_bool(6)
    vec[getattr(antibody,cdr)] = True

cdr_sel = py.rosetta.protocols.antibody.residue_selector.CDRResidueSelector(ab_info, vec)
where = np.where(cdr_sel.apply(pose))

u = mda.Universe(ref_top,ref_traj)
v = mda.Universe(mobile_top,mobile_traj)

u.transfer_to_memory(step=10,verbose=True)
v.transfer_to_memory(step=10,verbose=True)

if (cdr == 'h1') or (cdr=='h2') or (cdr=='h3') or (cdr=='heavy'):
    v_framework = np.concatenate((v.segments.resindices[0][5:11+1], v.segments.resindices[0][38:45+1],v.segments.resindices[0][80:90+1],v.segments.resindices[0][105:len(pose.chain_sequence(1))]))
    v_loop = v.segments.resindices[0][where]
elif (cdr == 'l1') or (cdr=='l2') or (cdr=='l3') or (cdr=='light'):
    v_framework = np.concatenate((v.segments.resindices[1][6:19+1], v.segments.resindices[1][37:43+1],v.segments.resindices[1][75:84+1],v.segments.resindices[1][100:len(pose.chain_sequence(2))]))
    v_loop = v.segments.resindices[1][where[0]-len(pose.chain_sequence(1))]
elif cdr == 'all':
    framework1 = np.concatenate((v.segments.resindices[0][5:11+1], v.segments.resindices[0][38:45+1],v.segments.resindices[0][80:90+1],v.segments.resindices[0][105:len(pose.chain_sequence(1))]))
    framework2 = np.concatenate((v.segments.resindices[1][6:19+1], v.segments.resindices[1][37:43+1],v.segments.resindices[1][75:84+1],v.segments.resindices[1][100:len(pose.chain_sequence(2))]))
    v_framework = np.concatenate((framework1,framework2))
    
resids = ''
for i in v_framework:
    resids += 'resindex ' + str(i) + ' or '
v_framework_selection = 'not name O and (backbone and (' + resids[:-4] +'))'

resids = ''
for i in v_loop:
    resids += 'resindex ' + str(i) + ' or '
v_loop_selection = 'not name O and (backbone and (' + resids[:-4] +'))'


if (cdr == 'h1') or (cdr=='h2') or (cdr=='h3') or (cdr=='heavy'):
    u_framework = np.concatenate((u.segments.resindices[0][5:11+1], u.segments.resindices[0][38:45+1],u.segments.resindices[0][80:90+1],u.segments.resindices[0][105:len(pose.chain_sequence(1))]))
    u_loop = u.segments.resindices[0][where]
elif (cdr == 'l1') or (cdr=='l2') or (cdr=='l3') or (cdr=='light'):
    u_framework = np.concatenate((u.segments.resindices[1][6:19+1], u.segments.resindices[1][37:43+1],u.segments.resindices[1][75:84+1],u.segments.resindices[1][100:len(pose.chain_sequence(2))]))
    u_loop = u.segments.resindices[1][where[0]-len(pose.chain_sequence(1))]
elif cdr == 'all':
    framework1 = np.concatenate((u.segments.resindices[0][5:11+1], u.segments.resindices[0][38:45+1],u.segments.resindices[0][80:90+1],u.segments.resindices[0][105:len(pose.chain_sequence(1))]))
    framework2 = np.concatenate((u.segments.resindices[1][6:19+1], u.segments.resindices[1][37:43+1],u.segments.resindices[1][75:84+1],u.segments.resindices[1][100:len(pose.chain_sequence(2))]))
    u_framework = np.concatenate((framework1,framework2))

resids = ''
for i in u_framework:
    resids += 'resindex ' + str(i) + ' or '
u_framework_selection = 'not name O and (backbone and (' + resids[:-4] +'))'

resids = ''
for i in u_loop:
    resids += 'resindex ' + str(i) + ' or '
u_loop_selection = 'not name O and (backbone and (' + resids[:-4] +'))'

z = align.AlignTraj(v,u,select=(v_framework_selection,u_framework_selection),verbose=True).run()
frmsd = z.rmsd # Framework RMSD

u_loop_position = u.select_atoms(u_loop_selection).positions
lrmsd = [] #loop rmsd
for _ in v.trajectory:
    lrmsd.append(rms.rmsd(v.select_atoms(v_loop_selection).positions, u_loop_position, center=False, superposition=False))

lrmsd = np.array(lrmsd)
print(name,index_ref,index_mob,cdr,rmsd_val,str(lrmsd.min())[:4],sep='\t')


plt.figure(figsize=(4,3))
plt.rcParams['font.size'] = '12'
plt.plot(lrmsd)
plt.plot(frmsd)
plt.xlabel('#Frame')
plt.ylabel('RMSD ($\AA$)')
plt.legend([cdr,'Framework'])
plt.savefig("{}_nativermsd_{}_{}_{}.png".format(name,index_ref,index_mob,cdr),dpi=300, facecolor='white',bbox_inches='tight')
