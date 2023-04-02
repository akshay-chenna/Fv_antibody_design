import numpy as np
import pyrosetta as py
from pyrosetta.rosetta.protocols import antibody
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
py.init()

pose = py.pose_from_pdb(sys.argv[1])
ab_info = antibody.AntibodyInfo(pose, antibody.Chothia_Scheme, antibody.North)

cdr_selector = CDRResidueSelector(ab_info)
sele = np.array(cdr_selector.apply(pose))
sele = np.where(sele)[0] + 1

np.savetxt(sys.argv[2],sele,fmt='%d', newline=" " , footer='\n', comments='')
