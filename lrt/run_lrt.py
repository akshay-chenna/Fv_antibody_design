import os
import sys
import numpy as np
from LRT import md_reu_lrt

lrt = md_reu_lrt('AB_CD')
os.system("bash relax_mdposes.sh")
lrt.compute_be_reu('(chainID A or chainID B)')
np.save('dG_5-25.npy',lrt.dG)
lrt.compute_forces()
lrt.generate_enm('mdframes/mdframe99_r.pdb',310,12)
lrt.apply_lrt()
np.save('d_lrt_5-25.npy',lrt.d_lrt)
lrt.top_modes()
np.save('O_lrt_5-25.npy',lrt.O_lrt)
for i in range(5,35,5):
    lrt.excite(i,'mdframe.pdb')
    lrt.save_gro('md_corrected.gro',i,cycle=sys.argv[1])
np.save('v_lrt_5-25.npy',lrt.v_extra.flatten())
