import sys
import numpy as np
import pyrosetta as py
import MDAnalysis as mda
import matplotlib.pyplot as plt
from LRT import md_reu_lrt

lrt = md_reu_lrt('AB_CD')
print("\t".join("%.2f" % x for x in lrt.cdr_bb_rmsd(sys.argv[1], sys.argv[2] , sys.argv[3] )))
