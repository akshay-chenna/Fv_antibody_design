import sys
import numpy as np
import pyrosetta as py
import MDAnalysis as mda
import matplotlib.pyplot as plt
from LRT import md_reu_lrt

lrt = md_reu_lrt('AB_CD')
with open('out.txt', 'a') as f:
    print((lrt.vecrmsd_lastframe(sys.argv[1], sys.argv[2] , sys.argv[3] ,sys.argv[4])), file=f)
