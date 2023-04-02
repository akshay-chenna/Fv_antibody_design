import sys
import MDAnalysis as mda
import numpy as np

u = mda.Universe(sys.argv[1],sys.argv[1])
x1 = np.where( ((u.segments[0].residues.resids >= 6) & (u.segments[0].residues.resids<=22))  | ((u.segments[0].residues.resids >= 39) & (u.segments[0].residues.resids<=46))  |  ((u.segments[0].residues.resids >= 81) & (u.segments[0].residues.resids<=91))| ((u.segments[0].residues.resids >= 106)))
x2 = np.where( ((u.segments[1].residues.resids >= 7) & (u.segments[1].residues.resids<=20))  | ((u.segments[1].residues.resids >= 38) & (u.segments[1].residues.resids<=44))  |  ((u.segments[1].residues.resids >= 76) & (u.segments[1].residues.resids<=85))| ((u.segments[1].residues.resids >= 101)))
fh = u.segments[0].residues.resindices[x1]
fl = u.segments[1].residues.resindices[x2]
f = np.concatenate((fh,fl))
np.savetxt(sys.argv[2],f,fmt='%d',newline='\n')
