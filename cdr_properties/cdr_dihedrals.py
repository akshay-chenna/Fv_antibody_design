# Get the torsions of the CDR regions. The purpose behind this code was to later use the pyigclassify based clustering method for our set (~1L poses from igfold).

import pyrosetta as py
import numpy as np
from pyrosetta.rosetta.protocols import antibody
from pyrosetta.rosetta.protocols.antibody.residue_selector import CDRResidueSelector
import multiprocessing
import threading
from joblib import Parallel, delayed

'''
Run anarci_.py before running this. To convert a given PDB to chothia numbering. Taken from igfold utilities
See: https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/12.01-RosettaAntibody-Framework-and-SimpleMetrics.ipynb
See: https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/02.05-Protein-Geometry.ipynb
'''

py.init()

class cdr_torsions:
    '''
    The purpose of this class is two fold:
        1) Saves the lengths of CDRs of a given pdb.
        2) Saves the phi and psi dihedral backones.
    '''
    __instance_map = {}
    __lock = multiprocessing.Lock()

    @classmethod
    def instance(cls,pdb):
        with cls.__lock:
            if pdb not in cls.__instance_map:
                cls.__instance_map[pdb] = cls(pdb)
            return cls.__instance_map[pdb]

    def __init__(self, pdb):
        self.pose = py.pose_from_pdb(pdb)
        self.ab_info = antibody.AntibodyInfo(self.pose, antibody.Chothia_Scheme, antibody.North)
        self.pdb = pdb

    def get_length(self):
        np.save(self.pdb+'_cdr_length', np.array([self.ab_info.get_CDR_length(i) for i in self.ab_info.get_all_cdrs()]))
    
    def get_dihedrals(self):
        for i in self.ab_info.get_all_cdrs():
            cdr = py.rosetta.utility.vector1_bool(6)
            cdr[i] = True
            cdr_selector = CDRResidueSelector(self.ab_info, cdr)
            sele = cdr_selector.apply(self.pose)
            phi = []
            psi = []
            for s in range(1,len(sele)):
                if sele[s]:
                    phi.append(self.pose.phi(s))
                    psi.append(self.pose.psi(s))
            np.save(self.pdb+'_'+self.ab_info.get_CDR_name(i)+'_phi',np.array(phi))
            np.save(self.pdb+'_'+self.ab_info.get_CDR_name(i)+'_psi',np.array(psi))

def execute(in_file):
    try:
        x = cdr_torsions(in_file)
        x.get_length()
        x.get_dihedrals()
    except:
        pass

with open('list.txt','r') as f:
    names = [line.strip('\n') for line in f]

Parallel(n_jobs=90, require='sharedmem')(delayed(execute)('renumbered/'+l) for l in names)
