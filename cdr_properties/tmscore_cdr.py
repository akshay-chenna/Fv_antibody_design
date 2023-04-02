# Computes tmscores between two CDRs.

import pyrosetta as py
import numpy as np
from tmtools import tm_align
from pyrosetta.rosetta.protocols import antibody
from pyrosetta.rosetta.protocols.antibody.residue_selector import CDRResidueSelector
import multiprocessing
import threading
from joblib import Parallel, delayed

'''
Run anarci_.py before running this. To convert a given PDB to chothia numbering.
See: https://pypi.org/project/tmtools/
See: https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/12.01-RosettaAntibody-Framework-and-SimpleMetrics.ipynb
'''

py.init()
class cdr_data:
    '''
    i) Finds the CDR sequence in a given numbered PDB.
    ii) Finds the CA coordinates of the CDR sequences from a PDB.
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

    def cdr_sequence(self):
        return ''.join([self.ab_info.get_CDR_sequence_with_stem(i) for i in self.ab_info.get_all_cdrs()])

    def cdr_CA_xyz(self):
        cdr_selector = CDRResidueSelector(self.ab_info)
        sele = cdr_selector.apply(self.pose)
        c=[]
        for i in range(1, len(sele)):
            if sele[i]:
                c.append([self.pose.residue(i).xyz('CA').x, self.pose.residue(i).xyz('CA').y, self.pose.residue(i).xyz('CA').z])
        return np.array(c)

def compute_tmscore(pdb1,pdb2):
    a = cdr_data(pdb1)
    b = cdr_data(pdb2)
    score = tm_align(a.cdr_CA_xyz(), b.cdr_CA_xyz(), a.cdr_sequence(), b.cdr_sequence())
    with open('cdr_tmscores_1k-1k.txt','a') as f:
        print(pdb1, pdb2, score.tm_norm_chain1, score.tm_norm_chain2, file=f)
    #print(pdb1, pdb2, score.tm_norm_chain1, score.tm_norm_chain2)

### Run the above functions for arbitrary inputs ####
with open('1000_names.txt') as f: # 1000 sequentially distinct pdbs from kmeans clustering on 100k OAS igfold predictions
    names_1000 = [line.strip('\n') for line in f]

with open('10k_random_names.txt') as f: # 10k random names
    names_random10k = [line.strip('\n') for line in f]

Parallel(n_jobs=90, require="sharedmem")(delayed(compute_tmscore)('renumbered/'+x,'renumbered/'+y) for x in names_1000 for y in names_1000)
