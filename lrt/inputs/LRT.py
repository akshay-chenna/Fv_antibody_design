import glob
import numpy as np
import pandas as pd
import pyrosetta as py
import MDAnalysis as mda
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.optimize import curve_fit
from MDAnalysis.analysis import align
from pyrosetta.rosetta.core import select
from scipy.linalg import fractional_matrix_power
from pyrosetta.rosetta.protocols import antibody
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, move_map_action

py.init('-out:level 0 -ex1 -ex2')
#py.init('-relax:default_repeats 1 -out:level 0 -ex1 -ex2')
#py.init('-relax:default_repeats 1')

class pylrt:
    '''
    Perform LRT excitation and relaxation purely on CPUs.
    '''

    def __init__(self, pdb, selection):
        '''
        Initializes the pylrt object with a PDB file and an atom selection string.
        '''
        self.name = pdb[:-4]
        self.selection = selection
        self.sfxn = py.get_score_function(True) 
        self.pose = py.pose_from_pdb(pdb)
        self.intf_mover = py.rosetta.protocols.analysis.InterfaceAnalyzerMover('LH_VW') ##Manual entry

    def checkpoint_pose(self):
        self.cpt_pose = self.pose.clone()
        
    def assign_pose(self):
        self.pose = self.cpt_pose.clone()
     
    def tor_relax(self,max_iters):
        '''
        Performs full relaxation on the protein.
        '''
        #self.start_pose = self.pose.clone()
        relax = FastRelax()
        relax.set_scorefxn(self.sfxn)
        relax.max_iter(max_iters)
        relax.apply(self.pose)
        self.relaxed_name = self.name+'_relax.pdb' 
        self.pose.dump_pdb(self.relaxed_name)
    
    def tor_intf_relax(self,max_iters):
        '''
        Performs full relaxation on the protein.
        '''
        #self.start_pose = self.pose.clone()
        relax = FastRelax(scorefxn_in=py.create_score_function("ref2015.wts"),standard_repeats=1,script_file="InterfaceRelax2019")
        relax.max_iter(max_iters)
        relax.apply(self.pose)
        self.relaxed_name = self.name+'_relax.pdb' 
        self.pose.dump_pdb(self.relaxed_name)
    
    def cart_intf_relax(self,max_iters):
        '''
        Relaxes in cartesian space
        '''
        relax = FastRelax(scorefxn_in=py.create_score_function("ref2015_cart.wts"),standard_repeats=1,script_file="InterfaceRelax2019")
        relax.cartesian(True) # Activate InterfaceRelax2019.dualspace
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
        relax.max_iter(max_iters)
        relax.apply(self.pose)
        self.relaxed_name = self.name+'_deformed_0.pdb' 
        self.pose.dump_pdb(self.relaxed_name)
    
    def tor_shortrelax(self,max_iters):
        '''
        Performs few iterations of relaxation on the protein in the torsional space.
        '''
        #self.start_pose = self.pose.clone()
        relax = FastRelax()
        relax.set_scorefxn(self.sfxn)
        relax.max_iter(max_iters)
        relax.apply(self.pose)
        self.relaxed_name = self.name+'_relax_temp.pdb'
        self.pose.dump_pdb(self.relaxed_name)
        
    def tor_intf_shortrelax(self,max_iters):
        '''
        Performs few iterations of relaxation on the protein using InterfaceRelax script.
        '''
        #self.start_pose = self.pose.clone()
        relax = FastRelax(scorefxn_in=py.create_score_function("ref2015.wts"),standard_repeats=1,script_file="InterfaceRelax2019")
        relax.max_iter(max_iters)
        relax.apply(self.pose)
        self.relaxed_name = self.name+'_relax_temp.pdb'
        self.pose.dump_pdb(self.relaxed_name)

    def resttor_intf_shortrelax(self, max_iters):
        '''
        Performs few iterations of restricted relaxation on the protein with partial packing and custom movemaps.
        '''
        #self.pose = self.cpt_pose.clone()
        relax = FastRelax(scorefxn_in=py.create_score_function("ref2015.wts"),standard_repeats=1,script_file="InterfaceRelax2019")
        
        ab_chains = select.residue_selector.OR_combine(select.residue_selector.ChainSelector('L'), select.residue_selector.ChainSelector('H'))
        nbr = select.residue_selector.NeighborhoodResidueSelector()
        nbr.set_focus_selector(ab_chains)
        nbr.set_include_focus_in_subset(True)

        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        tf.push_back(operation.RestrictToRepacking())
        tf.push_back(operation.IncludeCurrent())
        tf.push_back(operation.NoRepackDisulfides())
        prevent_repacking_rlt = operation.PreventRepackingRLT()
        prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr, flip_subset=True)
        tf.push_back(prevent_subset_repacking)

        mmf = MoveMapFactory()
        mmf.add_bb_action(move_map_action.mm_enable, nbr) #Enable whatever you want
        #mmf.add_chi_action(move_map_action.mm_enable, ab_chains)
        mm = mmf.create_movemap_from_pose(self.pose)

        relax.set_movemap(mm)
        relax.set_task_factory(tf)
        
        relax.max_iter(max_iters)
        relax.apply(self.pose)
        
        self.relaxed_name = self.name+'_relax_temp.pdb'
        self.pose.dump_pdb(self.relaxed_name)

    def restcart_intf_shortrelax(self, max_iters):
        '''
        Performs few iterations of restricted relaxation on the protein in cartesian space with partial packing and custom movemaps.
        '''
        #self.pose = self.cpt_pose.clone()
        relax = FastRelax(scorefxn_in=py.create_score_function("ref2015_cart.wts"),standard_repeats=1,script_file="CustomRelax2022")
        
        ab_chains = select.residue_selector.OR_combine(select.residue_selector.ChainSelector('L'), select.residue_selector.ChainSelector('H'))
        nbr = select.residue_selector.NeighborhoodResidueSelector()
        nbr.set_focus_selector(ab_chains)
        nbr.set_include_focus_in_subset(True)

        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        tf.push_back(operation.RestrictToRepacking())
        tf.push_back(operation.IncludeCurrent())
        tf.push_back(operation.NoRepackDisulfides())
        prevent_repacking_rlt = operation.PreventRepackingRLT()
        prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr, flip_subset=True)
        tf.push_back(prevent_subset_repacking)

        mmf = MoveMapFactory()
        '''
        #First disable everything, by default these are false.
        mmf.all_bb(setting=False) 
        mmf.all_bondangles(setting=False)
        mmf.all_bondlengths(setting=False)
        mmf.all_chi(setting=False)
        mmf.all_jumps(setting=False)
        mmf.set_cartesian(setting=False)
        '''
        mmf.add_bb_action(move_map_action.mm_enable, nbr) #Enable whatever you want
        #mmf.add_chi_action(move_map_action.mm_enable, ab_chains)
        mm = mmf.create_movemap_from_pose(self.pose)

        relax.set_movemap(mm)
        relax.set_task_factory(tf)
        
        relax.max_iter(max_iters)
        relax.apply(self.pose)
        
        self.relaxed_name = self.name+'_relax_temp.pdb'
        self.pose.dump_pdb(self.relaxed_name)
        
    def restcstcart_intf_shortrelax(self, max_iters):
        '''
        Performs few iterations of restricted relaxation on the protein in cartesian space with partial packing and custom movemaps and position restraints.
        No constraints applied on the first iteration, but applied from the second onwards.
        '''
        #self.pose = self.cpt_pose.clone()        
        relax = py.rosetta.protocols.relax.FastRelax(scorefxn_in=py.create_score_function("ref2015_cart.wts"), standard_repeats=1,script_file="ConstrainedRelax2023")
        
        sfxn = py.create_score_function("ref2015_cart.wts")
        sfxn.set_weight(py.rosetta.core.scoring.ScoreType.coordinate_constraint, 1) 
        relax.set_scorefxn(sfxn)
        
        ab_chains = select.residue_selector.OR_combine(select.residue_selector.ChainSelector('L'), select.residue_selector.ChainSelector('H'))
        nbr = select.residue_selector.NeighborhoodResidueSelector()
        nbr.set_focus_selector(ab_chains)
        nbr.set_include_focus_in_subset(True)

        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        tf.push_back(operation.RestrictToRepacking())
        tf.push_back(operation.IncludeCurrent())
        tf.push_back(operation.NoRepackDisulfides())
        prevent_repacking_rlt = operation.PreventRepackingRLT()
        prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr, flip_subset=True)
        tf.push_back(prevent_subset_repacking)

        mmf = MoveMapFactory()
        mmf.add_bb_action(move_map_action.mm_enable, nbr) #Enable whatever you want
        mm = mmf.create_movemap_from_pose(self.pose)

        relax.set_movemap(mm)
        relax.set_task_factory(tf)
        
        relax.max_iter(max_iters)
        relax.apply(self.pose)
        
        self.relaxed_name = self.name + '_relax_temp.pdb'
        self.pose.dump_pdb(self.relaxed_name)
        
    def binding_energy(self):
        '''
        Calculates the per-residue binding energy of the protein.
        '''
        self.intf_mover.apply(self.pose)
        return np.array(self.intf_mover.get_all_per_residue_data().dG)
    
    def compute_positions_mda(self, name):
        '''
        Computes the relative positions of atoms in the antibody with respect to the centre of mass of the antigen's interface, using MDanalysis.
        '''
        u = mda.Universe(name,name)
        antibody = u.select_atoms(self.selection)
        self.ab_len = antibody.n_residues
        ab_ca = u.select_atoms('(' + self.selection + ') and (name CA)')
        an_interface = u.select_atoms('(around 5 ((' + self.selection + ') and not (name H* or name [1-9]H*)))')

        relative_positions = ab_ca.positions - an_interface.center_of_mass()
        return relative_positions
    
    def compute_positions_py(self):
        '''
        Computes the relative positions of atoms in the antibody with respect to the centre of geometry of the antigen's interface, using Pyrosetta.
        '''
        ab_chains = select.residue_selector.OR_combine(select.residue_selector.ChainSelector('L'), select.residue_selector.ChainSelector('H'))
        ab = ab_chains.apply(self.pose)
        ab_ca = []
        for i in range(1,len(ab)+1):
            if ab[i]:
                ab_ca.append([self.pose.residue(i).xyz('CA').x, self.pose.residue(i).xyz('CA').y, self.pose.residue(i).xyz('CA').z])
        ab_ca = np.array(ab_ca)
        
        nbr = select.residue_selector.NeighborhoodResidueSelector()
        nbr.set_focus_selector(ab_chains)
        nbr.set_include_focus_in_subset(False)
        interface = nbr.apply(self.pose)
        centroid = []
        for i in range(1,len(interface)+1):
            if interface[i]:
                centroid.append([self.pose.residue(i).xyz('CA').x, self.pose.residue(i).xyz('CA').y, self.pose.residue(i).xyz('CA').z])
        centroid = np.mean(centroid,axis=0)
        return ab_ca - centroid
        
    
    def quadratic_hypersurface(self,P, a, b, c, d, e, f, g, h,i, j):
        '''
        Fitting function for calculating the gradients of binding energies
        '''
        x, y, z = P
        return  a + b*x + c*y + d*z + e*x**2 + f*y**2 + g*z**2 + h*x*y + i*y*z + j*z*x

    def func_forces(self, pos_arr, deltaG_arr, mean_pos_arr):
        '''
        Forces on each residue "i"
        '''
        coeff, _ = curve_fit(self.quadratic_hypersurface, pos_arr.T, deltaG_arr)
    
        def dgx(P):
            x, y, z = P
            return -1*(coeff[1]+2*coeff[4]*x + coeff[7]*y + coeff[9]*z)
        def dgy(P):
            x, y, z = P
            return -1*(coeff[2]+2*coeff[5]*y + coeff[7]*x + coeff[8]*z)
        def dgz(P):
            x, y, z = P
            return -1*(coeff[3]+2*coeff[6]*z + coeff[8]*y + coeff[9]*x)  
        return np.stack((dgx(mean_pos_arr),dgy(mean_pos_arr),dgz(mean_pos_arr)))

    def generate_trajectory(self, max_iters, frames):
        '''
        Generates a relaxation trajectory of the given number of steps.
        '''
        self.relative_pos = []
        self.dG = []
        for i in range(frames):
            self.restcstcart_intf_shortrelax(max_iters)
            self.relative_pos.append(self.compute_positions_mda(self.relaxed_name))
            self.dG.append(self.binding_energy())
            print('Frame '+str(i+1)+' generated')
            #self.report_energies()
        self.relative_pos = np.array(self.relative_pos)
        self.dG = np.array(self.dG)
        self.dG = self.dG * 4.184 # cal to Joule
        self.mean_pos = np.mean(self.relative_pos, axis=0)
        
    def compute_forces(self):
        '''
        Generates an array of forces for all residues
        '''
        self.forces = np.array([self.func_forces(self.relative_pos[:,i,:], self.dG[:,i], self.mean_pos[i]) for i in range(self.ab_len)]) 

    def generate_enm(self,temperature,cutoff):
        '''
        Generates an elastic network model including the environment.
        Calculates the eigenvectors and eigenvalues from the mass-weighted Hessian.
        Calculates the covariance matrix.
        '''
        u = mda.Universe(self.relaxed_name,self.relaxed_name)
        antibody_CA = u.select_atoms('(' + self.selection + ') and (name CA)')
        n_CA = self.ab_len
        complex_CA = u.select_atoms('name CA')
        n_res = complex_CA.n_residues
        positions = complex_CA.positions

        distance_map = distance.cdist(positions,positions)
        distance_map = np.repeat(np.repeat(distance_map,3,axis=1),3,axis=0)
        inv_distsq_map = distance_map**-2
        inv_distsq_map[inv_distsq_map == np.inf] = 0

        kirchoff_matrix = np.ones((n_res*3,n_res*3))
        kirchoff_matrix[distance_map>cutoff] = 0
        spring_constant = 1
        kirchoff_matrix = spring_constant*kirchoff_matrix

        xx_distance_map = distance.cdist(positions[:,0].reshape(-1,1),positions[:,0].reshape(-1,1), lambda u, v: u-v)
        yy_distance_map = distance.cdist(positions[:,1].reshape(-1,1),positions[:,1].reshape(-1,1), lambda u, v: u-v)
        zz_distance_map = distance.cdist(positions[:,2].reshape(-1,1),positions[:,2].reshape(-1,1), lambda u, v: u-v)

        H = np.zeros((n_res*3, n_res*3))
        H[0::3,0::3] = xx_distance_map * xx_distance_map
        H[1::3,1::3] = yy_distance_map * yy_distance_map
        H[2::3,2::3] = zz_distance_map * zz_distance_map
        H[0::3,1::3] = xx_distance_map * yy_distance_map
        H[0::3,2::3] = xx_distance_map * zz_distance_map
        H[1::3,2::3] = yy_distance_map * zz_distance_map
        H[1::3,0::3] = xx_distance_map * yy_distance_map
        H[2::3,0::3] = xx_distance_map * zz_distance_map
        H[2::3,1::3] = yy_distance_map * zz_distance_map
        H = -1 * kirchoff_matrix * inv_distsq_map * H

        for i in range(n_res):
            H[i*3+0,i*3+0] = -1*np.sum(H[i*3+0,0::3])
            H[i*3+1,i*3+1] = -1*np.sum(H[i*3+1,1::3])
            H[i*3+2,i*3+2] = -1*np.sum(H[i*3+2,2::3])
            H[i*3+0,i*3+1] = -1*np.sum(H[i*3+0,1::3])
            H[i*3+1,i*3+0] = H[i*3+0,i*3+1]
            H[i*3+0,i*3+2] = -1*np.sum(H[i*3+0,2::3])
            H[i*3+2,i*3+0] = H[i*3+0,i*3+2]
            H[i*3+1,i*3+2] = -1*np.sum(H[i*3+1,2::3])
            H[i*3+2,i*3+1] = H[i*3+1,i*3+2]

        H_PP = H[:n_CA*3,:n_CA*3]
        H_PL = H[:n_CA*3,n_CA*3:]
        H_LL = H[n_CA*3:,n_CA*3:]

        H = H_PP - H_PL @ np.linalg.inv(H_LL) @ H_PL.T
       
        residue_mass = antibody_CA.residues.masses/1000
        residue_mass = np.repeat(residue_mass,3)

        M = np.diag(residue_mass)
        M_f = fractional_matrix_power(M,-0.5)

        H_m = (M_f@H)@M_f

        [v, U_m] = np.linalg.eig(H_m)
        U_m = U_m[:,np.argsort(v)[6:]]
        v = v[np.argsort(v)[6:]]


        RT = temperature*8.314/1000
        C = np.zeros((np.shape(H_m)))
        for i in range(len(v)):
            C += (np.outer((M_f@U_m[:,i]),U_m[:,i])@M_f)/v[i]
        C = RT*C

        self.M = M
        self.v = v
        self.U_m = U_m
        self.residue_mass = residue_mass
        self.covariance = C
        
    def apply_lrt(self):
        '''
        Predicts the direction of change from linear response theory.
        Finds the overlap of the modes with the direction of change.
        '''
        d_lrt = self.covariance @ self.forces.flatten() # This computes the direction
        d_lrt = fractional_matrix_power(self.M,0.5) @ d_lrt
        self.d_lrt = d_lrt / np.linalg.norm(d_lrt)
        self.O_lrt = np.dot(self.d_lrt, self.U_m) # This computes the column-wise dot product

    def top_modes(self):
        '''
        Finds the required top 'N' modes (a combination) with 95% overlap with the predicted direction.
        '''
        O_lrt_index = np.flip(np.argsort(np.abs(self.O_lrt)))
        distance = []
        dp = []
        for i in range(len(self.O_lrt)):
            direction = np.zeros(U_m.shape[0])
            for j in range(i+1):
                direction += self.O_lrt[O_lrt_index[j]] * self.U_m[:,O_lrt_index[j]]
            direction = direction/np.linalg.norm(direction)
            distance.append(np.linalg.norm(self.d_lrt-direction))
            dp.append(np.dot(self.d_lrt,direction))
            #if dp[-1] >=0.95:
            if distance[-1] <= 0.2:
                break
        print('Take top {} modes'.format(i+1))

        self.O_lrt_top = self.O_lrt[O_lrt_index[:i+1]]
        self.U_m_top = self.U_m[:,O_lrt_index[:i+1]]
        self.v_top = self.v[O_lrt_index[:i+1]]

    def topN_deformation_vector(self, rmsd):
        '''
        Generates the deformation vector of the protein for a given RMSD, using topN modes
        '''
        deformation_direction = (self.O_lrt_top@(self.U_m_top.T)).reshape(-1,3)
        deformation_magnitude = rmsd * np.sqrt(len(self.residue_mass)/3) / (np.sqrt(np.sum((self.O_lrt_top**2)))) ## Check '/3*3'
        self.deformation_vec = deformation_magnitude * deformation_direction

    def tot_deformation_vector(self, rmsd):
        '''
        Generates the deformation vector of the protein for a given RMSD, using all modes
        '''
        deformation_direction = (self.O_lrt@(self.U_m.T)).reshape(-1,3)
        deformation_magnitude = rmsd * np.sqrt(len(self.residue_mass)/3) / (np.sqrt(np.sum((self.O_lrt**2)))) ## Check '/3*3'
        self.deformation_vec = deformation_magnitude * deformation_direction
    
    def deform_protein(self):
        '''
        Deforms the coordinates of the protein along the deformation vector
        '''
        for r in range(1,self.pose.chain_end(2)+1): ##Manual entry
            for a in range(1,self.pose.residue(r).natoms()+1):
                new_coords = np.array(self.pose.residue(r).xyz(a)) + self.deformation_vec[r-1]
                self.pose.residue(r).set_xyz(a, py.rosetta.numeric.xyzVector_double_t(new_coords[0], new_coords[1], new_coords[2]))
                
    def tor_relax_deform(self,cycle):
        '''
        Relaxes the deformed protein in torsional space.
        '''
        relax = py.rosetta.protocols.relax.FastRelax()
        relax.set_scorefxn(self.sfxn)
        relax.max_iter(200)
        relax.constrain_coords(True)
        relax.ramp_down_constraints(True)
        relax.apply(self.pose)
        self.relaxed_name = self.name+'_deformed_'+str(cycle)+'.pdb' 
        self.pose.dump_pdb(self.relaxed_name)
    
    def cart_intf_relax_deform(self,cycle):
        '''
        Relaxes the deformed protein, in cartesian space, using InterfaceRelax script.
        '''
        relax = py.rosetta.protocols.relax.FastRelax(scorefxn_in=py.create_score_function("ref2015_cart_cst.wts"),standard_repeats=1,script_file="CustomRelax2022")

        relax.cartesian(True)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
        relax.max_iter(200)
        relax.constrain_coords(True)        
        relax.ramp_down_constraints(True)
        relax.apply(self.pose)
        
        self.relaxed_name = self.name+'_deformed_'+str(cycle)+'.pdb' 
        self.pose.dump_pdb(self.relaxed_name)
        
    def cstcart_intf_relax_deform(self,cycle,cst_weight,stddev):
        '''
        Relaxes the deformed protein, in cartesian space, using constraints.
        '''
        relax = py.rosetta.protocols.relax.FastRelax(scorefxn_in=py.create_score_function("ref2015_cart.wts"), standard_repeats=1,script_file="ConstrainedRelax2023")
        
        sfxn = py.create_score_function("ref2015_cart.wts")
        sfxn.set_weight(py.rosetta.core.scoring.ScoreType.coordinate_constraint, cst_weight) 
        relax.set_scorefxn(sfxn)
        
        relax.cartesian(True)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
        relax.max_iter(200)
        
        true_selector = py.rosetta.core.select.residue_selector.TrueResidueSelector()
        # Apply a virtual root onto the pose to prevent large lever-arm effects while minimizing with coordinate constraints.
        virtual_root = py.rosetta.protocols.simple_moves.VirtualRootMover()
        virtual_root.set_removable(True)
        virtual_root.set_remove(False)
        virtual_root.apply(self.pose)

        # Construct the CoordinateConstraintGenerator
        coord_constraint_gen = py.rosetta.protocols.constraint_generator.CoordinateConstraintGenerator()
        coord_constraint_gen.set_id("constrain_all_backbone_atoms")
        coord_constraint_gen.set_ambiguous_hnq(False)
        coord_constraint_gen.set_bounded(False) # Uses harmonic constraints.
        coord_constraint_gen.set_sidechain(False) # This in combination with the bottom, set_ca_only, means that only the backbone atoms are constrained.
        coord_constraint_gen.set_sd(stddev) # Sets a standard deviation of contrained atoms to stddev angstroms RMSD.
        coord_constraint_gen.set_ca_only(False)
        coord_constraint_gen.set_residue_selector(true_selector) # applies to all residues
        self.coord_constraint_gen = coord_constraint_gen

        # Apply the CoordinateConstraintGenerator using the AddConstraints mover
        add_constraints = py.rosetta.protocols.constraint_generator.AddConstraints()
        add_constraints.add_generator(coord_constraint_gen)
        add_constraints.apply(self.pose)
        
        relax.ramp_down_constraints(False)
        relax.apply(self.pose)
        self.relaxed_name = self.name+'_deformed_'+str(cycle)+'.pdb' 
        self.pose.dump_pdb(self.relaxed_name)
        
    def remove_constraints(self):
        remove_constraints = py.rosetta.protocols.constraint_generator.RemoveConstraints()
        remove_constraints.add_generator(self.coord_constraint_gen)
        remove_constraints.apply(self.pose)
    
    def bb_rmsd(self,pdb_relaxed, pdb_deformed,selection='(chainID H or chainID L) and backbone'):
        u = mda.Universe(pdb_relaxed,pdb_relaxed)
        v = mda.Universe(pdb_deformed,pdb_deformed)
        return mda.analysis.rms.RMSD(u,v,select=selection).run().results.rmsd[0,2]

    def cdr_bb_rmsd(self,pose_chothia):
        '''
        Compute OCD; FR/Loop RMSDs of self.pose wrt to pose_chothia
        '''
        pose_renum = py.pose_from_pdb(pose_chothia)
        ab_info = antibody.AntibodyInfo(pose_renum, antibody.Chothia_Scheme, antibody.North)

        H_selector = select.residue_selector.ChainSelector("H")
        L_selector = select.residue_selector.ChainSelector("L")
        fv_selector = select.residue_selector.OrResidueSelector(selector1=H_selector, selector2=L_selector)

        pose_fv = py.Pose()
        py.rosetta.core.pose.pdbslice(pose_fv, self.pose, fv_selector.selection_positions(self.pose))

        return np.array(antibody.cdr_backbone_rmsds(pose_fv, pose_renum, ab_info, ab_info)) 

    def report_energies(self):
        '''
        Reports the complex energy and the binding energy
        '''
        self.intf_mover.apply(self.pose)
        self.complex_energy = self.intf_mover.get_complex_energy()
        self.be = self.intf_mover.get_interface_dG()
        print('Complex energy: ', self.complex_energy)
        print('dG: ', self.be)
    
    def dump(self, pdb_name):
        '''
        Dumps the structure
        '''
        self.pose.dump_pdb(pdb_name)

class md_reu_lrt:
    '''
    Perform LRT excitation and relaxation using MD and calculates eneriges from REU.
    '''
    
    def __init__(self, py_selection):
        '''
        Initializes the pylrt object with a PDB file and an atom selection string.
        '''
        self.intf_mover = py.rosetta.protocols.analysis.InterfaceAnalyzerMover(py_selection)
    
    def compute_positions_mda_traj(self, topology, traj, selection):
        '''
        Computes the relative positions of atoms in the antibody with respect to the centre of mass of the antigen's interface, using MDanalysis.
        '''
        u = mda.Universe(topology,traj)
        self.selection = selection
        antibody = u.select_atoms(self.selection)
        self.antibody = antibody
        self.ab_len = antibody.n_residues

        ab_cb = u.select_atoms(self.selection +' and (name CA)')
        an_interface = u.select_atoms('(around 5 ((' + self.selection + ') and not (name H* or name [1-9]H*)))')

        relative_positions = []
        for _ in u.trajectory[-101:]:
            relative_positions.append(ab_cb.positions - an_interface.center_of_mass())
        self.relative_pos = np.delete(relative_positions,0,0)
        self.mean_pos = np.mean(self.relative_pos, axis=0) 
    
    def relax(self,pose,max_iters=200,cst_weight=1,stddev=1):
        '''
        Relaxes the deformed protein, in cartesian space, using constraints.
        '''
        self.pose = py.pose_from_pdb(pose)
        relax = py.rosetta.protocols.relax.FastRelax(scorefxn_in=py.create_score_function("ref2015_cart.wts"), standard_repeats=1,script_file="CustomRelax2022")
        
        sfxn = py.create_score_function("ref2015_cart.wts")
        sfxn.set_weight(py.rosetta.core.scoring.ScoreType.coordinate_constraint, cst_weight) 
        relax.set_scorefxn(sfxn)
        
        relax.cartesian(True)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
        
        true_selector = py.rosetta.core.select.residue_selector.TrueResidueSelector()
        # Apply a virtual root onto the pose to prevent large lever-arm effects while minimizing with coordinate constraints.
        virtual_root = py.rosetta.protocols.simple_moves.VirtualRootMover()
        virtual_root.set_removable(True)
        virtual_root.set_remove(False)
        virtual_root.apply(self.pose)

        # Construct the CoordinateConstraintGenerator
        coord_constraint_gen = py.rosetta.protocols.constraint_generator.CoordinateConstraintGenerator()
        coord_constraint_gen.set_ambiguous_hnq(False)
        coord_constraint_gen.set_bounded(False) # Uses harmonic constraints.
        coord_constraint_gen.set_sidechain(False) # This in combination with the bottom, set_ca_only, means that only the backbone atoms are constrained.
        coord_constraint_gen.set_sd(stddev) # Sets a standard deviation of contrained atoms to stddev angstroms RMSD.
        coord_constraint_gen.set_ca_only(False)
        coord_constraint_gen.set_residue_selector(true_selector) # applies to all residues
        self.coord_constraint_gen = coord_constraint_gen

        # Apply the CoordinateConstraintGenerator using the AddConstraints mover
        add_constraints = py.rosetta.protocols.constraint_generator.AddConstraints()
        add_constraints.add_generator(coord_constraint_gen)
        add_constraints.apply(self.pose)
        
        # Controlling backbone and torsional DOFs
        ab_chains = select.residue_selector.OR_combine(select.residue_selector.ChainSelector('A'), select.residue_selector.ChainSelector('B'))
        nbr = select.residue_selector.NeighborhoodResidueSelector()
        nbr.set_focus_selector(ab_chains)
        nbr.set_include_focus_in_subset(True)

        ## Only repack in antibody and neighbour residues of antibody
        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        tf.push_back(operation.RestrictToRepacking())
        tf.push_back(operation.IncludeCurrent())
        tf.push_back(operation.NoRepackDisulfides())
        prevent_repacking_rlt = operation.PreventRepackingRLT()
        prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr, flip_subset=True)
        tf.push_back(prevent_subset_repacking)

        ## Only energy minimize the backbone,side of antibody and neighbour residues of antibody
        mmf = MoveMapFactory()
        mmf.set_cartesian(setting=True)
        mmf.add_bb_action(move_map_action.mm_enable, nbr) #Enable whatever you want
        mmf.add_chi_action(move_map_action.mm_enable, nbr)
        mm = mmf.create_movemap_from_pose(self.pose)

        relax.set_movemap(mm)
        relax.set_task_factory(tf)
        
        relax.max_iter(max_iters) 
        relax.ramp_down_constraints(True)
        relax.apply(self.pose)

        relaxed_name = pose[:-4]+'_r.pdb'
        self.pose.dump_pdb(relaxed_name)

    def get_dG(self,pose):
        p = py.pose_from_pdb(pose)
        sfxn = py.create_score_function("ref2015.wts")
        sfxn(p)
        self.intf_mover.apply(p)
        return np.array(self.intf_mover.get_all_per_residue_data().dG)

    def compute_positions_mda(self, name, selection):
        '''
        Computes the relative positions of atoms in the antibody with respect to the centre of mass of the antigen's interface, using MDanalysis.
        '''
        u = mda.Universe(name,name)
        self.selection = selection
        antibody = u.select_atoms(selection)
        self.antibody = antibody
        self.ab_len = antibody.n_residues
        ab_ca = u.select_atoms('(' + selection + ') and (name CA)')
        an_interface = u.select_atoms('(around 5 ((' + selection + ') and not (name H* or name [1-9]H*)))')

        relative_positions = ab_ca.positions - an_interface.center_of_mass()
        return relative_positions

    def compute_be_reu(self,selection):
        dG = []
        relative_pos = []
        
        for i in range(49,100):
            #self.relax('mdframes/mdframe{}.pdb'.format(i))
            dG.append(self.get_dG('mdframes/mdframe{}_r.pdb'.format(i)))
            relative_pos.append(self.compute_positions_mda('mdframes/mdframe{}_r.pdb'.format(i),selection))
        
        self.relative_pos = np.array(relative_pos)
        self.mean_pos = np.mean(self.relative_pos, axis=0)
        
        self.dG = np.array(dG)[:,:self.ab_len] * 4.2 # Convert to kJoule/mole

    def compute_be_gmmpbsa(self,file):
        def load_data(file_list):
            dg_files = glob.glob(file_list)
            data = np.concatenate([np.genfromtxt(file, skip_header=1) for file in dg_files])
            data = data[np.argsort(data[:,0])]
            return data

        mm = load_data(file+'/*mm*dat')
        pol = load_data(file+'/*_pol*dat')
        apol = load_data(file+'/*_apol*dat')

        mmpbsa = mm+pol+apol
        mmpbsa[:,0] = mmpbsa[:,0]/3
        self.dG = mmpbsa[:,1:self.ab_len+1]
        
    def compute_be_gmx(self, files):
        with open(files) as f:
            for i, line in enumerate(f):
                if line.startswith('Total'):
                    skiprows = i
        df = pd.read_csv(files,skiprows=skiprows-1).dropna()
        x = df[df['Residue'].str.startswith('R')]['TOTAL'].to_numpy() * 4.2 #Convert to kJoule/mol
        self.dG = x.reshape(-1,self.ab_len)

    def quadratic_hypersurface(self,P, a, b, c, d, e, f, g, h,i, j):
        '''
        Fitting function for calculating the gradients of binding energies
        '''
        x, y, z = P
        return  a + b*x + c*y + d*z + e*x**2 + f*y**2 + g*z**2 + h*x*y + i*y*z + j*z*x

    def func_forces(self, pos_arr, deltaG_arr, mean_pos_arr):
        '''
        Forces on each residue "i"
        '''
        coeff, _ = curve_fit(self.quadratic_hypersurface, pos_arr.T, deltaG_arr)
    
        def dgx(P):
            x, y, z = P
            return -1*(coeff[1]+2*coeff[4]*x + coeff[7]*y + coeff[9]*z)
        def dgy(P):
            x, y, z = P
            return -1*(coeff[2]+2*coeff[5]*y + coeff[7]*x + coeff[8]*z)
        def dgz(P):
            x, y, z = P
            return -1*(coeff[3]+2*coeff[6]*z + coeff[8]*y + coeff[9]*x)  
        return np.stack((dgx(mean_pos_arr),dgy(mean_pos_arr),dgz(mean_pos_arr)))

    def compute_forces(self):
        '''
        Generates an array of forces for all residues
        '''
        self.forces = np.array([self.func_forces(self.relative_pos[:,i,:], self.dG[:,i], self.mean_pos[i]) for i in range(self.ab_len)]) 

    def generate_enm(self,pdb,temperature,cutoff):
        '''
        Generates an elastic network model including the environment.
        Calculates the eigenvectors and eigenvalues from the mass-weighted Hessian.
        Calculates the covariance matrix.
        Units: Mass: kg; length: Angs;
        '''
        u = mda.Universe(pdb,pdb)
        self.antibody_CA = u.select_atoms('(' + self.selection + ') and (name CA)')
        n_CA = self.ab_len
        complex_CA = u.select_atoms('name CA')
        n_res = complex_CA.n_residues
        positions = complex_CA.positions

        distance_map = distance.cdist(positions,positions)
        distance_map = np.repeat(np.repeat(distance_map,3,axis=1),3,axis=0)
        inv_distsq_map = distance_map**-2
        inv_distsq_map[inv_distsq_map == np.inf] = 0

        kirchoff_matrix = np.ones((n_res*3,n_res*3))
        kirchoff_matrix[distance_map>cutoff] = 0
        spring_constant = 1
        kirchoff_matrix = spring_constant*kirchoff_matrix

        xx_distance_map = distance.cdist(positions[:,0].reshape(-1,1),positions[:,0].reshape(-1,1), lambda u, v: u-v)
        yy_distance_map = distance.cdist(positions[:,1].reshape(-1,1),positions[:,1].reshape(-1,1), lambda u, v: u-v)
        zz_distance_map = distance.cdist(positions[:,2].reshape(-1,1),positions[:,2].reshape(-1,1), lambda u, v: u-v)

        H = np.zeros((n_res*3, n_res*3))
        H[0::3,0::3] = xx_distance_map * xx_distance_map
        H[1::3,1::3] = yy_distance_map * yy_distance_map
        H[2::3,2::3] = zz_distance_map * zz_distance_map
        H[0::3,1::3] = xx_distance_map * yy_distance_map
        H[0::3,2::3] = xx_distance_map * zz_distance_map
        H[1::3,2::3] = yy_distance_map * zz_distance_map
        H[1::3,0::3] = xx_distance_map * yy_distance_map
        H[2::3,0::3] = xx_distance_map * zz_distance_map
        H[2::3,1::3] = yy_distance_map * zz_distance_map
        H = -1 * kirchoff_matrix * inv_distsq_map * H

        for i in range(n_res):
            H[i*3+0,i*3+0] = -1*np.sum(H[i*3+0,0::3])
            H[i*3+1,i*3+1] = -1*np.sum(H[i*3+1,1::3])
            H[i*3+2,i*3+2] = -1*np.sum(H[i*3+2,2::3])
            H[i*3+0,i*3+1] = -1*np.sum(H[i*3+0,1::3])
            H[i*3+1,i*3+0] = H[i*3+0,i*3+1]
            H[i*3+0,i*3+2] = -1*np.sum(H[i*3+0,2::3])
            H[i*3+2,i*3+0] = H[i*3+0,i*3+2]
            H[i*3+1,i*3+2] = -1*np.sum(H[i*3+1,2::3])
            H[i*3+2,i*3+1] = H[i*3+1,i*3+2]

        H_PP = H[:n_CA*3,:n_CA*3]
        H_PL = H[:n_CA*3,n_CA*3:]
        H_LL = H[n_CA*3:,n_CA*3:]

        H = H_PP - H_PL @ np.linalg.inv(H_LL) @ H_PL.T
       
        residue_mass = self.antibody_CA.residues.masses/1000
        residue_mass = np.repeat(residue_mass,3)

        M = np.diag(residue_mass)
        M_f = fractional_matrix_power(M,-0.5)

        H_m = (M_f@H)@M_f

        [v, U_m] = np.linalg.eig(H_m)
        U_m = U_m[:,np.argsort(v)[6:]]
        v = v[np.argsort(v)[6:]]


        RT = temperature*8.314/1000
        C = np.zeros((np.shape(H_m)))
        for i in range(len(v)):
            C += (np.outer((M_f@U_m[:,i]),U_m[:,i])@M_f)/v[i]
        C = RT*C

        self.M = M
        self.v = v
        self.U_m = U_m
        self.residue_mass = residue_mass
        self.covariance = C

    def apply_lrt(self):
        '''
        Predicts the direction of change from linear response theory.
        Finds the overlap of the modes with the direction of change.
        '''
        d_lrt = self.covariance @ self.forces.flatten() # This computes the direction
        d_lrt = fractional_matrix_power(self.M,0.5) @ d_lrt
        self.d_lrt = d_lrt / np.linalg.norm(d_lrt)
        self.O_lrt = np.dot(self.d_lrt, self.U_m) # This computes the column-wise dot product

    def top_modes(self):
        '''
        Finds the required top 'N' modes (a combination) with 80% overlap (0.2) with the predicted direction.
        '''
        O_lrt_index = np.flip(np.argsort(np.abs(self.O_lrt)))
        distance = []
        dp = []
        for i in range(len(self.O_lrt)):
            direction = np.zeros(self.U_m.shape[0])
            for j in range(i+1):
                direction += self.O_lrt[O_lrt_index[j]] * self.U_m[:,O_lrt_index[j]]
            direction = direction/np.linalg.norm(direction)
            distance.append(np.linalg.norm(self.d_lrt-direction))
            dp.append(np.dot(self.d_lrt,direction))
            if dp[-1] >= 0.90:
                break
            #if distance[-1] <= 0.2:
            #    break
        print('Take top {} modes'.format(i+1))

        self.O_lrt_top = self.O_lrt[O_lrt_index[:i+1]]
        self.U_m_top = self.U_m[:,O_lrt_index[:i+1]]
        self.v_top = self.v[O_lrt_index[:i+1]]
        
    def excite(self, delT, box):
        '''
        Computes the velocities of excitation for a given delta T.
        Converts all vectors to atom-basis from residue(CA) basis.
        '''
        v_extra = np.zeros(self.U_m_top.shape[0])
        for i in range(self.v_top.shape[0]):
            v_extra += (self.v_top[i]**0.5) * (self.O_lrt_top[i]) * (self.U_m_top[:,i])
        v_extra = v_extra/ np.linalg.norm(v_extra)
        
        M_exite = np.zeros((self.antibody.n_atoms*3, self.antibody.n_atoms*3))
        np.fill_diagonal(M_exite,(self.antibody.masses.repeat(3)/1000) ** 0.5)
        
        v_excite = np.zeros(self.antibody.n_atoms*3)
        for i in range(self.antibody.n_residues):
            for j in self.antibody.residues[i].atoms.indices:
                v_excite[3*j] = v_extra[3*i]
                v_excite[3*j+1] = v_extra[3*i+1]
                v_excite[3*j+2] =  v_extra[3*i+2]
        
        u = mda.Universe(box,box)
        alpha = (3 * u.atoms.n_atoms * 0.008314) / (np.linalg.norm(M_exite @ v_excite)**2)
        lmda = (alpha * delT)**0.5
        v_extra = 0.316 *lmda * v_extra # Convert to Angs/ps
        self.v_extra = v_extra.reshape(-1,3)

    def save_gro(self,gro_in,delT, cycle):
        '''
        Save a gro file with new velocities.
        '''
        gro = mda.Universe(gro_in,gro_in)
        frame = gro.trajectory[-1]
        for i in range(self.ab_len):
            for j in gro.residues[i].atoms.indices:
                frame.velocities[j] += self.v_extra[i] #Careful with units. MDanalysis converts nm to angs automatically when reading and back to nm when writing. Therefore velocity length units must be in angstroms.
        gro.atoms.write('excited_it{}_{}K.gro'.format(cycle,delT), reindex=False)

    def cdr_bb_rmsd(self,p1, p2, pose_chothia):
        '''
        Compute OCD; FR/Loop RMSDs of two poses
        '''
        self.p1 = py.pose_from_pdb(p1)
        self.p2 = py.pose_from_pdb(p2)

        pose_renum = py.pose_from_pdb(pose_chothia)
        ab_info = antibody.AntibodyInfo(pose_renum, antibody.Chothia_Scheme, antibody.North)

        H_selector = select.residue_selector.ChainSelector("H")
        L_selector = select.residue_selector.ChainSelector("L")
        fv_selector = select.residue_selector.OrResidueSelector(selector1=H_selector, selector2=L_selector)

        p1_fv = py.Pose()
        p2_fv = py.Pose()
        py.rosetta.core.pose.pdbslice(p1_fv, self.p1, fv_selector.selection_positions(self.p1))
        py.rosetta.core.pose.pdbslice(p2_fv, self.p2, fv_selector.selection_positions(self.p2))

        return np.array(antibody.cdr_backbone_rmsds(p1_fv, p2_fv, ab_info, ab_info)) 

    def vecrmsd_lastframe(self,pdb_mdenm,pdb_ref,vector,selection='((chainID A or chainID B) and name CA)'):
        '''
        Finds the rmsd of a pdb wrt to the reference frame along a vector.
        '''
        u = mda.Universe(pdb_mdenm,pdb_mdenm)
        v = mda.Universe(pdb_ref,pdb_ref)
        vector = np.load(vector)
        vector = vector / np.linalg.norm(vector)

        u_fvca = u.select_atoms(selection)
        v_fvca = v.select_atoms(selection)
        _ = align.alignto(u_fvca,v_fvca)
        deformation = (u_fvca.positions - v_fvca.positions).flatten()
        return (np.dot(deformation,vector)**2/u_fvca.n_residues)**0.5

    def vecrmsd_traj(self,pdb_mdenm,xtc_mdenm,pdb_ref,vector,selection='((chainID A or chainID B) and name CA)'):
        '''
        Finds the RMSD of a traj wrt to a reference frame along a vector.
        '''
        u= mda.Universe(pdb_mdenm,xtc_mdenm)
        v = mda.Universe(pdb_ref,pdb_ref)
        vector = np.load(vector)
        vector = vector / np.linalg.norm(vector)
    
        u_fvca = u.select_atoms(selection)
        v_fvca = v.select_atoms(selection)
    
        dp = []
        for _ in u.trajectory:
            _ = align.alignto(u_fvca,v_fvca)
            deformation = (u_fvca.positions - v_fvca.positions).flatten()
            dp.append((np.dot(deformation,vector)**2/u_fvca.n_residues)**0.5)
        return np.array(dp)
