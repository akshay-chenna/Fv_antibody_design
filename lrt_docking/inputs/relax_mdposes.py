import sys
import pyrosetta as py
from pyrosetta.rosetta.core import select
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, move_map_action
py.init('-out:level 0 -ex1 -ex2')

def relax(pose, reps, max_iters=200, cst_weight=1,stddev=1):
    '''
    Relax script
    '''
    p = pose
    pose = py.pose_from_pdb(pose)
    relax = py.rosetta.protocols.relax.FastRelax(scorefxn_in=py.create_score_function("ref2015_cart.wts"), standard_repeats=int(reps),script_file="CustomRelax2022")

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
    virtual_root.apply(pose)

    # Construct the CoordinateConstraintGenerator
    coord_constraint_gen = py.rosetta.protocols.constraint_generator.CoordinateConstraintGenerator()
    coord_constraint_gen.set_ambiguous_hnq(False)
    coord_constraint_gen.set_bounded(False) # Uses harmonic constraints.
    coord_constraint_gen.set_sidechain(False) # This in combination with the bottom, set_ca_only, means that only the backbone atoms are constrained.
    coord_constraint_gen.set_sd(stddev) # Sets a standard deviation of contrained atoms to stddev angstroms RMSD.
    coord_constraint_gen.set_ca_only(False)
    coord_constraint_gen.set_residue_selector(true_selector) # applies to all residues
    coord_constraint_gen = coord_constraint_gen

    # Apply the CoordinateConstraintGenerator using the AddConstraints mover
    add_constraints = py.rosetta.protocols.constraint_generator.AddConstraints()
    add_constraints.add_generator(coord_constraint_gen)
    add_constraints.apply(pose)

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
    mm = mmf.create_movemap_from_pose(pose)

    relax.set_movemap(mm)
    relax.set_task_factory(tf)

    relax.max_iter(max_iters) 
    relax.ramp_down_constraints(True)
    relax.apply(pose)

    relaxed_name = p[:-4]+'_r.pdb'
    pose.dump_pdb(relaxed_name)
    
relax(sys.argv[1], sys.argv[2])
