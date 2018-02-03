#!/usr/bin/env python2

import os

import pyrosetta
from pyrosetta import rosetta


def get_pose_scores(pose):
    '''Get the scores for a pose into a list of dictionaries.'''
    # Get the non-zero weight score types
    
    weights = pose.energies().weights()
    score_types = rosetta.core.scoring.ScoreType.__dict__
    non_zero_weight_types = []
 
    for k in score_types.keys():
        if k.startswith('__'): continue
        
        print score_types[k]
        if weights.get(score_types[k]) != 0:
            non_zero_weight_types.append(k)

    print non_zero_weight_types

    print pose.energies().total_energies().show_nonzero()

def minimize(pose, residues_bb_movable, residues_sc_movable):
    '''Minimize the pose.'''
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn.set_weight(rosetta.core.scoring.chainbreak, 5)

    mm = rosetta.core.kinematics.MoveMap()
    mm.set_chi(False)
    mm.set_bb(False)

    for i in residues_bb_movable:
        mm.set_bb(i, True)
        mm.set_chi(i, True)
    
    for i in residues_sc_movable:
        mm.set_chi(i, True)

    min_opts = rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, True )

    min_mover = rosetta.protocols.simple_moves.MinMover()
    min_mover.movemap(mm)
    min_mover.min_options(min_opts)
    min_mover.score_function(sfxn)

    min_mover.apply(pose)

def relax_structure(pdb_file, fold_tree, relax_fun):
    '''Relax the structure from a pdb file.
    Dump the relaxed structure.
    '''
    # Load pose

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)
    pose.fold_tree(fold_tree)
    rosetta.core.pose.correctly_add_cutpoint_variants(pose)

    # Score

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose)
    print 'score_before_relax =', pose.energies().total_energy()

    # Relax
   
    relax_fun(pose)

    # Output

    sfxn(pose)
    print 'score_after_relax =', pose.energies().total_energy()

    pose.dump_pdb(os.path.join('outputs', 'relaxed_' + os.path.basename(pdb_file)))

    get_pose_scores(pose)


if __name__ == '__main__':
    pyrosetta.init()

    loop_residues_4am3 = list(range(311,321))

    surrounding_residues_4am3 = [2, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 64, 66, 67, 68, 71, 87, 88, 89, 90, 91, 92, 93, 94, 98, 99, 100, 101, 102, 103, 219, 220, 221, 222, 223, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 274, 275, 276, 278, 285, 287, 291, 292, 293, 294, 295, 306, 307, 308, 309, 310, 322, 323, 324, 325]

    ft = rosetta.core.kinematics.FoldTree()
    ft.add_edge(1, 322, 1)
    ft.add_edge(1, 320, -1)
    ft.add_edge(322, 321, -1)
    ft.add_edge(322, 433, -1)

    min_relax = lambda pose : minimize(pose, loop_residues_4am3, surrounding_residues_4am3)

    relax_structure('inputs/4ma3_native.pdb', ft, min_relax)

