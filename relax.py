#!/usr/bin/env python2

import os

import pyrosetta
from pyrosetta import rosetta


def get_pose_scores(pose):
    '''Get the scores for a pose into a list of dictionaries.
    The first dict stores the scores of the whole pose. The
    rest of dicts store the scores of each residue.
    '''
    # Get the non-zero weight score types
    
    weights = pose.energies().weights()
    score_types = rosetta.core.scoring.ScoreType.__dict__
    non_zero_weight_types = []
 
    for k in score_types.keys():
        if k.startswith('__'): continue
        
        if weights[score_types[k]] != 0:
            non_zero_weight_types.append(k)

    # Get whole pose scores

    total_energies = pose.energies().total_energies() 
    d_pose = {}

    for t in non_zero_weight_types:
        d_pose[t] = total_energies[score_types[t]]
    d_pose['total_energy'] = pose.energies().total_energy()
    
    scores = [d_pose]

    # Get residue scores

    for i in range(1, pose.size() + 1):
        res_energies = pose.energies().residue_total_energies(i) 
        d_res = {}

        for t in non_zero_weight_types:
            d_res[t] = res_energies[score_types[t]]
        d_res['total_energy'] = pose.energies().residue_total_energy(i)

        scores.append(d_res)

    return scores

def sort_residues_by_abs_total_energy_diff(scores1, scores2):
    '''Sort the residues by the absolute total score differences
    in two list of scores. Then sort the energy terms for each
    of the residue.
    '''
    energy_diff = []

    for i in range(len(scores1)):
        energy_diff_d = {}
        for k in scores1[i].keys():
            energy_diff_d[k] = scores2[i][k] - scores1[i][k]
        
        energy_diff.append((i, energy_diff_d))

    energy_diff_sorted_dicts = sorted(energy_diff, key=lambda x : abs(x[1]['total_energy']), reverse=True)

    # Generate the sorted table

    energy_diff_sorted_table = []
    for i, d in energy_diff_sorted_dicts:
        pairs = sorted([(k, d[k]) for k in d.keys()], key=lambda x: abs(x[1]), reverse=True)
        energy_diff_sorted_table.append((i, pairs))

    return energy_diff_sorted_table

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
    scores_before_relax = get_pose_scores(pose)

    # Relax
   
    relax_fun(pose)

    # Output

    sfxn(pose)
    print 'score_after_relax =', pose.energies().total_energy()
    scores_after_relax = get_pose_scores(pose)

    energy_diff_sorted_table = sort_residues_by_abs_total_energy_diff(scores_before_relax, scores_after_relax)

    for i, v in energy_diff_sorted_table:
        if abs(v[0][1]) < 0.3: continue
        print i

        big_diffs = ['{0}'.format(x) for x in v if abs(x[1]) > 0.3]
        print '\t'.join(big_diffs)

    pose.dump_pdb(os.path.join('outputs', 'relaxed_' + os.path.basename(pdb_file)))



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
    #relax_structure('inputs/4ma3_lowest_rmsd.pdb', ft, min_relax)

