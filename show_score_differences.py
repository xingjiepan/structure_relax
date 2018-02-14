#!/usr/bin/env python2

import sys

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

def print_score_differences(pdb_file1, pdb_file2):
    '''Get the score differences between two structures in pdb files.
    '''
    # Load pose

    pose1 = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose1, pdb_file1)
    pose2 = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose2, pdb_file2)

    # Score

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn(pose1)
    scores1 = get_pose_scores(pose1)
    sfxn(pose2)
    scores2 = get_pose_scores(pose2)

    # Output

    energy_diff_sorted_table = sort_residues_by_abs_total_energy_diff(scores1, scores2)

    for i, v in energy_diff_sorted_table:
        if abs(v[0][1]) < 0.3: continue
        print i, 'pose' if i == 0 else pose1.residue(i).name3()

        big_diffs = ['{0}'.format(x) for x in v if abs(x[1]) > 0.3]
        print '\t'.join(big_diffs)


if __name__ == '__main__':
    pyrosetta.init()

    pdb_file1 = sys.argv[1]
    pdb_file2 = sys.argv[2]

    print_score_differences(pdb_file1, pdb_file2)

