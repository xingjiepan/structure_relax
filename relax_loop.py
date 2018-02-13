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

def fast_relax(pose, residues_bb_movable, residues_sc_movable):
    '''Fast relax the pose'''
    # Set up move map

    mm = rosetta.core.kinematics.MoveMap()
    mm.set_chi(False)
    mm.set_bb(False)

    for i in residues_bb_movable:
        mm.set_bb(i, True)
        mm.set_chi(i, True)
    
    for i in residues_sc_movable:
        mm.set_chi(i, True)

    # Set up task factory

    turn_off_packing = rosetta.core.pack.task.operation.PreventRepacking()

    repackable_residues = residues_bb_movable + residues_sc_movable
    for i in range(1, pose.size() + 1):
        if not i in repackable_residues:
            turn_off_packing.include_residue(i)

    extra_rotamers = rosetta.core.pack.task.operation.ExtraRotamersGeneric()
    extra_rotamers.ex1(True)
    extra_rotamers.ex2(True)
    extra_rotamers.extrachi_cutoff(0)

    task_factory = rosetta.core.pack.task.TaskFactory()
    task_factory.push_back(rosetta.core.pack.task.operation.RestrictToRepacking())
    task_factory.push_back(turn_off_packing)

    # Do fast relax

    fast_relax_rounds = 5
    sfxn = rosetta.core.scoring.get_score_function()
    fast_relax = rosetta.protocols.relax.FastRelax(sfxn, fast_relax_rounds)
    fast_relax.set_movemap(mm) 
    fast_relax.set_task_factory(task_factory)

    fast_relax.apply(pose)
 
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

    # Print score differences

    #for i, v in energy_diff_sorted_table:
    #    if abs(v[0][1]) < 0.3: continue
    #    print i

    #    big_diffs = ['{0}'.format(x) for x in v if abs(x[1]) > 0.3]
    #    print '\t'.join(big_diffs)

    pose.dump_pdb(os.path.join('outputs', 'relaxed_' + os.path.basename(pdb_file)))

    return pose

def relax_loops(native_input_pdb_file, lowest_rmsd_file, lowest_score_file, loops_file):
    '''Relax loops and the surrounding residues. The surrounding residues
    are determined by the loops in the native_input_pdb_file.
    '''
    # Load the native pose and loops
    
    native_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(native_pose, native_input_pdb_file)
    
    loops = rosetta.protocols.loops.Loops(loops_file)

    # Get the fold tree
    
    fold_tree = rosetta.core.kinematics.FoldTree()
    rosetta.protocols.loops.fold_tree_from_loops(native_pose, loops, fold_tree, True)

    # Get loop residues

    loop_residues = []
    for i in range(loops.size()):
        loop_residues += list(range(loops[i + 1].start(), loops[i + 1].stop() + 1))

    # Get surrounding residues
    
    is_packable = rosetta.utility.vector1_bool(native_pose.size())
    for i in range(1, native_pose.size() + 1):
        is_packable[i] = False

    native_pose.update_residue_neighbors()
    rosetta.protocols.loops.select_loop_residues(native_pose, loops, True, is_packable, 10.0)

    surrounding_residues = [i for i in range(1, native_pose.size() + 1) if is_packable[i]]

    # Relax the structures

    relax_fun = lambda pose : fast_relax(pose, loop_residues, surrounding_residues)

    relax_structure(native_input_pdb_file, fold_tree, relax_fun)
    relax_structure(lowest_rmsd_file, fold_tree, relax_fun)
    relax_structure(lowest_score_file, fold_tree, relax_fun)

if __name__ == '__main__':
    pyrosetta.init()

    pdb_id = '4m6m'

    relax_loops('inputs/ama_ii_with_crystal_mates/{0}.pdb'.format(pdb_id), 
            'inputs/ama_ii_with_crystal_mates/{0}_lowest_rmsd.pdb.gz'.format(pdb_id), 
            'inputs/ama_ii_with_crystal_mates/{0}_lowest_score.pdb.gz'.format(pdb_id), 
            'inputs/ama_ii_with_crystal_mates/{0}.loop'.format(pdb_id))
