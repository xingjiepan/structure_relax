#!/usr/bin/env python2

import sys

import pyrosetta
from pyrosetta import rosetta


if __name__ == '__main__':
    pyrosetta.init()

    pdb_file = sys.argv[1]

    pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(pose, pdb_file)

    sfxn = rosetta.core.scoring.get_score_function()
    sfxn.set_weight(rosetta.core.scoring.fa_dun, 0)
    sfxn.set_weight(rosetta.core.scoring.fa_rep, 0)
    sfxn(pose)

    print "Total score =", pose.energies().total_energy()
