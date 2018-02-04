#!/usr/bin/env python2
'''Copy the rotamers of a structure to the backbone of another structure.
Usage:
    ./copy_rotamers.py source_file target_file
'''

import sys
import os

import pyrosetta
from pyrosetta import rosetta


if __name__ == '__main__':
    pyrosetta.init()

    source_file = sys.argv[1]
    target_file = sys.argv[2]
    

    source_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(source_pose, source_file)
    target_pose = rosetta.core.pose.Pose()
    rosetta.core.import_pose.pose_from_file(target_pose, target_file)

    for seqpos in range(1, source_pose.size() + 1):
        for i in range(1, source_pose.residue(seqpos).nchi() + 1):
            target_pose.set_chi(i, seqpos, source_pose.chi(i, seqpos))


    target_pose.dump_pdb('rotamer_changed_' + os.path.basename(target_file))
