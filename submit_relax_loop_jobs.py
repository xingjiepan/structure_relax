#!/usr/bin/env python2
'''Submit the loop relax jobs to the SGE cluster.
Usage:
    ./submit_relax_loop_jobs.py input_path output_file num_jobs 
'''

import os
import sys
import subprocess


if __name__ == '__main__':
    input_path = sys.argv[1]
    output_file = sys.argv[2]
    num_jobs = int(sys.argv[3])

    for f in os.listdir('job_outputs'):
        os.remove(os.path.join('job_outputs', f))
    
    cmd = ['qsub',
          '-e', 'job_outputs',
          '-o', 'job_outputs',
          '-t', '1-{0}'.format(num_jobs),
          './relax_loop.py',
          input_path,
          output_file]

    subprocess.check_call(cmd)
