"""
Operations on hydra.

This is a simple wrapper for launching commands on hydra.
"""
import os
import random
import json
import time
import string
import subprocess
import fileOps
import pipeline


hydra_template = '''#!/usr/bin/env bash
#$ -N "{name}"
#$ -V
#$ -pe threads {cpu}
#$ -l mem_free={mem}
#$ -cwd
#$ -o {name}_stdout.txt
#$ -e {name}_stderr.txt
#$ -S "/usr/bin/env bash"

{cmd}
'''


def run_hydra(cmd, cpu=1, mem='6G', sleep_timer=10, name=None):
    """
    Runs cmd on hydra.
    :param cmd: string
    :param cpu: CPU to request
    :param mem: memory to request
    :param sleep_timer: Time in seconds to wait before checking on status. if None, will not wait.
    :param name: Job name. If not set, will be randomly generated
    :return: None
    """
    if name is None:
        name = ''.join([random.choice(string.ascii_letters) for _ in xrange(10)])
    with fileOps.TemporaryFilePath() as sub:
        hydra_str = hydra_template.format(cpu=cpu, mem=mem, cmd=cmd, name=name)
        with open(sub, 'w') as outf:
            outf.write(hydra_str)
        with open(sub) as fh:
            job_id = subprocess.check_output(['qsub'], stdin=fh).split()[2]
        if sleep_timer is None:
            return
        while True:
            exit_code = get_job_exit_code(job_id)
            if exit_code is None:
                time.sleep(sleep_timer)
            elif exit_code == 0:
                return
            else:
                raise pipeline.ProcException('Job failed with exit code {}'.format(exit_code))


def get_job_exit_code(job_id):
    # check if hydra succeeded and failed paths exist
    succeeded_path = os.path.join(os.getenv('HYDRA_ROOT'), 'archive', 'succeeded', job_id)
    failed_path = os.path.join(os.getenv('HYDRA_ROOT'), 'archive', 'failed', job_id)

    # if neither exist, job is still running
    if not os.path.exists(succeeded_path) and not os.path.exists(failed_path):
        return None

    if os.path.exists(succeeded_path):
        return 0

    complete_file = os.path.join(failed_path, 'complete')
    if not os.path.exists(complete_file):
        return 42  # lost host, exit value is 42

    results = json.load(open(complete_file))
    if 'exit_code' in results:
        exit_code = results['exit_code']
    else:
        exit_code = 1
    return exit_code