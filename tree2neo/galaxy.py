from __future__ import print_function
import sys
from time import sleep
import os.path
from bioblend.galaxy import GalaxyInstance
# from .db import get_galaxy_api_key

def get_gi(api_key):
    # api_key = get_galaxy_api_key(email)
    if api_key is not None:
        # 192.168.2.218 is ctbgx.sanbi.ac.za
        gi = GalaxyInstance(url='http://192.168.2.218', key=api_key)
    else:
        gi = None
    return gi


def submit_fasttree_job(api_key, fasta_filename, history_name='working_history'):
    gi = get_gi(api_key)
    if gi is not None:
        history = gi.histories.create_history(name=history_name)
        history_id = history['id']
        upload_result = gi.tools.upload_file(fasta_filename, history_id)
        if 'outputs' in upload_result and len(upload_result['outputs']) == 1:
            uploaded_dataset_id = upload_result['outputs'][0]['id']
            run_result = gi.tools.run_tool(history_id, 'fasttree', dict(
                input_alignment=dict(src='hda', id=uploaded_dataset_id)
            ))
            return (run_result, history_id)
    return None


def get_job_state(api_key, job_id):
    gi = get_gi(api_key)
    job_state = gi.jobs.get_state(job_id)
    return job_state


def wait_on_output(api_key, job_id):
    job_state = 'new'
    good_job_states = set(('new', 'queued', 'running'))
    while job_state != 'ok' and job_state != 'error' and job_state in good_job_states:
        print('galaxy job state:', job_state, file=sys.stderr)
        job_state = get_job_state(api_key, job_id)
        sleep(1)
    print("final Galaxy job state:", job_state, file=sys.stderr)
    if job_state == 'ok':
        print('Galaxy job completed successfully', file=sys.stderr)
        gi = get_gi(api_key)
        if gi is not None:
            job_result = gi.jobs.show_job(job_id)
            # fasttree generates a nhx (newick tree) format output called output_tree
            if 'exit_code' in job_result and job_result['exit_code'] == 0:
                print("job completed with exit code 0, all good")
                if 'outputs' in job_result and 'output_tree' in job_result['outputs']:
                    tree_output_id = job_result['outputs']['output_tree']['id']
                    return tree_output_id
    return None


def fetch_output(api_key, output_path, output_id):
    gi = get_gi(api_key)
    if gi is not None:
        output_filename = os.path.join(output_path, output_id + '.nhx')
        gi.datasets.download_dataset(output_id, file_path=output_filename, wait_for_completion=True,
                                     use_default_filename=False)
        return output_filename
    return None


def delete_history(api_key, history_id):
    gi = get_gi(api_key)
    if gi is not None:
        gi.delete_history(history_id, purge=True)
