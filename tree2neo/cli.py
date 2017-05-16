import os
import os.path
from shutil import rmtree
from tempfile import NamedTemporaryFile, tempdir
import click
from .db import build_relationships, variants_to_fasta
from .docker import Docker
from .treeproc import FastTree
from .galaxy import submit_fasttree_job, wait_on_output, fetch_output

@click.group()
def cli():
    """
    This script parses a TREE file and builds a Neo4j Graph database.
    """
    pass


@cli.command()
@click.argument('tree_dir', type=click.Path(exists=True, dir_okay=True), required=True)
@click.argument('history_id', type=str, required=True)
@click.argument('refdb_dir', type=click.Path(exists=True, dir_okay=True), required=False)
# When running tree2neo with Dockerfile/docker-compose, we don't want docker inside docker.
@click.option('-d/-D', default=True, help='Run Neo4j docker container.')
def init(tree_dir, d, history_id, refdb_dir=None):
    """
    Copy reference database and load TREE to Neo4j Graph database.
    :param history_id:
    :param tree_dir:
    :param refdb_dir:
    :param d:
    :return:
    """
    if d:
        docker = Docker(refdb_dir=refdb_dir)
        docker.run()
    tree = FastTree(history_id, tree_dir=tree_dir)
    tree.process()
    build_relationships()


def load_tree_from_vsets(email, history_ids, outputdir=None):
    dir_made = False
    if outputdir is None:
        outputdir = os.path.join(tempdir, 'ft_' + str(os.getpid()) + '_working')
        os.mkdir(outputdir, 0o600)
        dir_made = True
    with NamedTemporaryFile(delete=False) as tmpfile:
        snp_count = variants_to_fasta(history_ids=history_ids, fasta_file=tmpfile)
        if snp_count > 0:
            tmpfile.close()
            history_name = ','.join(vset_names)
            run_result = submit_fasttree_job(email=email, fasta_filename=tmpfile.name, history_name=history_name)
            if run_result is not None:
                if 'jobs' in run_result and len(run_result['jobs']) == 1:
                    job_id = run_result['jobs'][0]['id']
                    output_id = wait_on_output(email=email, job_id=job_id)
                    if output_id is not None:
                        output_filename = fetch_output(email, outputdir, output_id)
                        if output_filename is not None:
                            tree = FastTree(history_name, tree_dir=outputdir)
                            tree.process()
                            build_relationships()
            os.remove(tmpfile.name)
    if dir_made:
        rmtree(outputdir)


@cli.command()
@click.argument('email', type=str, required=True)
@click.argument('history_ids', type=str, nargs=-1)
@click.option('--outputdir', type=str)
def tree_from_vsets(email, history_ids, outputdir):
    load_tree_from_vsets(email, history_ids, outputdir)


if __name__ == '__main__':
    cli()
