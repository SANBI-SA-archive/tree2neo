from __future__ import print_function
import sys
import os
import os.path
from shutil import rmtree
from tempfile import NamedTemporaryFile, gettempdir
import click
from .db import GraphDb
from .docker import Docker
from .treeproc import FastTree
from .galaxy import submit_fasttree_job, wait_on_output, fetch_output,\
    delete_history


@click.group()
def cli():
    """
    This script parses a TREE file and builds a Neo4j Graph database.
    """
    pass


@cli.command()
@click.argument('tree_dir', type=click.Path(exists=True, dir_okay=True),
                required=True)
@click.argument('history_id', type=str, required=True)
@click.argument('refdb_dir', type=click.Path(exists=True, dir_okay=True),
                required=False)
# When running tree2neo with Dockerfile/docker-compose,
# we don't want docker inside docker.
@click.option('--db_host', type=str, default='localhost')
@click.option('-d/-D', default=False, help='Run Neo4j docker container.')
def init(tree_dir, d, history_id, db_host, refdb_dir=None):
    """
    Copy reference database and load TREE to Neo4j Graph database.
    :param history_id:
    :param tree_dir:
    :param refdb_dir:
    :param d:
    :return:
    """
    if d:
        if refdb_dir is None:
            exit("When running in Docker mode we need an output dir")
        docker = Docker(refdb_dir)
        docker.run()
        db_host = 'localhost'
        http_port = docker.http_port
        bolt_port = docker.bolt_port
    else:
        http_port = 7474
        bolt_port = 7687
    tree = FastTree(history_id, tree_dir=tree_dir)
    tree.process()
    db = GraphDb(host=db_host, password='', use_bolt=False,
                 http_port=http_port, bolt_port=bolt_port)
    db.build_relationships()
    d.stop()


def load_tree_from_vsets(api_key, history_ids, outputdir=None,
                         db_host='localhost',
                         galaxy_url='http://192.168.2.218'):
    db = GraphDb(host=db_host, password='', use_bolt=False,
                 http_port=7474, bolt_port=7687)
    dir_made = False
    if outputdir is None:
        outputdir = os.path.join(gettempdir(), 'ft_' + str(os.getpid()) +
                                 '_' + api_key + '_working')
        if os.path.isdir(outputdir):
            rmtree(outputdir)
        os.mkdir(outputdir, 700)
        dir_made = True
    with NamedTemporaryFile(delete=False) as tmpfile:
        print("tempfile:", tmpfile.name, file=sys.stderr)
        snp_count = db.variants_to_fasta(history_ids=history_ids,
                                         fasta_file=tmpfile)
        if snp_count > 0:
            tmpfile.close()
            history_name = ','.join(history_ids)
            print("submitting job to Galaxy", file=sys.stderr)
            result = submit_fasttree_job(api_key=api_key,
                                         fasta_filename=tmpfile.name,
                                         history_name=history_name,
                                         galaxy_url=galaxy_url)
            print("Galaxy job submitted, result:",
                  result[0], file=sys.stderr)
            if result is not None:
                (run_result, history_id) = result
                if 'jobs' in run_result and len(run_result['jobs']) == 1:
                    job_id = run_result['jobs'][0]['id']
                    output_id = wait_on_output(api_key=api_key,
                                               job_id=job_id,
                                               galaxy_url=galaxy_url)
                    if output_id is not None:
                        print("fetching Galaxy output", file=sys.stderr)
                        output_filename = fetch_output(api_key, outputdir,
                                                       output_id,
                                                       galaxy_url=galaxy_url)
                        if output_filename is not None:
                            print("loading FastTree node", file=sys.stderr)
                            tree = FastTree(history_name, tree_dir=outputdir)
                            tree.process()
                            db.build_relationships()
                            delete_history(api_key=api_key,
                                           history_id=history_id,
                                           galaxy_url=galaxy_url)
            os.remove(tmpfile.name)
    if dir_made:
        rmtree(outputdir)


@cli.command()
@click.argument('api_key', type=str, required=True)
@click.argument('history_ids', type=str, nargs=-1)
@click.option('--outputdir', type=str)
@click.option('--db_host', type=str, default='localhost')
@click.option('--galaxy_url', type=str, default='http://192.168.2.218')
def tree_from_vsets(api_key, history_ids, outputdir, db_host, galaxy_url):
    load_tree_from_vsets(api_key, history_ids, outputdir, db_host, galaxy_url)


if __name__ == '__main__':
    cli()
