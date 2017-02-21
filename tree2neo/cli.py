import click
from db import build_relationships
from docker import Docker
from treeproc import FastTree


@click.group()
def cli():
    """
    This script parses a TREE file and builds a Neo4j Graph database.
    """
    pass


@cli.command()
@click.argument('tree_dir', type=click.Path(exists=True, dir_okay=True), required=True)
@click.argument('refdb_dir', type=click.Path(exists=True, dir_okay=True), required=False)
@click.argument('history_id', type=str, required=True)
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


if __name__ == '__main__':
    cli()
