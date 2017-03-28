"""
Interface to the Neo4j Database
"""
from tree2neo.combat_tb_model.model.core import *
from tree2neo.combat_tb_model.model.vcfmodel import *
from tree2neo.combat_tb_model.model.galaxyuser import *
from tree2neo.combat_tb_model.model.fasttree import *

from py2neo import Graph, getenv, watch

graph = Graph(host=getenv("DB", "thoba.sanbi.ac.za"), http_port=7474, bolt=True, password=getenv("NEO4J_PASSWORD", ""))

watch("neo4j.bolt")


def create_tree_nodes(name, data, history_id):
    """
    Create VariantSet Nodes
    :return:
    """
    v_set = FastTree(name=str(name), data=str(
        data), history_id=str(history_id))
    graph.create(v_set)


def build_relationships():
    """
    Build Relationships
    :return:
    """
    t_sets = FastTree.select(graph)
    v_sets = VariantSet.select(graph)
    for v_set in v_sets:
        for t_set in t_sets:
            # TODO: Find a better way to handle this.
            if v_set.history_id == t_set.history_id:
                t_set.has_var.add(v_set)
                graph.push(t_set)
