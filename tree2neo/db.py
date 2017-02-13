"""
Interface to the Neo4j Database
"""
from combat_tb_model.model import VariantSet, CallSet, VariantSite, Call, Gene, Feature, Tree

from py2neo import Graph, getenv, watch

graph = Graph(host=getenv("DB", "localhost"), http_port=7474, bolt=True, password=getenv("NEO4J_PASSWORD", ""))
watch("neo4j.bolt")


def create_tree_nodes(name, data):
    """
    Create VariantSet Nodes
    :return:
    """
    v_set = Tree(name=str(name), data=str(data))
    graph.create(v_set)

def build_relationships():
    """
    Build Relationships
    :return:
    """
    t_sets = Tree.select(graph)
    v_sets = VariantSet.select(graph)
    for v_set in v_sets:
        for t_set in t_sets:
            # TODO: Find a better way to handle this.
            if v_set.name == t_set.name:
                t_set.has_calls_in.add(v_set)
                graph.push(t_set)