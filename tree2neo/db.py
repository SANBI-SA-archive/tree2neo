"""
Interface to the Neo4j Database
"""
from tree2neo.combat_tb_model.model.core import *
from tree2neo.combat_tb_model.model.vcfmodel import *
from tree2neo.combat_tb_model.model.galaxyuser import *


from py2neo import Graph, getenv, watch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

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
                t_set.from_variant_set.add(v_set)
                graph.push(t_set)


def get_galaxy_api_key(email):
    query = GalaxyUser.select(graph).where("_.email = '{}'".format(email))
    user = query.first()
    if user is None:
        return None
    else:
        return user.user_key


def variants_to_fasta(history_ids, fasta_file=open('output.fasta', 'w')):
    TB_LEN = 4410929  # length of H37Rv reference sequence. TODO: store Chromosome in DB and compute this
    ref_list = ['R'] * TB_LEN
    snp_positions = set()
    total_variant_count = 0
    snp_count = 0
    for history_id in history_ids:
        if history_id == 'refvcf':  # hack because refvcf has not history ID
            variant_set = VariantSet.select(graph).where("_.name = 'refvcf'").first()
        else:
            variant_set = VariantSet.select(graph).where("_.history_id = '{}'".format(history_id)).first()
        if variant_set is not None:
            for variant in variant_set.has_variant:
                total_variant_count += 1
                if len(variant.ref_allele) == 1:
                    for call in variant.has_call:
                        if len(call.alt_allele) != 1:
                            is_snp = False
                            break
                    else:
                        is_snp = True
                        ref_list[variant.pos - 1] = variant.ref_allele
                    if is_snp:
                        snp_positions.add(variant.pos)

        sequences = []
        snp_count = len(snp_positions)
        for callset in variant_set.has_callsets:
            seq_list = ref_list[:]
            for call in callset.has_call:
                if call.pos in snp_positions:
                    seq_list[call.pos - 1] = call.alt_allele[0]
            seq = ''.join([base for base in seq_list if base != 'R'])
            assert len(seq) == snp_count, "sequence length does not match SNP count for {} (length {} count {})".format(
                callset.name, len(seq), snp_count)
            seq_record = SeqRecord(id=callset.name, description=variant_set.history_id, seq=Seq(seq))
            sequences.append(seq_record)

    if snp_count > 0:
        ref_str = ''.join([base for base in ref_list if base != 'R'])
        assert len(ref_str) == snp_count, "reference sequence length does not match SNP count (length {} count {})".format(
                len(ref_str), snp_count)
        # print("total variants:", total_variant_count, "SNP count:", snp_count)
        ref_seq = SeqRecord(id='H37Rv', description='reference', seq=Seq(ref_str))
        sequences.insert(0, ref_seq)
        SeqIO.write(sequences, fasta_file, "fasta")
    return snp_count
