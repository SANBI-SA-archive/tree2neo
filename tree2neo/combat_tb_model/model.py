from neomodel import (StructuredNode, StringProperty,
                      IntegerProperty, FloatProperty,
                      RelationshipTo, RelationshipFrom,
                      Relationship, One, OneOrMore)


class Gene(StructuredNode):
    """
    Genes
    """
    # print 'Gene Nodes'
    gene_id = StringProperty(unique_index=True)
    uniprot_entry = StringProperty(index=True)
    dbxref = StringProperty(index=True)
    ncbi_id = StringProperty(index=True)
    ncbi_acc = StringProperty(index=True)
    ensembl_id = StringProperty(index=True)
    name = StringProperty(index=True)
    preffered_name = StringProperty()
    locus_tag = StringProperty(index=True)
    gene_synonym = StringProperty(index=True)
    coding = StringProperty()
    protein_id = StringProperty(index=True)
    pseudo = StringProperty()
    cdc_ortholog = StringProperty(index=True)
    sub_feature = StringProperty()
    start = IntegerProperty()
    end = IntegerProperty()
    strand = StringProperty()
    description = StringProperty()
    citation = StringProperty()
    location = StringProperty()
    length = IntegerProperty()
    fasta = StringProperty()

    # @property
    # def getLocation(self):
    #   location = '{}...{} ({})'.format(self.start, self.end, self.strand)
    #   return location
    transcribed = RelationshipTo('Transcript', 'TRANSCRIBED', One)
    translated = RelationshipTo('CDS', 'TRANSLATED', One)
    has_ortholog = RelationshipTo('Ortholog', 'ORTHOLOG', OneOrMore)
    has_go_terms = RelationshipTo('GoTerm', 'ASSOCIATED_WITH', OneOrMore)
    has_interpro_terms = RelationshipTo('InterPro', 'ASSOCIATED_WITH', OneOrMore)
    mentioned_in = RelationshipTo('Publication', 'MENTIONED_IN', OneOrMore)


class Pseudogene(StructuredNode):
    """
    Transcripts
    """
    # print 'Pseudogene Nodes'
    pseudogene_id = StringProperty(unique_index=True)
    name = StringProperty()
    gene_id = StringProperty(index=True)
    description = StringProperty(index=True)
    biotype = StringProperty()
    start = IntegerProperty()
    end = IntegerProperty()
    strand = IntegerProperty()
    location = StringProperty()
    sequence = StringProperty()

    has_ortholog = RelationshipTo('Ortholog', 'ORTHOLOG', OneOrMore)


class Transcript(StructuredNode):
    """
    Transcripts
    """
    # print 'Transcript Nodes'
    # 'transcript_id' will be the name
    transcript_id = StringProperty(Unique_Index=True, index=True, required=True)
    # print 'Transcript Nodes'
    # transcript_id = StringProperty(unique_index=True)
    name = StringProperty()
    gene = StringProperty(index=True)
    note = StringProperty(index=True)
    coding = StringProperty()
    start = IntegerProperty()
    end = IntegerProperty()
    strand = IntegerProperty()
    location = StringProperty()
    # TODO: @thobalose I am planning to group these types
    # [ncRNA, tRNA, rRNA, transcript]
    _type = StringProperty(index=True)
    product = StringProperty(index=True)
    parent = StringProperty()  # The Parent gene
    # gene = StringProperty()
    sequence = StringProperty()

    transcribed = RelationshipFrom('Gene', One)
    translated = RelationshipTo('Protein', 'TRANSLATED', One)
    encodes = RelationshipTo('CDS', 'PROCESSED_INTO', One)
    part_of = RelationshipFrom('Exon', OneOrMore)

    yields = Relationship('Trna', 'HAS', OneOrMore)
    __yields = Relationship('NCrna', 'HAS', OneOrMore)
    __yields_ = Relationship('Rrna', 'HAS', OneOrMore)


class Trna(StructuredNode):
    """
    tRNA
    """
    # print 'Trna Nodes'
    trna_id = StringProperty(unique_index=True)
    name = StringProperty()
    gene_id = StringProperty(index=True)
    note = StringProperty(index=True)
    biotype = StringProperty()
    start = IntegerProperty()
    end = IntegerProperty()
    strand = IntegerProperty()
    location = StringProperty()
    sequence = StringProperty()

    holds = RelationshipFrom('Transcript', OneOrMore)


class NCrna(StructuredNode):
    """
    ncRNA
    """
    # print 'NCrna Nodes'
    ncrna_id = StringProperty(unique_index=True)
    name = StringProperty()
    gene_id = StringProperty(index=True)
    note = StringProperty(index=True)
    biotype = StringProperty()
    start = IntegerProperty()
    end = IntegerProperty()
    strand = IntegerProperty()
    location = StringProperty()
    sequence = StringProperty()

    holds = RelationshipFrom('Transcript', OneOrMore)


class Rrna(StructuredNode):
    """
    rRNA
    """
    # print 'Rrna Nodes'
    rrna_id = StringProperty(unique_index=True)
    name = StringProperty()
    gene_id = StringProperty(index=True)
    note = StringProperty(index=True)
    biotype = StringProperty()
    start = IntegerProperty()
    end = IntegerProperty()
    strand = IntegerProperty()
    location = StringProperty()
    sequence = StringProperty()

    holds = RelationshipFrom('Transcript', OneOrMore)


class Exon(StructuredNode):
    """
    Exon
    """
    # print 'Exon Nodes'
    exon_id = StringProperty()
    name = StringProperty()
    location = StringProperty()
    transcript = StringProperty(index=True)

    part_of = RelationshipTo('Transcript', 'PART_OF', OneOrMore)


class CDS(StructuredNode):
    """
    CDS
    """
    # print 'CDS Nodes'
    cds_id = StringProperty(unique_index=True)
    name = StringProperty(index=True)
    transcript = StringProperty()
    product = StringProperty(index=True)
    protein_id = StringProperty(index=True)

    composed_of = RelationshipTo('Exon', 'COMPOSED_OF', OneOrMore)
    translated_ = RelationshipTo('Protein', 'TRANSLATED', OneOrMore)
    translated = RelationshipFrom('Gene', One)


class Protein(StructuredNode):
    """
    Proteins
    """
    # print 'Protein Nodes'
    protein_id = StringProperty(unique_index=True)
    dbxref = StringProperty(index=True)
    ncbi_id = StringProperty(index=True)
    ncbi_acc = StringProperty(index=True)
    uniprot_id = StringProperty(index=True)
    swissprot_id = StringProperty(index=True)
    pdb = StringProperty()
    ensembl_id = StringProperty(index=True)
    parent = StringProperty(index=True)
    name = StringProperty(index=True)
    recommended_name = StringProperty()
    sequence = StringProperty()
    domain = StringProperty()
    family = StringProperty()
    function = StringProperty()
    three_d = StringProperty()
    length = IntegerProperty()
    mass = StringProperty()
    pdb_id = StringProperty()
    transcript = StringProperty(index=True)
    start = IntegerProperty()
    end = IntegerProperty()
    strand = IntegerProperty()
    interactor = StringProperty(index=True)

    # TODO: clarify the meaning here. 'interacts' means
    # 'interacts with TB protein', 'interacts_' means
    # 'interacts with human protein'
    interacts = RelationshipTo('Protein', 'INTERACTS_WITH', OneOrMore)
    interacts_ = RelationshipTo('HumanProtein', 'INTERACTS_WITH', OneOrMore)
    associated_with = RelationshipTo('GoTerm', 'ASSOCIATED_WITH', OneOrMore)
    associated_ = RelationshipTo('InterPro', 'ASSOCIATED_WITH', OneOrMore)
    translated_ = RelationshipFrom('CDS', OneOrMore)


class HumanProtein(StructuredNode):
    """
    HumanProteins
    """
    # print 'Protein Nodes'
    protein_id = StringProperty(unique_index=True)
    tb_protein = StringProperty(index=True)
    dbxref = StringProperty(index=True)
    ncbi_id = StringProperty(index=True)
    ncbi_acc = StringProperty(index=True)
    uniprot_id = StringProperty(index=True)
    swissprot_id = StringProperty(index=True)
    pdb = StringProperty()
    ensembl_id = StringProperty(index=True)
    parent = StringProperty(index=True)
    name = StringProperty(index=True)
    sequence = StringProperty()
    length = IntegerProperty()
    transcript = StringProperty(index=True)
    start = IntegerProperty()
    end = IntegerProperty()
    strand = IntegerProperty()
    interactor = StringProperty(index=True)

    interacts_ = RelationshipFrom('Protein', OneOrMore)
    # associated_with = RelationshipTo('GoTerm', 'ASSOCIATED_WITH', OneOrMore)
    # associated_ = RelationshipTo('InterPro', 'ASSOCIATED_WITH', OneOrMore)


class Ortholog(StructuredNode):
    """
    Ortholog
    """
    # print 'Ortholog Nodes'
    locus_name = StringProperty(unique_index=True)
    uniprot_id = StringProperty(index=True)
    organism = StringProperty()

    has_ortholog = RelationshipFrom('Gene', OneOrMore)


class GoTerm(StructuredNode):
    """
    GO Terms
    """
    # print 'GO Nodes'
    go_id = StringProperty(unique_index=True)
    name = StringProperty(index=True)
    namespace = StringProperty(index=True)
    is_a = StringProperty()

    is_a_ = Relationship('GoTerm', 'IS_A', OneOrMore)
    _genes = RelationshipFrom('Gene', OneOrMore)


class InterPro(StructuredNode):
    """
    InterPro
    """
    # print 'InterPro Nodes'
    interpro_id = StringProperty(unique_index=True)
    name = StringProperty()
    _genes = RelationshipFrom('Gene', OneOrMore)


class Pfam(StructuredNode):
    """
    Pfam
    """
    # print 'Pfam Nodes'
    pfam_id = StringProperty(index=True)
    name = StringProperty()


class Domain(StructuredNode):
    # print 'Domain Nodes'
    name = StringProperty(index=True)


class Publication(StructuredNode):
    pubmed_id = StringProperty(unique_index=True)
    pubmed_ci = StringProperty()
    title = StringProperty(index=True)
    authors = StringProperty(index=True)
    source = StringProperty()
    journal = StringProperty(index=True)

    mentioned_in = RelationshipFrom('Gene', OneOrMore)


class VariantCollection(StructuredNode):
    collection_name = StringProperty(unique_index=True)
    contains = RelationshipTo('Variant', 'CONTAINED_IN', OneOrMore)

    def __repr__(self):
        return 'VariantCollection({name})'.format(name=self.collection_name)


class Variant(StructuredNode):
    chrom = StringProperty(required=True)
    pos = IntegerProperty(required=True)
    start = IntegerProperty()
    end = IntegerProperty()
    id = StringProperty()
    ref = StringProperty(required=True)
    alt = StringProperty(required=True)
    qual = FloatProperty()
    filter = StringProperty()
    info = StringProperty()
    collection = RelationshipFrom('VariantCollection', One)

    @property
    def __repr__(self):
        if self.qual is None:
            qual = '.'
        else:
            qual = self.qual
        return (
            'Variant({chr} {pos} {id} {ref} {alt} {qual} {filter} {info})'.
                format(chr=self.chrom, pos=self.pos, id=self.id,
                       ref=self.ref, alt=self.alt,
                       qual=qual, filter=self.filter, info=self.info
            )
        )
