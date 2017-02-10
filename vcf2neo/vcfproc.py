"""
Interface to handle VCF files
"""
import glob
import time

import vcf
from db import create_variant_set_nodes, create_call_set_nodes, create_variant_site_nodes


class Vcf(object):
    """
    Handling VCF processing.
    """

    def __init__(self, vcf_dir=None):
        self.vcf_dir = vcf_dir

    def process(self):
        print("We have the following VCF files in directory ({}):\n".format(self.vcf_dir))
        for vcf_file in glob.glob(self.vcf_dir + "/*.vcf"):
            # TODO: Remove the two files from data
            if 'Drug' not in str(vcf_file):
                print(vcf_file)
                print("Processing: {}!".format(vcf_file))
                start = time.time()
                vcf_reader = vcf.Reader(open(vcf_file, 'r'))
                # TODO: Have a standard way of identifying variant_set_names
                vcf_file_name = str(vcf_file).replace(str(self.vcf_dir) + "/", "")
                # TODO: Let's use the file name for now
                create_variant_set_nodes(set_name=vcf_file_name)
                create_call_set_nodes(set_name=vcf_file_name)
                self.get_variant_sites(vcf_reader, vcf_file_name)
                end = time.time()
                print("Processed {} in {}!".format(vcf_file_name.upper(), end - start))
                time.sleep(2)

    def get_variant_sites(self, vcf_reader=None, vcf_file_name=None):
        for record in vcf_reader:
            print("\n")
            print(record)
            annotation = self.get_variant_ann(record)
            create_variant_site_nodes(record, annotation, vcf_file_name)

    @staticmethod
    def get_variant_ann(record=None):
        ann = record.INFO['ANN'][0].split('|')
        print(ann)
        return ann
