"""
Interface to handle VCF files
"""
import glob
import time

from db import create_tree_nodes


class FastTree(object):
    """
    Handling Tree processing.
    """

    def __init__(self, tree_dir=None):
        self.tree_dir = tree_dir

    def process(self):
        print("We have the following TREE files in directory ({}):\n".format(self.tree_dir))
        for tree_file in glob.glob(self.tree_dir + "/*.nhx"):
            start = time.time()

            # read the file contents
            with open(tree_file) as tree_f:
                data_str = tree_f.read()

            # TODO: Pass Tree File Name in the tool arguments
            #tree_file_name = str(tree_file).replace(str(self.tree_dir) + "/", "")
            tree_file_name = str(tree_file).rsplit('/', 1)[-1]

            # TODO: Let's use the file name for now
            create_tree_nodes(name=tree_file_name, data=data_str)

            end = time.time()
            print("Processed {} in {}!".format(tree_file_name.upper(), end - start))
            time.sleep(2)