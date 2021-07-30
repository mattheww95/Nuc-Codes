"""
Convert genome coordinates based off of gff
Fasta is single line, and 1 based
"""
import os
import sys
import platform
from dataclasses import dataclass
import pickle

codons = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


class GffFinder:
    """
    This is the start up class to read in attributes which will be dumped and read back in
    """
    def __init__(self, ref_fasta, ref_gff):
        self.ref_fasta = self.get_fasta(ref_fasta)
        self.ref_gff = self.get_gff(ref_gff)

    @staticmethod
    def get_fasta(fasta_path):
        """
        Return only the reference sequence used with no header
        :param fasta_path:
        :return: Reference sequence
        """
        fastalines = []
        with open(fasta_path, "r") as fasta:
            fastalines = fasta.readlines()

        carat_count = 0
        for ii in fastalines:
            for iii in ii:
                if iii == '>':
                    carat_count += 1

        if carat_count > 1:
            raise AssertionError("Only one fasta can be specified as a reference.")
        else:
            fasta_header = fastalines.pop(0)
            print(f"Reference sequence used: {fasta_header}")
            fasta_seq = ''.join(fastalines)
            if platform.system() == "Windows":
                fasta_seq.replace("\r\n", "")
            else:
                fasta_seq.replace("\n", "")

        return fasta_seq

    @staticmethod
    def get_gff(gff_path):
        """
        Read in the gff path, and parse the needed information
        :param gff_path: Path to the gff file
        :return: A parsed gff with attributes easily accessible
        """

        with open(gff_path, "r") as gff:
            gff_obj = gff.readlines()

        # remove un-needed information from the gff, such as link to its source

        # The gff is split by tabs and semicolons it seems
        # gff always has 9 columns
        # seqid, source (where its from), type, start, end, score, strand (f or r), phase
        # (0 first base of the codon,1 second base of the codon , 2 third codon base), attributes semicolon seperated
        gff_attributes = []
        for i in gff_obj:
            try:
                gff_lines = i.split("\t")
                attributes = GFFFile.split_attributes(gff_lines[8])
                score = gff_lines[5].replace(".", "0.0")
                start = gff_lines[3].replace(".", "0")
                end = gff_lines[4].replace(".", "0")
                phase = gff_lines[7].replace(".", "0")
                #print(gff_lines)
                gff_attributes.append(GFFFile(gff_lines[0], gff_lines[1], gff_lines[2], int(start),
                                              int(end), float(score), gff_lines[6], int(phase),
                                              attributes))
            except IndexError:
                ...

        return gff_attributes

    def query_attributes(self, attribure_query, attribute_option,  *argv):
        """
        select criteria from gff
        :param attribure_query: Which Key in attributes will be used to filter your data
        :param attribute_option: The option chosen for the filter
        :param argv: convert to amminoacids and return that position, can take multiple sites
        :return:
        """
        for i in self.ref_gff:
            attribute_name = i.attibutes.get(attribure_query)
            if attribute_name == attribute_option:
                for arg_i in argv[0]:
                    gene_placement = 1
                    for ii in range(i.start-1, i.end, 3):  # deal with base 1 encoding
                        if int(arg_i) == gene_placement:
                            codon = self.ref_fasta[ii:ii+3]
                            print(f"AA: {codons[codon]} Codon: {codon} {attribute_option} Position: {gene_placement} "
                                  f"Genomic Position: {ii+1}")
                        gene_placement += 1

    def display_gff_keys(self):
        """
        Display gff keys that can be chosen
        :return: sys exit
        """

        gff_keys = []
        for i in self.ref_gff:
            att_keys = [ii for ii in i.attibutes.keys()]
            gff_keys.extend(att_keys)

        print("The possible keys for you to filter your sequence based off of are:")
        print(set(gff_keys))
        sys.exit()  # leave the program




@dataclass
class GFFFile:
    """
    This class will split up the attributes of the file
    """
    seqid: str
    source: str
    type_feature: str
    start: int
    end: int
    score: float
    strand: str
    phase: int
    attibutes: dict

    @staticmethod
    def split_attributes(the_attributes_portion):
        """
        Split up all the attributes as needed
        :param the_attributes_portion: Literally the last column of the gff3 file
        :return:
        """
        return_dict = {}
        attribute = the_attributes_portion.replace("attributes=", "")
        attribute = attribute.replace("'", "").replace('"', "")
        attribute = attribute.strip()
        attribute = attribute.split(";")
        for i in attribute:
            attribute_value = i.split("=")
            return_dict[attribute_value[0].strip()] = attribute_value[1].strip()  # make sure they are naked!

        return return_dict


#def dump_gff_finder(obj: GFFFile):
#    """
#    Pickle the classs so ref fasta does not need to be specified each time
#    :param obj:
#    :return:
#    """


def help_message():
    """

    ISSUE: These messages are coupled to tightly to the other methods used to call them.
    :return:
    """
    print("nuc-code help options")
    print("A program to help Cody and hopefully make him help me more.")
    print("")
    print("Initialize Program to save a gff and reference sequence.")
    print("nuc-code init <reference_fasta> <reference_gff> <project_name>")
    print("")
    print("Check possible attributes to use for gff extraction.")
    print("nuc-code <project_name> attributes")
    print("e.g. nuc-code <project_name> attributes")
    print("")
    print("Query a gene positions, you can query one or multiple conditions")
    print("nuc-code <project_name> <gff attribute key> <gff attribute key value> <Amino acid positions>")
    print("e.g. nuc-code SARS Name orf1ab 4 5 6")
    print("")
    print("Help note")
    print("If no attribute is displayed check that the values your are using exist in your gff.")
    print("Note of Importance:")
    print("The Phase attribute of GFF file is not included in the tabulation of nucleotides.")
    print("It is also assumed that the GFF being used is base 1 indexed.")
    return 0


def initialize_class_storage(fasta: os.path, gff: os.path, project_name):
    """
    Make the directory to store the pickled class to be read in later
    :param project_name: The directory to write the pickled class out too.
    :param fasta: The reference fasta to be used.
    :param gff: The reference gff file to be used.
    :return:
    """
    fasta = os.path.abspath(fasta)
    gff = os.path.abspath(gff)
    new_project = os.path.join(os.path.dirname(__file__), project_name)
    if os.path.isdir(new_project):
        sys.exit(f"project {project_name} already exits.")

    if "fa" not in os.path.basename(fasta) or "gff" not in os.path.basename(gff):
        print(f"Files entered are: {fasta} and {gff}")
        print("Could not identify one of the file types. The program can accept: ")
        print("Fasta files with the extensions: fasta, fa, fas")
        print("gff files with extensions of: gff, gff3")
        sys.exit("Check that your fasta and gff files are entered in the appropriate order.")
    try:
        test_fa = open(fasta, "r")
        test_gff = open(gff, "r")
        os.mkdir(new_project)
        test_fa.close()
        test_gff.close()
    except FileNotFoundError:
        sys.exit("Could not open fasta or gff files")

    new_gff = GffFinder(fasta, gff)
    new_proj_file = os.path.join(new_project, project_name)
    with open(new_proj_file, "wb") as project:  # folder and pickled file will have same name
        pickle.Pickler(project, pickle.DEFAULT_PROTOCOL).dump(new_gff)
    print(f"Created new project {project_name}")
    print(f"Query positions in the future enter nuc-code {project_name} <options>")
    sys.exit()


def read_initilized_class(project_name, attribure_query, attribute_option, *argv):
    """
    execute an exisiting query
    :param project_name: The name of an existing project
    :param attribure_query: The query to execute (all these are in the class description)
    :param attribute_option: The attribute key in the gff
    :param argv: The gene positions to search
    :return:
    """
    project_path = os.path.join(os.path.dirname(__file__), project_name, project_name)
    if not os.path.exists(project_path):
        print(f"Could not find {project_name}, are you sure it has been initialized?")
        sys.exit("To initialize a project enter: nuc-code init <ref fasta> <ref gff> <project name>")

    with open(project_path, "rb") as prev_class:
        gff_class = pickle.load(prev_class)
        gff_class.query_attributes(attribure_query, attribute_option, *argv)

    return 0


def get_class_attrs(project_name):
    """
    Probable could be part of read initialized class but made this seperate may refactor later
    :param project_name:
    :return:
    """
    project_path = os.path.join(os.path.dirname(__file__), project_name, project_name)
    if not os.path.exists(project_path):
        print(f"Could not find {project_name}, are you sure it has been initialized?")
        sys.exit("To initialize a project enter: nuc-code init <ref fasta> <ref gff> <project name>")

    with open(project_path, "rb") as prev_class:
        gff_class = pickle.load(prev_class)
        gff_class.display_gff_keys()

    return 0


def control_gff_finder(sys_args: list):
    help_msg = {"help", "-h", "--help"}
    file_loc = os.path.dirname(__file__)
    directory_info = os.listdir(file_loc)
    if len(sys_args) < 3:
        if sys_args[1] in help_msg:
            help_message()
    elif sys_args[1] == "init":
        # need to add check for if enough args
        initialize_class_storage(sys_args[2], sys_args[3], sys_args[4])
    elif sys_args[1] in directory_info:
        if sys_args[2] == "attributes":
            get_class_attrs(sys_args[1]) # needs proj name to display attributes
        else:
            read_initilized_class(sys_args[1], sys_args[2], sys_args[3], sys_args[4:])

    return 0


if __name__ == "__main__":
    control_gff_finder(sys.argv)
