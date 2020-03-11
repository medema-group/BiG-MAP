# ! /usr/bin/env python3

"""
--------------- MGC/BGC module ---------------
Author: Koen van den Berg
University: Wageningen University and Research
Department: Department of Bioinformatics
Date: 21/05/2019
----------------------------------------------
This code was inspired by the BiG-SCAPE project, which can be found
at: https://git.wageningenur.nl/medema-group/BiG-SCAPE/blob/master/. I
highly recommend someone that would like to alter this script to also
take a look there.
The general idea of this script is to create a GCFs database (.fasta)
to use as input for the pipeline. It parses the DNA and protein
sequences from the antiSMASH genbank files, calculates GCFs using
mash, and outputs the database.
Additional information:
This script takes gut- and antiSMASH output directories as input,
converts the .gbk files to fasta files and then calculates gene
cluster families based on the similarity treshold (default=0.9) using
fastANI. In addition, HMMer is used to find several relevant
housekeeping genes from the whole genome genbank files, which are also
present in the antiSMASH output. These housekeeping genes will be
essential downstream for comparing the resutls of the gene clusters
with resutls that are known a priori. The output consists of the
following items: GCFs clusters, GCF fasta file, Mash results, and
optionally bigscape GCF clusters. Dependencies: BioPython, awk
"""

# Import statements
import os.path
import subprocess
from sys import argv
import sys
import argparse
from Bio import SeqIO
import json
import re
from pathlib import Path
import numpy as numpy
from Bio.SeqFeature import FeatureLocation, ExactPosition, BeforePosition, AfterPosition
import pickle
import time
from subprocess import Popen, PIPE

import shutil


def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
                                     usage=''' 
______________________________________________________________________
  BiG-MAP Family: creates a redundancy filtered reference fna 
______________________________________________________________________
Generic command: python3 BiG-MAP.family.py [Options]* 
-D [input dir(s)] -O [output dir]
Create a redundancy filtered fasta reference file from multiple
anti/gutSMASH outputs. Use BiG-MAP_process conda environment.
Obligatory arguments: 
    -D   Specify the path to the directory containing the gut- or
         antiSMASH outputs here. This could be a singular directory,
         or a space seperated list of directories.
    -O   Put path to the folder where the MASH filtered gene cluster
         files should be located here. The folder should be an
         existing folder. Default = current folder (.)
Options:
    -tg  Fraction between 0 and 1; the similarity treshold that
         determines when the protein sequences of the gene clusters
         can be considered similar. Default = 0.6.
    -th  Fraction between 0 and 1; the similarity treshold that
         determines when the protein sequences of the housekeeping genes
         can be considered similar. Default = 0.8
    -f   Specify here the number of genes that are flanking the core
         genes of the gene cluster. 0 --> only the core, n --> n
         genes included that flank the core. default = 0
    -g   Output whole genome fasta files for the MASH filtered gene 
         clusters as well. This uses more disk space in the output 
         directory. 'True' | 'False'. Default = False
    -b   Name of the path to bigscape.py. Default = False
    -p   Number of used parallel threads in the MASH and BiG-SCAPE
         filtering step. Default = 6
______________________________________________________________________
''')

    parser.add_argument("-D", "--indir", help=argparse.SUPPRESS, nargs="+", required=True)
    parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-tg", "--treshold_GC", help=argparse.SUPPRESS,
                       required=False, default=0.8, type=float)
    parser.add_argument("-th", "--treshold_HG", help=argparse.SUPPRESS,
                       required=False, default=0.1, type=float)
    parser.add_argument("-f", "--flank_genes",
                       help=argparse.SUPPRESS, required=False, type=int, default=0)
    parser.add_argument("-g", "--genomefiles",
                       help=argparse.SUPPRESS, required=False, type=bool, default=False)
    parser.add_argument("-b", "--bigscape_path", help=argparse.SUPPRESS, default=False)
    parser.add_argument("-pf", "--pfam", help=argparse.SUPPRESS, default=False)
    parser.add_argument("-c", "--cut-off", help=argparse.SUPPRESS, 
                       type=float, default= 0.2)
    parser.add_argument("-p", "--threads", help=argparse.SUPPRESS,
                       type=int, required=False, default=6)
    return (parser.parse_args())


######################################################################
# Functions for extracting the gene clusters from anti/gutSMASH outputs
######################################################################
def retrieveclusterfiles(indir):
    """Retrieves the .gbk files from the input directory including their paths.
    parameters
    ----------
    indir
        string, path to the gut/antiSMASH dirs
    returns
    ----------
    """
    try:
        for dirpath, dirnames, files in os.walk(indir):
            for f in files:
                if f.endswith(".gbk") and "region" in f:
                    yield (os.path.join(dirpath, f))
    except:  # Include actual exception here
        # In the case that the user just already inputs fasta
        # files, I can make sure that an error is captured
        # here. In that way, the program will skip this step and
        # move straigth towards the fastANI computation.
        pass


def getgenomegbk(infile):
    """Yields the genome.gbk file from region input file
    parameters
    ----------
    infile
        string, name of the input .gbk file
    returns
    ----------
    genomegbkfile = genomefilename
    """
    genome = infile.split("/")[-2]
    path = "/".join(infile.split("/")[:-1])
    if genome.endswith(".gbff"):
        genome = genome[:-5]
        return ("{}/{}.gbk".format(path, genome))
    else:
        return ("{}/{}.gbk".format(path, genome))


def parsegbkcluster(infile, nflank):
    """Parses the genbank files for DNA, protein, cluster, organism
    parameters
    ----------
    infile
        string, name of the input .gbk file
    returns
    ----------
    DNA = DNA sequence
    proteins = protein sequence
    clustername = name of the cluster
    organism = name of the organism
    #cluster_enzymes = {loc:gene_kind}
    """
    DNA = ""
    core_DNA = ""
    proteins = []
    GCs = []
    organism = ""
    feature_count = 0
    CDS_index = []
    core_index = []
    core_relative_locs = []
    absolute_locs = []
    gbkcontents = SeqIO.parse(infile, "genbank")
    for record in gbkcontents:
        for feature in record.features:
            # Parsing the regionname
            if "region" in feature.type:
                if "product" in feature.qualifiers:
                    for product in feature.qualifiers["product"]:
                        GCs.append(product)
            # Parsing the protein sequence
            if feature.type == "CDS":
                CDS_index.append(feature_count)  # remembering every CDS gene index
                if "translation" in feature.qualifiers.keys():
                    proteins.append(feature.qualifiers['translation'][0])
                # Parsing the relative core locations
                if "gene_kind" in feature.qualifiers.keys():
                    kind = feature.qualifiers['gene_kind'][0]
                    if kind == "biosynthetic":
                        core_index.append(feature_count)
            feature_count += 1

        # VALIDATION
        ##############################
        absolute_loc_start = record.annotations["structured_comment"]["antiSMASH-Data"]["Orig. start"]
        absolute_loc_end = record.annotations["structured_comment"]["antiSMASH-Data"]["Orig. end"]
        absolute_locs = [absolute_loc_start, absolute_loc_end]
        ##############################

        if CDS_index.index(min(core_index)) - nflank < 0 or CDS_index.index(max(core_index)) + nflank + 1 > len(
                CDS_index):
            print(
                f"!!!flank_genes (-f) is higher than the number of flanking genes in the cluster of file: {infile}, using whole gene cluster instead!!!")
        if nflank == 0:
            core_relative_locs = [record.features[i].location for i in core_index]
        else:
            if CDS_index.index(min(core_index)) - nflank >= 0:
                core_region = CDS_index[
                              CDS_index.index(min(core_index)) - nflank:CDS_index.index(max(core_index)) + nflank + 1]
            else:
                core_region = CDS_index[0:CDS_index.index(max(core_index)) + nflank + 1]
            core_relative_locs = [record.features[i].location for i in core_region]

        # Parsing the DNA sequence
        DNA = record.seq
        # organism = record.annotations['definition'].replace(" ", "_")
        organism = record.description
        organism = "_".join(organism.split(",")[0].split()[:-1])
        organism = organism.replace("(", "")
        organism = organism.replace(")", "")
        organism = organism.replace("/", "-")
    return (DNA, "".join(proteins), ":".join(GCs), organism, core_relative_locs, absolute_locs)


def writefasta(sequences, seqstype, cluster, organism, infile, outdir):
    """Writes the fasta file for each sequence
    parameters
    ----------
    sequences
        string, nucleotide/AA sequence
    seqstype
        string, either DNA or PROT
    region
        string, the gene cluster name
    organism
        string, name of the organism
    infile
        string, name of the input .gbk file
    outdir
        string, path to output directory
    returns
    ----------
    Outfile = the name of the written fasta file
    orgID = ID for the input organism genome
    fasta_header = header of the cluster dna sequence
    """
    regionno = infile.split("/")[-1].split(".")[-2]
    orgID = infile.split("/")[-1].split(".")[0]
    versionno = infile.split("/")[-1].split(".")[1]
    orgID = orgID[8:] if 'PROT' in orgID else orgID  # For housekeeping orgIDs
    Outfile = f"{outdir}{seqstype + cluster if 'HG' in seqstype else seqstype}-{orgID}.{organism}.{regionno}.fasta"
    fasta_header = f">gb|{orgID}.{versionno}.{regionno}|{seqstype}--Entryname={cluster}--OS={organism}--SMASHregion={regionno}"
    seq = "\n".join(str(sequences)[i:i + 80] for i in range(0, len(str(sequences)), 80))
    if not os.path.exists(Outfile):
        with open(Outfile, "w") as o:
            o.write(f"{fasta_header}\n")
            o.write(f"{seq}\n")
    return (Outfile, orgID, fasta_header[1:])


def locs2bedfile(indict, outfile):
    """saves the enzymes dictionary to a bed file format
    parameters
    ----------
    indict
        dictionary, enzyme location dictionary
    GCF_dict
        dict, dictionary of the gene cluster families
    returns
    ----------
    None
    """
    bedfile = outfile
    with open(bedfile, "w") as w:
        for clust, locs in indict.items():
            for loc in locs:
                start = str(loc.start)
                start = start.replace("<", "")
                end = str(loc.end)
                end = end.replace(">", "")
                w.write(f"{clust}\t{start}\t{end}\n")
    return(bedfile)


######################################################################
# similarity between gene clusters using MASH
######################################################################
def make_sketch(outdir, threads, option="GC"):
    """
    Calculates the distance between the query fasta files
    stored in the sketch file by using mash.
    Parameters
    ----------
    outdir
        string, the path to output directory
    threads
        int, number of threads
    option
        string, either 'GC' for the gene clusters or 'HG' for the housekeeping genes
    returns
    ----------
    """
    with open (f"{outdir}log.file", "wb") as log_file:
        try:
            if option == "GC":
                cmd_mash = f"mash sketch -o {outdir}mash_sketch -k 16 -p {threads} -s 5000 -a {outdir}GC_PROT*"
                p = Popen(cmd_mash, shell=True, stdout=PIPE, stderr=PIPE)
                stdout, stderr = p.communicate()
                log_file.write(stderr)

            elif option == "HG":
                cmd_mash = f"mash sketch -o {outdir}mash_sketch -k 4 -p {threads} -s 1000 -a {outdir}HG_PROT*"
                p = Popen(cmd_mash, shell=True, stdout=PIPE, stderr=PIPE)
                stdout, stderr = p.communicate()
                log_file.write(stderr)

        except(subprocess.CalledProcessError):
            # Raise error here for error table
            pass

def calculate_distance(outdir, output_file="mash_output_GC.tab"):
    """
    Calculates the distance between the query fasta files
    stored in the sketch file by using mash.
    Parameters
    ----------
    outdir
        string, the path to output directory
    returns
    ----------
    """
    infile = f"{outdir}mash_sketch.msh"
    outfile = f"{outdir}{output_file}"
    try:
        cmd_mash = f"mash dist {infile} {infile} > {outfile}"
        res_download = subprocess.check_output(cmd_mash, shell=True)

    except(subprocess.CalledProcessError):
        # Raise error here for error table
        pass
    return ()

######################################################################
# Extracting LSH clusters (buckets) using cut-off & calculate medoid
######################################################################
def calculate_medoid(outdir, cut_off, med={}, input_file="mash_output_GC.tab"):
    """
    calculates the GCFs based on similarity treshold
    parameters and calculates the medoid of that GCF
    ----------
    outdir
        string, the path to output directory
    cut_off
        float, between 0 and 1
    returns
    ----------
    dict_medoids = {fasta file of medoid: similar fasta files}
    """
    # Parse the input into a dictionary of gene families
    family_by_gene = {}
    family_by_gene_filtered = {}
    family_members = {}
    family_distance_matrices = {}
    dict_medoids = med
    infile = f"{outdir}{input_file}"
    with open(infile, "r") as input:
        for line in input:
            # Skip lines starting with '#'
            if line.startswith('#'):
                continue
            # Split into tab-separated elements
            gene1, gene2, distance, number2, overlap = line.strip().split('\t')
            no1, no2 = overlap.split("/")
            overlap = 1-(float(no1)/float(no2))
            # Look up the family of the first gene
            if gene2 in family_by_gene.keys():
                family_name = family_by_gene[gene2]
            else:
                family_name = gene2
                family_by_gene[gene2] = family_name
                family_members[family_name] = []
                family_distance_matrices[family_name] = []
            # filter genes that are already part of another family
            if gene2 == family_by_gene[gene2]:
                family_by_gene_filtered[gene2] = family_by_gene[gene2]
            # If gene1 and gene2 don't overlap, then there are two options:
            # 1. gene1 and gene2 do belong to the same family, and the "not overlap" is "odd"
            #    in this case we accept gene2 into our family
            # 2. gene1 and gene2 actually belong to different families, and we might have to create that family
            if not overlap <= cut_off:
                # gene1 doesn't overlap at all, so put that one into a separate family
                if gene1 in family_by_gene.keys():
                    if family_by_gene[gene1] == family_name:
                        # gene1 is in our family, so record the distance
                        add_to_distance_matrix(family_distance_matrices[family_name], family_members[family_name], gene2,
                                            gene1, float(distance))
                    else:
                        # gene is above cut off or doesn't belong to the family
                        pass
                else:
                    # gene1 doesn't have a family yet, make one
                    gene1_family_name = gene1
                    family_by_gene[gene1] = gene1_family_name
                    family_members[gene1_family_name] = []
                    family_distance_matrices[gene1_family_name] = []
                    # insert gene1 into that family as only member, with a distance of 0
                    add_to_distance_matrix(family_distance_matrices[gene1_family_name], family_members[gene1_family_name],
                                        gene1, gene1, float(distance))
            else:
                # There is some overlap, and we want gene1 in this family
                family_by_gene[gene1] = family_name
                add_to_distance_matrix(family_distance_matrices[family_name], family_members[family_name], gene2, gene1,
                                    float(distance))
    # For each family: Build a distance matrix, and then work out the mediod
    for family_name in family_by_gene_filtered.keys():
        # Calculate the mediod from the distances
        np_array = numpy.asarray(family_distance_matrices[family_name])
        medoid_index = numpy.argmin(np_array.sum(axis=0))
        # Create a dictionary using the medoid as key and the family_members as values
        dict_medoids[family_members[family_name][medoid_index]] = family_members[family_name]
    return(dict_medoids, family_distance_matrices)

def add_new_gene(distance_matrix, gene_list, gene):
    """
    Adds a distance matrix to
    ----------
    distance_matrix
        {fasta file name: distance matrix}
    gene_list
        list, fasta file name
    gene
        fasta file name
    returns
    ----------
    gene_list.index(gene)
        float, index number of the gene
    """
    if gene not in gene_list:
        # Add to list
        gene_list.append(gene)
        gene_list.index(gene)
        # Extend distance matrix: One new row, one new column
        for row in distance_matrix:
            row.append(0)
        distance_matrix.append([0] * len(gene_list))
    return(gene_list.index(gene))

def add_to_distance_matrix(distance_matrix, gene_list, gene1, gene2, distance):
    """
    Adds the distance of gene1-gene2 to the distance matrix and index
    ----------
    distance_matrix
        {fasta file name: distance matrix}
    gene_list
        list, fasta file name
    gene1
        fasta file name
    gene2
        fasta file name
    distance
        float, between 0 and 1
    returns
    ----------
    """
    index1 = add_new_gene(distance_matrix, gene_list, gene1)
    index2 = add_new_gene(distance_matrix, gene_list, gene2)
    distance_matrix[index1][index2] = distance
    distance_matrix[index2][index1] = distance
    return()

def writeGCFfasta(sim_dict, outdir, outfile):
    """Writes the GCFs reference fasta file and adds n_repr
    parameters
    ----------
    sim_dict
        dict, similarity dictionary for GCs
    outdir
        string, the path to output directory
    outfile
        string, name of the outfile
    returns
    ----------
    outfile = name of the GCF fasta file in outdir
    """
    infiles = sim_dict.keys()
    outfile = f"{outdir}{outfile}"
    with open(outfile, "w") as fout:
        for fkey in infiles:
            n_repr = len(sim_dict[fkey])
            fkey_contents = fkey.split("/")
            fkey_contents[-1] = fkey_contents[-1].replace("GC_PROT", "GC_DNA")
            fkey_contents[-1] = fkey_contents[-1].replace("HG_PROT", "HG_DNA")
            fname = "/".join(fkey_contents)
            with open(fname, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        line = line.strip()
                        if n_repr < 1:
                            fout.write(f"{line}--NR=1\n")
                        else:
                            fout.write(f"{line}--NR={n_repr}\n")
                    else:
                        fout.write(line)
    return (outfile)


def makefastaheadersim(sim_dict):
    """converts the file sim_dict to a fastaheader a similarity dictionary
    parameters
    ----------
    sim_dict
        dict, similarity dictionary for GCs
    returns
    ----------
    ret = {fastaheader: [fastaheaders that are similar]}
    """
    ret = {}
    infiles = sim_dict.keys()
    for fname in infiles:
        sim_fnames = sim_dict[fname]
        n_repr = len(sim_fnames)
        fkey_contents = fname.split("/")
        fkey_contents[-1] = fkey_contents[-1].replace("GC_PROT", "GC_DNA")
        fkey_contents[-1] = fkey_contents[-1].replace("HG_PROT", "HG_DNA")
        fname = "/".join(fkey_contents)
        with open(fname, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    fastaheader = f"{line[1:]}--NR={1 if n_repr < 1 else n_repr}"  # stripping '>'
                    ret[fastaheader] = []
        for sf in sim_fnames:
            fkey_contents = sf.split("/")
            fkey_contents[-1] = fkey_contents[-1].replace("GC_PROT", "GC_DNA")
            fkey_contents[-1] = fkey_contents[-1].replace("HG_PROT", "HG_DNA")
            sf = "/".join(fkey_contents)

            # sf = sf.replace("PROT", "DNA")
            with open(sf, "r") as s:
                for line in s:
                    line = line.strip()
                    if line.startswith(">"):
                        sim_fastaheader = line[1:]  # stripping '>'
                        ret[fastaheader].append(sim_fastaheader)
    return (ret)


def writejson(dictionary, outdir, outfile_name):
    """writes results in a dict to json format
    parameters
    ----------
    dictionary
        dict, dicionary containing some results (here GCFs results)
    outdir
        string, the path to output directory
    outfile_name
        string, name for the outfile
    returns
    ----------
    outfile
        name and path of the output file
    """
    outfile = f"{outdir}{outfile_name}.json"
    with open(outfile, "w") as w:
        w.write(json.dumps(dictionary, indent=4))
    return(outfile)


def applyfiltering(enzyme_locs, fastasimdict):
    """Corrects the enzyme locs for fastani filtering
    """
    ret = {}
    for cluster_NR in fastasimdict.keys():
        NR_index = cluster_NR.find("--NR")
        cluster = cluster_NR[:NR_index]
        """
        if cluster == "gb|QRHI01000014|DNA--Entryname=D-R_unassigned_test--OS=Fusobacterium_mortiferum_strain_AM25-18LB--SMASHregion=region001--NR=2":
            print(enzyme_locs[cluster])
            print("________________________________________________________________________________")
            ret[cluster_NR] = enzyme_locs[cluster]
            print(ret)
            sys.exit()
        """
        ret[cluster_NR] = enzyme_locs[cluster]
    return (ret)


######################################################################
# Writing and purging files
######################################################################
def purge(d, pattern):
    """removes files matching a pattern
    parameters
    ----------
    d
        string, directory path
    pattern
        string, regex
    returns
    ----------
    """
    for f in os.listdir(d):
        if re.search(pattern, f):
            os.remove(os.path.join(d, f))


def append_fasta(main_file, append_file):
    """appends append_file to main_file
    parameters
    ----------
    main_file
        string, name of the target file
    append_file
        string, name of the query file
    returns
    ----------
    """
    with open(main_file, "a+") as main:
        with open(append_file, "r") as append:
            for line in append:
                main.write(line)


def gbktofasta(gbkfile, fastafile, outdir):
    """copies a file to outdir
    parameters
    ----------
    gbkfile
        string, name of the genbank file (including path)
    fastafile
        string, name of the fasta file (including path)
    outdir
        string, the path to output directory
    returns
    ----------
    None
    """
    fastalocation = Path(f"{sys.path[0]}")
    SeqIO.convert(gbkfile, "genbank", fastafile, "fasta")
    try:
        os.mkdir(f"{outdir}genome_files")
    except(FileExistsError):
        pass
    # Move files into new directory
    shutil.move(fastafile, os.path.join(outdir, "genome_files"))


######################################################################
# Housekeeping genes: HMMer
#
# Here I look through every protein sequence present in the whole
# genome sequence file. I first parse every protein sequence from the
# whole genome genbank file with the DNA coordination. Then, hmmsearch
# is used search this database for housekeeping sequences using
# housekeeping genes hmms. The best hits are taken out, and since the
# coordinates were stored in the fastaheaders, these are again used to
# take out the housekeeping DNA sequences from the whole genome
# genbank files.
######################################################################
def prepareseqdb(genome_gbkfile, outdir):
    """prepares the seqdb.faa for HMMsearch
    parameters
    ----------
    genome_gbkfile
        string, name of the whole genome gbkfile
    outdir
        string, the path to output directory
    returns
    ----------
    seqdb = filename of protein sequence database
    """
    # parsing all protein sequences: Note that the genomic genbank
    # files are still high quality drafts, meaning that they contain
    # the scaffolds still and are not fully combined already.
    # recordnumber = 0
    scaff_number = -1
    gbkcontents = SeqIO.parse(genome_gbkfile, "genbank")
    with open(f"{outdir}seqdb.faa", "w") as db:
        for record in gbkcontents:
            # recordnumber += 1
            scaff_number = int(record.name[-3:])
            for feature in record.features:
                if feature.type == "CDS":
                    if "translation" in feature.qualifiers.keys():
                        protseq = feature.qualifiers['translation'][0]
                        loc = feature.location
                        db.write(f">{scaff_number},{loc.start},{loc.end},{loc.strand}\n{protseq}\n")
    return (f"{outdir}seqdb.faa")


def hmmsearch(seqdb, hmmfiles, hmmerpath, outdir):
    """searches the seqdb for housekeeping genes profiles
    seqdb
        string, name of the seqdb file
    returns
    ----------
    hmmoutput = filename of parsable hmmsearch output
    """
    hmmoutput = f"{outdir}hmmoutput.result.txt"
    try:
        cmd_hmmsearch = f"{hmmerpath if hmmerpath else ''}hmmsearch --tblout {hmmoutput} {hmmfiles} {seqdb}"
        res_hmmsearch = subprocess.check_output(cmd_hmmsearch, shell=True)
    except(subprocess.CalledProcessError):
        pass  # raise error here
    return (hmmoutput)


def parsehmmoutput(hmmresult):
    """parses the hmm output for housekeeping gene locations
    parameters
    ----------
    hmmresult
        filename of parsable hmmsearch output
    returns
    ----------
    genelocs = dict, {gene:location}
    """
    genelocs = {}
    with open(hmmresult, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                els = line.split()
                loc = els[0]
                gene = els[2]
                if not gene in genelocs:
                    genelocs[gene] = loc
    return (genelocs)


def getprotseqfromdb(seqdb, locheader):
    """Parse the HMM hit protein sequence from the seqdb
    parameters
    ----------
    seqdb
        filename, name of the seqdb file used for HMMering
    locheader
        >{recordnumber}_{loc}, header for each protein sequence
    returns
    ----------
    ret = protein sequence of housekeeping gene
    """
    ret = ""
    hit = -10
    with open(seqdb, "r") as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if locheader in line:
                hit = idx
            if idx == hit + 1:
                ret = line
    return (ret)


def getgenefromgbk(gbkfile, location):  # change to work with locations
    """parses a genesequence from a gbk file using the gene location
    parameters
    ----------
    gbkfile
        string, path to gbk file + file
    location
        string of coordinates, example: "[start:end>](+)"
    returns
    ----------
    ret = DNA sequence of housekeepinggene from featurelocation
          coordinates
    abs_loc = validation, contains the location of HG on specific
              scaffold. [scaffold, start, end]
    """
    ret = ""
    scaff_number, start, end, strand = location.split(",")
    scaff_number = int(scaff_number)

    # Making the FeatureLocation
    f_start = BeforePosition(start.strip("<")) if "<" in start else ExactPosition(start)
    f_end = AfterPosition(end.strip(">")) if ">" in end else ExactPosition(end)
    f = FeatureLocation(f_start, f_end, int(strand))

    gbkcontents = SeqIO.parse(gbkfile, "genbank")
    for record in gbkcontents:
        scaff_check = int(record.name[-3:])  # = scaffold number
        if scaff_check == scaff_number:
            DNA = record.seq
    ret = f.extract(DNA)  # The DNA sequence of the housekeepinggene

    # VALIDATION
    start = start.replace(">", "")
    start = start.replace("<", "")
    start = int(start)
    end = end.replace(">", "")
    end = end.replace("<", "")
    end = int(end)
    abs_loc = [scaff_number, start, end]
    return (ret, abs_loc)

######################################################################
# Preparing BiG-SCAPE
######################################################################

def movegbk(outdir, gbk_file, listoffiles):
    """Retrieves the .gbk files from the input directory including their paths.
    ----------
    outdir
        string, path to the output directory
    gbk_file
        string, path to the gut/antiSMASH dir
    listoffiles
        list, list of filenames to move
    returns
    ----------
    """
    # Make .gbk folder
    try:
        os.mkdir(f"{outdir}gbk_files")
    except(FileExistsError):
        pass
    # Copy .gbk into new directory
    path = os.path.join(f"{outdir}gbk_files")
    genome = gbk_file.split("/")[-1]
    for filename in listoffiles:
        if filename == genome:
            try:
                shutil.copy(gbk_file, os.path.join(path, genome))
            except(FileExistsError):
                pass
        else:
            pass
    return


def list_representatives(GCF_dict):
    """retrieves the represenatives
    ----------
    GCF_dict
        dictionary, {Representative: family members}
    returns
    ----------
    outlist: list, list of filenames
    """
    outlist = []

   # with open(json_file, "r") as jfile:
   #     GCF_dict = json.load(jfile)

    for key in GCF_dict.keys():
        GC_or_HG = key.split("|")[2].split("_")[0]
        if GC_or_HG == "GC":
            orgID = key.split("|")[1]
            filename = f"{orgID}.gbk"
            outlist.append(filename)
        else:
            pass
    return (outlist)


######################################################################
# Running BiG-SCAPE
######################################################################
def run_bigscape(path_to_bigscape, pfamfiles, outdir, cores, cut_off):
    """Runs BiG-SCAPE
    Parameters
    ----------
    path_to_bigscape
        string, path to bigscape
    location_of_pfam
        string, path to pfam files
    outdir
        string, path to output directory
    cores
        int, number of threads
    returns
    ----------
    """
    output_dir = f"{outdir}bigscape_output"
    gbk_files = os.path.join(f"{outdir}gbk_files")
    try:
        cmd_bigscape = f"python3 {path_to_bigscape}bigscape.py -i {gbk_files} -o {output_dir} -c {cores} \
        --pfam_dir {pfamfiles} --cutoffs {cut_off} --clans-off --hybrids-off --mibig"
        subprocess.check_output(cmd_bigscape, shell=True)
    except(subprocess.CalledProcessError):
        pass  # raise error here
    return


######################################################################
# Parsing BiG-SCAPE results
######################################################################
def retrieve_bigscapefiles(outdir):
    """Retrieves the .tsv files from the output directory including their paths.
    parameters
    ----------
    outdir
        string, path to the BiG-SCAPE output directory
    returns
    ----------
    """
    network_dir = F"{outdir}bigscape_output/network_files"
    try:
        for dirpath, dirnames, files in os.walk(network_dir):
            for filename in files:
                if filename.endswith(".tsv") and "clustering" in filename:
                    yield (os.path.join(dirpath, filename))
    except:
        pass

def parse_bigscape_result(bigscape_result):
    """parses the BiG-SCAPE results to a dictionary
    parameters
    ----------
    bigscape_result
        string, path and name of the .tsv file
    returns
    ----------
    GCFs = dict, {family orgID: family members}
    """
    GCFs = {}
    family_to_key = {}
    with open(bigscape_result, "r") as bigscape_file:
        for line in bigscape_file:
            if not line.startswith("#"):
                gene_lists = line.split()
                gene = gene_lists[0]
                familynr = gene_lists[1]
                if familynr not in family_to_key.keys():
                    family_to_key[familynr] = gene
                    key = gene
                    GCFs[key] = []
                else:
                    key = family_to_key[familynr]
                GCFs[key].append(gene)
    return (GCFs)

def make_jsondict(new_GCFs, old_GCFs):
    """
    creates a dictionary containing the BiG-SCAPE families with the full organism names
    ----------
    new_GCFs
        dictionary, {bigscape family ID number: bigscape family ID numbers}
    old_GCFs
        dictionary, {full organism name: full organism names}
    returns
    ----------
    json_dict: dict, {bigscape family full ID: bigscape family full ID}
    """
    json_dict = {}
    full_names = {}

    for key_old in old_GCFs.keys():
        GC_or_HG = key_old.split("|")[2].split("_")[0]
        if GC_or_HG == "GC":
            old_orgID = key_old.split("|")[1]
            # make dict with the full org names as values and the new orgIDs as keys
            for gene_names in new_GCFs.values():
                for new_orgID in gene_names:
                    if old_orgID == new_orgID:
                        full_names[new_orgID] = key_old
            # make a dict with the full org names as keys
            if old_orgID in new_GCFs.keys():
                json_dict[key_old] = []
        if GC_or_HG == "HG":
            json_dict[key_old] = [key_old]

    # append the full org names as values to the corresponding keys
    for keys in json_dict.keys():
        valueID = keys.split("|")[1]
        for values in new_GCFs.values():
            if values[0] == valueID:
                for value_item in values:
                    for key_names in full_names.keys():
                        if value_item == key_names:
                            json_dict[keys].append(full_names[key_names])
    return (json_dict)

def pickle_files(GCF_dict, fasta_file, bed_file, outdir, BGCF_dict = ""):
    """Pickle the output files for later use
    ----------
    GCF_dict
        dictionary, {family name: family members}
    fasta_file
        >fasta header, fasta sequence
    bed_file
        contains the cluster ID and the start and stop position
    outdir
        string, path to the output directory
    returns
    ----------
    """
    fasta_dict = {}
    bed_list = []
    outfile = open(f"{outdir}BiG-MAP.pickle", "wb")
    with open (fasta_file, "r") as fasta_file:
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            fasta_dict[seq.id] = seq.seq
    with open (bed_file, "r") as bed_files:
        for line in bed_files:
            bed_list.append(line)

    pickle.dump(fasta_dict, outfile)
    pickle.dump(GCF_dict, outfile)
    pickle.dump(BGCF_dict, outfile)
    pickle.dump(bed_list, outfile)

    outfile.close()
    return

######################################################################
# MAIN
######################################################################
def main():
    """
    STEPS
    ----------
    1) Parsing all the gene clusters (GCs) protein and DNA sequences to fasta
    2) Redundancy filtering the GCs on the protein level to obtain GC representatives
    3) Finding corresponding housekeeping genes (HGs) and parsing to fasta for prot and DNA
    4) Redundancy filtering the HGs on the protein level to obtain HG representatives
    5) Using the representative GCs and HGs to create one big reference.fna
    6) Optionally run second round of redundancy filtering
    """
    args = get_arguments()

    print("___________Extracting fasta files___________________")

    ################################
    # parsing data to fasta files
    ################################
    GC_enzyme_locs = {}  # storing the core + flanking genes locations
    genomedict = {}  # will be used to extract housekeeping genes
    GC_prot_files = []  # will be used as fastANI input
    GC_DNA_files = []  # will be used as fastANI input
    absolute_locs = {}  # VALIDATION
    for d in args.indir:
        for f in retrieveclusterfiles(d):
            # Parsing each .gbk file
            DNAseq, AAseq, GC, organism, core_locs, abs_locs = parsegbkcluster(f, args.flank_genes)
            # writing the full gene clusters to .fasta for fastANI
            prot_file, orgID, fasta_header = writefasta(AAseq, "GC_PROT", GC, organism, f, args.outdir)
            dna_file, orgID, fasta_header_DNA = writefasta(DNAseq, "GC_DNA", GC, organism, f, args.outdir)
            GC_prot_files.append(prot_file)
            GC_DNA_files.append(dna_file)
            # remembering the whole genome gbk files
            genomegbk = getgenomegbk(f)
            genomedict[orgID] = genomegbk
            # remembering the enzymatic genes locations in dictionary
            GC_enzyme_locs[fasta_header_DNA] = core_locs
            # Remembring the absolute cluster locations for VALIDATION:
            if not orgID in absolute_locs:
                absolute_locs[orgID] = []
            absolute_locs[orgID].append({fasta_header_DNA: abs_locs})

    ################################
    # Mash: similarity
    ################################

    make_sketch(args.outdir, args.threads)

    #checks the output of the mash sketch
    reruns = 0
    total_reruns = 25
    for i in range(total_reruns):
        with open (f"{args.outdir}log.file", "r") as log_file:
            if "ERROR" in log_file.read():
                make_sketch(args.outdir, args.threads)
                reruns += 1
                print("Encountered error in making sketch file. Retry attempt:", reruns)
                if reruns == total_reruns:
                    sys.exit("Maximum number of reruns is reached. Try increasing the number of reruns, or decrease the amount of input data.")
            else:
                break

    calculate_distance(args.outdir)


    ###############################
    # LSH buckets clustering
    ###############################
    GCFs, distance_matrix = calculate_medoid(args.outdir, args.treshold_GC)

    print("___________Adding housekeeping genes________________")

    ################################
    # housekeeping genes: Obtainging
    ################################
    hmmfile = Path(f"{sys.path[0]}") / "pfamdomains/combined.hmm"
    processed = []
    HG_prot_files = []
    HG_DNA_files = []
    absolute_locs["hgenes"] = []  # VALIDATION
    for fname in GCFs.keys():
        orgID = fname.split("/")[-1].split(".")[0][8:]  # was [5:] <-- the last one
        organism = ".".join(fname.split("/")[-1].split(".region")[0].split(".")[1:])

        # find housekeeping genes for org using HMMer
        if organism not in processed:
            seqdb = prepareseqdb(genomedict[orgID], args.outdir)
            hmmresult = hmmsearch(seqdb, hmmfile, "", args.outdir)
            genelocs = parsehmmoutput(hmmresult)
            for gene, loc in genelocs.items():
                prot_seq = getprotseqfromdb(seqdb, loc)
                DNA_seq, abs_hloc = getgenefromgbk(genomedict[orgID], loc)
                f_prot, ID_prot, housekeepingheader = writefasta(prot_seq, "HG_PROT", gene, organism, fname,
                                                                 args.outdir)
                f_DNA, ID_DNA, housekeepingheader = writefasta(DNA_seq, "HG_DNA", gene, organism, fname, args.outdir)
                HG_prot_files.append(f_prot)
                HG_DNA_files.append(f_DNA)

                # Remembering the locations
                GC_enzyme_locs[housekeepingheader] = [FeatureLocation(0, len(DNA_seq))]
                absolute_locs["hgenes"].append({housekeepingheader: abs_hloc})
            print(f"__________DONE: {organism}")

            # Convert genome.gbk --> genome.fasta and store in outdir
            if args.genomefiles:
                gbktofasta(genomedict[orgID], f"{args.outdir}{organism}.fasta", args.outdir)
        processed.append(organism)

    ################################
    # housekeeping genes: Redundancy filtering
    ################################

    make_sketch(args.outdir, args.threads, "HG")
    reruns = 0
    total_reruns = 25
    for i in range(total_reruns):
        with open (f"{args.outdir}log.file", "r") as log_file:
            if "ERROR" in log_file.read():
                make_sketch(args.outdir, args.threads, "HG")
                reruns += 1
                print("Encountered error in making sketch file. Retry attempt:", reruns)
                if reruns == total_reruns:
                    sys.exit("Maximum number of reruns is reached. Try increasing the number of reruns, or decrease the amount of input data.")
            else:
                break
    calculate_distance(args.outdir, "mash_output_HG.tab")

    GCFs_ALL, distance_matrix_ALL = calculate_medoid(args.outdir, args.treshold_HG, GCFs, "mash_output_HG.tab")

    fastadict_ALL = makefastaheadersim(GCFs_ALL)

    # Writing results to outdir
    #writejson(GCFs_ALL, args.outdir, "BiG-MAP.GCF_HGF")
    writejson(fastadict_ALL, args.outdir, "BiG-MAP.GCF_HGF")
    writejson(distance_matrix, args.outdir, "BiG-MAP.dist_GC")
    writejson(distance_matrix_ALL, args.outdir, "BiG-MAP.dist_HG")
    fasta_file = writeGCFfasta(GCFs_ALL, args.outdir, "BiG-MAP.GCF_HGF.fna")

    # Remembering the enzymatic gene locations
    GCF_enzyme_locs = applyfiltering(GC_enzyme_locs, fastadict_ALL)
    bed_file = locs2bedfile(GCF_enzyme_locs, f"{args.outdir}BiG-MAP.GCF_HGF.bed")
    writejson(absolute_locs, args.outdir, "absolute_cluster_locs")  # VALIDATION

    ###################################
    # Run second redundancy filtering
    ###################################
    if args.bigscape_path:
        list_gbkfiles = list_representatives(fastadict_ALL)
        print("________Preparing BiG-SCAPE input__________")
        for files in args.indir:
            for gbk_file in retrieveclusterfiles(files):
                movegbk(args.outdir, gbk_file, list_gbkfiles)
        print("__________Running BiG-SCAPE________________")
        run_bigscape(args.bigscape_path, args.pfam, args.outdir, args.threads, args.cut_off)

        for tsv_file in retrieve_bigscapefiles(args.outdir):
            parsed = parse_bigscape_result(tsv_file)

        dict_json = make_jsondict(parsed, fastadict_ALL)
        writejson(dict_json, args.outdir, "BiG-MAP.GCF")

    #############################
    # Pickling
    #############################
    if args.bigscape_path:
        pickle_files(fastadict_ALL, fasta_file, bed_file, args.outdir, dict_json)
    else:
        pickle_files(fastadict_ALL, fasta_file, bed_file, args.outdir)

    #############################
    # Cleaning output dir
    #############################
    purge(args.outdir, ".fasta")
    purge(args.outdir, ".txt")
    purge(args.outdir, ".faa")
    purge(args.outdir, "log.file")


if __name__ == "__main__":
    main()
