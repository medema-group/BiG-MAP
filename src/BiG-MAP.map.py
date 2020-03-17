#! /usr/bin/env python3

"""
--------------- Mapping module ---------------
Author: Koen van den Berg & Hannah Augustijn
University: Wageningen University and Research
Department: Department of Bioinformatics
Date: 21/01/2019
----------------------------------------------

The purpose of this script is to map the metagenomic and
metatranscriptomic samples to the fasta database that has been created
by module 2. This will allow to find the abundance and expression for
the found metabolic gene clusters and biosynthetic gene clusters by
gutSMASH and antiSMASH respecively. The core of this part of the
pipeline will consist of bowtie2 which, according to my BSC thesis,
performs most optimal using the sensitive-local setting. 
"""

# Import statements:
import os
import subprocess
from sys import argv
import sys
import argparse
from pathlib import Path
import json
import pandas as pd
import shutil
import re
import textwrap
import pickle

# Functions:
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
    usage='''
______________________________________________________________________
     BiG-MAP map: maps the paired reads to the predicted MGCs
______________________________________________________________________
Generic command: python3 BiG-MAP.map.py {-I1 [mate-1s] -I2 [mate-2s] | -U [samples]} -O [outdir] -F [family] [Options*]
Maps the metagenomic/metatranscriptomic reads to the fasta reference
file and outputs RPKM read counts in .csv and BIOM format. Use
BiG-MAP_process conda environment.
Data inputs: either paired or unpaired
    -I1   Provide the mate 1s of the paired metagenomic and/or
          metatranscriptomic samples here. These samples should be
          provided in fastq-format (.fastq, .fq, .fq.gz). Also, this 
          can be a space seperated list from the command line.
    -I2   Provide the mate 2s of the paired metagenomic and/or
          metatranscriptomic samples here. These samples should be
          provided in fastq-format (.fastq, .fq, .fq.gz). Also, this 
          can be a space seperated list from the command line.
    -U    Provide the unpaired metagenomic/metatranscriptomic samples
          here. These samples should be provided in fastq-format
          (.fastq, .fq, .fq.gz). Also, this can be a space seperated
          list from the command line.
File inputs: either separated or pickled:
    -F    Directory with all the output files from the family module 
    -P    Input files are in pickled format (named: BiG-MAP.[name].pickle). 
          The format of the pickled file: fasta file, GCF json file, and 
          optionally a bed file and/or BiG-SCAPE GCF dictionary.

Obligatory arguments:
    -O    Name of the output directory for where the output files are going 
          to be written. Default = current folder (.)
Options:
    -b    Outputs the resulting read counts in biom format (v1.0) as
          well. This will be useful to analyze the results in
          BiG-MAP.analyse. Therefore, it  is important to include
          the metadata here as well: this metagenomical data should
          be in the same format as the example metadata
    -f    Input files are in fasta format (.fna, .fa, .fasta): True/False. 
          Default = False.
    -s    Bowtie2 setting: 
          END-TO-END mode: very-fast, fast, sensitive, very-sensitive
          LOCAL mode: very-fast-local, fast-local, sensitive-local, 
          very-sensitive-local. DEFAULT = fast
    -th   Number of used threads in the bowtie2 mapping step. Default = 6
______________________________________________________________________
''')
    parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-I1","--fastq1", nargs='+',help=argparse.SUPPRESS, required=False)
    parser.add_argument("-I2","--fastq2",nargs='+',help=argparse.SUPPRESS, required = False)
    parser.add_argument("-U","--U_fastq",nargs='+',help=argparse.SUPPRESS, required = False)
    parser.add_argument("-F", "--family", help=argparse.SUPPRESS, required=False)
    parser.add_argument("-P", "--pickle_file", help=argparse.SUPPRESS, required=False)
    parser.add_argument( "-b", "--biom_output",
                         help=argparse.SUPPRESS, type=str, required = False)
    parser.add_argument( "-f", "--fasta", help=argparse.SUPPRESS,
                         type=str, required = False, default=False)
    parser.add_argument( "-s", "--bowtie2_setting", help=argparse.SUPPRESS,
                         type=str, required = False, default="fast")
    parser.add_argument( "-th", "--threads", help=argparse.SUPPRESS,
                         type=int, required = False, default=6)
    return(parser, parser.parse_args())

######################################################################
# Functions for mapping the reads against GCFs and % aligned
######################################################################
def bowtie2_index(reference, outdir):
    """indexes the fasta reference file
    parameters
    ----------
    reference
        string, the name of the reference fasta file (GCFs)
    outdir
        string, the path of the output directory
    returns
    ----------
    index name = the name of the built bowtie2 index
    """
    try:
        stem = Path(reference).stem
        index_name = os.path.join(outdir, stem)
        if not os.path.exists(index_name + ".1.bt2"):
            cmd_bowtie2_index = f"bowtie2-build {reference} {index_name}"
            res_index = subprocess.check_output(cmd_bowtie2_index, shell=True)
    except(subprocess.CalledProcessError):
        print("Error-code M3:001, check error table")
        # Proper error here, also exit code
    return (index_name)


def bowtie2_map(outdir, mate1, mate2, index, fasta, bowtie2_setting, threads):
    """Maps the .fq file to the reference (fasta)
    parameters
    ----------
    outdir
        string, the path of the output directory
    mate1
    mate2
    index
        string, the stemname of the bowtie2 index
    fasta
        Boolean, is the input fasta?
    bowtie2_setting
        string, the build-in setting for bowtie2
    threads
        int, number of threads used in the alignment
    returns
    ----------
    samfile = the .sam filename that contains all the results
    writes the mapping percentage to bowtie2_log.txt
    """
    stem = Path(mate1).stem
    sample = stem.split("_")[0]
    samfile = os.path.join(outdir, sample + ".sam")
    # In the case of unpaired, m1 and m2 are identical. Thus the following works:
    sample_command = f"-U {mate1}" if mate1 == mate2 else f"-1 {mate1} -2 {mate2}"
    if fasta == True:
        cmd_bowtie2_map = f"bowtie2 --{bowtie2_setting} --no-unal --threads {threads} -x {index} {sample_command} -S {samfile} -f"
    else:
        cmd_bowtie2_map = f"bowtie2 --{bowtie2_setting} --no-unal --threads {threads} -x {index} {sample_command} -S {samfile}"
    try:
        if not os.path.exists(samfile):

 #           print(f"the following command will be executed by bowtie2:\n\
#_____________________________________________________\n\
#{cmd_bowtie2_map}\n\
#_____________________________________________________\n")
            print(f"  Dealing with sample {sample}")
            res_map = subprocess.check_output(cmd_bowtie2_map, shell=True, stderr=subprocess.STDOUT)
            # Saving mapping percentage:
            with open(os.path.join(outdir, "bowtie2_log.txt"), "a+") as f:
                f.write(f"#{sample}\n{res_map.decode('utf-8')}")
    #except(subprocess.CalledProcessError):
    except:
        print('Unable to run bowtie')
    return (samfile)


def parse_perc(outdir):
    """parses the percentage from the bowtie2 stdout
    parameters
    ----------
    outdir
        string, the path of the output directory
    returns
    ----------
    ret = dictionary containing the mapping % for each sample
    """
    ret = {}
    sample = ""
    infile = os.path.join(outdir, "bowtie2_log.txt")
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                sample = line[1:]
            if "overall" in line:
                perc = line.split(" ")[0][:-1]
                perc = float(perc) / 100
                ret[sample] = [perc]
    return (ret)

def unpickle_files(pickled_file, outdir):
    """Extracts the files from the pickled_file
    parameters
    ----------
    pickled_files
        pickled file, contains the reference in
        fasta format, GCF family dict, optionally BGCF dict
        and bed file
    Returns
    ----------
    fasta_file, reference in fasta format
    GCF_dict, {family name: family members}
    BGCF_dict, {family name: family members}
    bed_file, orgID, loc1, loc2
    """
    f = open(pickled_file, "rb")
    bed_file = os.path.join(outdir, "BiG-MAP.GCF_HGF.bed")
    fasta_file = os.path.join(outdir, "BiG-MAP.GCF_HGF.fna")

    fasta = pickle.load(f)
    GCF_dict = pickle.load(f)
    BGCF_dict = pickle.load(f)
    bed = pickle.load(f)
    f.close()

    with open (bed_file, "w") as bed_file:
        for line in bed:
            bed_file.write(line)
    with open (fasta_file, "w") as fasta_file:
        for key in fasta.keys():
            sequence = fasta[key]
            fasta_seq = f">{key}\n{sequence}\n"
            fasta_file.write(fasta_seq)

    return(fasta_file.name, GCF_dict, BGCF_dict, bed_file.name)


######################################################################
# Functions for reading SAM and BAM files
######################################################################
def samtobam(sam, outdir):
    """converts .sam to .bam using samtools view
    parameters:
    ----------
    sam
        string, name of the outputted bowtie2 mapping
    outdir
        string, the path of the output directory
    returns
    ----------
    bamfile = the name of the .bam file
    """
    stem = Path(sam).stem
    bamfile = os.path.join(outdir, stem + ".bam")
    try:
        cmd_samtobam = f"samtools view\
        -b {sam}\
        > {bamfile}"
        res_samtobam = subprocess.check_output(cmd_samtobam, shell=True)
    except(subprocess.CalledProcessError):
        print("Unable to convert SAM file to BAM")
    return (bamfile)


def sortbam(bam, outdir):
    """sorts the bam file
    parameters
    ----------
    bam
        string, the name of the accession bamfile, ".bam"-file
    outdir
        string, the path of the output directory
    returns
    ----------
    sortedbam = name of the sorted bam file
    """
    stem = Path(bam).stem
    sortedbam = os.path.join(outdir, stem + ".sorted.bam")
    try:
        cmd_sortbam = f"samtools sort {bam} > {sortedbam}"
        res_sortbam = subprocess.check_output(cmd_sortbam, shell=True)
    except(subprocess.CalledProcessError):
        print('Unable to sort BAM file')
    return (sortedbam)


def indexbam(sortedbam, outdir):
    """Builds a bam index
    parameters
    ----------
    sortedbam
        string, the name of the sorted bam file
    outdir
        string, the path of the output directory
    returns
    ----------
    none
    """
    try:
        cmd_bam_index = f"samtools index {sortedbam}"
        res_index = subprocess.check_output(cmd_bam_index, shell=True)
    except(subprocess.CalledProcessError):
        print("Unable to build index file")
    return ()


def countbam(sortedbam, outdir):
    """calculates the raw counts from a BAM index
    parameters
    ----------
    sortedbam
        string, the name of the sorted bam file
    outdir
        string, the path of the output directory
    returns
    ----------
    counts_file = file containing the counts
    """
    counts_file = f"{sortedbam[:-3]}count"
    try:
        cmd_count = f"samtools idxstats {sortedbam} > {counts_file}"
        res_count = subprocess.check_output(cmd_count, shell=True)
    except(subprocess.CalledProcessError):
        print('Unable to calculate raw counts from BAM')
    return (counts_file)

def correct_counts(countsfile, family):
    """Corrects the number of counts for the BiG-SCAPE families
    ----------
    countsfile
        string, the name of the sorted counts file
    family
        json, {HGF/GCF representative: HGF/GCF members}
    returns
    ----------
    corrected_countsfile = file containing the counts
    """
    corrected_countsfile = f"{countsfile[:-12]}corrected.count"

    with open(corrected_countsfile, "w") as counts_adj:
        with open(countsfile, "r") as counts:
            for line in counts:
                cluster, length, nreads, nnoreads = line.strip().split("\t")
                for key in family.keys():

                    # If the BiG-SCAPE family is larger than 1, adjust the number of family members
                    if cluster == key and len(family[key]) > 1:
                        adj_key, number = correct_family_size(cluster, key, nreads, family)
                        counts_adj.write(f"{adj_key}\t{length}\t{number}\t{nnoreads}\n")

                    # The BiG-SCAPE family size is equal to 1, no correction is needed
                    elif cluster == key and len(family[key]) == 1:
                        counts_adj.write(f"{key}\t{length}\t{nreads}\t{nnoreads}\n")

                    # The cluster is already in another family
                    else:
                        pass
    return(corrected_countsfile)

def correct_family_size(cluster, key, number_reads, GCF_dict):
    """
    calculate the total sum of counts and number of family members of a BiG-SCAPE GCF family
    ----------
    cluster
        string, the name of the cluster
    key
        string, the name of the key
    number_reads
        int, the number of reads
    GCF_dict
        dictionary, {HGF/GCF representative: HGF/GCF members}
    returns
    ----------
    adj_key = adjusted key for total number of family members
    read_total = total sum of reads
    """
    fam_size = []
    read_total = 0
    for name in GCF_dict[key]:
        fam_size.append(float(name.split("--")[-1].split("=")[-1]))
        total_fam_size = int(sum(fam_size))
        if cluster in GCF_dict[key]:
            read_total += int(number_reads)
    adj_key = name.split("NR=")[0]
    adj_key = f"{adj_key}NR={total_fam_size}"

    return(adj_key, read_total)



def extractcorefrombam(bam, outdir, bedfile):
    """extracts regions in bedfile format from bam file
    parameters:
    ----------
    bam
        string, the name of the accession bamfile, ".bam"-file
    outdir
        string, the path of the output directory
    bedfile
        the name of the bedfile with core coordinates
    returns
    ----------
    bamfile = the name of the .bam file
    """
    #  samtools view -b -L BiG-MAP.enzymes.bedfile SRR5947807.bam > bedfile.bam
    bamstem = Path(bam).stem
    bamfile = os.path.join(outdir, "core_" + bamstem + ".bam")
    if os.path.exists(bedfile):
        try:
            cmd_extractcore = f"samtools view\
            -b {bam}\
            -L {bedfile}\
            > {bamfile}"
            res_extractcore = subprocess.check_output(cmd_extractcore, shell=True)
        except(subprocess.CalledProcessError):
            print('Unable to extract core locations frmom bedfile')
    else:
        # raise bedfile error here!!!
        pass
    return (bamfile)


######################################################################
# RPKM and TPM counting
######################################################################
def calculateTPM(countsfile):
    """Calculates the TPM values for a sample
    TPM = rate/sum(rate) * 10^6
    rate = nreads/cluster_length (kb)
    parameters
    ----------
    counts_file
        file containing the counts
    core
        bool, skip housekeeping genes
    returns
    ----------
    TPM = dictionary containing TPM counts per cluster
    """
    rates = {}
    ratesum = 0
    with open(countsfile, "r") as f:
        for line in f:
            line = line.strip()
            cluster, length, nreads, nnoreads = line.split("\t")
            if "NR" in cluster:
                NR = float(cluster.split("--")[-1].split("=")[-1])
                nreads = float(nreads) / NR
            try:
                rate = float(nreads) / float(length)
                rates[cluster] = rate
                ratesum += rate
            except(ZeroDivisionError):
                print('Unable to calculate TPM counts')
                pass
    TPM = {}
    for key in rates:
        try:
            TPM[key] = rates[key] / ratesum
        except(ZeroDivisionError):
            TPM[key] = 0
    return (TPM)


def calculateRPKM(countsfile):
    """Calculates the RPKM values for a sample
    RPKM = read_counts/(cluster_length * sum(read_counts)) * 10^9
    parameters
    ----------
    counts_file
        file containing the counts
    returns
    ----------
    RPKM = dictionary containing RPKM counts per cluster
    """
    sum_reads = 0
    read_counts = {}
    cluster_lengths = {}
    with open(countsfile, "r") as f:
        for line in f:
            if "*" not in line:
                line = line.strip()
                cluster, length, nreads, nnoreads = line.split("\t")
                if "NR" in cluster:
                    NR = float(cluster.split("--")[-1].split("=")[-1])
                    nreads = float(nreads) / NR
                    read_counts[cluster] = nreads
                else:
                    read_counts[cluster] = float(nreads)
                cluster_lengths[cluster] = float(length)
                sum_reads += float(nreads)

    RPKM = {}
    for key in read_counts:
        try:
            RPKM[key] = read_counts[key] / (sum_reads * cluster_lengths[key]) * 1000000000

        except(ZeroDivisionError):
            RPKM[key] = 0
    return (RPKM)


def parserawcounts(countsfile):
    """parses the raw counts from a countsfile
    parameters
    ----------
    counts_file
        file containing the counts
    returns
    ----------
    raw_counts = dictionary containing raw counts per cluster
    """
    raw_counts = {}
    with open(countsfile, "r") as f:
        for line in f:
            if "*" not in line:
                line = line.strip()
                cluster, length, nreads, nnoreads = line.split("\t")
                if "NR" in cluster:
                    NR = float(cluster.split("--")[-1].split("=")[-1])
                    raw_counts[cluster] = float(nreads) / NR
                else:
                    raw_counts[cluster] = float(nreads)
    return (raw_counts)


######################################################################
# Functions for analysing coverage with Bedtools genomecov
######################################################################
def preparebedtools(outdir, reference):
    """makes the -g genome.file for bedtools from a fasta file. This
    file is formatted as follows:
    fastaheader, length
    hseq1, len(seq1)
     :       :
    hseqn, len(seqn)
    parameters
    ----------
    outdir
        string, the path of the output directory
    reference
        string, the name of the reference fasta file (GCFs)
    returns
    ----------
    genome_file = name of the by bedtools required genome file
    """
    c = ""
    genome_file = os.path.join(outdir, "genome.file")
    with open(genome_file, "w") as w:
        with open(reference, "r") as f:
            for line in f:
                name, length, no1, no2 = line.strip().split("\t")
                w.write(f"{name}\t{length}\n")

    return (genome_file)


def bedtoolscoverage(gfile, outdir, sortedbam):
    """computes the coverage for each mapped region
    parameters
    ----------
    bedtoolspath
        path to the bedtools installation
    gfile
        genome file build by preparebedtools()
    outdir
        string, the path of the output directory
    sortedbam
        name of the sorted bam file
    returns
    ----------
    bg_file = the name of the bedgraph file
    """
    stem = Path(sortedbam).stem
    bg_file = os.path.join(outdir, stem.split('.')[0] + ".bg")

    try:
        cmd_bedtools = f"bedtools genomecov -bga -ibam {sortedbam} > {bg_file}"
        res_bedtools = subprocess.check_output(cmd_bedtools, shell=True)
    except(subprocess.CalledProcessError):
        print('Unable to compute the coverage')
    return (bg_file)


def computetotalcoverage(bgfile):
    """computes the total coverage of a gene cluster from a .bg file
    parameters
    ----------
    bgfile
        name of the bedgraph file
    returns
    ----------
    total_coverage = {cluster: totalcov}
    """
    nocov = {}
    clusterlen = {}
    with open(bgfile, "r") as f:
        for line in f:
            line = line.strip()
            cluster, start, end, cov = line.split("\t")
            clusterlen[cluster] = float(end)  # last encounter is length
            if not cluster in nocov:  # make entry
                nocov[cluster] = 0
            if float(cov) == 0:  # enter no coverage values
                nocov[cluster] += (float(end) - float(start))
    total_coverage = {}
    for key in nocov.keys():
        perc = (clusterlen[key] - nocov[key]) / clusterlen[key]
        # Set treshold here!!!
        total_coverage[key] = perc
    return (total_coverage)

def correct_coverage(coverage, reference):
    """
    """
    outdict = {}
    with open (reference, "r") as ref:
        for line in ref:
            name, no1, no2, no3 = line.strip().split("\t")
            for key in coverage.keys():
                orgname, number = key.split("NR=")
                if orgname in name:
                    outdict[name] = coverage[key]
    return(outdict)


def computecorecoverage(bedgraph, bedfile):
    """computes the core "enzymatic" coverage for gene clusters
    EXPLANATION:
    This computation is based on the bedfile that contains the
    coordinates for the enzymatic genes. What happens is that the
    algorithm finds entries that have 0 coverage, and compares if the
    regions of these entries are within the enzymatic core
    regions. Based on this, certain computations will be made as
    defined in the local function. In the end, the amount of bases for
    which 0 coverage was found are added up, and then substracted from
    the total length of the core:
    core_coverage = (length_core - bases_not_covered)/length_core
    EXAMPLE:
    say that of a core of length 3000 340 bases are not covered, then:
    core_coverage = (3000-340)/3000 = 0.887
    parameters
    ----------
    bedgraph
        name of the bedgraph file
    bedfile
        the name of the bedfile with core coordinates
    returns
    ----------
    core_coverage = dict, {cluster: corecov}
    """

    def local_computecov(start_list, end_list, local_entry):
        """Computes the local cov value
        Ls = local entry start
        Le = local entry end
        Ts = true start of core gene
        Te = true end of core gene
        parameters
        ----------
        start_list
            list, all the start coords for the enzymatic core genes
        end_list
            list, all the start coords for the enzymatic core genes
        local_entry
            list, [start, end] of local entry
        returns
        ----------
        local_cov = int, value for not covered bases for entry
        """
        ret_cov = 0
        Ls = local_entry[0]
        Le = local_entry[1]
        for Ts, Te in zip(start_list, end_list):
            if Ls > Ts and Ls < Te and Le > Ts and Le < Te:
                # in
                ret_cov += Le - Ls
            if Ls < Ts and Ls < Te and Le > Ts and Le < Te:
                # start out
                ret_cov += Le - Ts
            if Ls > Ts and Ls < Te and Le > Ts and Le > Te:
                # end out
                ret_cov += Te - Ls
            if Ls < Ts and Ls < Te and Le > Ts and Le > Te:
                # start&end out
                ret_cov += Te - Ts
        return (ret_cov)

    # Parsing bedfile for Ts, Te, and core lengths:
    core_starts = {}
    core_ends = {}
    core_lengths = {}
    with open(bedfile, "r") as bf:
        for line in bf:
            line = line.strip()
            clust, start, end = line.split("\t")
            if not clust in core_lengths.keys():
                core_lengths[clust] = 0
            core_lengths[clust] += int(end) - int(start)
            if not clust in core_starts.keys():
                core_starts[clust] = []
            if not clust in core_ends.keys():
                core_ends[clust] = []
            core_starts[clust].append(int(start))
            core_ends[clust].append(int(end))
    # Parsing bedgraph for entries
    nocov = {}
    with open(bedgraph, "r") as f:
        for line in f:
            line = line.strip()
            cluster, start, end, cov = line.split("\t")
            # if "--NR" in cluster_NR:
            #    NR_index = cluster_NR.find("--NR")
            #    cluster = cluster_NR[:NR_index]
            # else:
            #    cluster = cluster_NR
            if not cluster in nocov:  # make entry
                nocov[cluster] = 0
            if float(cov) == 0:  # enter no coverage values
                not_covered = local_computecov(core_starts[cluster], core_ends[cluster], [int(start), int(end)])
                nocov[cluster] += not_covered
    # Final coverage calculation:
    total_coverage = {}
    for key in nocov.keys():
        perc = (core_lengths[key] - nocov[key]) / core_lengths[key]
        total_coverage[key] = perc
    return (total_coverage)


def familycorrect(c_dict, family):
    """uses family to change values
    parameters
    ----------
    c_dict
        dictionary, {clustername:value}
    family
        json, {HGF representative: HGF members}
    """
    ret = {}

    for GC, v in c_dict.items():
        if "HG_DNA" in GC:
            for HGF_member in family[GC]:
                key_NR = GC[GC.index("--NR"):]
                ret[f"{HGF_member}{key_NR}"] = v
        else:
            ret[GC] = v

    return (ret)


######################################################################
# Functions for writing results and cleaning output directory
######################################################################
def writejson(dictionary, outdir, outfile_name):
    """writes results in a dict to json format
    parameters
    ----------
    dictionary
        dict, dicionary containing some results (here mapping results)
    outdir
        string, the path to output directory
    outfile_name
        string, name for the outfile
    returns
    ----------
    outfile = name of the .json outfile
    """
    outfile = os.path.join(outdir, outfile_name)
    with open(outfile, "w") as w:
        w.write(json.dumps(dictionary, indent=4))
    return (outfile)


def export2biom(outdir, core="",):
    """writes the results to biom format for easy loading into metagenomeSeq
    parameters
    ----------
    outdir
        string, the path to output directory
    returns
    ----------
    biom_file = the created biom-format file (without metadata)
    """
    biom_file = os.path.join(outdir, "BiG-MAP.map" + core + ".biom")
    inp_counts = os.path.join(outdir, "BiG-MAP.map.results." + core + 'RPKM_filtered.txt')
    cmd_export2biom = f"biom convert -i {inp_counts} -o {biom_file} --table-type='Pathway table' --to-json"
    res_export = subprocess.check_output(cmd_export2biom, shell=True)
    return (biom_file)


def decoratebiom(biom_file, outdir, metadata, core=""):
    """inserts rows and column data
    """
    out_biom = '.'.join(biom_file.split('.')[0:-1]) + '.meta.biom'
    cmd_sample = f"biom add-metadata -i {biom_file} -o {out_biom} -m {metadata}"
    res_add = subprocess.check_output(cmd_sample, shell=True)
    if core == "core":
        metadata_f = os.path.join(outdir, 'BiG-MAP.map.core.coverage.txt')
        cmd_feature = f"biom add-metadata --observation-metadata-fp {metadata_f} -i {biom_file} -o {out_biom}"
        res_feature = subprocess.check_output(cmd_feature, shell=True)
    return (out_biom)

def convert2json(out_biom):
    """writes the results to biom format for easy loading into metagenomeSeq
    parameters
    ----------
    outdir
        string, the path to output directory
    returns
    ----------
    biom_file = the created biom-format file (without metadata)
    """
    json_file = '.'.join(out_biom.split('.')[0:-1]) + 'json.biom'
    cmd_convert2json = f"biom convert -i {out_biom} -o {json_file} --to-json"
    res_export = subprocess.check_output(cmd_convert2json, shell=True)
    return (json_file)


def decode_biom(json_file):
    outf = '.'.join(json_file.split('.')[0:-1]) + '.dec.biom'
    with open(json_file, 'rb') as f:
        data = f.read()
    jobj = json.loads(data)
    with open(outf, 'w') as out:
        out.write(json.dumps(jobj, indent=4))

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


def movetodir(outdir, dirname, pattern):
    """moves files matching a patterd into new directory
    parameters
    ----------
    outdir
        string, the path to output directory
    dirname
        string, name of the new direcory
    pattern
        string, regex
    returns
    ----------
    None
    """
    # Make directory
    try:
        os.mkdir(os.path.join(outdir, dirname))
    except:
        pass
    # Move files into new directory
    for f in os.listdir(outdir):
        if re.search(pattern, f):
            try:
                shutil.move(os.path.join(outdir, f), os.path.join(outdir, dirname))
            except:
                pass

######################################################################
# MAIN
######################################################################
def main():
    """
    The following steps are performed for each sample
    1) preparation of mapping (=alignment)
    2) bowtie2 for mapping, counting reads
    3) bedtools for computing coverage for each cluster
    4) using the fastANI result for adding HGF values
    5) saving all the results in dictionary (=memory)
    5) writing the results to .json (=BIOM) and .csv
    6) cleaning output directory
    """
    parser, args = get_arguments()
    
    #get the output results from the family module
    family_results = os.listdir(args.family)
    bf = [file for file in family_results if 'BiG-MAP.GCF_HGF.bed' in file][0]
    bed_file = os.path.join(args.family, bf)
    ref = [file for file in family_results if 'BiG-MAP.GCF_HGF.fna' in file][0]
    reference = os.path.join(args.family, ref)
    js = [file for file in family_results if 'BiG-MAP.GCF_HGF.json' in file][0]
    json_file = os.path.join(args.family, js)
    bjs = [file for file in family_results if 'BiG-MAP.GCF.json' in file][0]
    bjson_file = os.path.join(args.family, bjs)


    if args.fastq1 and args.fastq2 and not args.U_fastq:
        print("__________Fastq-files_________________________________")
        print("\n".join(args.fastq1))
        print("______________________________________________________")
        print("\n".join(args.fastq2))
        fastq_files = zip(args.fastq1, args.fastq2)
    elif not args.fastq1 and not args.fastq2 and args.U_fastq:
        print("__________Fastq-files_________________________________")
        print("\n".join(args.U_fastq))
        fastq_files = zip(args.U_fastq, args.U_fastq)  # In the case of unpaired, m1 and m2 are identical
    else:
        parser.print_help()
        print("ERROR: -I1/-I2 and -U are mutually exclusive")
        sys.exit()

    if json_file and bjson_file and not args.pickle_file:       
        with open(json_file, "r") as jfile:
            family = json.load(jfile)
        with open(bjson_file, "r") as bjfile:
            BGCF = json.load(bjfile)
    elif not reference and not family and args.pickle_file:
        reference, family, BGCF, bed_file = unpickle_files(args.pickle_file, args.outdir + os.sep)
    else:
        parser.print_help()
        print("ERROR: -R/-F and -P are mutually exclusive")
        sys.exit()
    
    try:
        os.mkdir(args.outdir)
    except:
        pass

    results = {}  # Will be filled with TPM,RPKM,coverage for each sample
    results_core = {}
    mapping_percentages = {}  # Mappping percs for each sample

    ##############################
    # Preparing mapping
    ##############################
    i = bowtie2_index(reference, args.outdir + os.sep)

    ##############################
    # Whole cluster calculation
    ##############################
    print('Mapping reads using bowtie')
    for m1, m2 in fastq_files:
        s = bowtie2_map(args.outdir + os.sep, m1, m2, i, args.fasta, args.bowtie2_setting, args.threads)
        b = samtobam(s, args.outdir + os.sep) 
        sortb = sortbam(b, args.outdir + os.sep)
        indexbam(sortb, args.outdir + os.sep)
        countsfile = countbam(sortb, args.outdir + os.sep)
        if BGCF:
            countsfile = correct_counts(countsfile, BGCF)

        TPM = calculateTPM(countsfile)
        RPKM = calculateRPKM(countsfile)
        raw = parserawcounts(countsfile)

        ##############################
        # bedtools: coverage
        ##############################
        bedtools_gfile = preparebedtools(args.outdir + os.sep, countsfile)
        bedgraph = bedtoolscoverage(bedtools_gfile, args.outdir + os.sep, sortb)
        coverage = computetotalcoverage(bedgraph)

        if BGCF:
            coverage = correct_coverage(coverage, countsfile)
            # GCF and HGF consideration:
            TPM = familycorrect(TPM, BGCF)
            RPKM = familycorrect(RPKM, BGCF)
            raw = familycorrect(raw, BGCF)
            coverage = familycorrect(coverage, BGCF)
        else:
            TPM = familycorrect(TPM, family)
            RPKM = familycorrect(RPKM, family)
            raw = familycorrect(raw, family)
            coverage = familycorrect(coverage, family)

        ##############################
        # saving results in one dictionary
        ##############################
        sample = Path(b).stem
        results[f"{sample}.TPM"] = [TPM[k] for k in RPKM.keys()]
        results[f"{sample}.RPKM"] = [RPKM[k] for k in RPKM.keys()]
        results[f"{sample}.RAW"] = [raw[k] for k in RPKM.keys()]
        results[f"{sample}.cov"] = [coverage[k] for k in RPKM.keys()]
        results["gene_clusters"] = list(RPKM.keys())  # add gene clusters as well

        ##############################
        # Core calculation
        ##############################
        if bed_file:
            sortb = extractcorefrombam(sortb, args.outdir + os.sep, bed_file)
            indexbam(sortb, args.outdir + os.sep)
            countsfile = countbam(sortb, args.outdir + os.sep)
            if BGCF:
                countsfile = correct_counts(countsfile, BGCF)

            core_TPM = calculateTPM(countsfile)
            core_RPKM = calculateRPKM(countsfile)
            core_raw = parserawcounts(countsfile)
            # Coverage
            core_bedgraph = bedtoolscoverage(bedtools_gfile, args.outdir + os.sep, sortb)
            core_coverage = computecorecoverage(core_bedgraph, bed_file)
            if BGCF:
                core_coverage = correct_coverage(core_coverage, countsfile)
            # GCF and HGF consideration:
            if BGCF:
                core_TPM = familycorrect(core_TPM, BGCF)
                core_RPKM = familycorrect(core_RPKM, BGCF)
                core_raw = familycorrect(core_raw, BGCF)
                core_coverage = familycorrect(core_coverage, BGCF)
            else:
                core_TPM = familycorrect(core_TPM, family)
                core_RPKM = familycorrect(core_RPKM, family)
                core_raw = familycorrect(core_raw, family)
                core_coverage = familycorrect(core_coverage, family)

            # core_coverage = computetotalcoverage(core_bedgraph)
            results[f"{sample}.coreTPM"] = [core_TPM[k] for k in core_RPKM.keys()]
            results[f"{sample}.coreRPKM"] = [core_RPKM[k] for k in core_RPKM.keys()]
            results[f"{sample}.coreRAW"] = [core_raw[k] for k in core_RPKM.keys()]
            results[f"{sample}.corecov"] = [core_coverage[k] if "GC_DNA--" in k else 0 for k in core_RPKM.keys()]

    ##############################
    # writing results file: pandas
    ##############################
    # writing all the results to csv
    df = pd.DataFrame(results)
    df.set_index("gene_clusters", inplace=True)
    df.to_csv(os.path.join(args.outdir, "BiG-MAP.map.results.ALL.csv"))

    # writing RPKM (core) filtered results
    headers_RPKM = [rpkmkey for rpkmkey in results.keys() if ".RPKM" in rpkmkey]
    df_RPKM = df[headers_RPKM]
    df_RPKM.columns = [h[:-5] for h in headers_RPKM]
    df_RPKM.to_csv(os.path.join(args.outdir, "BiG-MAP.map.results.RPKM_filtered.csv"))
    df_RPKM.to_csv(os.path.join(args.outdir, "BiG-MAP.map.results.RPKM_filtered.txt"), sep="\t")

    headers_coreRPKM = [rpkmkey for rpkmkey in results.keys() if ".coreRPKM" in rpkmkey]
    df_coreRPKM = df[headers_coreRPKM]
    df_coreRPKM.columns = [h[:-9] for h in headers_coreRPKM]
    df_coreRPKM.to_csv(os.path.join(args.outdir, "BiG-MAP.map.results.coreRPKM_filtered.csv"))
    df_coreRPKM.to_csv(os.path.join(args.outdir, "BiG-MAP.map.results.coreRPKM_filtered.txt"), sep="\t")

    # Writing row coverages:
    headers_cov = [corekey for corekey in results.keys() if ".corecov" in corekey]
    df_cov = df[headers_cov]
    df_cov.columns = [h[:-8] for h in headers_cov if ".corecov" in h]
    df_cov.index.names = ['#gene_clusters']
    df_cov.to_csv(os.path.join(args.outdir, "BiG-MAP.map.core.coverage.txt"), sep="\t")

    # writing the results to biom format:
    print('Adding metadeta to biom and converting files into json format')
    if args.biom_output:
        try:
            if bed_file:
                biomfile = export2biom(args.outdir)
                biom_out1 = decoratebiom(biomfile, args.outdir, args.biom_output, "core")
                json_out1 = convert2json(biom_out1)
                decode_biom(json_out1)
                biomfile2 = export2biom(args.outdir, "core")
                biom_out2 = decoratebiom(biomfile2, args.outdir, args.biom_output, "core")
                json_out2 = convert2json(biom_out2)
                decode_biom(json_out2)
            else:
                biomfile = export2biom(args.outdir)
                biom_out = decoratebiom(biomfile, args.outdir, args.biom_output)
                json_out = convert2json(biom_out)
                decode_biom(json_out)
        except(EOFError):
            biomfile = export2biom(args.outdir)

    # writing mapping percentages for each sample to csv
    mapping_percentages = parse_perc(args.outdir)
    df_perc = pd.DataFrame(mapping_percentages)
    df_perc.to_csv(os.path.join(args.outdir, "BiG-MAP.percentages.csv"))

    ##############################
    # Moving and purging files
    ##############################
    print('Reorganizing output directory')
    movetodir(args.outdir + os.sep, "bowtie2-index", ".bt2")
    movetodir(args.outdir + os.sep, "bedtools-results", ".bg")
    movetodir(args.outdir + os.sep, "bedtools-results", ".file")
    movetodir(args.outdir + os.sep, "bowtie2-map-results", ".bam")
    movetodir(args.outdir + os.sep, "bowtie2-map-results", ".sam")
    movetodir(args.outdir + os.sep, "bowtie2-map-results", ".bai")
    movetodir(args.outdir + os.sep, "bowtie2-raw-counts", ".count")
    movetodir(args.outdir + os.sep, "csv-results", ".csv")
    movetodir(args.outdir + os.sep, "csv-results", ".txt")
    movetodir(args.outdir + os.sep, "biom-results", ".biom")
    all_biom = os.listdir(args.outdir + os.sep + "biom-results")
    to_delete = [file for file in all_biom if not 'dec' in file]
    for f in to_delete:
        os.remove(args.outdir + os.sep + "biom-results/" + f)

if __name__ == "__main__":
    main()
