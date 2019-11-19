#! usr/bin/env python3

"""
--------------- Mapping module ---------------
Author: Koen van den Berg
University: Wageningen University and Research
Department: Department of Bioinformatics
Date: 01/07/2019
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
import os.path
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

# Functions:
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
    usage='''
______________________________________________________________________

     BiG-MAP map: maps the paired reads to the predicted MGCs
______________________________________________________________________

Generic command: python3 BiG-MAP.map.py {-I1 [mate-1s] -I2 [mate-2s] | -U [samples]} -R [reference] -O [outdir] -F [family] [Options*]

Maps the metagenomic/metatranscriptomic reads to the fasta reference
file and outputs RPKM read counts in .csv and BIOM format

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

Obligatory arguments:
    -R    Provide the reference fasta file in .fasta or .fna format
    -O    Put path to the output folder where the results should be
          deposited. Default = current folder (.)
    -F    Input the by fastANI defined GCFs and HGFs file named:
          'fastani.GCFheaders.json' here. 

Options:
    -cc   Also calculate the RPKM and coverage values for the core of
          the cluster present in the bedfile. Specify the bedfile
          here. Bedfiles are outputted by BiG-MAP.genecluster.py
          automatically and are named: BiG-MAP.GCF_HGF.bed
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
                      very-sensitive-local
          DEFAULT = fast
    -th   Number of used threads in the bowtie2 mapping step. Default = 6
______________________________________________________________________

''')
    parser.add_argument("-R", "--reference", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-I1","--fastq1", nargs='+',help=argparse.SUPPRESS, required=False)
    parser.add_argument("-I2","--fastq2",nargs='+',help=argparse.SUPPRESS, required = False)
    parser.add_argument("-U","--U_fastq",nargs='+',help=argparse.SUPPRESS, required = False)
    parser.add_argument("-F", "--family", help=argparse.SUPPRESS, required=True)
    parser.add_argument( "-cc", "--corecalculation",
                         help=argparse.SUPPRESS, required = False)
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
        index_name = f"{outdir}{stem}"
        if not os.path.exists(index_name + ".1.bt2"):
            cmd_bowtie2_index = f"bowtie2-build {reference} {index_name}"
            res_index = subprocess.check_output(cmd_bowtie2_index, shell=True)
    except(subprocess.CalledProcessError):
        print("Error-code M3:001, check error table")
        # Proper error here, also exit code
    return(index_name)

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
    samfile = f"{outdir}{sample}.sam"
    # In the case of unpaired, m1 and m2 are identical. Thus the following works:
    sample_command = f"-U {mate1}" if mate1 == mate2 else f"-1 {mate1} -2 {mate2}"
    try:
        if not os.path.exists(samfile):
            cmd_bowtie2_map = f"bowtie2\
            --{bowtie2_setting}\
            --no-unal\
            --threads {threads} \
            {'-f' if fasta else ''} \
            -x {index} \
            {sample_command}\
            -S {samfile}" #  The .sam file will contain only the map results for 1 sample
            
            print(f"the following command will be executed by bowtie2:\n\
_____________________________________________________\n\
{cmd_bowtie2_map}\n\
_____________________________________________________\n")
            res_map = subprocess.check_output(cmd_bowtie2_map, shell=True, stderr = subprocess.STDOUT)
            # Saving mapping percentage:
            with open(f"{outdir}bowtie2_log.txt", "a+") as f:
                f.write(f"#{sample}\n{res_map.decode('utf-8')}")
    except(subprocess.CalledProcessError):
        pass # raise error here
    return(samfile)

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
    infile = f"{outdir}bowtie2_log.txt"
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                sample = line[1:]
            if "overall" in line:
                perc = line.split(" ")[0][:-1]
                perc = float(perc)/100
                ret[sample] = [perc]
    return(ret)

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
    bamfile = f"{outdir}{stem}.bam"
    try:
        cmd_samtobam = f"samtools view\
        -b {sam}\
        > {bamfile}"
        res_samtobam = subprocess.check_output(cmd_samtobam, shell=True)
    except(subprocess.CalledProcessError):
        pass # raise error here
    return(bamfile)

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
    sortedbam = f"{outdir}{stem}.sorted.bam"
    try:
        cmd_sortbam = f"samtools sort {bam} > {sortedbam}"
        res_sortbam = subprocess.check_output(cmd_sortbam, shell=True)
    except(subprocess.CalledProcessError):
        pass # raise error here
    return(sortedbam)

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
        pass # raise error here
    return()

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
        pass # raise error here
    return(counts_file)

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
    bamfile = f"{outdir}core_{bamstem}.bam"
    if os.path.exists(bedfile):
        try:
            cmd_extractcore = f"samtools view\
            -b {bam}\
            -L {bedfile}\
            > {bamfile}"
            res_extractcore = subprocess.check_output(cmd_extractcore, shell=True)
            print(cmd_extractcore)
        except(subprocess.CalledProcessError):
            pass # raise error here
    else:
        # raise bedfile error here!!!
        pass
    return(bamfile)

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
                    nreads = float(nreads)/NR
            try:
                rate = float(nreads)/float(length)
                rates[cluster] = rate
                ratesum += rate
            except(ZeroDivisionError):
                pass
    TPM = {}
    for key in rates:
        try:
            TPM[key] = rates[key]/ratesum
        except(ZeroDivisionError):
            TPM[key] = 0
    return(TPM)

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
                    nreads = float(nreads)/NR
                    read_counts[cluster] = nreads
                else:
                    read_counts[cluster] = float(nreads)
                cluster_lengths[cluster] = float(length)
                sum_reads += float(nreads)

    RPKM = {}
    for key in read_counts:
        try:
            RPKM[key] = read_counts[key]/(sum_reads*cluster_lengths[key]) * 1000000000                    
                
        except(ZeroDivisionError):
            RPKM[key] = 0
    return(RPKM)

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
                    raw_counts[cluster] = float(nreads)/NR
                else:
                    raw_counts[cluster] = float(nreads)
    return(raw_counts)
    
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
    genome_file = f"{outdir}genome.file"
    with open(genome_file, "w") as w:
        with open(reference, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    line = line.replace("core_", "")
                    w.write(f"{c}\n")
                    c = 0
                    w.write(f"{line[1:]}\t")
                else:
                    c += len(line)
            w.write(f"{c}")

    #try:
    #    # This could also be done using samtools idxstats and extracting $1 & $2
    #    cmd_prepare = """awk \
    #    '$0 ~ ">" {\
    #         print c; \
    #         c=0;\
    #         printf substr($0,2,150) "\t";\
    #         } \
    #    $0 !~ ">" {\
    #        c+=length($0);\
    #    } \
    #    END { \
    #    print c; \
    #    }\
    #    ' %s > %s"""%(reference, genome_file) 
    #    res_prepare = subprocess.check_output(cmd_prepare, shell=True)
    #except(subprocess.CalledProcessError):
    #    pass # raise error here

    return(genome_file)
    
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
    bg_file = f"{outdir}{stem.split('.')[0]}.bg"
    
    try:
        cmd_bedtools = f"bedtools genomecov -bga -ibam {sortedbam} -g {gfile} > {bg_file}"
        res_bedtools = subprocess.check_output(cmd_bedtools, shell=True)
    except(subprocess.CalledProcessError):
        pass # raise error here
    return(bg_file)

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
            clusterlen[cluster] = float(end) # last encounter is length
            if not cluster in nocov: # make entry
                nocov[cluster] = 0
            if float(cov) == 0: # enter no coverage values
                    nocov[cluster] += (float(end)-float(start))
    total_coverage = {}
    for key in nocov.keys():
        perc = (clusterlen[key] - nocov[key])/clusterlen[key]
        # Set treshold here!!!
        total_coverage[key] = perc
    return(total_coverage)

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
            if Ls>Ts and Ls<Te and Le>Ts and Le<Te:
                # in
                ret_cov += Le-Ls
            if Ls<Ts and Ls<Te and Le>Ts and Le<Te:
                # start out
                ret_cov += Le-Ts
            if Ls>Ts and Ls<Te and Le>Ts and Le>Te:
                # end out
                ret_cov += Te-Ls  
            if Ls<Ts and Ls<Te and Le>Ts and Le>Te:
                # start&end out
                ret_cov += Te-Ts
        return(ret_cov)
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
            #if "--NR" in cluster_NR:
            #    NR_index = cluster_NR.find("--NR")
            #    cluster = cluster_NR[:NR_index]
            #else:
            #    cluster = cluster_NR
            if not cluster in nocov: # make entry
                nocov[cluster] = 0
            if float(cov) == 0: # enter no coverage values
                not_covered = local_computecov(core_starts[cluster], core_ends[cluster], [int(start), int(end)])
                nocov[cluster] += not_covered
    # Final coverage calculation:
    total_coverage = {}
    for key in nocov.keys():
        perc = (core_lengths[key] - nocov[key])/core_lengths[key]
        total_coverage[key] = perc
    return(total_coverage)



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
    with open(family, "r") as j:
        family = json.load(j)
    
    for GC, v in c_dict.items():
        if "HG_DNA"in GC:
            for HGF_member in family[GC]:
                key_NR = GC[GC.index("--NR"):]
                ret[f"{HGF_member}{key_NR}"] = v
        else:
            ret[GC] = v
    return(ret)




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
    outfile = f"{outdir}{outfile_name}"
    with open(outfile, "w") as w:
        w.write(json.dumps(dictionary, indent=4))
    return(outfile)

def export2biom(outdir, core = ""):
    """writes the results to biom format for easy loading into metagenomeSeq
    parameters
    ----------
    outdir
        string, the path to output directory
    returns
    ----------
    biom_file = the created biom-format file (without metadata)
    """    
    biom_file = f"{outdir}BiG-MAP.map{core}.biom"
    cmd_export2biom = f"biom convert -i {outdir}BiG-MAP.map.results.{core}RPKM_filtered.txt -o {biom_file} --table-type='Pathway table' --to-json"
    res_export = subprocess.check_output(cmd_export2biom, shell=True)
    return(biom_file)

def decoratebiom(biom_file, outdir, metadata, core = ""):
    """inserts rows and column data
    """
    cmd_sample = f"biom add-metadata -i {biom_file} -o {biom_file} -m {metadata}"
    res_add = subprocess.check_output(cmd_sample, shell=True)
    if core == "core":
        cmd_feature = f"biom add-metadata --observation-metadata-fp {outdir}BiG-MAP.map.core.coverage.txt -i {biom_file} -o {biom_file}"
        res_feature = subprocess.check_output(cmd_feature, shell=True)
    with open(biom_file, "r") as f:
        biom_dict = json.load(f)
    with open(biom_file, "w") as w:
        w.write(json.dumps(biom_dict, indent=4))
    return(biom_file)

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
        os.mkdir(f"{outdir}{dirname}")
        print(f"Directory {outdir}{dirname} created")
    except(FileExistsError):
        print(f"Directory {outdir}{dirname} already exists")
    # Move files into new directory
    for f in os.listdir(outdir):
        if re.search(pattern, f):
            shutil.move(os.path.join(outdir,f), os.path.join(outdir, dirname))

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

    if args.fastq1 and args.fastq2 and not args.U_fastq:
        print("__________Fastq-files_________________________________")
        print("\n".join(args.fastq1))
        print("______________________________________________________")
        print("\n".join(args.fastq2))
        fastq_files = zip(args.fastq1, args.fastq2)
    elif not args.fastq1 and not args.fastq2 and args.U_fastq:
        print("__________Fastq-files_________________________________")
        print("\n".join(args.U_fastq))
        fastq_files = zip(args.U_fastq, args.U_fastq) # In the case of unpaired, m1 and m2 are identical
    else:
        parser.print_help()
        print("ERROR: -I1/-I2 and -U are mutually exclusive")
        sys.exit()

    

    results = {} #Will be filled with TPM,RPKM,coverage for each sample
    results_core = {}
    mapping_percentages = {} #Mappping percs for each sample
    
    ##############################
    # Preparing mapping
    ##############################
    i = bowtie2_index( args.reference, args.outdir)

    ##############################
    # Whole cluster calculation
    ##############################
    for m1, m2 in fastq_files:
        s = bowtie2_map(args.outdir, m1, m2, i, args.fasta, args.bowtie2_setting, args.threads)
        b = samtobam(s, args.outdir)
        sortb = sortbam(b, args.outdir)
        indexbam(sortb, args.outdir)
        countsfile = countbam(sortb, args.outdir)
        TPM =  calculateTPM(countsfile)
        RPKM = calculateRPKM(countsfile)
        raw = parserawcounts(countsfile)

        ##############################
        # bedtools: coverage
        ##############################
        bedtools_gfile = preparebedtools(args.outdir, args.reference)
        bedgraph = bedtoolscoverage(bedtools_gfile, args.outdir, sortb)
        coverage = computetotalcoverage(bedgraph)

        # GCF and HGF consideration:
        TPM = familycorrect(TPM, args.family)
        RPKM = familycorrect(RPKM, args.family)
        raw = familycorrect(raw, args.family)
        coverage = familycorrect(coverage, args.family)

        ##############################
        # saving results in one dictionary
        ##############################
        sample = Path(b).stem
        results[f"{sample}.TPM"] = [TPM[k] for k in RPKM.keys()]
        results[f"{sample}.RPKM"] = [RPKM[k] for k in RPKM.keys()]
        results[f"{sample}.RAW"] = [raw[k] for k in RPKM.keys()]
        results[f"{sample}.cov"] = [coverage[k] for k in RPKM.keys()]
        results["gene_clusters"] = list(RPKM.keys()) # add gene clusters as well

        ##############################
        # Core calculation
        ##############################
        if args.corecalculation:
            sortb = extractcorefrombam( sortb, args.outdir, args.corecalculation)
            indexbam( sortb, args.outdir)
            countsfile = countbam(sortb, args.outdir)
            core_TPM =  calculateTPM(countsfile)
            core_RPKM = calculateRPKM(countsfile)
            core_raw = parserawcounts(countsfile)
            # Coverage
            core_bedgraph = bedtoolscoverage(bedtools_gfile, args.outdir, sortb)
            core_coverage = computecorecoverage(core_bedgraph, args.corecalculation)
            # GCF and HGF consideration:
            core_TPM = familycorrect(core_TPM, args.family)
            core_RPKM = familycorrect(core_RPKM, args.family)
            core_raw = familycorrect(core_raw, args.family)
            core_coverage = familycorrect(core_coverage, args.family)
            #core_coverage = computetotalcoverage(core_bedgraph)
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
    df.to_csv(f"{args.outdir}BiG-MAP.map.results.ALL.csv")
    
    # writing RPKM (core) filtered results
    headers_RPKM = [rpkmkey for rpkmkey in results.keys() if ".RPKM" in rpkmkey]
    df_RPKM = df[headers_RPKM]
    df_RPKM.columns = [h[:-5] for h in headers_RPKM]
    df_RPKM.to_csv(f"{args.outdir}BiG-MAP.map.results.RPKM_filtered.csv")
    df_RPKM.to_csv(f"{args.outdir}BiG-MAP.map.results.RPKM_filtered.txt", sep="\t")

    headers_coreRPKM = [rpkmkey for rpkmkey in results.keys() if ".coreRPKM" in rpkmkey]
    df_coreRPKM = df[headers_coreRPKM]
    df_coreRPKM.columns = [h[:-9] for h in headers_coreRPKM]
    df_coreRPKM.to_csv(f"{args.outdir}BiG-MAP.map.results.coreRPKM_filtered.csv")
    df_coreRPKM.to_csv(f"{args.outdir}BiG-MAP.map.results.coreRPKM_filtered.txt", sep="\t")

    # Writing row coverages:
    headers_cov = [corekey for corekey in results.keys() if ".corecov" in corekey]
    df_cov = df[headers_cov]
    df_cov.columns = [h[:-8] for h in headers_cov if ".corecov" in h]
    df_cov.index.names = ['#gene_clusters']
    df_cov.to_csv(f"{args.outdir}BiG-MAP.map.core.coverage.txt", sep="\t")

    # writing the results to biom format:
    if args.biom_output:
        try:
            if args.corecalculation:
                biomfile = export2biom(args.outdir)
                biomdict = decoratebiom(biomfile, args.outdir, args.biom_output, "core")
                biomfile = export2biom(args.outdir, "core")
                biomdict = decoratebiom(biomfile, args.outdir, args.biom_output, "core")
            else:
                biomfile = export2biom(args.outdir)
                biomdict = decoratebiom(biomfile, args.outdir, args.biom_output)    
        except(EOFError):
            biomfile = export2biom(args.outdir)
    
    # writing mapping percentages for each sample to csv
    mapping_percentages = parse_perc(args.outdir)
    df_perc = pd.DataFrame(mapping_percentages)
    df_perc.to_csv(f"{args.outdir}BiG-MAP.percentages.csv")

    ##############################
    # Moving and purging files
    ##############################
    movetodir(args.outdir, "bowtie2-index", ".bt2")
    movetodir(args.outdir, "bedtools-results", ".bg")
    movetodir(args.outdir, "bedtools-results", ".file")
    movetodir(args.outdir, "bowtie2-map-results", ".bam")
    movetodir(args.outdir, "bowtie2-map-results", ".sam")
    movetodir(args.outdir, "bowtie2-map-results", ".bai")
    movetodir(args.outdir, "bowtie2-raw-counts", ".count")
    movetodir(args.outdir, "csv-results", ".csv")
    movetodir(args.outdir, "csv-results", ".txt")
    movetodir(args.outdir, "biom-results", ".biom")
    
if __name__ == "__main__":
    main()
