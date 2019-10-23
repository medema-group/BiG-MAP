#! usr/bin/env python3

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
fastani, and outputs the database.

Additional information: 
This script takes gut- and antiSMASH output directories as input,
converts the .gbk files to fasta files and then calculates gene
cluster families based on the similarity treshold (default=0.9) using
fastANI. In addition, HMMer is used to find several relevant
housekeeping genes from the whole genome genbank files, which are also
present in the antiSMASH output. These housekeeping genes will be
essential downstream for comparing the resutls of the gene clusters
with resutls that are known a priori. The output consists of the
following items: GCFs clusters, GCF fasta file, fastANI results, and a
fastANI histogram. Dependencies: BioPython, awk

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
from Bio.SeqFeature import FeatureLocation
import shutil

def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
    usage=''' 
______________________________________________________________________

  Metaclust genecluster: creates a redundancy filtered reference fna 
______________________________________________________________________

Generic command: python3 metaclust.genecluster.py [Options]* 
-D [input dir(s)] -O [output dir]

Create a redundancy filtered fasta reference file from multiple
anti/gutSMASH outputs.  


Obligatory arguments: 
    -D   Specify the path to the directory containing the gut- or
         antiSMASH outputs here. This could be a singular directory,
         or a space seperated list of directories.
    -O   Put path to the folder where the fastANI filtered gene cluster
         files should be located here. The folder should be an
         existing folder. Default = current folder (.)

Options:
    -t   Fraction between 0 and 1; the similarity treshold that
         determines when the protein sequences of the gene clusters
         can be considered similar (>0.80 is assumed to have the same
         function). Default = 0.9.
    -f   Specify here the number of genes that are flanking the core
         genes of the gene cluster. 0 --> only the core, n --> n
         genes included that flank the core. defualt = 0
    -g   Output whole genome fasta files for the fastANI filtered gene 
         clusters as well. This uses more disk space in the output 
         directory. 'True' | 'False'. Default = False
    -th  Number of used threads in the fastANI filtering step. Default = 6
______________________________________________________________________
''')
    
    parser.add_argument( "-D", "--indir",help=argparse.SUPPRESS, nargs = "+", required = True)
    parser.add_argument( "-O", "--outdir",help=argparse.SUPPRESS, required = True)
    parser.add_argument( "-t", "--treshold", help=argparse.SUPPRESS,
                         required = False, default=0.9, type=float)
    parser.add_argument( "-f", "--flank_genes",
                         help=argparse.SUPPRESS, required=False, type=int, default=0)
    parser.add_argument( "-g", "--genomefiles",
                         help=argparse.SUPPRESS, required=False, type=bool, default=False)
    parser.add_argument( "-th", "--threads", help=argparse.SUPPRESS,
                         type=int, required = False, default=6)
    return(parser.parse_args())


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
                    yield(os.path.join(dirpath,f))
    except: # Include actual exception here
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
        return("{}/{}.gbk".format(path, genome))
    else:
        return("{}/{}.gbk".format(path, genome))
    

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
    gbkcontents = SeqIO.parse(infile,"genbank")
    for record in gbkcontents:
        for feature in record.features:
            # Parsing the regionname
            if "region" in feature.type:
                if "product" in feature.qualifiers:
                    for product in feature.qualifiers["product"]:
                        GCs.append(product)
            # Parsing the protein sequence
            if feature.type == "CDS":
                CDS_index.append(feature_count) # remembering every CDS gene index
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

        if CDS_index.index(min(core_index))-nflank < 0 or CDS_index.index(max(core_index))+nflank+1 > len(CDS_index):
            print(f"!!!flank_genes (-f) is higher than the number of flanking genes in the cluster of file: {infile}, using whole gene cluster instead!!!")
        if nflank == 0:
            core_relative_locs = [record.features[i].location for i in core_index]
        else:
            if CDS_index.index(min(core_index)) - nflank >= 0:
                core_region = CDS_index[CDS_index.index(min(core_index))-nflank:CDS_index.index(max(core_index))+nflank+1]
            else:
                core_region = CDS_index[0:CDS_index.index(max(core_index))+nflank+1]
            core_relative_locs = [record.features[i].location for i in core_region]
        
        # Parsing the DNA sequence
        DNA = record.seq
        #organism = record.annotations['definition'].replace(" ", "_")
        organism = record.description
        organism = "_".join(organism.split(",")[0].split()[:-1])
        organism = organism.replace("(", "")
        organism = organism.replace(")","")
    return(DNA, "".join(proteins), ":".join(GCs), organism, core_relative_locs, absolute_locs)

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
    orgID =  infile.split("/")[-1].split(".")[0]
    orgID = orgID[5:] if 'PROT' in orgID else orgID # For housekeeping orgIDs
    Outfile = f"{outdir}{cluster if seqstype=='HousekeepingGene' else seqstype}-{orgID}.{organism}.{regionno}.fasta"
    fasta_header = f">gb|{orgID}|{seqstype}--Entryname={cluster}--OS={organism}--SMASHregion={regionno}"
    seq = "\n".join(str(sequences)[i:i+80] for i in range(0,len(str(sequences)),80))
    if not os.path.exists(Outfile):
        with open(Outfile, "w") as o:
            o.write(f"{fasta_header}\n")
            o.write(f"{seq}\n")
    return(Outfile, orgID, fasta_header[1:])

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
                start = start.replace("<","")
                end = str(loc.end)
                end = end.replace(">","")
                w.write(f"{clust}\t{start}\t{end}\n")

######################################################################
# similarity between gene clusters using fastANI
######################################################################
def preparefastANI(infile_list, outdir):
    """Makes the query and reference files for fastANI
    parameters
    ----------
    infilelist
        list, contains the .fasta files written by "dirtofasta"
    outdir
        string, the path to output directory
    To be more precise here, fastANI requires as input a file
    containing query .fasta files and reference .fasta files for a
    many-to-many comparison. Therefore, this function will be used
    to create those files. In this case, the query and reference
    files will have identical contents, since a many-to-many
    comparison is performed.
    returns
    ----------
    """
    query = "{}fastani.query.txt".format(outdir)
    reference = "{}fastani.reference.txt".format(outdir)
    with open(query, "w") as wqry, open(reference, "w") as wref:
        for f in infile_list:
            wqry.write("{}\n".format(f))
            wref.write("{}\n".format(f))

def computesimilarity(outdir, threads):
    """computes the similarity between sequences using fastANI
    parameters
    ----------
    outdir
        string, the path to output directory
    pathtofastANI
        string, path to fastANI location
    returns
    ----------
    """
    try:
        cmd_fastani = f"fastANI --ql {outdir}fastani.query.txt --rl {outdir}fastani.reference.txt --fragLen 30 --minFrag 1 -k 16 -t {threads} -o {outdir}fastani.results"
        res_download = subprocess.check_output(cmd_fastani, shell=True)
    except(subprocess.CalledProcessError):
        # Raise error here for error table
        pass
    return()

def fastanihistogram(outdir):
    """visualizes the fastani results using R
    parameters
    ----------
    fastani_out
        str, name of the fastani results file
    returns
    ----------
    """
    # Obtaining directory name
    abspath = os.path.abspath("metaclust.genecluster.py")
    dname = os.path.dirname(abspath)
    try:
        # Reducing file size first using awk
        cmd_awk = """awk -F '\t' '{print $4/$5}' %sfastani.results > %sfastani.Rinput"""
        cmd_awk = cmd_awk %(outdir, outdir)
        cmd_R = "Rscript {}/GCF.visualization.R {}fastani.Rinput {}fastani.histogram.png".format(dname, outdir, outdir)
        res_reduce = subprocess.check_output(cmd_awk, shell=True)
        res_R = subprocess.check_output(cmd_R, shell=True)
    except(subprocess.CalledProcessError):
        # Raise error here for error table
        pass
    return()

######################################################################
# Extracting LSH clusters (buckets) using cut-off
######################################################################
def makeGCF(tsim, outdir, dict_to_add = {}):
    """calculates the GCFs based on similarity treshold
    parameters
    ----------
    tsim
        float, between 0 and 1
    outdir
        string, the path to output directory        
    returns
    ----------
    dsim = {fasta file: similar fasta files}
    """
    # Filtering:
    # This could also be done in python during the clustering step
    cmd_filter = """awk -F '\t' '{if ($4/$5>%f) print $0}' %sfastani.results > %sfastani.results.filtered"""
    cmd_filter = cmd_filter %(tsim-0.05, outdir, outdir)
    res_filter = subprocess.check_output(cmd_filter, shell=True)
    # Clustering:
    dsim = dict_to_add
    infile = f"{outdir}fastani.results.filtered"
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            els = line.split("\t")
            sim = float(els[3])/float(els[4])
            x = [els[0],els[1],sim]
            if x[0] not in dsim:
                if dsim.values():
                    B = any(x[0] in e for e in list(dsim.values()))
                    if not B:
                        dsim[x[0]] = []
                        if float(x[2]) >= tsim:
                            dsim[x[0]].append(x[1])
                else:
                    dsim[x[0]] = []
                    if float(x[2]) >= tsim:
                        dsim[x[0]].append(x[1])
            elif float(x[2]) >= tsim:
                dsim[x[0]].append(x[1])
    return(dsim)

def writeGCFfasta(sim_dict, outdir, outfile, tocombine = "DNA-"):
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
            fkey_contents[-1] = fkey_contents[-1].replace("PROT-", tocombine)
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
    return(outfile)

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
        fkey_contents[-1] = fkey_contents[-1].replace("PROT-", "DNA-")
        fname = "/".join(fkey_contents)
        #fname = fkey.replace("PROT-", tocombine)
        #fname = fname.replace("PROT", "DNA")
        with open(fname, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    fastaheader = f"{line[1:]}--NR={1 if n_repr < 1 else n_repr}" # stripping '>'
                    ret[fastaheader] = []
        for sf in sim_fnames:
            fkey_contents = sf.split("/")
            fkey_contents[-1] = fkey_contents[-1].replace("PROT-", "DNA-")
            sf = "/".join(fkey_contents)

            #sf = sf.replace("PROT", "DNA")
            with open(sf, "r") as s:
                for line in s:
                    line = line.strip()
                    if line.startswith(">"):
                        sim_fastaheader = line[1:] # stripping '>'
                        ret[fastaheader].append(sim_fastaheader)
    return(ret)

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
    """
    outfile = f"{outdir}{outfile_name}.json"
    with open(outfile, "w") as w:
        w.write(json.dumps(dictionary, indent=4))


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
    return(ret)
        
        
    
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
    recordnumber = 0
    gbkcontents = SeqIO.parse(genome_gbkfile,"genbank")
    with open(f"{outdir}seqdb.faa", "w") as db:
        for record in gbkcontents:
            recordnumber += 1
            for feature in record.features:
                if feature.type == "CDS":
                    if "translation" in feature.qualifiers.keys():
                        protseq = feature.qualifiers['translation'][0]
                        loc = feature.location
                        db.write(f">{recordnumber}_{loc}\n{protseq}\n")
    return(f"{outdir}seqdb.faa")
                
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
        pass # raise error here
    return(hmmoutput)

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
    return(genelocs)

def getgenefromgbk(gbkfile, location): # change to work with locations
    """parses a genesequence from a gbk file using the gene location
    parameters
    ----------
    gbkfile
        string, path to gbk file + file
    location
        string of coordinates, example: "[start:end>](+)"
    returns
    ----------
    DNA = dna sequence of the gene
    [scaffoldnr, start, end] of housekeeping gene locations
    """
    recordnumber = int(location.split("_")[0])
    location = location.split("_")[1]
    start = location.split(":")[0][1:]
    start = start.replace(">","") # rare occurence in feature.location
    start = start.replace("<","") # rare occurence in feature.location
    start = int(start)
    end = location.split(":")[1].split("]")[0]
    end = end.replace(">","") # rare occurence in feature.location
    end = end.replace("<","") # rare occurence in feature.location
    end = int(end)
    sign = 1 if "(+)" in location else 0
    sign = -1 if "(-)" in location else sign
    loc = FeatureLocation(start,end,sign)
    abs_loc = [recordnumber, start, end]
    gbkcontents = SeqIO.parse(gbkfile,"genbank")
    recordnumbercheck = 0
    for record in gbkcontents:
        recordnumbercheck += 1
        if recordnumber == recordnumbercheck:
            DNA = record.seq
    return(loc.extract(DNA), abs_loc)

######################################################################
# MAIN
######################################################################
def main():
    """
    """
    args = get_arguments()

    print("___________Extracting fasta files___________________")
    
    ##############################
    # parsing data to fasta files
    ##############################
    GC_enzyme_locs = {} # storing the core + flanking genes locations
    genomedict = {} # will be used to extract housekeeping genes
    prot_fasta_files = [] # will be used as fastANI input
    DNA_fasta_files = [] # will be used as fastANI input
    absolute_locs = {} # VALIDATION
    for d in args.indir:
        for f in retrieveclusterfiles(d):
            # Parsing each .gbk file
            DNAseq, AAseq, GC, organism, core_locs, abs_locs = parsegbkcluster(f, args.flank_genes)
            # writing the full gene clusters to .fasta for fastANI
            prot_file, orgID, fasta_header = writefasta(AAseq, "PROT", GC, organism, f, args.outdir)
            dna_file, orgID, fasta_header_DNA = writefasta(DNAseq, "DNA", GC, organism, f, args.outdir)
            prot_fasta_files.append(prot_file)
            DNA_fasta_files.append(dna_file)
            # remembering the whole genome gbk files
            genomegbk = getgenomegbk(f)
            genomedict[orgID] = genomegbk
            # remembering the enzymatic genes locations in dictionary
            GC_enzyme_locs[fasta_header_DNA] = core_locs
            # Remembring the absolute cluster locations for VALIDATION:
            if not orgID in absolute_locs:
                absolute_locs[orgID] = []
            absolute_locs[orgID].append({fasta_header_DNA: abs_locs})

    ##############################
    # fastani: similarity
    ##############################
    preparefastANI(prot_fasta_files, args.outdir)
    computesimilarity(args.outdir, args.threads)
    #fastanihistogram(args.outdir)

    ##############################
    # LSH buckets clustering
    ##############################
    GCFs = makeGCF(args.treshold, args.outdir)
    #fastadict = makefastaheadersim(GCFs)
    #GCF_enzyme_locs = applyfiltering(GC_enzyme_locs, fastadict)
    # Writing results
    #writejson(GCFs, args.outdir, "fastani.fam")
    #writejson(fastadict, args.outdir, "fastani.GCFheaders")
    


    print("___________Adding housekeeping genes________________")
    
    ##############################
    # housekeeping genes: Obtainging
    ##############################
    hmmfile = Path(f"{sys.path[0]}") / "pfamdomains/combined.hmm"
    processed = []
    HG_fasta_files = []
    absolute_locs["hgenes"] = [] # VALIDATION
    for fname in GCFs.keys():
        orgID = fname.split("/")[-1].split(".")[0][5:] # was [5:] <-- the last one 
        organism = ".".join(fname.split("/")[-1].split(".region")[0].split(".")[1:])
        # find housekeeping genes for org using HMMer
        if organism not in processed: 
            seqdb = prepareseqdb(genomedict[orgID], args.outdir) 
            hmmresult = hmmsearch(seqdb, hmmfile, "", args.outdir)
            genelocs = parsehmmoutput(hmmresult)
            for gene, loc in genelocs.items():
                seq, abs_hloc = getgenefromgbk(genomedict[orgID], loc)
                f, ID, housekeepingheader = writefasta(seq, "HousekeepingGene", gene, organism, fname, args.outdir)
                HG_fasta_files.append(f)
                #append_fasta(GCFs_fasta, f)
                GC_enzyme_locs[housekeepingheader] = [FeatureLocation(0,len(seq))]
                absolute_locs["hgenes"].append({housekeepingheader: abs_hloc})
            print(f"__________DONE: {organism}")
            # Convert genome.gbk --> genome.fasta and store in outdir
            if args.genomefiles:
                gbktofasta(genomedict[orgID], f"{args.outdir}{organism}.fasta", args.outdir)
        processed.append(organism)

    ##############################
    # housekeeping genes: Redundancy filtering
    ##############################
    preparefastANI(HG_fasta_files, args.outdir)
    computesimilarity(args.outdir, args.threads)
    
    GCFs_ALL = makeGCF(0.7, args.outdir, GCFs)
    fastadict_ALL = makefastaheadersim(GCFs_ALL)

    # Writing results to outdir
    writejson(GCFs_ALL, args.outdir, "fastani.fam")
    writejson(fastadict_ALL, args.outdir, "fastani.GCFheaders")
    GCFs_fasta = writeGCFfasta(GCFs_ALL, args.outdir, "metaclust.GCFs_DNA_reference.fna")

    # Remembering the enzymatic gene locations
    GCF_enzyme_locs = applyfiltering(GC_enzyme_locs, fastadict_ALL)
    bedfile = locs2bedfile(GCF_enzyme_locs, f"{args.outdir}metaclust.enzymes.bed") 
    writejson(absolute_locs, args.outdir, "absolute_cluster_locs") # VALIDATION
    
    ##############################
    # Cleaning output dir (10)
    ##############################    
    purge(args.outdir, ".fasta")
    purge(args.outdir, ".txt")
    purge(args.outdir, ".faa")


if __name__ == "__main__":
    main()

    
