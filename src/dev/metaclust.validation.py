#! usr/bin/env python3

"""
--------------- Validation ------------------
Author: Koen van den Berg
University: Wageningen University and Research
Department: Department of Bioinformatics
Date: 16/09/2019
----------------------------------------------

Calculates sensitivity and specificity of the pipeline.
"""

# Import statements:
import os.path
import subprocess
from sys import argv
import sys
import argparse
from pathlib import Path
import json

# Defenitions
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
    usage='''
____________________________________________________________________
                                                                    
            Metaclust Validate: validates the pipeline              
____________________________________________________________________

Generic command: 
python3 metaclust.map.py [Options]* -C [command] -R [sim_reads_file] 
-J [bed_file]

Validates the pipeline by splitting the sim_reads_file and finding the
expected raw read count for the clusters in the simulated reads. 


Obligatory arguments:
    -C    Command: 'validate' or 'split'
    -R    Provide the simulated reads file
    -J    Json absolute locations file from MODULE 2
    -S    Sam file of bowtie2 predictions
    -F    GCFs of fastaheaders
    -O    Outdir
______________________________________________________________________
''')
    parser.add_argument("-R", "--sim_reads", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-J", "--json_locs", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-C", "--command", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-S", "--sam_file", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-F", "--GCFfasta", help=argparse.SUPPRESS, required=True)
    return(parser.parse_args())

######################################################################
# Functions for split
######################################################################
def split_paired_file(paired_file):
    """Splits the paired.fasta file into two seperate files
    parameters
    ----------
    paired_file
        fasta file with both mate pairs in them
    returns
    ----------
    None
    """
    DNA1 = False
    DNA2 = False
    outfile_1 = f"{Path(paired_file).stem}_1.fasta"
    outfile_2 = f"{Path(paired_file).stem}_2.fasta"
    with open(paired_file, "r") as f:
        with open(outfile_1, "w") as w1, open(outfile_2, "w") as w2:
            for line in f:
                #line = line.strip()
                if DNA1 and not "/2" in line:
                    w1.write(line)
                if "/1" in line:
                    w1.write(line)
                    DNA1 = True
                    DNA2 = False
                if DNA2 and not "/1" in line:
                    w2.write(line)
                if "/2" in line:
                    w2.write(line)
                    DNA1 = False
                    DNA2 = True

######################################################################
# Functions for validate
######################################################################
def find_ground_truth(sim_reads, json_locs):
    """Counts the raw read counts present in the simulated reads
    parameters
    ----------
    sim_reads
        fasta file, contains the simulated reads
    json_locs
        bedfile, contains the cluster gene locations
    returns
    ----------
    ret = {fastaheader:true_reads_that_belong_there}
    """
    ret = {}
    with open(json_locs, "r") as j:
        abs_locs_dict = json.load(j) # locations
    # Preparing the return dict:
    for orgID, l in abs_locs_dict.items():
        for d in l:
            for fh in d.keys():
                ret[fh] = []
    # Predicting raw read counts
    with open(sim_reads, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                contents = line.split()
                position = contents[2][20:-1] if "complement" in contents[2] else contents[2][9:]
                read_start = int(position.split("..")[0])
                read_end = int(position.split("..")[1])
                orgID = contents[1].split("=")[1][:-2]
                scaf = int(orgID[-4:])
                description = "_".join(contents[3:][:-5])
                read_organism = description[13:]
                read_ID = contents[0][1:]
                try:
                    for d in abs_locs_dict["hgenes"]:
                        for fh, coords in d.items():
                            scaf_nr, hgene_start, hgene_end = int(coords[0]), int(coords[1]), int(coords[2])
                            org_start = fh.find("--OS=") + 5
                            org_end = fh.find("--SMASH")
                            hgene_organism = fh[org_start:org_end]
                            if scaf == scaf_nr and read_organism == hgene_organism:
                                if read_start >= hgene_start and read_end <= hgene_end:
                                    ret[fh].append(read_ID)
                    for d in abs_locs_dict[orgID]:
                        for fh, coords in d.items():
                            clust_start, clust_end = int(coords[0]), int(coords[1])
                            if read_start >= clust_start and read_end <= clust_end:
                                ret[fh].append(read_ID)
                except(KeyError):
                    pass
    return(ret)

def fastani_validate(GCFheaders, ground_truth):
    """Extend ground truth using fastANI filtered fastaheaders
    parameters
    ----------
    GCFheaders
    ground_truth
    returns
    ----------
    ret = {fastaheader : read_IDs with fastani filtering}
    """
    with open(GCFheaders, "r") as j:
        GCFfasta = json.load(j) # fastani headers
    for headerkey in ground_truth.keys():
        if headerkey in GCFfasta.keys():
            for sim_header in GCFfasta[headerkey]:
                if sim_header != headerkey:
                    for read_ID in ground_truth[sim_header]:
                        ground_truth[headerkey].append(read_ID)
    return(ground_truth)

def find_bowtie2_maps(sam_file):
    """Retrieves the bowtie2 prediction from the sam file
    parameters
    ----------
    sam_file
       sam-file, output of bowtie2 mapping module
    returns
    ----------
    ret = {fastaheader : predicted_read_IDs}
    """
    ret = {}
    with open(sam_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("@"):
                contents = line.split()
                read_ID = contents[0]
                fastaheader = contents[2]
                if "NR" in fastaheader:
                    NR_index = fastaheader.find("--NR")
                    fastaheader = fastaheader[:NR_index]
                if not fastaheader in ret:
                    ret[fastaheader] = []
                ret[fastaheader].append(read_ID)
    return(ret)
    
def validation_metrics(ground_truth, bowtie2_prediction):
    """Calculates the recall and specificity
    Find list instersection!!
    parameters
    ----------
    ground_truth
        dict, contains the ground truth as obtained from the sim reads
    bowtie2_prediction
        dict, contains the predicted read mappings
    returns
    ----------
    metrics = {fastaheader : [truth, prediction, intersection, recall, precision]}
    """
    metrics = {}
    for fh, predicted_IDs in bowtie2_prediction.items():
        true_IDs = ground_truth[fh]
        set_predicted_IDs = set(predicted_IDs)
        set_true_IDs = set(true_IDs)
        intersection = set_predicted_IDs.intersection(set_true_IDs)
        true_count = len(set_true_IDs)
        predicted_count = len(set_predicted_IDs)
        intersection_count = len(intersection)
        try:    
            specificity = intersection_count/predicted_count
        except(ZeroDivisionError):
            specificity = "0 (ZeroDivision)"
        try:
            recall = intersection_count/true_count
        except(ZeroDivisionError):
            recall = "0 (ZeroDivision)"
        if not recall == "0 (ZeroDivision)":
            metrics[fh] = [true_count, predicted_count, intersection_count, recall, specificity]
            
        
    return(metrics)


def make_csv(pcd, outdir):
    """makes csv for predicted counts dict
    parameters
    ----------
    pcd
        dict, predicted counts dict
    outdir
        string, path to outdir
    returns
    ----------
    None
    """
    outfile = f"{outdir}metrics.csv"
    with open(outfile, "w") as w:
        w.write("gene_cluster, true_count, predicted_count, intersection, recall, precision\n")
        for fh, metrics in pcd.items():
            if not "HousekeepingGene" in fh:
                w.write(f"{fh},{metrics[0]},{metrics[1]},{metrics[2]},{metrics[3]},{metrics[4]}\n")



def main():
    """main of validation script
    """
    args = get_arguments()

    if args.command == "validate":
        print("__________Finding ground truth__________")
        ground_truth = find_ground_truth(args.sim_reads, args.json_locs)
        ground_truth = fastani_validate(args.GCFfasta, ground_truth)
        print("__________Finding bowtie2 prediction__________")
        bowtie2_mapped_truth = find_bowtie2_maps(args.sam_file)
        print("__________Calculating validation metrics__________")
        metrics = validation_metrics(ground_truth, bowtie2_mapped_truth)
        make_csv(metrics, args.outdir)
        print("__________Finished calculating__________")

    if args.command == "split":
        split_paired_file(args.sim_reads)
    
if __name__ == "__main__":
    main()
