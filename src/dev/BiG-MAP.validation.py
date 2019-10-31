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
import pandas as pd

# Defenitions
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
    usage='''
____________________________________________________________________
                                                                    
            BiG-MAP Validate: validates the pipeline              
____________________________________________________________________

Generic command: 
python3 BiG-MAP.validate.py [Options]* -C [command] -R [sim_reads_file] 
-J [bed_file]

Validates the pipeline by splitting the sim_reads_file and finding the
expected raw read count for the clusters in the simulated reads. 


Obligatory arguments:
    -C    Command: 'validate' | 'split' | 'results'
    -R    Provide the simulated reads file
    -J    Json absolute locations file from MODULE 2
    -S    Sam file of bowtie2 predictions
    -F    GCFs of fastaheaders
    -O    Outdir

Options
    -g    ground truth json file
    -r    results .csv files
    -s    fastANI treshold setting
    -n    name of output files
______________________________________________________________________
''')
    parser.add_argument("-R", "--sim_reads", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-J", "--json_locs", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-C", "--command", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-S", "--sam_file", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-F", "--GCFfasta", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-r", "--results", help=argparse.SUPPRESS, required=False, nargs = "+")
    parser.add_argument("-g", "--ground_truth", help=argparse.SUPPRESS, required=False)
    parser.add_argument("-s", "--fastANI_setting", help=argparse.SUPPRESS, required=False)
    parser.add_argument("-n", "--outname", help=argparse.SUPPRESS, required=False)
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
                                    """
                                    print("__________________________________________________")
                                    print(f"from HG header:{hgene_organism}")
                                    print(f"from read:{read_organism}")
                                    print(f"{fh}")
                                    """
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

    ret = {}
    for GCFheaderkey in GCFfasta.keys():
        ret[GCFheaderkey] = []
        for sim_header in GCFfasta[GCFheaderkey]:
            for read_ID in ground_truth[sim_header]:
                ret[GCFheaderkey].append(read_ID)
    return(ret)

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
                #if "NR" in fastaheader:
                #    NR_index = fastaheader.find("--NR")
                #    fastaheader = fastaheader[:NR_index]
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
        # if ground_truth, prediction and intersection is 0, recall and precision are 1
        if true_count == 0 and predicted_count == 0: # Implies intersection is also 0
            true_count = 1e-99
            predicted_count = 1e-99
            intersection = 1e-99
        try:    
            specificity = intersection_count/predicted_count
        except(ZeroDivisionError):
            specificity = "0 (ZeroDivision)"
        try:
            recall = intersection_count/true_count
        except(ZeroDivisionError):
            recall = "0 (ZeroDivision)"
        #if not recall == "0 (ZeroDivision)":
        metrics[fh] = [true_count, predicted_count, intersection_count, recall, specificity]
    return(metrics)

def validation_metrics_v2(ground_truth, bowtie2_prediction):
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
    for fh, true_IDs in ground_truth.items():
        if fh in bowtie2_prediction.keys():
            predicted_IDs = bowtie2_prediction[fh]
        else:
            predicted_IDs = []
        
        set_predicted_IDs = set(predicted_IDs)
        set_true_IDs = set(true_IDs)
        intersection = set_predicted_IDs.intersection(set_true_IDs)
        true_count = len(set_true_IDs)
        predicted_count = len(set_predicted_IDs)
        intersection_count = len(intersection)
        if true_count == 0 and predicted_count == 0: # Implies intersection is also 0
            true_count = 1e-99
            predicted_count = 1e-99
            intersection_count = 1e-99
        try:    
            specificity = intersection_count/predicted_count
        except(ZeroDivisionError):
            specificity = "0 (ZeroDivision)"
        try:
            recall = intersection_count/true_count
        except(ZeroDivisionError):
            recall = "0 (ZeroDivision)"
        #if not recall == "0 (ZeroDivision)":
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
    outfile_genecluster = f"{outdir}genecluster.metrics.csv"
    outfile_housekeepinggene = f"{outdir}housekeepinggene.metrics.csv"
    with open(outfile_genecluster, "w") as wg, open(outfile_housekeepinggene, "w") as wh:
        wg.write("gene_cluster, true_count, predicted_count, intersection, recall, precision\n")
        wh.write("housekeepinggene, true_count, predicted_count, intersection, recall, precision\n")
        for fh, metrics in pcd.items():
            if "HG_" in fh:
                wh.write(f"{fh},{metrics[0]},{metrics[1]},{metrics[2]},{metrics[3]},{metrics[4]}\n")
            else:
                wg.write(f"{fh},{metrics[0]},{metrics[1]},{metrics[2]},{metrics[3]},{metrics[4]}\n")

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

######################################################################
# Functions for result
######################################################################
def process_results(csvfile):
    """Processes the csv results files into summary
    """
    avg_recall = 0
    avg_precision = 0
    csvstem = Path(csvfile).stem
    fastani_t = ".".join(csvstem.split(".")[0:2])
    bowtie2_s = csvstem.split(".")[2]
    
    #fastani_t = csvstem.split(".")[0]
    #bowtie2_s = csvstem.split(".")[1]
    with open(csvfile, "r") as f:
        for idx, line in enumerate(f):
            if idx > 0:
                line = line.strip()
                recall = 0 if line.split(",")[4] == "0 (ZeroDivision)" else line.split(",")[4]
                precision = 0 if line.split(",")[5] == "0 (ZeroDivision)" else line.split(",")[5]
                avg_recall += float(recall)
                avg_precision += float(precision)
        avg_recall = avg_recall/(idx)
        avg_precision = avg_precision/(idx)
    return([bowtie2_s, fastani_t, avg_recall, avg_precision])

def makepivottable(sum_list, summary_type, outdir, name):
    """Creates a pandas dataframe from summary list
    parameters
    ----------
    sum_list
        list, [[entry1], [entry2], ...]
    summary_type
        string, geneclusters | housekeepinggenes
    outdir
        string, the pathname to the output directory
    returns
    ----------
    df_recall = dataframe containing recall values
    df_precision = dataframe containing precision values
    """
    df = pd.DataFrame(sum_list)
    df.columns = ["bowtie_setting", "fastANI_setting", "recall", "precision"]
    df_recall = df[['bowtie_setting', 'fastANI_setting', 'recall']]
    df_precision = df[['bowtie_setting', 'fastANI_setting', 'precision']]
    df_recall = pd.pivot_table(df_recall, values = 'recall', index = ['bowtie_setting'], columns = 'fastANI_setting')
    df_precision = pd.pivot_table(df_precision, values = 'precision', index = ['bowtie_setting'], columns = 'fastANI_setting')

    df_recall.to_csv(f"{outdir}{name}.recall.{summary_type}.csv")
    df_precision.to_csv(f"{outdir}{name}.precision.{summary_type}.csv")

    
    
    return(df_recall, df_precision)
    
    


def makegradienttable(dataframe, outdir):
    """Makes a gradient table for all the settings
    parameters
    ----------
    dataframe
        pandas dataframe
    outdir
        output directory
    returns
    ----------
    tablepng = .png of table, filename
    """
    dataframe.style.background_gradient(cmap=blues)
    pass
# https://stackoverflow.com/questions/12286607/making-heatmap-from-pandas-dataframe


def main():
    """main of validation script
    """
    args = get_arguments()

    if args.command == "validate":
        print("__________Finding ground truth")
        if os.path.exists(f"{args.ground_truth}.{args.fastANI_setting}"):
            print(f"__________opening:{args.ground_truth}.{args.fastANI_setting}")
            with open(f"{args.ground_truth}.{args.fastANI_setting}", "r") as j:
                ground_truth = json.load(j) # locations
        else:
            ground_truth = find_ground_truth(args.sim_reads, args.json_locs)
            writejson(ground_truth, args.outdir, f"ground_truth.{args.fastANI_setting}")
        ground_truth = fastani_validate(args.GCFfasta, ground_truth)
        
        print("__________Finding bowtie2 prediction")
        bowtie2_mapped_truth = find_bowtie2_maps(args.sam_file)
        #writejson(bowtie2_mapped_truth, args.outdir, f"bowtie2_prediction.{args.fastANI_setting}")

        print("__________Calculating validation metrics")
        metrics = validation_metrics_v2(ground_truth, bowtie2_mapped_truth)
        make_csv(metrics, args.outdir)
        print("__________Finished calculating")

    if args.command == "split":
        split_paired_file(args.sim_reads)

    if args.command == "results":
        summary_geneclusters = []
        summary_housekeepinggenes = []
        for csvf in args.results:
            if "genecluster" in csvf:
                s = process_results(csvf)
                summary_geneclusters.append(s)
            else:
                s = process_results(csvf)
                summary_housekeepinggenes.append(s)
        df_rg, df_pg = makepivottable(summary_geneclusters, "genecluster", args.outdir, args.outname)
        df_rh, df_ph = makepivottable(summary_housekeepinggenes, "housekeepinggene", args.outdir, args.outname)
        
        print("__________RECALL__________")
        print(df_rh)
        print("__________PRECISION_______")
        print(df_ph)

        #sns.heatmap(df_recall, annot=True, fmt=".1f")
        #plt.show()
        #makegradienttable(df_recall, args.outdir)
        pass
    
if __name__ == "__main__":
    main()
