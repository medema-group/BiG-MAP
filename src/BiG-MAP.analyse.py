#! usr/bin/env python3

"""
--------------- Analyse module ---------------
Author: Koen van den Berg
University: Wageningen University and Research
Department: Department of Bioinformatics
Date: 01/07/2019
----------------------------------------------

This is a wrapper script for the analysis of the
metagenomic/metatranscriptomic samples. This script first normalizes
and filters the data. Next either a ZIG model or Kruskall-Wallis model
will be implemented to compute differentially abundant/expressed gene
clusters. Benjamini-Hochber FDR is implemented to compensate for
multiple hypothesis testing. Output will be heatmaps of significantly
DA genes.

"""
# Import statements:
import os.path
import subprocess
from sys import argv
import sys
import argparse
from pathlib import Path
import json
import re
import shutil

######################################################################
# Argparse functionality
######################################################################
class Arguments(object):
    def __init__(self):
        parser = argparse.ArgumentParser(description="", usage='''

______________________________________________________________________

   BiG-MAP analyse: analyse the biom-outputs (ZIG/Kruskall-Wallis)
______________________________________________________________________

generic command: python3 BiG-MAP.analyse.py <command> [<args>]

Analyse the .BIOM output from BiG-MAP.map.py 


Available commands
    inspect    show available comparison options
    test       analyse a BIOM file
______________________________________________________________________
''')
        parser.add_argument("command", help="Subcommand to run")
        self.args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, self.args.command):
            print("unrecognized argument")
            parser.print_help()
            exit(1)
        getattr(self,self.args.command)()
        
    def inspect(self):
        parser = argparse.ArgumentParser(description="",
        usage='''
______________________________________________________________________

  BiG-MAP analyse: analyse the biom-outputs (ZIG/Kruskall-Wallis)
______________________________________________________________________

Generic command: python3 BiG-MAP.analyse.py inspect -B [biom_file]

Display the avalaible comparisons that are present in the biom-file
file here


Obligatory arguments:
    -B    Provide the Biom file here 

Optional arguments:
    -e    Explore and create a heatmap of the biom-file. This function
          filters the values for highest core coverage and outputs a 
          heatmap. Specify outdir here. Requires -s and -m
    -s    Sampletype: metagenomic | metatranscriptomic
    -m    provide the metagroup here. This is the first column in the 
          inspect output. Examples: DiseaseStatus, Longitude, etc...

______________________________________________________________________
''')
        parser.add_argument("-B", "--biom_file", help=argparse.SUPPRESS, required=True)
        parser.add_argument("-e", "--explore", help=argparse.SUPPRESS, required=False)
        parser.add_argument("-s", "--sampletype", help=argparse.SUPPRESS, required=False)
        parser.add_argument("-m", "--metagroup", help=argparse.SUPPRESS, required=False)
        self.inspect = parser.parse_args(sys.argv[2:])
        #return(self.inspect)

    def test(self):
        parser = argparse.ArgumentParser(description="", usage ='''


______________________________________________________________________

BiG-MAP analyse: analyse the biom-outputs (ZIG/Kruskall-Wallis)
______________________________________________________________________

Generic command: python3 BiG-MAP.analyse.py test -B [biom_file] 
-T [SampleType] -M [meta_group] -G [[groups]] -O [outdir]



Tests the present biom file using either a fitZIG model or a
Kruskall-Wallis model. Note that it is also possible to work in R
studio with the R script: meteclust.norm.R

Obligatory arguments:
    -B    Provide the Biom file here
    -T    metagenomic/metatranscriptomic
    -M    provide the metagroup here. This is the first column in the 
          options output. Examples: DiseaseStatus, Longitude, etc...
    -G    Space separated list of 2 groups that are to be compared. 
          Example: UC and non-IBD --> UC non-IBD
    -O    Put path to the output folder where the results should be
          deposited. Default = current folder (.)
______________________________________________________________________
''')
        parser.add_argument("-B", "--biom_file",
                            help=argparse.SUPPRESS, required=True)
        parser.add_argument("-T", "--sample_type",
                            help=argparse.SUPPRESS, type=str,
                            required=True)
        parser.add_argument("-M", "--metagroup",
                            help=argparse.SUPPRESS, type=str,
                            required=True)
        parser.add_argument("-G", "--groups", help=argparse.SUPPRESS,
                            type=str, nargs='+', required=True)
        parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS,
                            required=True)
        self.test = parser.parse_args(sys.argv[2:])

######################################################################
# Inspect functionality
######################################################################
def extractoptions(biom_file):
    """extracts the meta and group options from the biom file
    parameters
    ----------
    biom_file
        biom-format, filename of biom input file
    returns
    ----------
    ret = {metagroup:groups}
    """
    with open(biom_file, "r") as f:
        biom_dict = json.load(f)
    results = {}
    for key in biom_dict["columns"]:
        for k,meta in key["metadata"].items():
            if k not in results:
                results[k] = []
            else:
                results[k].append(meta)
    ret = {key:set(results[key]) for key in results}
    return(ret)

def pprint(d, indent=0):
    """prints a dictionary to the terminal in readable format
    parameters
    ----------
    d
        dict, a dictionary
    indent
        int, the level of indentation
    returns
    ----------
    None
    """
    for key, value in d.items():
        print('\t' * indent + str(key))
        if isinstance(value, dict):
            pretty(value, indent+1)
        else:
            print('\t' * (indent+1) + str(value))
        print("______________________________________________________________________")


def makeexploreheatmap(biom_file, sampletype, outdir, MT):
    """wrapper for R script to analyse biom formats
    parameters
    ----------
    biom_file
        biom-format, filename of biom input file
    outdir
        string, path to output dir
    MT
        string, name of the metagroup
    returns
    ----------
    None
    """
    # Obtaining directory name
    abspath = os.path.abspath("BiG-MAP.analyse.py")
    dname = os.path.dirname(abspath)
    Rloc = "/mnt/scratch/berg266/programs/R-3.6.1/bin/"
    s_type = "METAGENOMIC" if "genomic" in sampletype else "METATRANSCRIPTOMIC"
    try:
        cmd_R = f"{Rloc}Rscript {dname}/BiG-MAP.norm.R {biom_file} {s_type} {outdir} '{MT}' 'xxx' 'xxx' TRUE"
        res_R = subprocess.check_output(cmd_R, shell=True)
    except(subprocess.CalledProcessError):
        print("__________ERROR________________________________________")
        print(f"An error occured, did you input the correct values? Current values:\nbiom file: {biom_file}\nSampletype: {sampletype} -->{s_type}\nMetagroup: {MT if MT else 'metagroup is not specified, use the inspect function solely without -e tag'}")
        sys.exit()
    return()



        
######################################################################
# Analyse functionality
######################################################################
def analysebiom(biom_file, sampletype, outdir, MT, groups):
    """wrapper for R script to analyse biom formats
    parameters
    ----------
    biom_file
        biom-format, filename of biom input file
    outdir
        string, path to output dir
    MT
        string, name of the metagroup
    groups
        list, 2 strings of groups to be compared
    returns
    ----------
    None
    """
    # Obtaining directory name
    abspath = os.path.abspath("BiG-MAP.analyse.py")
    dname = os.path.dirname(abspath)
    Rloc = "/mnt/scratch/berg266/programs/R-3.6.1/bin/"
    s_type = "METAGENOMIC" if "genomic" in sampletype else "METATRANSCRIPTOMIC"
    try:
        # Reducing file size first using awk
        cmd_R = f"{Rloc}Rscript {dname}/BiG-MAP.norm.R {biom_file} {s_type} {outdir} '{MT}' '{groups[0]}' '{groups[1]}' FALSE"
        res_R = subprocess.check_output(cmd_R, shell=True)
    except(subprocess.CalledProcessError):
        print("########## ERROR ####################################")
        print(f"An error occured, did you input the correct values? Current values:\nbiom file: {biom_file}\nSampletype: {sampletype} -->{s_type}\nMetagroup: {MT}\nGroup_1: {groups[0]}\nGroup_2: {groups[1]}")
        print("Visit the available options through the 'inspect' command")
        sys.exit()
    return()

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
    abspath = os.path.abspath("BiG-MAP.analyse.py")
    dname = os.path.dirname(abspath)
    # Make directory
    """
    try:
        os.mkdir(f"{outdir}{dirname}")
        print(f"Directory {outdir}{dirname} created")
    except(FileExistsError):
        print(f"Directory {outdir}{dirname} already exists")
    """
    # Move files into new directory
    for f in os.listdir(dname):
        if re.search(pattern, f):
            shutil.move(os.path.join(dname,f), os.path.join(outdir))

######################################################################
# MAIN
######################################################################
def main():
    main_args = Arguments()

    if main_args.args.command == "inspect":
        print("__________Extracting options________________________")        
        d = extractoptions(main_args.inspect.biom_file)
        pprint(d)
        if main_args.inspect.explore and main_args.inspect.metagroup:
            print("__________Creating Heatmap__________________________")        
            makeexploreheatmap(main_args.inspect.biom_file, \
                               main_args.inspect.sampletype, \
                               main_args.inspect.explore, \
                               main_args.inspect.metagroup)
            print("__________Creation Successful_______________________")        
        else:
            print(f"Specify the -s and -m tags as well! Current values:\n-s\t{main_args.inspect.sampletype}\n-m\t{main_args.inspect.metagroup}\nThe correct values for the -s and -m tag can be found using only the inspect function without the -e tag.")


        print("___________Moving result files to outdir_____________")
        movetodir(main_args.inspect.explore, "analyse_results", ".png")


    if main_args.args.command == "test":
        print("__________Loading R Functions_______________________")
        analysebiom(main_args.test.biom_file,\
                    main_args.test.sample_type,\
                    main_args.test.outdir,\
                    main_args.test.metagroup,\
                    main_args.test.groups)
        print("__________Testing Successful________________________")        

        print("___________Moving result files to outdir_____________")
        movetodir(main_args.test.outdir, "analyse_results", ".png")
        movetodir(main_args.test.outdir, "analyse_results", ".csv")
    
if __name__ == "__main__":
    main()
    
