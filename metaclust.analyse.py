#! usr/bin/env python3

"""Author: Koen van den Berg
University: Wageningen University and Research
Department: Department of Bioinformatics
Date: 01/07/2019

The main purpose of this script is to analyse the biom files that are
outputted by module 3. This script will function as a Python wrapper
script.

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
        parser = argparse.ArgumentParser(description="This is a\
        wrapper script for the analysis of the\
        metagenomic/metatranscriptomic samples. This script first\
        normalizes and filters the data. Next either a ZIG model or\
        Kruskall-Wallis model will be implemented to compute\
        differentially abundant/expressed gene\
        clusters. Benjamini-Hochber FDR is implemented to compensate\
        for multiple hypothesis testing. Output will be heatmaps of\
        significantly DA genes.", usage='''python3 metaclust.analyse.py <command> [<args>]

The available commands are:
    inspect    show available comparison options
    test       analyse a BIOM file''')
        parser.add_argument("command", help="Subcommand to run")
        self.args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, self.args.command):
            print("unrecognized argument")
            parser.print_help()
            exit(1)
        getattr(self,self.args.command)()
        
    def inspect(self):
        parser = argparse.ArgumentParser(description="Display the\
        avalaible comparisons that are present in the biom-file file\
        here", usage='''python3 metaclust.analyse.py inspect -B <biom_file> ''')
        parser.add_argument("-B", "--biom_file", help="Provide the\
        biom file here", required=False)
        self.inspect = parser.parse_args(sys.argv[2:])
        #return(self.inspect)

    def test(self):
        parser = argparse.ArgumentParser(description="Tests the\
        present biom file using either a fitZIG model or a\
        Kruskall-Wallis model.", usage = "python3 metaclust.analyse.py\
 test -B <biom_file> -M <meta_group> -G <[groups]> -O <outdir>")
        parser.add_argument("-B", "--biom_file", help="Provide the\
        biom file here", required=True)
        parser.add_argument("-M", "--metagroup", help="provide the\
        metagroup here. This is the first column in the options\
        output. Examples: DiseaseStatus, SampleType etc...", type=str,
                            required=True)
        parser.add_argument("-G", "--groups", help="Space separated\
        list of 2 groups that are to be compared. Example UC non-IBD",
        type=str, nargs='+', required=True)
        parser.add_argument("-O", "--outdir", help="Put the path to\
        the output folder for the results here. The folder should be\
        an existing folder. Default = current folder (.)",
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
        print("--------------------------------------------------")


######################################################################
# Analyse functionality
######################################################################
def analysebiom(biom_file, outdir, MT, groups):
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
    abspath = os.path.abspath("metaclust.genecluster.py")
    dname = os.path.dirname(abspath)
    Rloc = "/mnt/scratch/berg266/programs/R-3.6.1/bin/"
    try:
        # Reducing file size first using awk
        cmd_R = f"{Rloc}Rscript {dname}/metaclust.norm.R {biom_file} {outdir} {MT} {groups[0]} {groups[1]}"
        print(cmd_R)
        res_R = subprocess.check_output(cmd_R, shell=True)
    except(subprocess.CalledProcessError):
        # Raise error here for error table
        pass
    return()

def movetooutdir(outdir, pattern):
    """moves files matching a patterd into output directory
    parameters
    ----------
    outdir
        string, the path to output directory
    pattern
        string, regex
    returns
    ----------
    None
    """
    abspath = os.path.abspath("metaclust.genecluster.py")
    dname = os.path.dirname(abspath)
    # Move files into outdir
    for f in os.listdir(dname):
        if re.search(pattern, f):
            shutil.move(os.path.join(dname,f), outdir)

######################################################################
# MAIN
######################################################################
def main():
    main_args = Arguments()

    if main_args.args.command == "inspect":
        d = extractoptions(main_args.inspect.biom_file)
        pprint(d)
        sys.exit()

    try:
        analysebiom(main_args.test.biom_file,\
                    main_args.test.outdir,\
                    main_args.test.metagroup,\
                    main_args.test.groups)
    except:
        # IF there is no hit (error) then there is no significance
        pass

    movetooutdir(main_args.test.outdir, ".png")
    




    #analysebiom(args.biom_file, args.outdir, args.meta, args.group_1, args.group_2)

    
if __name__ == "__main__":
    main()
    
