#!/usr/bin/env python3

"""
--------------- Download module ---------------
Author: Koen van den Berg, Hannah Augustijn, Victoria Pascal Andreu
University: Wageningen University and Research
Department: Department of Bioinformatics
Date: 21/05/2019
-----------------------------------------------
"""

# Import statements
import os
import subprocess
from sys import argv
import sys
import argparse


# Definitions
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
    usage='''


______________________________________________________________________

    BiG-MAP Download: correctly download using fastq-dump program    
______________________________________________________________________

Generic command: python3 BiG-MAP.download.py -A [accession_list_file] -O [path_to_outdir]

Downloads user specified WGS data from the NCBI. USE BiG-MAP_process
conda environment.

Obligatory arguments:
    -A    Put the filepath to the accession list file here. This file
          should consist of one accession on each line only. 
          EXAMPLE:

          SRR5947807
          SRR5947812
          SRR5947811
          SRR5947810
          SRR5947809
          SRR5947808

          Such accession lists can be downloaded from the SRA run
          selector: https://www.ncbi.nlm.nih.gov/Traces/study/?
    -O    Put path to the output folder where the results should be
          deposited. Default = current folder (.)

Options:
    -f    Boolean: Utilize fasterq-dump for faster file conversion.
          Will not result in gzipped fastq files, but fastq conversion
          can be up to ten times faster. Otherwise uses fastq-dump.
______________________________________________________________________
''')
    parser.add_argument( "-O", "--outdir", help=argparse.SUPPRESS, required = True)
    parser.add_argument( "-A", "--acclist", help=argparse.SUPPRESS, required = True)
    parser.add_argument("-f", "--faster", help=argparse.SUPPRESS, action="store_true", required = False)
    return(parser.parse_args())

######################################################################
# Functions for parsing the user input
######################################################################
def parseacclist(infile):
    """
    parses the user specified accession list file
    parameters
    ----------
    infile
        the name of the accession list input file
    returns
    ----------
    acclist = a list containing each accession number from file
    """
    acclist = []
    try:
        with open(infile, "r") as f:
            for line in f:
                line = line.strip()
                acclist.append(line)
    except(FileNotFoundError):
        # raise error
        print(f"file: {infile} not found!")
    return(acclist)
    
def parsemetadata(number):
    """parses the recommended dataset specified by user
    parameters
    ----------
    number
        the number from the datatable in the documentation
    returns
    ----------
    acclist = a list containing each accession number from file
    """
    acclist = []
    script_path = sys.path[0]
    try:
        with open("{}/Metadata/{}_metadata.txt".format(script_path, number), "r") as f:
            f.readline()
            for line in f:
                line = line.strip()
                accs.append(line.split(",")[0])
    except:
        # Raise error here
        pass
    return(accs)


######################################################################
# Functions for downloading the sra files
######################################################################
def downloadSRA(acc, outdir):
    """Retrieves the SRA file from the online NCBI repository
    parameters
    ----------
    acc 
        str, name of the processed accession    
    outdir 
        str, path to the output directory
    returns
    ----------
   """
    try:
        print(f'Downloading sample {acc}')
        cmd_download = f"prefetch {acc} -O {outdir}"
        res_download = subprocess.check_output(cmd_download, shell=True)
    except:
        print(f'Unable to download sample {acc}')
        pass
    return()

def convertSRAtofastq(acc, outdir, faster):
    """Extracts the .fastq files from the downloaded .sra files
    parameters
    ----------
    outdir
        str, path to the output directory    
    acc 
        str, name of the processed accession
    faster
        bool, use fasterq-dump instead of fastq-dump
    returns
    ----------
    """
    out_fastq = os.path.join(outdir, acc + "_pass_1.fastq")
    if os.path.exists(out_fastq) or os.path.exists(out_fastq + ".gz"):
        print("The .fastq file for {} already exists in {}".format(acc, outdir))
        return

    print(f'Converting SRA {acc} into fastq format')
    sra_folder = os.path.join(outdir, acc)

    if faster:
        # correct output file format and names to be consistent with fastq-dump
        # commands are identical as added args to fastq-dump are default in fasterq-dump
        base_fastq = os.path.join(outdir, acc + "_pass.fastq")
        cmd = f"fasterq-dump --seq-defline '@$ac.$si.$ri $sn length=$rl'\
            --qual-defline '+$ac.$si.$ri $sn length=$rl' --outfile {base_fastq} {sra_folder}"
    else:
        cmd = f"fastq-dump --gzip --skip-technical\
        --readids --read-filter pass --dumpbase --split-3 --clip --outdir {outdir} {sra_folder}"
    res_download = subprocess.check_output(cmd, shell=True)
    return()
        
def sendnotification(acc, outdir, lst):
    """Sends a notification to the terminal
    parameters
    ----------
    acc
        str, name of the processed accession
    outdir
        str, path to the output directory
    lst 
        list, list that contains all the accessions
    returns
    ----------
    """
    print("{}.fastq file has been succesfully downloaded and extracted to {}".format(acc, outdir))
    print("Progress: {}/{}".format(lst.index(acc)+1, len(lst)))

######################################################################
# MAIN
######################################################################

def main():
    """
    The following steps are performed to download the fastq files:
    1) check if user uses recommended set, or SRA run table
    2) parse accession numbers
    3) download SRA data to directory
    4) convert to fastq and compress, remove .sra file
    5) send notification for each download success
    """
    args = get_arguments()
    accessions = []
    
    #make output dir if it doesn't exist
    try:
        os.mkdir(args.outdir)
    except:
        pass
    
    if args.acclist:
        accessions = parseacclist(args.acclist)
    else:
        accessions = parsemetadata(args.datanumber)
    for acc in accessions:
        downloadSRA(acc, args.outdir)
        convertSRAtofastq(acc, args.outdir, args.faster)
        sendnotification(acc, args.outdir, accessions)

# Main code
if __name__ == "__main__":
    main()
