 #! usr/bin/env python3

"""
--------------- Download module ---------------
Author: Koen van den Berg
University: Wageningen University and Research
Department: Department of Bioinformatics
Date: 21/05/2019
-----------------------------------------------
"""

# Import statements
import os.path
import subprocess
from sys import argv
import sys
import argparse


# Definitions
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
    usage='''


######################################################################
# Metaclust Download: correctly download using fastq-dump program    #
######################################################################
Generic command: python3 metaclust.download.py [Options]* -A [accession_list_file] -O [path_to_outdir]



Downloads user specified WGS data from the NCBI. Dependencies: fastq-dump

Obligatory arguments:
    -A    Put the filepath to the accession list file here. This file
          should consist of one accession on each line only. 
    -O    Put path to the output folder where the results should be
          deposited. Default = current folder (.)

''')
    parser.add_argument( "-n", "--datanumber", help="Input the WGS\
    data, as mentioned in the table with data recommendations (see\
    https://github.com/KoenvdBerg/MSc_thesis). Input number\
    here.", required = False)
    parser.add_argument( "-O", "--outdir", help=argparse.SUPPRESS, required = True)
    parser.add_argument( "-A", "--acclist", help=argparse.SUPPRESS, required = True)
    parser.add_argument( "-p", "--fastqdump", help="Specify the full\
    path to the fastq-dump program location here. default =\
    /bin/fastq-dump. In the case that the program is not installed,\
    download the binary from:\
    https://ncbi.github.io/sra-tools/fastq-dump.html. Then\
    unzip the package and use the path to that folder here.", required
    = False)
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
    if os.path.exists("{}/{}.sra".format(outdir, acc)):
        print("The srafile: {}.sra already exists in {}".format(acc,outdir))
    else:
        Path = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{}/{}/{}/{}".format(acc[:3], acc[:6], acc, acc+".sra")
        cmd_download = "wget {} -P {}".format(Path, outdir)
        res_download = subprocess.check_output(cmd_download, shell=True)
    return()

def convertSRAtofastq(acc, outdir, pathtofastqdump):
    """Extracts the .fastq files from the downloaded .sra files
    parameters
    ----------
    outdir
        str, path to the output directory    
    acc 
        str, name of the processed accession
    returns
    ----------
    """
    if os.path.exists("{}/{}_pass_1.fastq.gz".format(outdir, acc)):
        print("The .fastq file for {}*.fastq.gz already exists in {}".format(acc, outdir))
    else:
        cmd = f"{pathtofastqdump if pathtofastqdump else ''}fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --outdir {outdir} {outdir}/{acc}.sra"
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
    print(args)
    
    accessions = []
    
    if args.acclist:
        accessions = parseacclist(args.acclist)
    else:
        accessions = parsemetadata(args.datanumber)
    for acc in accessions:
        downloadSRA(acc, args.outdir)
        convertSRAtofastq(acc, args.outdir, args.fastqdump)
        sendnotification(acc, args.outdir, accessions)

# Main code
if __name__ == "__main__":
    main()
