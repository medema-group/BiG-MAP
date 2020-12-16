#!/usr/bin/env python3

"""
--------------- Analyse module ---------------
Authors: Hannah Augustijn, Koen van den Berg, Victoria Pascal Andreu
University: Wageningen University and Research
Department: Department of Bioinformatics
Date: 09/03/2020
----------------------------------------------
This script performs a statistical analysis of the
metagenomic/metatranscriptomic samples. First, the script
normalizes and filters the data. Next, the best covered gene
clusters can be observed and the Kruskal Wallis and fitZIG model
will be used to compute differentially abundant/expressed gene clusters.
Benjamini-Hochberg FDR compensates for multiple hypothesis testing.
The output of the script are heatmaps in pdf format.
"""
# Import statements:
import os.path
import subprocess
import sys
import argparse
import re
import shutil
import json
import pandas as pd
import numpy as np
from scipy.stats import mstats
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from rpy2.robjects.packages import importr, STAP
import rpy2.robjects.packages as rpackages

######################################################################
# Argparse functionality
######################################################################
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
                                     usage='''
______________________________________________________________________
BiG-MAP analyse: analyse the biom-outputs (ZIG/Kruskal-Wallis)
______________________________________________________________________
Generic command: python3 BiG-MAP.analyse.py --explore --compare
-B [biom_file] -T [SampleType] -M [meta_group] -O [outdir] [options*]
Tests the present biom file using either a fitZIG model or a
Kruskal-Wallis model. 
Obligatory arguments:
    -B          Provide the Biom file here
    -T          Metagenomic/metatranscriptomic
    -M          provide the metagroup here. This is the first column in the
                options output. Examples: DiseaseStatus, Longitude, etc...
    -O          Put path to the output folder where the results should be
                deposited. Default = current folder (.)
    --explore   Explore the best covered gene clusters.
    -t          Optional argument for -explore: adjust the
                amount of displayed gene clusters. Default = 20.
    -fe         File name for explore heatmap. Default = explore_heatmap.
    --compare   Make a comparison between two groups using fitZIG and Kruskal-Wallis.
    -g          Space separated list of 2 groups that are to be compared.
                Example: UC and Control --> UC Control.
    -af         Alpha value to determine significance of the fitZIG model. Default=0.05.
    -ak         Alpha value to determine significance of the Kruskal Wallis model. Default=0.05.
    -fc         Output file names for Kruskal Wallis and fitZIG heatmaps. Input first the 
                name for Kruskal Wallis, then for fitZIG. Example: map1_kruskall map2_fitZIG.
                Default = [group1]vs[group2]_[kw/fz].
_____________________________________________________________
''')

    parser.add_argument("-B", "--biom_file",
                        help=argparse.SUPPRESS, required=True)
    parser.add_argument("-T", "--sample_type",
                        help=argparse.SUPPRESS, type=str.upper,
                        required=True, choices=['METAGENOMIC', 'METATRANSCRIPTOMIC'])
    parser.add_argument("-M", "--metagroup",
                        help=argparse.SUPPRESS, type=str, required=True)
    parser.add_argument("-c", "--compare", action='store_true', help=argparse.SUPPRESS,
                        required=False)
    parser.add_argument("-e", "--explore", action='store_true', help=argparse.SUPPRESS,
                        required=False)
    parser.add_argument("-g", "--groups", help=argparse.SUPPRESS,
                        type=str, nargs='+', required=False)
    parser.add_argument("-t", "--threshold", help=argparse.SUPPRESS,
                        type=int, default=20, required=False)
    parser.add_argument("-af", "--alpha_fitzig", help=argparse.SUPPRESS,
                        type=float, default=0.05, required=False)
    parser.add_argument("-ak", "--alpha_kruskal", help=argparse.SUPPRESS,
                        type=float, default=0.05, required=False)
    parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS,
                        required=True)
    parser.add_argument("-fe", "--file_name_explore", help=argparse.SUPPRESS,
                        required=False)
    parser.add_argument("-fc", "--file_names_compare", help=argparse.SUPPRESS,
                        type=str, nargs='+', required=False)
    return parser.parse_args()

######################################################################
# CONVERT BIOM
######################################################################

def export2biom(biom_file, outdir):
    """Converts biom input file to json dictionary file
    parameters
    ----------
    biom_file
        string, path to biom file
    outdir
        string, path to output directory
    returns
    ----------
    json_file = the created json-format file
    """
    json_file = os.path.join(outdir, "BiG-MAP.table.txt")
    cmd_export2biom = f"biom convert -i {biom_file} -o {json_file} \
    --table-type='Pathway table' --to-json"
    res_export = subprocess.check_output(cmd_export2biom, shell=True)
    return json_file

def get_sample_type(sample_type, json_file, metagroup):
    """Parses the biom file to extract the sample information belonging
    to the inputted sample type (metagenomic/metatranscriptomic).
    ----------
    sample_type
        string, name of the sample type: metagenomic/metatranscriptomic.
    json_file
        dict, biom file converted to a dictionary format
    metagroup
        string, name of the metagroup
    returns
    ----------
    index_samples = dict, {index of sample: sample id}
    list_samples = list, all of sample names belonging to the metagroup
    metadata = dict, {sample id: metadata}
    """
    index_samples = {}
    metadata = {}
    list_samples = []

    for sample in json_file["columns"]:
        if sample_type == (sample["metadata"]["SampleType"]).upper():
            list_samples.append(sample["id"])
            index = json_file["columns"].index(sample)
            index_samples[index] = sample
            metadata[sample["id"]] = sample["metadata"][metagroup]
    return(index_samples, list_samples, metadata)

def filter_rows(json_file, sample_index):
    """
    Removes rows which don't have at least the amount of positives of
    the cutoff
    ----------
    json_file
        dict, biom file converted to a dictionary format
    sample_index
        dict, index of samples and sample ids
    returns
    ----------
    out_list = list, indexes which are above the cut-off
    """
    counts = 0
    index_number = 0
    com = []
    dict_counts = {}
    tot_hits = []
    out_list = []
    gc_ids = []

    # Cut-off such that 25% of the samples have a hit
    cutoff = int((len(sample_index.keys()))*0.25)
    # Determines the amount of hits for a GC
    for data_index in json_file["data"]:
        for key in sample_index.keys():
            if data_index[1] == key:
                if data_index[0] not in com:
                    com.append(data_index[0])
                    counts = 0
                    dict_counts[data_index[0]] = 0
                if data_index[0] in com:
                    counts += 1
                    dict_counts[data_index[0]] = counts
    # Compares the total amount of hits for a GC with the cut-off
    for ids in json_file["rows"]:
        if ids["id"] not in gc_ids:
            gc_ids.append(ids["id"])
    for ids2 in gc_ids:
        for key in dict_counts.keys():
            if "HG_" in ids2 and key not in tot_hits:
                if key == index_number:
                    tot_hits.append(key)
            else:
                if dict_counts[key] >= cutoff and key not in tot_hits:
                    tot_hits.append(key)
        index_number += 1

    # Create list containing the data above the cut-off
    for data_index in json_file["data"]:
        for index_no in tot_hits:
            if index_no == data_index[0]:
                out_list.append(data_index)
    return out_list

def get_gc_ids(json_file, index_list, index_samples):
    """ Get the GC names and RPKM values from the filtered indexes
    ----------
    json_file
        dict, biom file converted to a dictionary format
    index_list
        list, list of the indexes above the threshold
    index_samples:
        dict, all the indexes with corresponding GC names
    returns
    ----------
    out_dict = dict, {GC_ids: RPKM values}
    """
    out_dict = {}
    rpkm_values = []
    sample_index = []
    gc_ids = []
    index_number = 0

    for ids in json_file["rows"]:
        if ids["id"] not in gc_ids:
            gc_ids.append(ids["id"])
    # filtered index numbers of the samples
    s_ids = index_samples.keys()
    for ids in json_file["rows"]:
        for index in index_list:
            # filter index lists by selecting present sample index numbers
            if index[0] == index_number and index[1] in s_ids:
                rpkm_values.append(index[2])
                sample_index.append(index[1])
                out_dict[gc_ids[index[0]]] = rpkm_values
        # add zero's for missing values
        if sample_index == []:
            pass
        else:
            mis_index = set(list(index_samples.keys())).difference(sample_index)
            for indx in mis_index:
                rpkm_values.insert((indx-list(index_samples.keys())[0]), 0.0)
        rpkm_values = []
        sample_index = []
        index_number += 1
    return out_dict

def make_dataframe(gc_ids, sample_ids):
    """ makes a pandas dataframe
    ----------
    GC_ids
        dict, GC names as keys and ints as values
    sample_ids
        list, names of the sample ids
    returns
    ----------
    df = pandas dataframe
    """
    df = (pd.DataFrame(gc_ids, index=sample_ids)).T
    return df

def norm_log2_data(df):
    """Normalization and log2 convertion as performed by metagenomeseq
    ----------
    df
        dataframe, rows with samples, columns with GC and RPKM values
    returns
    ----------
    norm_df = dataframe with normalized and log2 converted RPKM values
    """
    norm_df = pd.DataFrame()

    df_adj = df.replace(0, np.nan)
    quantile_dict = (df_adj.quantile(axis=0)).to_dict() # calculate the quantile
    # determine the numeric values in the dataframe
    numeric_cols = [col for col in df_adj if df_adj[col].dtype.kind != 'O']
    df_adj[numeric_cols] -= np.finfo(float).eps # substract the machine epsilon
    # calculate the normalization factor by taking the sum of the values below the quantile \
    # normalize the data by dividing the counts by the normalization factor
    for sample in quantile_dict.keys():
        normfac = df_adj[sample][df_adj[sample] <= quantile_dict[sample]].sum()
        norm_df = norm_df.append(df_adj[sample]/(normfac/1000))
    norm_df = norm_df.T
    norm_df[numeric_cols] += 1 # add 1 as correction for log2 calc
    norm_df = (np.log2(norm_df)).replace(np.nan, 0.0) # calculate log2
    return norm_df

def get_coverage(biom_dict, sample_list, gc_names):
    """ get coverage scores from the biom file
    ----------
    biom_dict
        dict, biom file converted to a dictionary format
    sample_list
        list, ids of all the relevant samples
    gc_names
        list, names of all the relevant GCs
    returns
    ---------
    out_dict = dict {GC name: [coverage scores]}
    """
    out_dict = {}
    list_samples = []
    list_values = []
    for ids in biom_dict["rows"]:
        if ids["id"] in gc_names:
            for sample_name in ids["metadata"].keys():
                if sample_name in sample_list and sample_name not in list_samples:
                    list_samples.append(sample_name)
                    list_values.append(float(ids["metadata"][sample_name]))
            out_dict[ids["id"]] = list_values
            list_samples = []
            list_values = []
    return out_dict

def best_cov(cov_df, df, display):
    """ determines the best covered GCs
    ---------
    cov_df
        dataframe, GC names as index and coverage scores as values
    df
        dataframe, dataframe with RPKM scores as values
    display
        int, the amount of displayable GC. Default = 20
    returns
    ---------
    sign_GC = dataframe, filtered on the highest cov scores
    sign_cov = dataframe, filtered coverage dataframe on the display amount
    """
    cov_mean = cov_df.mean(axis=1)
    cov_mean = cov_mean.nlargest(display)
    sign_gc = df[df.index.isin(list(cov_mean.index))]
    sign_cov = cov_df[cov_df.index.isin(list(cov_mean.index))]
    return sign_gc, sign_cov

def get_relevant_hg(sign_gc, df):
    """ extract the relevant housekeeping genes
    ---------
    sign_GC
        dataframe, filtered dataframe of GCs and RPKM values
    df
        dataframe, GC as index, sample IDs as columns and RPKM values
    returns
    ---------
    function_names = list, function desciptions of organisms
    sign_hg = list, full names of the relevant HGs
    id_list = list, GC and HG ids
    """
    id_list = []
    name_list = []
    sign_hg = []
    function_names = []

    index_list = list(sign_gc.index)
    index_full = list(df.index)
    for gc_name in index_list:
        name_list.append(gc_name.split("--OS=")[1].split("--SMASHregion=")[0])
        id_list.append(gc_name.split("|")[1].split(".")[0])

    # make a list of all the relevant HGs
    for id_name in name_list:
        for full_names in index_full:
            if id_name in full_names and "HG_" in full_names:
                sign_hg.append(full_names)

    # get the function from the HG name
    for hg_name in sign_hg:
        function = hg_name.split("Entryname=")[1].split("--OS=")[0]
        if function not in function_names:
            function_names.append(function)
    return(function_names, sign_hg, name_list)

def make_hg_df(function_names, sign_hg, id_list, df, groups="", metadata=""):
    """ create a dataframe of HG as index and functions as colums with values as means
    ---------
    function_names
        list, function desciptions of organisms
    sign_hg
        list, full names of the relevant HGs
    id_list
        list, GC and HG ids
    df
        dataframe, GC as index, sample IDs as columns and RPKM values
    groups
        list, inputted groups, for instance: [UC, CD]
    metadata
        dict, {sample id: metadata}
    returns
    --------
    relevant_hg = dataframe, index HG names, values RPKMs
    mean_df = dataframe, mean RPKM values of the HGs
    """
    if groups != "":
        row = pd.Series(metadata, name="Metadata")
        df = df.append(row).sort_values(by=["Metadata"], axis=1)
        df = df.replace(0, float(0.0)).T
        df = df.loc[df["Metadata"].isin(groups)].T
        df = df.drop(["Metadata"], axis=0)
    means_dict = ((df[df.index.isin(sign_hg)]).mean(axis=1)).to_dict()
    for gc_name in means_dict.keys():
        for function in function_names:
            for ids in id_list:
                if function in gc_name and ids in gc_name:
                    if type(means_dict[gc_name]) == type(id_list):
                        means_dict[gc_name] = means_dict[gc_name]
                    else:
                        means_dict[gc_name] = [means_dict[gc_name]]
                        means_dict[gc_name].append(function)
                        means_dict[gc_name].append(ids)
    mean_df = (pd.DataFrame(means_dict, index=["Mean", "Function", "ID"])).T
    mean_df = mean_df.pivot(index='ID', columns='Function', values='Mean')
    return mean_df

def sort_coverage(gc_df, cov, metadata):
    """ sort the GCs on coverage from high to low and adds metadata
    --------
    gc_df
        dataframe, filtered GC dataframe with RPKM values
    coverage
        dataframe, filtered coverage scores
    metadata
        dict, {sample id: metadata}
    returns
    --------
    GC = dataframe, sorted on coverage score
    meta_types = list, names of the metadata groups for each column
    """
    cov_mean = cov.mean(axis=1)
    gc_df = pd.concat([gc_df, cov_mean], axis=1, sort=True)
    gc_df = (gc_df.sort_values(by=[0], ascending=False)).drop([0], axis=1)
    # add metadata and sort on the groups
    row = pd.Series(metadata, name="Metadata")
    gc_df = gc_df.append(row).sort_values(by=["Metadata"], axis=1)
    meta_types = list(gc_df.loc["Metadata", :].values)
    gc_df = gc_df.drop(["Metadata"], axis=0)
    gc_df = gc_df.replace(0, float(0.0))
    return gc_df, meta_types

def sort_housekeeping(gc_df, hg_df, cov):
    """ Sorts HG on coverage and appends zeros for missing values
    --------
    gc_df
        dataframe, filtered GC dataframe with RPKM values
    hg_df
        dataframe, filtered HG dataframe with RPKM values
    coverage
        dataframe, filtered coverage scores
    returns
    --------
    HG_new = dataframe, sorted on coverage score
    """
    index_list = []
    sort_id = {}

    # add zero's for non existing HG
    index_names = list(gc_df.index)
    for index in index_names:
        index_list.append(index.split("--OS=")[1].split("--SMASHregion=")[0])
    hg_df = hg_df.replace(np.nan, 0.0)
    index_hg = list(hg_df.index)
    column_hg = len(list(hg_df.columns))
    cov_index = list(cov.index)
    for id_gc in index_list:
        if id_gc not in index_hg:
            hg_df.loc[id_gc] = [0] * column_hg
    # sort the HGs on coverage from high to low
    for cov_id in cov_index:
        for id_gc in index_list:
            if id_gc in cov_id:
                sort_id[id_gc] = cov_id
    df_new = hg_df.rename(index=sort_id)
    cov_mean = cov.mean(axis=1)
    hg_new = pd.concat([df_new, cov_mean], axis=1, sort=True)
    hg_new = (hg_new.sort_values(by=[0], ascending=False)).drop([0], axis=1)
    return hg_new

######################################################################
# FITZIG MODEL
######################################################################

def run_fit_zig(biomfile, groups, data_group, biom_dict):
    """ run metagenomeSeq fitZIG model using an R function and Rpy2
    --------
    biomfile
        string, path to the biomfile
    groups
        list, names of the inputted groups
    data_type
        string, either metagenomic or metatranscriptomic
    data_group
        string, name of of the groups. For instance: DiseaseStatus
    returns
    --------
    fitzig_out = R dataframe containing the fitzig results
    """
    # import and install R packages
    base = importr('base')
    utils = importr('utils')
    utils = rpackages.importr('utils')
    packnames = ('biomformat', 'metagenomeSeq')
    for sample in biom_dict["columns"]:
        data_type = sample["metadata"]["SampleType"]

    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        for i in names_to_install:
            utils.install_packages(i, repos="http://cran.rstudio.com/")
    metagenomeseq = importr('metagenomeSeq')
    biomformat = importr('biomformat')

    r_command = ('''
         run_fitZIG <- function(path, group1, group2, datatype, datagroup){
             ### load and filter the biom file for the inputted groups and data types
             MR <- loadBiom(path)
             MR_sample <- MR[,which(pData(MR)$SampleType==datatype)]
             cut_off <- floor(length(colnames(MRcounts(MR_sample, norm=F, log=T)))*0.25)
             MR_sample <- filterData(MR_sample, present = cut_off)
             MR_mod <- MR_sample[,which(pData(MR_sample)[datagroup]==group1 | \
             pData(MR_sample)[datagroup] == group2)]
             d1 <- pData(MR_mod)[datagroup][,1]
             ### calculate the normalization factor
             normFactor <- calcNormFactors(obj=MR_mod,p=0.5)
             normFactor <- log2(normFactor/median(normFactor$normFactors) + 1)
             mod <- model.matrix(~d1 + normFactor$normFactors)
             pData(MR_mod@expSummary$expSummary)$normFactors = normFactor
             ### run the fitZIG model
             fit <- fitZig(obj = MR_mod, mod = mod, useCSSoffset = F)
             ### perform FDR correction
             MR_coefs <- MRcoefs(fit, by=2, number = length(d1), group = 2, \
             adjustMethod = "BH", alpha = 0,01, taxa = fit@taxa)
             MR_coefs$clust_name = rownames(MR_coefs)
             return(MR_coefs)
         }
         ''')
    # call the R function
    r_pkg = STAP(r_command, "r_pkg")
    try:
        fitzig_out = r_pkg.run_fitZIG(biomfile, groups[0], groups[1], data_type, data_group)
    except subprocess.CalledProcessError:
        print("There has been an error in the R-script, please check the input (meta)data")
        sys.exit()
    return fitzig_out

######################################################################
# KRUSKAL WALLIS
######################################################################
def kruskal_wallis(norm_df, metadata, groups):
    """ performes the kruskal wallis test and corrects the obtained p-values
    using Benjamini Hochberg FDR correction
    --------
    norm_df
        dataframe, normalized GCs
    metadata
        dict, {sample id: metadata}
    groups
        list, names of the inputted groups
    returns
    --------
    fdr_df = dataframe, contains the adjusted P-values for the GCs
    """
    gc_groups = {}
    p_values = []

    group1 = groups[0]
    group2 = groups[1]
    row = pd.Series(metadata, name="Metadata")
    df = norm_df.append(row).sort_values(by=["Metadata"], axis=1)
    df = df.replace(0, float(0.0)).T
    df = df.loc[df["Metadata"].isin(groups)]

    for gc_name in df.columns:
        if "GC_DNA--" in gc_name: # filter out the housekeeping genes
            gc_groups[gc_name] = {}
            for grp in df['Metadata'].unique():
                # make arrays of the groups per GC {GC: {group1: array, group2: array}}
                gc_groups[gc_name][grp] = df[gc_name][df['Metadata'] == grp].values
    # perform Kruskal Wallis test
    for gc in gc_groups.keys():
        no, pval = mstats.kruskalwallis(gc_groups[gc][group1], gc_groups[gc][group2])
        p_values.append(pval)
    fdr = fdrcorrection(p_values, alpha=0.05, method="i")
    fdr_df = pd.DataFrame(data=fdr, columns=gc_groups.keys(), index=["T/F", "pval"]).T
    return fdr_df

######################################################################
# STRUCTURE RESULTS
######################################################################

def parse_results(results, all_counts, groups, metadata, alpha, res_type):
    """ parse the output results of the fitZIG or Kruskal Wallis model
    --------
    results
        dataframe, R or pandas structured dataframe
    all_counts
        dict, normalized GCs and RPKM values
    groups
        list, names of the groups
    metadata
        dict, {sample id: metadata}
    alpha
        float, value to determine significance
    res_type
        string, either fitzig or kruskal
    returns
    --------
    counts = dataframe containing the significant counts according to fitzig
    meta_types = list of metagroups
    """
    meta_types = []
    counts = []
    if res_type == "fitzig":
        gc_names = list(results[3])
        adj_pvals = list(results[2])
    else:
        gc_names = results.index
        adj_pvals = list(results["pval"])

    df = pd.DataFrame(adj_pvals, index=gc_names, columns=["Adj p-values"])
    df = df[df.index.str.contains("GC_DNA--")]
    df = df[df["Adj p-values"] < alpha]
    sign_gc = list(df.index)
    if sign_gc == [] and res_type == "fitzig":
        print("There are no significant gene clusters found. Try increasing the alpha (-af flag)")
        sys.exit()
    elif sign_gc == [] and res_type == "kruskal":
        print("There are no significant gene clusters found. Try increasing the alpha (-ak flag)")
    else:
        row = pd.Series(metadata, name="Metadata")
        counts = all_counts.append(row).sort_values(by=["Metadata"], axis=1)
        counts = counts.replace(0, float(0.0)).T
        counts = counts.loc[counts["Metadata"].isin(groups)].T
        meta_types = list(counts.loc["Metadata", :].values)
        counts = counts.drop(["Metadata"], axis=0)
        counts = counts[all_counts.index.isin(sign_gc)]
    return counts, meta_types

def order_coverage(metadata, cov, groups):
    """ order the coverage scores on the metadata
    --------
    metadata
        dict, {sample id: metadata}
    cov
        dataframe, relavant coverage scores
    returns
    --------
    cov = dataframe, filtered coverage dict on the metadata
    """
    group_meta = {}
    for sample_id in metadata.keys():
        if metadata[sample_id] in groups:
            group_meta[sample_id] = metadata[sample_id]
    row = pd.Series(group_meta, name="Metadata")
    cov = cov.append(row).sort_values(by=["Metadata"], axis=1).T
    cov = cov.replace(0, float(0.0))
    cov = cov.groupby(["Metadata"]).mean().T
    return cov, group_meta

def get_log2fold(gc_df, groups, group_meta):
    """ calculate the log2 fold change
    --------
    gc_df
        dataframe, relevant GCs
    groups
        list, inputted groups
    returns
    --------
    log_fold = dataframe, log2 fold change of the relevant GCs
    """
    groups = sorted(groups)
    group1 = groups[0]
    group2 = groups[1]
    row = pd.Series(group_meta, name="Metadata")
    sign_fit_adj = gc_df.append(row).sort_values(by=["Metadata"], axis=1).T
    sign_fit_adj = sign_fit_adj.replace(0, float(0.0))
    sign_fit_adj = sign_fit_adj.groupby(["Metadata"]).mean().T
    sign_fit_adj["log2_fold"] = sign_fit_adj[group1] - sign_fit_adj[group2]
    log_fold = sign_fit_adj.drop([group1, group2], axis=1)
    return log_fold

def structure_data(gc_df, cov, log_fold):
    """ sort the data for the fitzig model
    --------
    hg_df
        dataframe, relevant HGs
    gc_df
        dataframe, relevant GCs
    cov
        dataframe, relevant coverage
    log_fold
        dataframe, log2 fold change
    returns
    --------
    hg_df = dataframe with the HG names as index and HG RPKMs as values
    sorted_cov = dataframe with GC names as index and coverage scores as values
    log_fold = dataframe with the log2 fold change
    sorted_gc = dataframe with the GC names as index and GC RPKMs as values
    """
    index_list = []
    index_names = list(gc_df.index)
    for index in index_names:
        index_list.append(index.split("|")[1].split(".")[0])
    gc_df = gc_df.replace(0, float(0.0))

    # sort the dataframe on the mean coverage from high to low
    cov_mean = cov.mean(axis=1)
    sorted_gc = pd.concat([gc_df, cov_mean], axis=1, sort=True)
    sorted_gc = (sorted_gc.sort_values(by=[0], ascending=False)).drop([0], axis=1)
    # sort the coverage on the mean coverage
    sorted_cov = pd.concat([cov, cov_mean], axis=1, sort=True)
    sorted_cov = (sorted_cov.sort_values(by=[0], ascending=False)).drop([0], axis=1)
    sorted_cov = pd.melt(sorted_cov.reset_index(), \
    id_vars="index", var_name="group", value_name="cov")
    # sort the log2 fold change on the mean coverage
    sorted_log_fold = pd.concat([log_fold, cov_mean], axis=1, sort=True)
    sorted_log_fold = (sorted_log_fold.sort_values(by=[0], ascending=False)).drop([0], axis=1)
    return sorted_cov, sorted_log_fold, sorted_gc

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
# VISUALIZATION
######################################################################

def make_explore(sign_gc, cov, metadata, metagroup, sample_type, outdir, \
                file_name="explore_heatmap", sign_hg=""):
    """ Creates the explore heatmap in pdf format.
    ----------
    sign_GC
        dataframe, rows contain samples, columns significant GCs, \
        values are normalized RPKM
    cov
        dataframe, rows contain samples, columns GCs, values are coverage scores
    sign_hg
        dataframe, rows contain housekeeping gene functions, columns GCs IDs, \
        values are RPKM
    metadata
        dict, {sample id: metadata}
    metagroup
        string, metagroup name
    sample_type
        string, metagenomic/metatranscriptomic
    outdir
        string, path to output directory
    returns
    ----------
    """
    y_axis_labels = []
    percentage = []

    pdf_file = os.path.join(outdir, f"{file_name}.pdf")
    eps_file = os.path.join(outdir, f"{file_name}.eps")
    index_names = list(sign_gc.index)
    for name in index_names:
        name_function = name.split("Entryname=")[1].split("--SMASH")[0]
        gc_id = name.split("|")[1]
        full_name=f"{name_function}--ID={gc_id}"
        if len(full_name) > 80:
            function, species = name_function.split("--OS=")
            full_name=f"{function}\nOS={species}--ID={gc_id}"
        y_axis_labels.append(full_name)

    with PdfPages(pdf_file) as pdf:
        fig = plt.figure(figsize=(90, 40))
        fig.subplots_adjust(wspace=0.0, hspace=0.0) # space between the plots

        # bar to display the groups
        ax0 = plt.subplot2grid((15, 14), (0, 2), colspan=9, rowspan=1)
        values, counts = np.unique(list(metadata.values()), return_counts=True)
        for count in counts:
            percentage.append(int(count)/sum(counts) * 100)
        df = pd.DataFrame(list(zip(values, percentage)), columns=[metagroup, 'val'])
        group_bar = df.set_index(metagroup).T.plot(kind='barh', stacked=True, \
        ax=ax0, colormap='summer')
        group_bar.axis('off')
        plt.xlim([0, 100])
        plt.title(f"Explore heatmap: {sample_type}", fontsize=60)

        # places legends on different coordinates for metagenomic or metatranscriptomic data
        if sample_type == "METAGENOMIC":
            ax0.text(120, -4, 'Abundance (DNA)', fontsize=40)
            legend = plt.legend(bbox_to_anchor=(1.31, -7), frameon=False, \
            prop={'size': 35}, title=metagroup)
            legend.get_title().set_fontsize('40')
            cbar_ax = fig.add_axes([0.85, .52, .01, .1])
        else:
            ax0.text(125.4, -2.5, 'Expression (RNA)', fontsize=40)
            legend = plt.legend(bbox_to_anchor=(1.37, -8.7), frameon=False, \
            prop={'size': 35}, title=metagroup)
            legend.get_title().set_fontsize('40')
            cbar_ax = fig.add_axes([0.87, .60, .01, .1])

        # heatmap of RPKM values on the first position of the grid
        ax1 = plt.subplot2grid((14, 14), (1, 2), colspan=9, rowspan=14)
        heatmap = sns.heatmap(sign_gc, ax=ax1, cbar_ax=cbar_ax, cmap='viridis', \
        xticklabels=False, linewidths=.05, linecolor='black')
        lower_half, upper_half = plt.ylim() # discover the values for bottom and top
        lower_half += 0.5 # Add 0.5 to the bottom
        upper_half -= 0.5 # Subtract 0.5 from the top
        plt.ylim(lower_half, upper_half)
        ax1.set_yticklabels(labels=y_axis_labels, rotation=0, fontsize=35)

        # coverage barplot in the second position of the grid
        ax2 = plt.subplot2grid((14, 14), (1, 11), colspan=1, rowspan=14)
        coverage = sns.barplot(x=cov.mean(axis=1).sort_values(ascending=False), \
        y=cov.index, ax=ax2, label="Coverage", color="grey")
        plt.xlabel("Coverage", fontsize=40)
        plt.tick_params(labelsize=30)
        coverage.set(yticks=[])

        if sample_type == "METATRANSCRIPTOMIC":
            # housekeeping genes heatmap in the thirth position of the grid
            ax3 = plt.subplot2grid((14, 14), (1, 12), colspan=1, rowspan=14)
            cbar_ax = fig.add_axes([0.87, .40, .01, .1])
            heatmap2 = sns.heatmap(sign_hg, xticklabels=True, \
            linewidths=.01, linecolor='black', yticklabels=False, ax=ax3, \
            cbar_ax=cbar_ax)
            ax3.text(6.2, 8.4, 'Housekeeping \n genes', fontsize=40)
            ax3.set_xticklabels(sign_hg.columns, rotation=90, fontsize=30)
            ax3.set_ylabel('')
        plt.savefig(eps_file, format='eps')
        pdf.savefig()
    return

def make_compare(sign_gc, log_fold, cov, metadata, metagroup, sample_type,\
                outdir, groups, file_name, plot_type, alpha=0.15, sign_hg=""):
    """ Creates the explore heatmap in pdf format
    ----------
    sign_GC
        dataframe, rows contain samples, columns significant GCs, \
        values are normalized RPKM
    cov
        dataframe, rows contain samples, columns GCs, values are coverage scores
    sign_hg
        dataframe, rows contain housekeeping gene functions, columns GCs IDs, \
        values are RPKM
    metadata
        dict, {sample id: metadata}
    metagroup
        string, metagroup name
    sample_type
        string, metagenomic/metatranscriptomic
    outdir,
        string, path to output directory
    groups
        list, names of the groups
    alpha
        float, cut-off value of fitZIG/Kruskal Wallis
    returns
    ----------
    """
    y_axis_labels = []
    percentage = []
    group1 = groups[0]
    group2 = groups[1]

    eps_file = os.path.join(outdir, f"{file_name}.eps")
    pdf_file = os.path.join(outdir, f"{file_name}.pdf")
    index_names = list(sign_gc.index)
    for name in index_names:
        name_function = name.split("Entryname=")[1].split("--SMASH")[0]
        gc_id = name.split("|")[1]
        full_name=f"{name_function}--ID={gc_id}"
        if len(full_name) > 90:
            function, species = name_function.split("--OS=")
            full_name=f"{function}\nOS={species}--ID={gc_id}"
        y_axis_labels.append(full_name)
    with PdfPages(pdf_file) as pdf:
        fig = plt.figure(figsize=(90, 40))
        fig.subplots_adjust(wspace=0.0, hspace=0.0) # space between the plots

        # bar to display the groups
        ax0 = plt.subplot2grid((15, 15), (0, 3), colspan=8, rowspan=1)
        values, counts = np.unique(list(metadata.values()), return_counts=True)
        for count in counts:
            percentage.append(int(count)/sum(counts) * 100)
        df = pd.DataFrame(list(zip(values, percentage)), columns=[metagroup, 'val'])
        group_bar = df.set_index(metagroup).T.plot(kind='barh', stacked=True, \
        ax=ax0, colormap='summer')
        group_bar.axis('off')
        plt.xlim([0, 100])
        if plot_type == "fitzig":
            plt.title(\
            f"fitZIG model: {sample_type} \n p<{alpha} -- {group1} vs {group2}", fontsize=60)
        else:
            plt.title(\
            f"Kruskal Wallis model: {sample_type} \n p<{alpha} -- {group1} vs {group2}", fontsize=60)
        # places legends on different coordinates for metagenomic or metatranscriptomic data
        if sample_type == "METAGENOMIC":
            ax0.text(126.0, -4.3, 'Abundance (DNA)', fontsize=40)
            legend = plt.legend(bbox_to_anchor=(1.40, -6.9), frameon=False, \
            prop={'size': 35}, title=metagroup)
            legend.get_title().set_fontsize('40')
            cbar_ax = fig.add_axes([0.82, .52, .01, .1])
        else:
            ax0.text(140.4, -2.5, 'Expression (RNA)', fontsize=40)
            legend = plt.legend(bbox_to_anchor=(1.54, -10.5), frameon=False, \
            prop={'size': 35}, title=metagroup)
            legend.get_title().set_fontsize('40')
            cbar_ax = fig.add_axes([0.87, .60, .01, .1])

        # heatmap of RPKM values on the first position of the grid
        ax1 = plt.subplot2grid((14, 15), (1, 3), colspan=8, rowspan=14)
        heatmap = sns.heatmap(sign_gc, ax=ax1, cbar_ax=cbar_ax, cmap='viridis', \
        xticklabels=False, linewidths=.05, linecolor='black')
        lower_half, upper_half = plt.ylim() # discover the values for bottom and top
        lower_half += 0.5 # Add 0.5 to the bottom
        upper_half -= 0.5 # Subtract 0.5 from the top
        plt.ylim(lower_half, upper_half)
        ax1.set_yticklabels(labels=y_axis_labels, rotation=0, fontsize=35)

        # log2 fold-change barplot in the second position of the grid
        ax2 = plt.subplot2grid((14, 15), (1, 11), colspan=1, rowspan=14)
        log2_fc = sns.barplot(x=log_fold.mean(axis=1), \
        y=log_fold.index, ax=ax2, label="LFC", color="grey")
        plt.xlabel("Log2 fc", fontsize=40)
        plt.tick_params(labelsize=30)
        log2_fc.set(yticks=[])

        # coverage stripplot in the thirth position of the grid
        ax3 = plt.subplot2grid((14, 15), (1, 12), colspan=1, rowspan=14)
        coverage = sns.stripplot(data=cov, y="index", hue="group", x="cov", ax=ax3, size=20)
        plt.xlabel("cov", fontsize=40)
        plt.ylabel(None)
        plt.tick_params(labelsize=30)
        coverage.set(yticks=[])

        if sample_type == "METATRANSCRIPTOMIC":
            legend2 = plt.legend(loc='center left', bbox_to_anchor=(2.3, 0.32), frameon=False, \
            prop={'size': 35}, title="Coverage", markerscale=3)
            legend2.get_title().set_fontsize('40')
            # housekeeping genes heatmap in the fourth position of the grid
            ax4 = plt.subplot2grid((14, 15), (1, 13), colspan=1, rowspan=14)
            cbar_ax = fig.add_axes([0.87, .40, .01, .1])
            heatmap2 = sns.heatmap(sign_hg, xticklabels=True, \
            linewidths=.01, linecolor='black', yticklabels=False, ax=ax4, \
            cbar_ax=cbar_ax)
            ax4.text(1.7, 0.6, 'Housekeeping \n genes', horizontalalignment='center', \
            verticalalignment='center', transform=ax4.transAxes, fontsize=40)
            ax4.set_xticklabels(sign_hg.columns, rotation=90, fontsize=30)
            ax4.set_ylabel('')
            ax4.set_xlabel('')
        else:
            legend2 = plt.legend(loc='center left', bbox_to_anchor=(1.2, 0.34), frameon=False, \
            prop={'size': 35}, title="Coverage", markerscale=3)
            legend2.get_title().set_fontsize('40')
        plt.savefig(eps_file, format='eps')
        pdf.savefig()
    return

######################################################################
# MAIN
######################################################################
def main():
    """
    Steps
    --------
    1) Load the biom file and convert to a dictionary
    2) Filter the biom dict on the inputted metatype and convert to a pandas dataframe
    3) Normalize the data based on the normalization of metagenomeSeq
    4) Structure the data and create the explore heatmap if requested
    5) Run the Kruskall-Wallis test if requested
    6) Run the fit-ZIG model if requested
    """
    args = get_arguments()

    if not args.compare and not args.explore:
        print("Please use the --compare or --explore flag")
        sys.exit()
    sample_type = (args.sample_type).upper()
    #create output dir if it does not exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Creating json formatted file from biom
    print("__________Loading biom file_______________________")
    json_file = export2biom(args.biom_file, args.outdir)

    with open(json_file, "r") as jfile:
        biom_dict = json.load(jfile)
    # Filtering input data
    sample_index, sample_list, metadata = get_sample_type(sample_type, \
    biom_dict, args.metagroup)
    index_list = filter_rows(biom_dict, sample_index)
    gc_values = get_gc_ids(biom_dict, index_list, sample_index)
    # Making a pandas dataframe and normalize the data
    gc_df = make_dataframe(gc_values, sample_list)
    norm_df = norm_log2_data(gc_df)
    gc_df.to_csv((os.path.join(args.outdir, 'all_RPKMs.tsv')), sep='\t')
    norm_df.to_csv((os.path.join(args.outdir, 'all_RPKMs_norm.tsv')), sep='\t')

    if args.explore:
        print("__________Making explore heatmap__________________")
        # Extract the coverage scores from the BIOM file
        dict_cov = get_coverage(biom_dict, sample_list, list(gc_values.keys()))
        cov = make_dataframe(dict_cov, sample_list)
        cov.to_csv((os.path.join(args.outdir, 'coverage.tsv')), sep='\t')
        sign_gc, sign_cov = best_cov(cov, norm_df, args.threshold)
        # Make the explore heatmap
        if sample_type == "METATRANSCRIPTOMIC":
            gc_names, sign_hg, id_list = get_relevant_hg(sign_gc, norm_df)
            pandas_hg = make_hg_df(gc_names, sign_hg, id_list, norm_df)
            gc_sorted, meta_types = sort_coverage(sign_gc, sign_cov, metadata)
            hg_sorted = sort_housekeeping(sign_gc, pandas_hg, sign_cov)
            hg_sorted.to_csv((os.path.join(args.outdir, 'explore_HGs.tsv')), sep='\t')
            gc_sorted.to_csv((os.path.join(args.outdir, 'explore_GCs.tsv')), sep='\t')
            if args.file_name_explore:
                make_explore(gc_sorted, sign_cov, metadata, args.metagroup, \
                sample_type, args.outdir, args.file_name_explore, hg_sorted)
            else:
                make_explore(gc_sorted, sign_cov, metadata, args.metagroup, \
                sample_type, args.outdir, "explore_heatmap", hg_sorted)
        else:
            gc_sorted, meta_types = sort_coverage(sign_gc, sign_cov, metadata)
            gc_sorted.to_csv((os.path.join(args.outdir, 'explore_GCs.tsv')), sep='\t')
            if args.file_name_explore:
                make_explore(gc_sorted, sign_cov, metadata, args.metagroup, \
                sample_type, args.outdir, args.file_name_explore)
            else:
                make_explore(gc_sorted, sign_cov, metadata, args.metagroup, \
                sample_type, args.outdir)
        if not args.compare:
            os.remove(os.path.join(args.outdir + "BiG-MAP.table.txt"))
        movetodir(args.outdir + os.sep, "tsv-results", ".tsv")

    if args.compare and args.groups:
        group1 = args.groups[0]
        group2 = args.groups[1]
        print("__________Kruskal-wallis model____________________")
        # run the Kruskal Wallis model
        kruskal_results = kruskal_wallis(norm_df, metadata, args.groups)
        # parse the output and get the significant GCs
        sign_kw, meta_types = parse_results(kruskal_results, norm_df, args.groups, \
        metadata, args.alpha_kruskal, "kruskal")
        if args.file_names_compare and len(args.file_names_compare) == 2:
            kruskal_file = args.file_names_compare[0]
            fitzig_file = args.file_names_compare[1]
        elif args.file_names_compare and len(args.file_names_compare) != 2:
            print("Please input two output file names for the Kruskal-Wallis \
            and fitZIG heatmaps. For instance: UCvsCD_kw UCvsCD_fz")
        else:
            kruskal_file = f"{group1}vs{group2}_kw"
            fitzig_file = f"{group1}vs{group2}_fz"
        if meta_types == []:
            pass
        else:
            # get coverage
            dict_cov = get_coverage(biom_dict, list(sign_kw), sign_kw.index)
            cov = make_dataframe(dict_cov, list(sign_kw))
            cov, group_meta = order_coverage(metadata, cov, args.groups)
            # get log2 fold change
            log_fold = get_log2fold(sign_kw, args.groups, group_meta)
            if sample_type == "METATRANSCRIPTOMIC":
                # get the relevant housekeeping genes
                gc_names, sign_hg, id_list = get_relevant_hg(sign_kw, norm_df)
                kw_hg_pandas = make_hg_df(gc_names, sign_hg, id_list, \
                norm_df, args.groups, metadata)
                sorted_cov, log_fold, sorted_gc = structure_data(sign_kw, cov, log_fold)
                kw_hg_pandas = sort_housekeeping(sign_kw, kw_hg_pandas, cov)
                kw_hg_pandas = kw_hg_pandas.replace(np.nan, 0.0)
                kw_hg_pandas.to_csv((os.path.join(args.outdir, \
                f'{group1}vs{group2}_HG_kw.tsv')), sep='\t')
                sorted_gc.to_csv((os.path.join(args.outdir, \
                f'{group1}vs{group2}_GC_kw.tsv')), sep='\t')
                # make the kruskal wallis figure
                make_compare(sorted_gc, log_fold, sorted_cov, group_meta, args.metagroup, \
                sample_type, args.outdir, args.groups, kruskal_file, "kruskal", \
                args.alpha_kruskal, kw_hg_pandas)
            else:
                sorted_cov, log_fold, sorted_gc = structure_data(sign_kw, cov, log_fold)
                sorted_gc.to_csv((os.path.join(args.outdir, \
                f'{group1}vs{group2}_GC_kw.tsv')), sep='\t')
                # make the kruskal wallis figure
                make_compare(sorted_gc, log_fold, sorted_cov, group_meta, args.metagroup, \
                sample_type, args.outdir, args.groups, kruskal_file, "kruskal", args.alpha_kruskal)

        print("__________Fit-ZIG model____________________")
        # run the fitzig model in R
        fitzig_results = run_fit_zig(args.biom_file, args.groups, args.metagroup, biom_dict)
        # parse the output and get the significant GCs
        sign_fit, meta_types = parse_results(fitzig_results, norm_df, args.groups, metadata, \
        args.alpha_fitzig, "fitzig")
        # get coverage
        dict_cov = get_coverage(biom_dict, list(sign_fit), sign_fit.index)
        cov = make_dataframe(dict_cov, list(sign_fit))
        cov, group_meta = order_coverage(metadata, cov, args.groups)
        # get log2 fold change
        log_fold = get_log2fold(sign_fit, args.groups, group_meta)

        if sample_type == "METATRANSCRIPTOMIC":
            # get the relevant housekeeping genes
            gc_names, sign_hg, id_list = get_relevant_hg(sign_fit, norm_df)
            fit_hg_pandas = make_hg_df(gc_names, sign_hg, id_list, norm_df, args.groups, metadata)
            sorted_cov, log_fold, sorted_gc= structure_data(sign_fit, cov, log_fold)
            fit_hg_pandas = sort_housekeeping(sign_fit, fit_hg_pandas, cov)
            fit_hg_pandas = fit_hg_pandas.replace(np.nan, 0.0)
            fit_hg_pandas.to_csv((os.path.join(args.outdir, \
            f'{group1}vs{group2}_HG_fz.tsv')), sep='\t')
            sorted_gc.to_csv((os.path.join(args.outdir, \
            f'{group1}vs{group2}_GC_fz.tsv')), sep='\t')
            # make the fitzig figure
            make_compare(sorted_gc, log_fold, sorted_cov, group_meta, args.metagroup, \
            sample_type, args.outdir, args.groups, fitzig_file, "fitzig", args.alpha_fitzig, fit_hg_pandas)
        else:
            sorted_cov, log_fold, sorted_gc = structure_data(sign_fit, cov, log_fold)
            sorted_gc.to_csv((os.path.join(args.outdir, \
            f'{group1}vs{group2}_GC_fz.tsv')), sep='\t')
            # make the fitzig figure
            make_compare(sorted_gc, log_fold, sorted_cov, group_meta, args.metagroup, \
            sample_type, args.outdir, args.groups, fitzig_file, "fitzig", args.alpha_fitzig)
        os.remove(os.path.join(args.outdir + "/BiG-MAP.table.txt"))
        movetodir(args.outdir + os.sep, "tsv-results", ".tsv")

    if args.compare and not args.groups:
        print("Please profide the group information")
        sys.exit()

if __name__ == "__main__":
    main()
