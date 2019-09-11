#!/usr/bin/env Rscript

#Author: Koen van den Berg
#University: Wageningen University and Research
#Department: Department of Bioinformatics
#Date: 01/07/2019

# Good metagenomeSeq visualisations:
# https://pdfs.semanticscholar.org/0cda/aa11820521bb6cb57fc576eb06c40ddb4650.pdf
# https://jokergoo.github.io/ComplexHeatmap-reference/book/more-examples.html

# On methods
# https://ieeexplore.ieee.org/document/6741402 (10.1109/BIGCOMP.2014.6741402)
# https://rachaellappan.github.io/16S-analysis/differential-abundance-with-metagenomeseqs-fitzig.html#fitzig-models

#source("https://bioconductor.org/biocLite.R")

####################################################
# import statements
####################################################
packages = c("metagenomeSeq","biomformat","ComplexHeatmap","viridisLite", "RColorBrewer","dplyr", "stringr")
package.check <- lapply(packages, FUN = function(x) {
    #if (!require(x, character.only = TRUE)) {
    #    print("not installed") #install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
    #suppressWarnings(suppressMessages(library(x, character.only = TRUE, quietly = T)))
    #}
})

####################################################
# Functions
####################################################
writestdout <- function(...) cat(sprintf(...), sep='', file=stderr())

CorrectHousegenes <- function(hgenes.df){
  # Duplicates housekeeping gene rowvalues for identical organisms
  # -parameters-
  # hgenes.df = dataframe containing housekeeping genes 
  # -outputs-
  # Corrected dataframe containing housekeeping genes
  ret <- hgenes.df
  organisms <- substr(rownames(hgenes.df),
                      regexpr("--OS=", rownames(hgenes.df))+5,
                      regexpr("--SMASH",rownames(hgenes.df))-1)
  for (row in rownames(hgenes.df)){
    avg.check <- rowMeans(hgenes.df[which(rownames(hgenes.df) %in% row),])
    if (avg.check > 0){
      organism <- substr(row,
                         regexpr("--OS=", row)+5,
                         regexpr("--SMASH",row)-1)
      alike.rows <- which(organisms %in% organism)
      ret[alike.rows,] <- hgenes.df[which(rownames(hgenes.df) %in% row),]
    }    
  }
  return(ret)
}

FindHouseGenes <- function(significant_hits, countsdf) {
  # Finds the relevant housekeeping genes for significant hits
  # -parameters-
  # significant_hits: sig hits from fitZIG or kruskall model
  # countsdf: dataframe containing all the counts
  # -outputs-
  # dataframe containing avg counts for each housekeeping gene per gene cluster

  # Finding the significant hits first, based on their unique fasta IDs
  significant_hits <- data.frame(significant_hits)
  countsdf <- data.frame(countsdf)
  hits_IDs <- substr(rownames(significant_hits),
                     0,
                     regexpr("DNA--", rownames(significant_hits))-1)
  relevant_hgenes <- countsdf[which(substr(rownames(countsdf), 0, regexpr("DNA--", rownames(significant_hits))-1) %in% hits_IDs),] # HousekeepingGenes
  relevant_hgenes <- relevant_hgenes %>%
    filter(str_detect(relevant_hgenes$rowname, "HousekeepingGene")) # Housekeeping genes later
  rownames(relevant_hgenes) <- relevant_hgenes$rowname
  relevant_hgenes$rowname <- NULL


# Reshaping the relevant hits to rows
  if (nrow(relevant_hgenes) > 0){
    h.avg <- rowMeans(relevant_hgenes) # Taking the average expression values
    hname <- substr(rownames(relevant_hgenes), 
                    regexpr("Entryname=",rownames(relevant_hgenes))+10,
                    regexpr("--OS",rownames(relevant_hgenes))-1)
    ID <- substr(rownames(relevant_hgenes),
                 0,
                 regexpr("HousekeepingGene", rownames(relevant_hgenes))-1) # later HousekeepingGene
    h.avg <- data.frame(cbind(h.avg, ID, hname))
    h.avg.hits <- reshape(h.avg, timevar = "hname", idvar = "ID", direction = "wide") # Creating a dataframe containing the clusters as rows and the hgene names as columns
    names_significant_hits <- data.frame(rownames(significant_hits))
    names_significant_hits$ID <- hits_IDs
    
    # Processing result
    result <- merge(names_significant_hits, h.avg.hits, by="ID", all = T)
    rownames(result) <- result$rownames.significant_hits. 
    result$ID <- NULL
    result$rownames.significant_hits. <- NULL
    ret <- mutate_all(result, function(x) as.numeric(as.character(x)))
    ret[is.na(result)] <- 0
    rownames(ret) <- rownames(result)
    ret <- ret[match(rownames(significant_hits), rownames(ret)),]
    return(ret)
  } else {
    return(data.frame())
  }
}

makeZIGmodel <- function(MRobj, meta, groups){
  # Makes a fitzig model and find DA genes
  # -parameters-
  # MRobj: MR experiment object from metagenomeSeq
  # meta: name of the metadata group condition (like DiseaseStatus)
  # groups: vector of two group conditions c("G1","G2")
  # -outputs-
  # result$data: the gaussian modelled counts dataframe
  # result$meta: the metadata groups
  # result$log2fold: logfold ratios between groups
  # result$cov: core cluster coverage values 
  # result$p-value: the adjusted pvalues 
  
  # Making the fitZig model:
  MR_mod <- MRobj[,which(pData(MRobj)[meta]==groups[1] | pData(MRobj)[meta] == groups[2])]
  d1 <- pData(MR_mod)[meta][,1]
  normFactor <- normFactors(MR_mod)
  normFactor <- log2(normFactor/median(normFactor) + 1) # Error value in model
  mod <- model.matrix(~d1 + normFactor ) # model matrix
  fit <- fitZig(obj = MR_mod, mod = mod, useCSSoffset = F)
  
  # Parsing DA genes:
  MR_coefs <- MRcoefs(fit, by=2, number = 30, group = 3)
  MR_coefs <- MR_coefs[order(MR_coefs$adjPvalues),]      
  countsdf <- data.frame(MRcounts(MR_mod, norm=T, log=T))
  countsdf$rowname <- rownames(countsdf)
  
  # Filtering housekeeping genes out:
  filtered_countsdf <- countsdf %>%
    filter(str_detect(countsdf$rowname, "DNA--"))
  rownames(filtered_countsdf) <- filtered_countsdf$rowname
  filtered_countsdf$rowname <- NULL
  countsmat <- as.matrix(filtered_countsdf)
  sig_hits <- countsmat[which(MR_coefs$adjPvalues<0.05),]


  # Finding the relevant housekeeping genes counts:
  if (length(sig_hits)>0){
    hgenes <- FindHouseGenes(sig_hits, countsdf)
    hgenes <- CorrectHousegenes(hgenes)
    # Filtering coverage values:
    cov = data.frame(fData(MR_mod))
    covd <- data.frame(sapply(cov, function(x) as.numeric(as.character(x))))
    rownames(covd) <- rownames(cov)
    covd = covd[which(rownames(covd) %in% rownames(sig_hits)),]
    covd <- covd[,which(colnames(covd) %in% colnames(sig_hits))]
    avg_cov <- cbind(rowMeans(covd[,which(d1==groups[1])]), rowMeans(covd[,which(d1==groups[2])]))
    
    # Calculating fold-change:
    logfold <- rowMeans(sig_hits[,which(d1==groups[1])]) - rowMeans(sig_hits[,which(d1==groups[2])])
    
    result <- list("data"=sig_hits, "meta"=d1, "log2fold"=logfold, "coverage"=avg_cov, "p-value"=MR_coefs, "hgenes"=hgenes)
    return(result)
  } else {
    return(list("data"=data.frame(), "p-value"=MR_coefs))
  }
}

makekruskalltest <- function(MRobj, meta, groups, alpha){
  # Makes a Kruskal non-parametric model and find DA genes
  # -parameters-
  # MRobj: MR experiment object from metagenomeSeq
  # meta: name of the metadata group condition (like DiseaseStatus)
  # groups: vector of two group conditions c("G1","G2")
  # alpha: significance rate for type I error
  # -outputs-
  # result$data: the non-parametric modelled counts dataframe
  # result$meta: the metadata groups
  # result$log2fold: logfold ratios between groups
  # result$cov: core cluster coverage values 
  # result$p-value: the adjusted pvalues 
  
  # Making the Kruskall-Wallis non-parametric model:
  MR_mod <- MRobj[,which(pData(MRobj)[meta]==groups[1] | pData(MRobj)[meta] == groups[2])]
  countsdf <- data.frame(MRcounts(MR_mod, norm=T, log=T))
  countsdf$rowname <- rownames(countsdf)

  # Filtering housekeeping genes out:
  filtered_countsdf <- countsdf %>%
    filter(str_detect(countsdf$rowname, "DNA--"))
  rownames(filtered_countsdf) <- filtered_countsdf$rowname
  filtered_countsdf$rowname <- NULL
  countsmat <- data.frame(filtered_countsdf)

  # statistical testing:
  d1 <- as.factor(pData(MR_mod)[meta][,1])
  p.value <- apply(countsmat, 1, function(x){kruskal.test(x, d1)$p.value})
  adj.pvalue <- p.adjust(p.value, method = "BH")
  sig_hits <- countsmat[which(adj.pvalue<alpha),]

  # Finding the relevant housekeeping genes counts:
  if (nrow(sig_hits)>0) {
    hgenes <- FindHouseGenes(sig_hits, countsdf)
    hgenes <- CorrectHousegenes(hgenes)
    # Filtering coverage values:
    cov = data.frame(fData(MR_mod))
    covd <- data.frame(sapply(cov, function(x) as.numeric(as.character(x))))
    rownames(covd) <- rownames(cov)
    covd <- covd[which(rownames(covd) %in% rownames(sig_hits)),]
    covd <- covd[,which(colnames(covd) %in% colnames(sig_hits))]
    avg_cov <- cbind(rowMeans(covd[,which(d1==groups[1])]), rowMeans(covd[,which(d1==groups[2])]))

  # Calculating fold-change:
    logfold <- rowMeans(sig_hits[,which(d1==groups[1])]) - rowMeans(sig_hits[,which(d1==groups[2])])
    result <- list("data"=sig_hits, meta=d1, "log2fold"=logfold, "coverage"=avg_cov, "p-value"=adj.pvalue, "hgenes"=hgenes)
    return(result)
  } else {
  return(list("data"=data.frame(), "p-value"=adj.pvalue))
  }
}

makecomplexheatmap <- function(test_result, title, samplename){
  # Makes a complex heatmap from test output of makeZIGmodel
  # or makekruskalltest
  # -parameters-
  # test_result: output from tests as input here
  # title: string, the title of the plot
  # -output-
  # complex heatmap
  ha = HeatmapAnnotation(Condition=test_result$meta)
                         #annotation_name_side = "left")
  ha_row = rowAnnotation(lfc=row_anno_barplot(test_result$log2fold),
                         cov=row_anno_points(matrix(test_result$coverage,
                                                    nc = 2),
                                             pch = 16:17,
                                             gp = gpar(col=2:3)))
  cluster.name <- substr(rownames(test_result$data),
                         regexpr("Entryname=",rownames(test_result$data))+10, 
                         regexpr("SMASH",rownames(test_result$data))-3)					     
  ha_row_left = rowAnnotation(row_names=anno_text(cluster.name, location = 1, just = "right"))
  ht_main = Heatmap(test_result$data,
                    name = samplename, 
                    column_title = title,
                    #column_km = 2,
                    #row_km = 3, 
                    #col = colorRamp2(c(0, 15), c("white", "deepskyblue")),
                    col = viridis(9),
                    top_annotation = ha, 
                    right_annotation = ha_row,
		    left_annotation = ha_row_left,
                    #row_names_gp = gpar(fontsize = 7),
                    show_column_names = FALSE,
		    show_row_names = FALSE,
                    row_title = NULL, 
                    show_row_dend = F,
                    column_order = order(test_result$meta),
                    #row_names_side = "left",
                    row_dend_side = "right",
                    cluster_rows = T)
  lgd_main = Legend(labels = c(unique(test_result$meta[1]),unique(test_result$meta)[2]), 
                    title = "cov",
                    type = "points",
                    pch = 16:17,
                    legend_gp = gpar(col=2:3))
  if (nrow(test_result$hgenes) > 0) {
     hg.names <- substr(colnames(test_result$hgene), 7, 100)
     ha_hg <- HeatmapAnnotation(col_names=anno_text(hg.names, location = 1, just = "right"), 
                               gp = gpar(fontsize=6)) 			       

     ht_house = Heatmap(test_result$hgenes,
                     name = "hgenes",
                     show_row_dend = F,
		     bottom_annotation = ha_hg,
		     show_column_names = F,
                     col = magma(9),
                     column_title = "HGs",
                     show_row_names = F,
                     column_names_gp = gpar(fontsize=7),
		     width = unit(2,"cm"))
     ht_list = ht_main + ht_house
     draw(ht_list, main_heatmap=samplename, ht_gap = unit(1, "mm"),  annotation_legend_list=lgd_main)
  } else {
  draw(ht_main, main_heatmap=samplename, ht_gap = unit(1, "mm"),  annotation_legend_list=lgd_main)
  }
}

####################################################
# MAIN
####################################################
# The expected input args are:
# $1 the input biom file name
# $2 Metagroup
# $3 sampletype 
# $4 Group 1
# $5 Group 2
# $6 Kruskall/fitZIG
args <- commandArgs(trailingOnly = T)

biom_file <- args[1]
sampletype <- args[2]
outdir = args[3]
MT <- args[4]
group_1 <- args[5]
group_2 <- args[6]

# If not using the command_line version, insert your own data above
# instead of the command line arguments.
#biom_file <- "gut.metaclust.map.biom" 
#MT <- "DiseaseStatus"
#sampletype <- "METAGENOMIC"
#group_1 <- "CD"
#group_2 <- "non-IBD"
#test <- "fitzig"

if (sampletype == "METAGENOMIC"){
   sample.name <- "Abundance (DNA)"
} else {
  sample.name <- "Expression (RNA)"
}

##################################
# Loading in a biom format:
##################################
# MR
MR <- loadBiom(biom_file)
MR_sample <- MR[,which(pData(MR)$SampleType==sampletype)]


##################################
# Normalization (CSS)
##################################
# filters the data so that >25% of the samples do have a hit. 
cut_off <- floor(length(colnames(MRcounts(MR_sample, norm=F, log=T)))*0.25)

MR_sample <- filterData(MR_sample, present = cut_off)
p <- cumNormStat(MR_sample)
MR_sample <- cumNorm(MR_sample, p=p)

##################################
# Differential expresssion analysis
##################################
# Moderated t-tests on gaussian modelled data using Fitzig:
ZIGtest_sample <- makeZIGmodel(MR_sample, MT, c(group_1, group_2))

# Showcasing normality after Fitzig
#par(mfrow=c(2,1))
#ggqqplot(data = Kg_CD_nonIBD$data[1,], title = "Kruskall Model")
#ggqqplot(Zg_CD_nonIBD$data[1,], title = "fitZIG model")

# Kruskall-Wallis testing
Kruskalltest_sample <- makekruskalltest(MR_sample, MT, c(group_1, group_2), 0.05)

# Post-Hoc testing with: 
# Dunns test, Conover-Iman test, Dwass-Steel-Citchlow-Flinger test


##################################
# Exploratory analysis
##################################
# MR gen
if (nrow(ZIGtest_sample$data) > 0){
  png(sprintf("ZIGmodel_%s_%s_%s.png",group_1, group_2, sampletype), width = 1400)
  makecomplexheatmap(ZIGtest_sample, sprintf("%s samples fitZIG model\np<0.05 -- %s vs %s", sampletype, group_1, group_2), sample.name)
  dev.off()
} else{
  writestdout("There are no significant differentially abundant gene clusters found! Therefore no heatmap could be produced. Consult the p-value csv output.\nSampletype: %s\nModel: fitZIG model\nMetagroup: %s\nGroups: %s and %s\n", sampletype, MT, group_1, group_2)
}

if (nrow(Kruskalltest_sample$data) > 0){
  png(sprintf("Kruskallmodel_%s_%s_%s.png",group_1,group_2,sampletype), width = 1400)
  makecomplexheatmap(Kruskalltest_sample, sprintf("%s samples Kruskall-Wallis model\np<0.05 -- %s vs %s", sampletype, group_1, group_2), sample.name)
  dev.off()
} else{
  writestdout("There are no significant differentially abundant gene clusters found! Therefore no heatmap could be produced. Consult the p-value csv output.\nSampletype: %s\nModel: Kruskall-Wallis model\nMetagroup: %s\nGroups: %s and %s\n", sampletype, MT, group_1, group_2)
}


##################################
# Writing pvalue results
##################################
write.csv(ZIGtest_sample$`p-value`, file = sprintf("ZIGtest_pvals_%s_%s_%s.csv", group_1, group_2, sampletype))
write.csv(Kruskalltest_sample$`p-value`, file = sprintf("Kruskalltest_pvals_%s_%s_%s.csv", group_1, group_2, sampletype))