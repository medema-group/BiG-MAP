#Author: Koen van den Berg
#University: Wageningen University and Research
#Department: Department of Bioinformatics
#Date: 01/07/2019

# Good metagenomeSeq visualisations:
# https://pdfs.semanticscholar.org/0cda/aa11820521bb6cb57fc576eb06c40ddb4650.pdf
# https://jokergoo.github.io/ComplexHeatmap-reference/book/more-examples.html

# import statements
library(metagenomeSeq)
library(biomformat)
library(pheatmap)
library(ggplot2)
library(reshape)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(viridisLite)

# Loading in a biom format:
MR <- load_biom("metaclust.map.biom")
MR_gen <- MR[,which(pData(MR)$SampleType=="METAGENOMIC")]
MR_trans <- MR[,which(pData(MR)$SampleType=="METATRANSCRIPTOMIC")]
disease_gen <- pData(MR_gen)$DiseaseStatus
disease_trans <- pData(MR_trans)$DiseaseStatus

# Normalization (CSS)
## metagenomical samples:
p <- cumNormStat(MR_gen)
MR_gen <- cumNorm(MR_gen, p=p)
## metatranscriptomical samples:
p <- cumNormStat(MR_trans)
MR_trans <- cumNorm(MR_trans, p=p)

# Exporting matrix
## metagenomical samples:
mat_gen <- MRcounts(MR_gen, norm=T, log=T)
mat_rownames <- substr(rownames(mat_gen), regexpr("Entryname=",rownames(mat_gen))+10,regexpr("OS",rownames(mat_gen))-3)
rownames(mat_gen) <- mat_rownames
mat_gen <- mat_gen[,order(disease_gen)]
disease_gen <- disease_gen[order(disease_gen)]

## metatranscriptomical samples:
mat_trans <- MRcounts(MR_trans, norm=T, log=T)
rownames(mat_trans) <- mat_rownames
mat_trans <- mat_trans[,order(disease_trans)]
disease_trans <- disease_trans[order(disease_trans)]

# Statistical testing
## DA testing
### Preparation
#controls = grep("Extraction.Control", pData(MR)$SampleType)
#objtrim = MR[,-controls]
sparseFeatures = which(rowSums(MRcounts(MR_gen)>0)<10)
MRtrim = MR_gen[-sparseFeatures,]
MRp = cumNormStat(MRtrim,pFlag=TRUE,main="trimmed obj data")

### Testing
normFactor <- normFactors(MR_gen)
normFactor <- log2(normFactor/median(normFactor) + 1)
mod <- model.matrix(~disease_gen + normFactor )
fit <- fitZig(obj = MR_gen, mod = mod, useCSSoffset = F)

## (OK) presence/absence testing
res <- fitPA(MR_gen, cl=disease_gen)

## (OK) Kruskal
p.value <- apply(mat_gen, 1, function(x){kruskal.test(x, disease_gen)$p.value})
matp_gen <- mat_gen[which(p.value<0.05),]

p.value <- apply(mat_trans, 1, function(x){kruskal.test(x, disease_trans)$p.value})
matp_trans <- mat_trans[which(p.value<0.05),]




# Exploratory analysis
## (OK) metagenomeSeq heatmap
## Metagenomical samples
heatmapColcolours = brewer.pal(9, "Set1")[as.integer(factor(disease_gen))]
heatmapCols = magma(10)#colorRampPalette(brewer.pal(9,"RdBu"))(50)
plotMRheatmap(mat_gen,
              n=94,
              cexRow=1,
              cexCol=1,
              margins =c(8,38),
              trace = "none",
              col = heatmapCols,
              ColSideColors=heatmapColcolours,
              scale="none",
              Colv=F)

## Metatranscriptomical samples
heatmapColcolours = brewer.pal(9, "Set1")[as.integer(factor(disease_trans))]
heatmapCols = magma(10)#colorRampPalette(brewer.pal(9,"RdBu"))(50)
plotMRheatmap(matp_trans,
              n=94,
              cexRow=1,
              cexCol=1,
              margins =c(8,38),
              trace = "none",
              col = heatmapCols,
              ColSideColors=heatmapColcolours,
              scale="none",
              Colv=F)

# Complex Heatmap
## Metagenomical samples
ha = HeatmapAnnotation(Condition=disease_gen,
                       annotation_name_side = "left")
ht_list = Heatmap(matp_gen[1:50,], 
                  name = "expression", 
                  column_title = 'Metagenomic samples\nSchirmer et al.',
                  #column_km = 3,
                  #row_km = 5, 
                  #col = colorRamp2(c(-2, 0, 2), c("green", "black", "red")),
                  col = magma(9),
                  top_annotation = ha, 
                  row_names_gp = gpar(fontsize = 10),
                  show_column_names = FALSE, 
                  row_title = NULL, 
                  show_row_dend = T,
                  column_order = order(disease_gen) )
draw(ht_list, ht_gap = unit(5, "mm"))

ha = HeatmapAnnotation(Condition=disease_trans, annotation_name_side = "left")
ht_list = Heatmap(matp_trans, 
                  name = "expression", 
                  column_title = 'Metatranscriptomic samples\nSchirmer et al.',
                  #column_km = 3,
                  #row_km = 5, 
                  #col = colorRamp2(c(-2, 0, 2), c("green", "black", "red")),
                  col = magma(9),
                  top_annotation = ha, 
                  row_names_gp = gpar(fontsize = 10),
                  show_column_names = FALSE, 
                  row_title = NULL, 
                  show_row_dend = T,
                  column_order = order(disease_gen) )
draw(ht_list, ht_gap = unit(5, "mm"))






## PlotORD
plotOrd(MR_gen, tran=F, usePCA = T, useDist = T, bg=factor(disease_gen), pch=21)
plotOrd(MR_trans, tran=F, usePCA = T, useDist = T, bg=factor(disease_trans), pch=21)

## Plotcorr
plotCorr(MR_gen, n=50, cexRow=1, cexCol=1, trace="none", dendrogram="none", col=heatmapCols)
plotCorr(MR_trans, n=50, cexRow=1, cexCol=1, trace="none", dendrogram="none", col=heatmapCols)
plotRare(MR_gen, cl=factor(disease_gen), pch=21)
plotRare(MR_trans, cl=factor(disease_trans), pch=21)

## PlotOTU
classindex <- list("UC" = which(disease_gen == "UC"))
classindex$nonIBD <- which(disease_gen == "non-IBD")
classindex$CD <- which(disease_gen == "CD")
otu = 2
plotOTU(MR_gen, otu = otu, classindex)

## Plotgenus
x = fData(MR)$DiseaseStatus[otu]
otulist <- grep(x, fData(MR)$DiseaseStatus)
plotGenus(MR, )

