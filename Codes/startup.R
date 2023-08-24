#### This file is a script that is used for setting up the workspace for the analysis in this project

rm(list=ls(all=TRUE))
Sys.setlocale("LC_MESSAGES", "en_US.utf8")
gc()
library("igraph")
library(org.Hs.eg.db)
library(tidyverse)
library(Biobase)
library(edgeR)
library(DESeq2)
library(limma)
library(here)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(VennDiagram)
library(gplots)
library(annotationTools)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(homologene)
library(reshape2)
library(ggvenn)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(hrbrthemes)
library(viridis)
library(stringr)
library(devtools)
library(enrichR)
library(STRINGdb)
library(scales)
library(UpSetR)
library(apeglm)
library(nlme)
library(lme4)
library("BiocParallel")
library(GeneNetworkBuilder)
library(fgsea)
library(ggpubr)
library(EDASeq)
library(sva)
library(RColorBrewer)
library(rcompanion)
library(readr)
library(WGCNA)
library(umap)
library(reticulate)
library(multiClust)
library(NbClust)
library(factoextra)
library(survcomp)
library(recount)
library(globaltest)
library(GEOmetadb)
library(biobroom)
library(compositions)
library(MLSeq)
library("factoextra")
library(SASxport)
library(BioNERO)
library(haven)
library(ontologyIndex)
library(readxl)
library(openxlsx)
library(ggtext)

'%!in%' <- function(x,y)!('%in%'(x,y))

round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
#head(brewer.pal(8, "RdYlGn"),n=1) # Red = Up #D73027
#tail(brewer.pal(8, "RdYlGn"),n=1) # Green = down #1A9850
## Up partners = #FFD320
## Down partners = #0072B2


##### General Notes

### Raw expression data (from Brandt) ###
# Dimensions #
# 4543 subjects and 41 metadata cols. 
# 60672 genes

## Then get the pheno data from preprocess (Wenjia recoded phenotypes)
# Please get vars "Anyflag" and "batch" from the raw files
# check to see if the rownames of pheno data from Wenjia matches the raw
# check to see if the colnames from the exprs from raw matches rownames of pheno

## Then add other pheno data based on preference
