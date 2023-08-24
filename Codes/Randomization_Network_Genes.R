#### This file is a script that is responsible for performing the randomization of the upstream regulators

#### Prerequisites: 1. Network_Genes.R

## Setup env
rm(list=ls(all=TRUE))
library(here)
#library(mediation)
library(stringr)
#library(purrr)
library(foreach)
library(parallel)
library(doParallel)
library(bigmemory)
library(STRINGdb)


## Load expression data
dat = readRDS(str_c(here("data/"), "/Filtered_ExpressionSet_Clustering_SelSubData.rds"))
genes = rownames(Biobase::exprs(dat))
rm(dat)

## String Database for TF annotations
string <- STRINGdb$new(version="11.5", species=9606,score_threshold=700)
string_net <- STRINGdb$new(version="11.5", species=9606,score_threshold=0)
methods = as.data.frame(STRINGdb$methods())

## TF-Gene database
tf = read.table("Z:/Data/GRN/FANTOM5_individual_networks/394_individual_networks/whole_blood_ribopure.txt.gz")
colnames(tf) = c("TF","Gene", "Score")
tf_filt = tf[tf$Score>0.4,] #Medium confidence
rm(tf)

## Extract true Upregulated GRN genes
name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes")
grn_act = list()
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Up/Genes"),"/",name[[i]], ".txt"), header = F, sep = "\t")
  grn_act[[i]] = as.character(df$V1)
}
sets = c("Set A", "Set B", "Set C", "Set D")
grn_act = setNames(grn_act, sets)


## Upregulated randomization
iter = 1000

random_fun = function(genes, grn_act, tf_filt) {
  DE = sample(genes, size = length(grn_act[[1]]), replace = F)
  
  grn = list()
  grn[[1]] = DE
  
  de_tfs = tf_filt[tf_filt$TF %in% DE,]
  tfs_degene = tf_filt[tf_filt$Gene %in% DE,]
  
  map = string$map(de_tfs, "TF", removeUnmappedRows = T)
  grn[[2]] = unique(map$TF)
  
  neighbors = string$get_neighbors(map$STRING_id)

  if (length(neighbors)>0) {
    int = string_net$get_interactions(union(neighbors,map$STRING_id))
    int = int[int$from %in% map$STRING_id,]
    int = int[order(-int$combined_score),]
    int = int[!duplicated(int[,c(1,2)]),]
    int = int[int$to %in% neighbors,]
    int = int[int$combined_score>700,]
    anno = string_net$add_proteins_description(data.frame(STRING_id = unique(int$to)))
    grn[[3]] = unique(anno$preferred_name)
  } else {
    grn[[3]] = NULL
  }
  
  map = string$map(tfs_degene, "TF", removeUnmappedRows = T)
  grn[[4]] = unique(map$TF)
  
  total = c()
  int = c()
  for (i in 1:length(grn)) {
    total[i] = length(unique(grn[[i]]))
    int[i] = length(intersect(grn[[i]], grn_act[[i]]))
  }
  res = cbind(total,int)
  
  fullres = list()
  fullres[[1]] = grn
  fullres[[2]] = res
  return(fullres)
}

test = random_fun(genes, grn_act, tf_filt)

strt<-Sys.time()
cl <- parallel::makeCluster(24)
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(STRINGdb))

#clusterExport(cl=cl, c('exprs', 'subData','treat','controls','meds'))
fullout = foreach(i=1:iter ,.packages = "foreach",.errorhandling = c("pass")) %dopar% {
  set.seed(i)
  out = random_fun(genes, grn_act, tf_filt)
  out
}
parallel::stopCluster(cl)
print(Sys.time()-strt)
saveRDS(fullout, str_c(here("Res/"), "/Randomization/Results_Up.rds"))

## Extract true Downegulated GRN genes
name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes")
grn_act = list()
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Down/Genes"),"/",name[[i]], ".txt"), header = F, sep = "\t")
  grn_act[[i]] = as.character(df$V1)
}
sets = c("Set A", "Set B", "Set C", "Set D")
grn_act = setNames(grn_act, sets)


## Downregulated randomization
iter = 1000

random_fun = function(genes, grn_act, tf_filt) {
  DE = sample(genes, size = length(grn_act[[1]]), replace = F)
  
  grn = list()
  grn[[1]] = DE
  
  de_tfs = tf_filt[tf_filt$TF %in% DE,]
  tfs_degene = tf_filt[tf_filt$Gene %in% DE,]
  
  map = string$map(de_tfs, "TF", removeUnmappedRows = T)
  grn[[2]] = unique(map$TF)
  
  neighbors = string$get_neighbors(map$STRING_id)
  
  if (length(neighbors)>0) {
    int = string_net$get_interactions(union(neighbors,map$STRING_id))
    int = int[int$from %in% map$STRING_id,]
    int = int[order(-int$combined_score),]
    int = int[!duplicated(int[,c(1,2)]),]
    int = int[int$to %in% neighbors,]
    int = int[int$combined_score>700,]
    anno = string_net$add_proteins_description(data.frame(STRING_id = unique(int$to)))
    grn[[3]] = unique(anno$preferred_name)
  } else {
    grn[[3]] = NULL
  }
  
  map = string$map(tfs_degene, "TF", removeUnmappedRows = T)
  grn[[4]] = unique(map$TF)
  
  total = c()
  int = c()
  for (i in 1:length(grn)) {
    total[i] = length(unique(grn[[i]]))
    int[i] = length(intersect(grn[[i]], grn_act[[i]]))
  }
  res = cbind(total,int)
  fullres = list()
  fullres[[1]] = grn
  fullres[[2]] = res
  return(fullres)
}

test = random_fun(genes, grn_act, tf_filt)

strt<-Sys.time()
cl <- parallel::makeCluster(24)
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(STRINGdb))

#clusterExport(cl=cl, c('exprs', 'subData','treat','controls','meds'))
fullout = foreach(i=1:iter ,.packages = "foreach",.errorhandling = c("pass")) %dopar% {
  set.seed(i)
  out = random_fun(genes, grn_act, tf_filt)
  out
}
parallel::stopCluster(cl)
print(Sys.time()-strt)
saveRDS(fullout, str_c(here("Res/"), "/Randomization/Results_Down.rds"))
