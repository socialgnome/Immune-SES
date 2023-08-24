#### This file is a script that is responsible for performing the randomization analysis of the upstream regulators

#### Prerequisites: 1. Randomization_Analysis.R

## Setup env
library(here)
dev.off()
source(here("R/startup.R"))

## Load expression data
dat = readRDS(str_c(here("data/"), "/Filtered_ExpressionSet_Clustering_SelSubData.rds"))
genes = rownames(Biobase::exprs(dat))
rm(dat)

## String Database for TF annotations
string <- STRINGdb$new(version="11.5", species=9606,score_threshold=700)
string_net <- STRINGdb$new(version="11.5", species=9606,score_threshold=0)
methods = as.data.frame(STRINGdb$methods())

## TF-Gene database
tf = read.table("~/Data/GRN/FANTOM5_individual_networks/394_individual_networks/whole_blood_ribopure.txt.gz")
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

## Extract Randomization results
library(metap)
iter = 1000
res = readRDS(str_c(here("Res/"), "/Randomization/Results_Up.rds"))
sumP = c()
for (i in 1:length(grn_act)) {
  test = lapply(res, `[[`, 1) 
  test = lapply(test, `[[`, i)
  test = table(unlist(test))
  P = c()
  for (j in 1:length(grn_act[[i]])) {
    t = grn_act[[i]][[j]]
    n = test[which(names(test)==t)]
    p = as.numeric((n+1)/(iter+1))
    P = c(P,p)
  }
  sumP = c(sumP, sumlog(P)$p)
}

sumP_Up = sumP


## Extract true Downregulated GRN genes
name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes")
grn_act = list()
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Down/Genes"),"/",name[[i]], ".txt"), header = F, sep = "\t")
  grn_act[[i]] = as.character(df$V1)
}
sets = c("Set A", "Set B", "Set C", "Set D")
grn_act = setNames(grn_act, sets)

## Extract Randomization results
library(metap)
iter = 1000
res = readRDS(str_c(here("Res/"), "/Randomization/Results_Down.rds"))
sumP = c()
for (i in 1:length(grn_act)) {
  test = lapply(res, `[[`, 1) 
  test = lapply(test, `[[`, i)
  test = table(unlist(test))
  P = c()
  for (j in 1:length(grn_act[[i]])) {
    t = grn_act[[i]][[j]]
    n = test[which(names(test)==t)]
    p = as.numeric((n+1)/(iter+1))
    P = c(P,p)
  }
  sumP = c(sumP, sumlog(P)$p)
}

sumP_Down = sumP
