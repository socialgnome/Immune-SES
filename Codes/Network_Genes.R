#### This file is a script that is responsible for identifying the upstream regulators of SES

#### Prerequisites: 1. TF_Gene.R

## Setup env
library(here)
dev.off()
source(here("R/startup.R"))


## Biomart (Downloaded 11.12.22)
sapiens_ensembl <- readRDS(str_c(here("data/"), "sapiens_ensembl.rds"))


## String Database for TF annotations
string <- STRINGdb$new(version="11.5", species=9606,score_threshold=700, network_type = "full") #Medium confidence
string_net <- STRINGdb$new(version="11.5", species=9606,score_threshold=0)
methods = as.data.frame(STRINGdb$methods())


## Load signatures
treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")

## DE 
DE = list()

for (i in 1:length(treatment)) {
  df = readRDS(str_c(here("Res/DE/"), treatment[[i]],".rds"))
  DE[[i]] = df[df$adj.P.Val<0.05,]
}

## TF - DE
tfde = list()
for (i in 1:length(treatment)) {
  df = read.table(str_c(here("Res/TF_Gene/"), treatment[[i]],".txt"), header = T)
  tfde[[i]] = df
}

treatment = c("SES Composite","Subjective Social Status","Occupation","Education","Income")
DE = setNames(DE, treatment)
tfde = setNames(tfde, treatment)


### Upregulated SES GRN
treat = treatment[[1]]
de_data = DE[[1]]
de_data = de_data[de_data$logFC>0,]
tfde_data = tfde[[1]]

grn = list()
grn[[1]] = de_data$gene

de_tfs = rbind(tfde_data[which(!is.na(tfde_data$TF_LFC) & tfde_data$TF_LFC>0),], ## TF's that are DE and Up
                tfde_data[which(tfde_data$TF_LFC<0 & tfde_data$Gene_LFC>0),]) ## TF's that are Down but control Up genes
tfs_degene = tfde_data[!is.na(tfde_data$Gene_LFC) & tfde_data$Gene_LFC>0 & is.na(tfde_data$TF_LFC),]  ## TF's connected to up genes and are not DE themselves

map = string$map(de_tfs, "TF", removeUnmappedRows = T)
grn[[2]] = unique(map$TF)

neighbors = string$get_neighbors(map$STRING_id)
int = string_net$get_interactions(union(neighbors,map$STRING_id))
int = int[int$from %in% map$STRING_id,]
int = int[order(-int$combined_score),]
int = int[!duplicated(int[,c(1,2)]),]
int = int[int$to %in% neighbors,]
int = int[int$combined_score>700,]

anno = string_net$add_proteins_description(data.frame(STRING_id = unique(int$to)))
grn[[3]] = setdiff(setdiff(unique(anno$preferred_name), grn[[1]]), grn[[2]])

grn[[4]] = setdiff(unique(tfs_degene$TF), grn[[3]])

name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes")
grn = setNames(grn, name)
for (i in 1:length(grn)) {
  write.table(grn[[i]],str_c(here("Res/Network/SES/Up/Genes/"),name[[i]], ".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
}
write.table(unique(unlist(grn)),here("Res/Network/SES/Up/Genes/All.txt"), row.names = F, col.names = F, sep = "\t", quote = F)


### Downregulated SES GRN
treat = treatment[[1]]
de_data = DE[[1]]
de_data = de_data[de_data$logFC<0,]
tfde_data = tfde[[1]]

grn = list()
grn[[1]] = de_data$gene

de_tfs = rbind(tfde_data[which(!is.na(tfde_data$TF_LFC) & tfde_data$TF_LFC<0),], ## TF's that are DE and Down
               tfde_data[which(tfde_data$TF_LFC>0 & tfde_data$Gene_LFC<0),]) ## TF's that are Up but control Down genes
tfs_degene = tfde_data[!is.na(tfde_data$Gene_LFC) & tfde_data$Gene_LFC<0 & is.na(tfde_data$TF_LFC),]  ## TF's connected to up genes and are not DE themselves

map = string$map(de_tfs, "TF", removeUnmappedRows = T)
grn[[2]] = unique(map$TF)

neighbors = string$get_neighbors(map$STRING_id)
int = string_net$get_interactions(union(neighbors,map$STRING_id))
int = int[int$from %in% map$STRING_id,]
int = int[order(-int$combined_score),]
int = int[!duplicated(int[,c(1,2)]),]
int = int[int$to %in% neighbors,]
int = int[int$combined_score>700,]

anno = string_net$add_proteins_description(data.frame(STRING_id = unique(int$to)))
grn[[3]] = setdiff(setdiff(unique(anno$preferred_name), grn[[1]]), grn[[2]])

grn[[4]] = setdiff(unique(tfs_degene$TF), grn[[3]])

name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes")
grn = setNames(grn, name)
for (i in 1:length(grn)) {
  write.table(grn[[i]],str_c(here("Res/Network/SES/Down/Genes/"),name[[i]], ".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
}
write.table(unique(unlist(grn)),here("Res/Network/SES/Down/Genes/All.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
