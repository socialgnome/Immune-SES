#### This file is a script that is responsible for performing the TF-Gene Connections

#### Prerequisites: 1. DE.R

## Setup env

library(here)
dev.off()

source(here("R/startup.R"))
df = read.table("~/Data/GRN/FANTOM5_individual_networks/394_individual_networks/whole_blood_ribopure.txt.gz")
colnames(df) = c("TF","Gene", "Score")
df_filt = df[df$Score>0.4,] #Medium confidence


## DE 
treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")

for (i in 1:length(treatment)) {
  df = readRDS(str_c(here("Res/DE/"), treatment[[i]],".rds"))
  df = df[df$adj.P.Val<0.05,]
  
  ## TF - DE
  df_filt$TF_LFC = df$logFC[match(df_filt$TF, df$gene)]
  
  ## Gene - DE
  df_filt$Gene_LFC = df$logFC[match(df_filt$Gene, df$gene)]
  #df_filt[is.na(df_filt)] <- 0
  write.table(df_filt, str_c(here("Res/TF_Gene/"), treatment[[i]],".txt"), col.names = T, row.names = F, sep = "\t", quote = F)
}