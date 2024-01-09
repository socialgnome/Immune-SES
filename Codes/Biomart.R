#### This file is a script that is responsible for creating an ensembl dataset that is used for gene name conversions

#### Run Order = 1

rm(list=ls(all=TRUE))
library(biomaRt)
library(stringr)
library(here)
sapiens_ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
saveRDS(sapiens_ensembl, str_c(here("data/"), "sapiens_ensembl.rds"))
