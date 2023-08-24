#### This file is a script that is responsible for performing the Gene Network Analysis for SES

#### Prerequisites: 1. Network_Genes.R

## Setup env
library(here)
dev.off()
source(here("R/startup.R"))


## Biomart (Downloaded 11.12.22)
sapiens_ensembl <- readRDS(str_c(here("data/"), "sapiens_ensembl.rds"))


## Extract Upregulated GRN genes
name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes")
grn = list()
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Up/Genes/"),name[[i]], ".txt"), header = F, sep = "\t")
  grn[[i]] = as.character(df$V1)
}
list = list()
list[[1]] = grn[[1]]
list[[2]] = c(grn[[2]],grn[[3]], grn[[4]])
write.table(unique(list[[2]]),here("Res/Network/SES/Up/Genes/Upstream.txt"), row.names = F, col.names = F, sep = "\t", quote = F)

sets = c("DE Up", "Upstream Up")
list = setNames(list, sets)
lists = list

## Extract Downregulated GRN genes
name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes")
grn = list()
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Down/Genes/"),name[[i]], ".txt"), header = F, sep = "\t")
  grn[[i]] = as.character(df$V1)
}
list = list()
list[[1]] = grn[[1]]
list[[2]] = c(grn[[2]],grn[[3]], grn[[4]])
write.table(unique(list[[2]]),here("Res/Network/SES/Down/Genes/Upstream.txt"), row.names = F, col.names = F, sep = "\t", quote = F)

sets = c("DE Down", "Upstream Down")
list = setNames(list, sets)
lists = c(lists, list)

list = fromList(lists)
elements <- unique(unlist(lists))

rownames(list) = elements


## Attach DE results
df = readRDS(str_c(here("Res/DE/"), "ses_sss_composite" , ".rds"))
## Merge results
list$LFC = df$logFC[match(rownames(list), df$gene)]
list$P = df$adj.P.Val[match(rownames(list), df$gene)]

list = list %>% rownames_to_column(var = "gene")
write.table(list,here("Res/Network/SES/All/Cytoscape/Node.txt"), row.names = F, col.names = T, sep = "\t", quote = F)


### Enrichment of the gene classes
grn = lists
#grn = list()
#grn[[1]] = c(list[[1]], list[[2]])
#grn[[2]] = c(list[[1]], list[[2]])
#grn = setNames(grn, c("Up", "Down"))


## Perform enrichment on the upregulated cluster
grn_ent = list()
conv = list()
for (i in 1:length(grn)) {
  df = biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values =grn[[i]], mart = sapiens_ensembl, useCache = F)
  conv[[i]] = df
  df = as.character(df$entrezgene_id[!is.na(df$entrezgene_id)])
  grn_ent[[i]] = df
}

grn_ent = setNames(grn_ent, names(grn))
eres = compareCluster(grn_ent, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH")
eres = eres@compareClusterResult

sets = names(grn)

eres$genes = 0
for (i in 1:length(eres$Cluster)) {
  eres$genes[[i]] = paste(conv[[which(sets==eres$Cluster[[i]])]]$hgnc_symbol[match(unlist(strsplit(eres$geneID[[i]], "/")),conv[[which(sets==eres$Cluster[[i]])]]$entrezgene_id)],collapse = ",")
}

relation = read.table("~/Projects/Helper files/Reactome/ReactomePathwaysRelation.txt", sep = "\t", stringsAsFactors = F)
names = read.table("~/Projects/Helper files/Reactome/ReactomePathways.txt", sep = "\t", stringsAsFactors = F, quote = "")
names = names[names$V3=="Homo sapiens",]
test = eres
Des = data.frame(Result = test$Description, stringsAsFactors = F)
Des$Level1 = 0
Des$Level2 = 0
Des$Level3 = 0

for (i in 1:length(test$ID)) {
  
  k3 = tryCatch(which(relation$V2 == test$ID[i]),error = function(e) NULL)
  if (length(k3)>0) {
    
    k2 = tryCatch(relation$V1[k3], error = function(e) NULL)
    k3 = tryCatch(relation$V2[k3], error = function(e) NULL)
    
    k1 = tryCatch(which(relation$V2 == k2[1]), error = function(e) NULL)
    k1 = tryCatch(relation$V1[k1], error = function(e) NULL)
    
    k0 = tryCatch(which(relation$V2 == k1[1]), error = function(e) NULL)
    while (length(k0)>0) {
      k3= k2
      k2 = k1
      k1 = tryCatch(relation$V1[k0], error = function(e) NULL)
      k0 = tryCatch(which(relation$V2==k1[1]), error = function(e) NULL)
    }  
    k = c(k1,k2,k3)
    k = k[complete.cases(k)]
    
    
    if (length(names$V2[which(names$V1==k[1])])>0) {
      Des$Level1[i] = names$V2[which(names$V1==k[1])]
    }
    if (length(names$V2[which(names$V1==k[2])])>0) {
      Des$Level2[i] = names$V2[which(names$V1==k[2])]
    }
    if (length(names$V2[which(names$V1==k[3])])>0) {
      Des$Level3[i] = names$V2[which(names$V1==k[3])]
    }
  } else {
    kp = which(names$V1 == test$ID[i])
    if (length(kp)>0) {
      Des$Level1[i] = names$V2[kp]
    }
  }
  
}

Des[Des==0] <- NA
Des$Show = 0
Des$ID = test$ID
test = test[complete.cases(Des[,c(2)]),]
Des = Des[complete.cases(Des[,c(2)]),]
for (i in 1:length(Des$Result)) {
  k = as.character(Des[i,-c(1,5,6)])
  k = k[complete.cases(k)]
  Des$Show[i] = k[length(k)]
}


Final = test
Final$Show = Des$Show
Final$Level1 = Des$Level1
Final$Level2 = Des$Level2
Final$Result = Des$Result
Final = Final[order(Final$Level1, Final$p.adjust),]

FinalT = Final[Final$Level1=="Immune System",]
sel = c()
for (i in 1:length(FinalT$Cluster)) {
  t = unlist(strsplit(FinalT$genes[[i]], ","))
  sel = c(sel,t)
}
sel = unique(sel)
sel = list[list$gene %in% sel,]
sel$Instance = 0
for (i in 1:length(sel$gene)) {
  sel$Instance[[i]] = length(grep(sel$gene[[i]], FinalT$genes))
}
sel$Class = 1
sel$Class[sel$Instance>=5 & sel$Instance<10] = 2
sel$Class[sel$Instance<5] = 3
sel$Class[sel$Instance<2] = 4


write.table(sel,here("Res/Network/SES/All/Cytoscape/Selected_Node.txt"), row.names = F, col.names = T, sep = "\t", quote = F)


### String
## String Database for TF annotations
rownames(sel) = seq(1, length(sel$gene))
string <- STRINGdb$new(version="11.5", species=9606,score_threshold=700, network_type = "full") #Medium confidence
methods = as.data.frame(STRINGdb$methods())

map = string$map(sel, "gene", removeUnmappedRows = T)
map = map[!duplicated(map[,c(1)]),] 

int = string$get_interactions(map$STRING_id)
int = int[int$from %in% map$STRING_id,]
int = int[order(-int$combined_score),]
int = int[!duplicated(int[,c(1,2)]),]
int = int[int$to %in% map$STRING_id,]
int = int[int$combined_score>700,]
int = int[!duplicated(apply(int[,1:2], 1, function(row) paste(sort(row), collapse=""))),] 

int$display = paste0(int$from, " (pp) ", int$to)
int$combined_score = int$combined_score/1000

int$weight = 0 
for (i in 1:length(int$from)) {
  f = map$Instance[which(map$STRING_id == int$from[[i]])]
  t = map$Instance[which(map$STRING_id == int$to[[i]])]
  int$weight[[i]] = mean(f,t)
}

write.table(int,here("Res/Network/SES/All/Cytoscape/Selected_Node_cEdge.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
