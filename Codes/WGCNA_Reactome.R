#### This file is a script that is responsible for the analysing the functional 
#### roles of each cluster that is identified by WGCNA

#### This script also creates Supplementary Figure S1, S4 and Dataset S2


#### Prerequisites: 1. WGCNA_Analysis.R

## Setup env

library(here)
dev.off()
source(here("R/startup.R"))

## Load clustering results
load(str_c(here("Res/WGCNA/"), "WGCNA_Res_Enrichment.RData"))

data = data.frame(Cluster = rep(mod_nos$Cluster,3), Value = c(mod_nos$Up,mod_nos$Down,rep(0, length(mod_nos$Cluster))),
                  Reg = c(rep("Upregulated", length(mod_nos$Cluster)), rep("Downregulated", length(mod_nos$Cluster)),rep("Non-DE", length(mod_nos$Cluster))))
data$Value[data$Reg=="Non-DE"] = mod_nos$Total-(mod_nos$Down+mod_nos$Up)
data$Reg = factor(data$Reg, levels = c("Downregulated", "Upregulated","Non-DE"))
## Mod plot
tiff(str_c(here("Res/WGCNA/"), "WGCNA_Numbers.tiff"), units="px", width=(6*800), height=(6*500), res=600)
p = ggplot(data, aes(x = Cluster, y = Value, fill = Reg)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.9) +
  theme_bw(base_size = 12) +
  theme(panel.spacing =unit(0.05, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth = 1,color='#476b6b')) +
  scale_fill_manual(values = c(brewer.pal(10, "RdYlGn")[10],brewer.pal(10, "RdYlGn")[1],"#669999"), name = "Regulation") +
  xlab("WGCNA Cluster Label") + ylab("Number of Genes") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri")) +
  theme(legend.text=element_text(size=12, family = "Calibri")) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position="bottom") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4000)) 
p
dev.off()

## Find enriched pathways
glists = split(moduleInfo$Gene, moduleInfo$Module)

## Biomart (Downloaded 11.12.22)
sapiens_ensembl <- readRDS(str_c(here("data/"), "sapiens_ensembl.rds"))

glists_ent = list()
conv = list()
for (i in 1:length(glists)) {
  df = biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values =glists[[i]], mart = sapiens_ensembl)
  conv[[i]] = df
  df = as.character(df$entrezgene_id[!is.na(df$entrezgene_id)])
  glists_ent[[i]] = df
}

clusternames = names(glists)

glists_ent = setNames(glists_ent, names(glists))
eres = compareCluster(glists_ent, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH")
eres = eres@compareClusterResult

eres$genes = 0
for (i in 1:length(eres$Cluster)) {
  eres$genes[[i]] = paste(conv[[which(clusternames==eres$Cluster[[i]])]]$hgnc_symbol[match(unlist(strsplit(eres$geneID[[i]], "/")),conv[[which(clusternames==eres$Cluster[[i]])]]$entrezgene_id)],collapse = ",")
}

eres$Cramer = 0
for (i in 1:length(eres$Cluster)) {
  n1 = as.numeric(strsplit(eres$GeneRatio[[i]], "/")[[1]][[1]])
  n2 = as.numeric(strsplit(eres$GeneRatio[[i]], "/")[[1]][[2]])
  n3 = as.numeric(strsplit(eres$BgRatio[[i]], "/")[[1]][[1]])
  n4 = as.numeric(strsplit(eres$BgRatio[[i]], "/")[[1]][[2]])
  dd = matrix(c(n1,n3-n1,n2-n1,n4-n2-n3+n1),nrow = 2,ncol = 2)
  eres$Cramer[i] = cramerV(dd)
}

lists = split(eres$ID,eres$Cluster)

## Contexualize reactome enrichment results
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
#Des = Des[,-c(4)]
Des$Show = 0
Des$ID = test$ID
test = test[complete.cases(Des[,c(2)]),]
Des = Des[complete.cases(Des[,c(2)]),]
for (i in 1:length(Des$Result)) {
  k = as.character(Des[i,-c(1,5,6)])
  k = k[complete.cases(k)]
  Des$Show[i] = k[length(k)]
}

s1 = Des
Final = test
Final$Show = Des$Level1
Final$Level1 = Des$Level1
Final$Level2 = Des$Level2
Final$Result = Des$Result
Final = Final[order(Final$Level1, Final$p.adjust),]
Final$FinalCount = 0
Final$Terms = 0
Final$Genes = 0
Final$MedC = 0
Final$MaxC = 0

for (i in 1:length(Final$ID)) {
  m = which(Final$Show==Final$Show[i] & Final$Cluster==Final$Cluster[i])
  k = strsplit(Final$geneID[m],"/")
  l = unique(unlist(k))
  Final$FinalCount[i] = length(l)
  Final$Terms[i] = length(m)
  Final$Genes[i] = paste(l, collapse = "/")
  Final$MedC[i] = median(Final$Cramer[m])
  Final$MaxC[i] = max(Final$Cramer[m])
}



# Write tabular results
tab = data.frame(Treatment = Final$Cluster, ReactomeID = Final$ID, Pathway = Final$Description,
                 ParentNode = Final$Level1, ChildNode =  Final$Level2, GeneRatio = Final$GeneRatio,
                 BgRatio = Final$BgRatio, Adj.P.Val = Final$p.adjust, CramerV = Final$Cramer,
                 CollapsedCramerV = Final$MedC,CollapsedCount = Final$FinalCount, Genes = Final$genes)
tab = tab[order(tab$Treatment),]
write.xlsx(tab, str_c(here("Res/Tables/"), "Reactome_WGCNA_All.xlsx"), overwrite = T)
tab = tab[,c(1,2,3,4,8,9,12)]
tab$Adj.P.Val  = sprintf("%.3f x 10^(%d)", tab$Adj.P.Val/10^floor(log10(abs(tab$Adj.P.Val))), floor(log10(abs(tab$Adj.P.Val))))
tab$CramerV  = sprintf("%.3f x 10^(%d)", tab$CramerV/10^floor(log10(abs(tab$CramerV))), floor(log10(abs(tab$CramerV))))
write.xlsx(tab, str_c(here("Res/Tables/"), "Reactome_WGCNA_Concise.xlsx"), overwrite = T)

## Plot reactome contextualized results


Final = Final[!duplicated(Final[,c(1,13)]),] #Cluster, Show
Final$FinalSize = 7

FinalT = Final
FinalT = droplevels(FinalT)
FinalT$Neg = -log(FinalT$p.adjust)
try = as.character(unique(FinalT$Show))
FinalT$Show = factor(FinalT$Show, levels = unique(FinalT$Show))
FinalT$Show = str_wrap(FinalT$Show, width = 60)
FinalT$Level1 = str_wrap(FinalT$Level1, width = 25)

FinalT$Cluster  = factor(FinalT$Cluster, levels =  clusternames)


tiff(str_c(here("Res/WGCNA/"), "WGCNA_Reactome.tiff"), units="px", width=(6*900), height=(6*700), res=600)

p <- ggplot(FinalT, aes(x = Cluster, y = Show, size = FinalCount, color = MedC, fill = MedC)) +
  geom_point(shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  scale_color_viridis(alpha = 0.9,
                      breaks = c(min(FinalT$MedC),(min(FinalT$MedC)+ max(FinalT$MedC))/2, max(FinalT$MedC)),
                      labels = c("0.02", "0.2", "0.4"),
                      discrete = F, name = "Cramer's V\n")+
  scale_fill_viridis(alpha = 0.9,
                     breaks = c(min(FinalT$MedC),(min(FinalT$MedC)+ max(FinalT$MedC))/2, max(FinalT$MedC)),
                     labels = c("0.02", "0.2", "0.4"),
                     discrete = F, name = "Cramer's V\n")+
  scale_size_continuous(range = c(3,7), name = "Gene Count", breaks = c(min(FinalT$FinalCount),10-0,max(FinalT$FinalCount)), labels = c("3","100", "200")) + 
  #facet_grid(Level1~.,scales="free",space="free")+
  theme(strip.text.y = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.05, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth = 1,color='#476b6b')) +
  scale_y_discrete(limits = rev) +
  xlab("WGCNA Cluster Labels") + ylab("Reactome Pathways") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri")) +
  theme(legend.text=element_text(size=12, family = "Calibri")) +
  scale_x_discrete(drop=FALSE)+
  theme(legend.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(legend.position="right") 

p
dev.off()


Des = s1
Final = test
Final$Show = Des$Show
Final$Level1 = Des$Level1
Final$Level2 = Des$Level2
Final$Result = Des$Result
Final = Final[order(Final$Level1, Final$p.adjust),]
Final$FinalCount = 0
Final$Terms = 0
Final$Genes = 0
Final$MedC = 0
Final$MaxC = 0

for (i in 1:length(Final$ID)) {
  m = which(Final$Show==Final$Show[i] & Final$Cluster==Final$Cluster[i])
  k = strsplit(Final$geneID[m],"/")
  l = unique(unlist(k))
  Final$FinalCount[i] = length(l)
  Final$Terms[i] = length(m)
  Final$Genes[i] = paste(l, collapse = "/")
  Final$MedC[i] = median(Final$Cramer[m])
  Final$MaxC[i] = max(Final$Cramer[m])
}
Final = Final[!duplicated(Final[,c(1,13)]),] #Cluster, Show

Final$Modules = 0

for (i in 1:length(Final$Cluster)) {
  m = which(Final$Show==Final$Show[[i]])
  Final$Modules[[i]] = paste0(sort(as.numeric(as.character(Final$Cluster[m]))), collapse = ",")
}

write.table(Final, here("Res/WGCNA/Reactome_level3.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
