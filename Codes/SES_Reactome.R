#### This file is a script that is responsible for performing the DE enrichment for SES

#### This script also creates Figure 2, Supplementary Figure S3, Dataset S1

#### Prerequisites: 1. DE.R

## Setup env
library(here)
dev.off()
source(here("R/startup.R"))

## Biomart (Downloaded 11.12.22)
sapiens_ensembl <- readRDS(str_c(here("data/"), "sapiens_ensembl.rds"))

## Load signatures
treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")

## DE genes
all = list()
up = list()
down = list()
for (i in 1:length(treatment)) {
  df = readRDS(str_c(here("Res/DE/"), treatment[[i]],".rds"))
  all[[i]] = df$gene[df$adj.P.Val<0.05]
  up[[i]] = df$gene[df$adj.P.Val<0.05 & df$logFC>0]
  down[[i]] = df$gene[df$adj.P.Val<0.05 & df$logFC<0]
}
treatment = c("SES Composite","Subjective Social Status","Occupation","Education","Income")

all = setNames(all, treatment)
down = setNames(down, treatment)
up = setNames(up, treatment)


## SES genes extraction
enrich = list()
enrich[[1]] = up[[1]]
enrich[[2]] = down[[1]]
reg = c("Upregulated", "Downregulated")
enrich = setNames(enrich, reg)

## SES enrichment  
enrich_ent = list()
conv = list()
for (i in 1:length(enrich)) {
  df = biomaRt::getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values =enrich[[i]], mart = sapiens_ensembl, useCache = F)
  conv[[i]] = df
  df = as.character(df$entrezgene_id[!is.na(df$entrezgene_id)])
  enrich_ent[[i]] = df
}

enrich_ent = setNames(enrich_ent, reg)

eres = compareCluster(enrich_ent, fun = "enrichPathway", organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH")
eres = eres@compareClusterResult
eres = eres[eres$Count>2,]

eres$genes = 0
for (i in 1:length(eres$Cluster)) {
  eres$genes[[i]] = paste(conv[[which(reg==eres$Cluster[[i]])]]$hgnc_symbol[match(unlist(strsplit(eres$geneID[[i]], "/")),conv[[which(reg==eres$Cluster[[i]])]]$entrezgene_id)],collapse = ",")
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

#tiff(str_c(here("Res/Figures/"), "DE_Pathway_SES_UpsetR.tiff"), units="px", width=(3*850), height=(3*675), res=300)
t1 = upset(fromList(lists),
           mainbar.y.label = "Intersecting enriched pathways",
           sets.x.label = "Total enriched pathways",
           keep.order = T,
           main.bar.color = "#555555",
           matrix.color  = "#326a97",
           set_size.show = T,
           point.size = 2.8,
           shade.alpha = 0.2,
           matrix.dot.alpha = 0.7,
           line.size = 0.6,
           mb.ratio = c(0.6,0.4), text.scale = 1.7, 
           set_size.numbers_size = 6, set_size.scale_max = 100)
t1
#dev.off()


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
Des$Show = 0
Des$ID = test$ID
test = test[complete.cases(Des[,c(2)]),]
Des = Des[complete.cases(Des[,c(2)]),]
for (i in 1:length(Des$Result)) {
  k = as.character(Des[i,-c(1,4,5,6)]) ### Choose the level that you want to contexualize the data
  k = k[complete.cases(k)]
  Des$Show[i] = k[length(k)]
}


## Combine similar results 
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

# Write tabular results
tab = data.frame(Regulation = Final$Cluster, ReactomeID = Final$ID, Pathway = Final$Description,
                 ParentNode = Final$Level1, ChildNode =  Final$Level2, GeneRatio = Final$GeneRatio,
                 BgRatio = Final$BgRatio, Adj.P.Val = Final$p.adjust, CramerV = Final$Cramer,
                 CollapsedCramerV = Final$MedC,CollapsedCount = Final$FinalCount, Genes = Final$genes)
tab = tab[order(tab$Regulation),]
write.xlsx(tab, str_c(here("Res/Tables/"), "Reactome_SES_All.xlsx"), overwrite = T)
tab = tab[,c(1,2,3,4,8,9,12)]
tab$Adj.P.Val  = sprintf("%.3f x 10^(%d)", tab$Adj.P.Val/10^floor(log10(abs(tab$Adj.P.Val))), floor(log10(abs(tab$Adj.P.Val))))
tab$CramerV  = sprintf("%.3f x 10^(%d)", tab$CramerV/10^floor(log10(abs(tab$CramerV))), floor(log10(abs(tab$CramerV))))
write.xlsx(tab, str_c(here("Res/Tables/"), "Reactome_SES_Concise.xlsx"), overwrite = T)

## Plot reactome contextualized results

sav = Final
Final = Final[!duplicated(Final[,c(1,13)]),] #Cluster, Show
Final$FinalSize = 7

FinalT = Final
FinalT = droplevels(FinalT)
FinalT$Show[which(FinalT$Show=="Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell")] = "Immunoregulatory interactions in Lymphoid cell"
FinalT$Show[which(FinalT$Show=="Nucleotide-binding domain, leucine rich repeat containing receptor (NLR) signaling pathways")] = "NLR Signaling pathways"
FinalT$Show[which(FinalT$Show=="Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.")] = "Oxidative phosphorylation"
FinalT$Show[which(FinalT$Show=="The citric acid (TCA) cycle and respiratory electron transport")] = "TCA cycle"
FinalT$Show[which(FinalT$Show=="Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)")] = "Independent Nonsense Mediated Decay"
FinalT$Show[which(FinalT$Show=="Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)")] = "Dependent Nonsense Mediated Decay"
FinalT$Level1[which(FinalT$Level1=="Organelle biogenesis and maintenance")] = "Organelle biogenesis"

FinalT$Neg = -log(FinalT$p.adjust)
try = as.character(unique(FinalT$Show))
FinalT$Show = factor(FinalT$Show, levels = unique(FinalT$Show))
FinalT$Show = str_wrap(FinalT$Show, width = 60)
FinalT$Level1 = str_wrap(FinalT$Level1, width = 25)

FinalT$Cluster  = factor(FinalT$Cluster, levels =  c("Upregulated", "Downregulated"))
breaks_col = c(min(FinalT$MedC),(min(FinalT$MedC)+ max(FinalT$MedC))/2, max(FinalT$MedC))
breaks_size = c(min(FinalT$FinalCount),35,max(FinalT$FinalCount))

tiff(str_c(here("Res/Figures/"), "Reactome_Enrichment_SES.tiff"), units="px", width=(6*750), height=(6*850), res=600)

p <- ggplot(FinalT, aes(x = Cluster, y = Show, size = FinalCount, color = MedC, fill = MedC)) +
  geom_point(shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  scale_color_viridis(alpha = 0.9,
                      breaks = c(min(FinalT$MedC),(min(FinalT$MedC)+ max(FinalT$MedC))/2, max(FinalT$MedC)),
                      labels = c("0.03", "0.19", "0.35"),
                      discrete = F, name = "Cramer's V\n")+
  scale_fill_viridis(alpha = 0.9,
                     breaks = c(min(FinalT$MedC),(min(FinalT$MedC)+ max(FinalT$MedC))/2, max(FinalT$MedC)),
                     labels = c("0.03", "0.19", "0.35"),
                     discrete = F, name = "Cramer's V\n")+
  scale_size_continuous(range = c(3,7), name = "Gene Count", breaks = c(min(FinalT$FinalCount),35,max(FinalT$FinalCount)), labels = c("3","35", "70")) + 
  facet_grid(Level1~.,scales="free",space="free")+
  theme(strip.text.y = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth =1,color='#476b6b')) +
  xlab("") + ylab("Reactome Pathways") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 14, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=13, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=13,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  scale_x_discrete(drop=FALSE)+
  theme(legend.text=element_text(size=12, family = "Calibri")) +
  theme(legend.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(legend.position="bottom") 

p
dev.off()


### Only Immune

Final = sav
Final = Final[Final$Level1=="Immune System",]
Final = Final[!duplicated(Final[,c(1,16)]),] #Cluster, Result
Final$FinalSize = 7

FinalT = Final
FinalT = droplevels(FinalT)
FinalT$Neg = -log(FinalT$p.adjust)
try = as.character(unique(FinalT$Result))
#FinalT$Result = substring(FinalT$Result, 16)

FinalT$Result[which(FinalT$Result=="Nucleotide-binding domain, leucine rich repeat containing receptor (NLR) signaling pathways")] = "NLR Signaling pathways"


FinalT = FinalT[order(FinalT$p.adjust),]
FinalT$Result = factor(FinalT$Result, levels = unique(FinalT$Result))
FinalT$Result = str_wrap(FinalT$Result, width = 40)
FinalT$Level2 = str_wrap(FinalT$Level2, width = 25)

FinalT$Cluster  = factor(FinalT$Cluster, levels =  c("Upregulated", "Downregulated"))
tiff(str_c(here("Res/Figures/"), "Reactome_Immune_Enrichment_SES.tiff"), units="px", width=(6*750), height=(6*810), res=600)

p <- ggplot(FinalT, aes(x = Cluster, y = Result, size = Count, color = Cramer, fill = Cramer)) +
  geom_point(shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  scale_color_viridis(alpha = 0.9,
                      breaks = breaks_col,
                      labels = c("0.03", "0.19", "0.35"),
                      discrete = F, name = "Cramer's V\n")+
  scale_fill_viridis(alpha = 0.9,
                     breaks = breaks_col,
                     labels = c("0.03", "0.19", "0.35"),
                     discrete = F, name = "Cramer's V\n")+
  scale_size_continuous(range = c(2,6), name = "Gene Count", breaks = breaks_size, labels = c("3","35", "70")) +
  facet_grid(Level2~.,scales="free",space="free")+
  theme(strip.text.y = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth=1,color='#476b6b')) +
  xlab("") + ylab("") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  scale_x_discrete(drop=FALSE)+
  theme(legend.text=element_text(size=11, family = "Calibri")) +
  theme(legend.title = element_text(size = 12, face = "bold", family = "Calibri")) +
  theme(legend.position="none") 

p
dev.off()

tiff(str_c(here("Res/Figures/"), "Reactome_Immune_Enrichment_SES_Full.tiff"), units="px", width=(6*800), height=(6*900), res=600)

p <- ggplot(FinalT, aes(x = Cluster, y = Result, size = Count, color = Cramer, fill = Cramer)) +
  geom_point(shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  scale_color_viridis(alpha = 0.9,
                      breaks = c(min(FinalT$Cramer),(min(FinalT$Cramer)+ max(FinalT$Cramer))/2, max(FinalT$Cramer)),
                      labels = c("0.03", "0.06", "0.09"),
                      discrete = F, name = "Cramer's V\n")+
  scale_fill_viridis(alpha = 0.9,
                     breaks = c(min(FinalT$Cramer),(min(FinalT$Cramer)+ max(FinalT$Cramer))/2, max(FinalT$Cramer)),
                     labels = c("0.03", "0.06", "0.09"),
                     discrete = F, name = "Cramer's V\n")+
  scale_size_continuous(range = c(4,8), name = "Gene Count", breaks = c(min(FinalT$Count),19,max(FinalT$Count)), labels = c("3","19", "35")) +
  facet_grid(Level2~.,scales="free",space="free")+
  theme(strip.text.y = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth=1,color='#476b6b')) +
  xlab("") + ylab("Immune Pathways") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title = element_text(size = 14, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_text(size=13, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=13,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  scale_x_discrete(drop=FALSE)+
  theme(legend.text=element_text(size=12, family = "Calibri")) +
  theme(legend.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(legend.position="bottom") 

p
dev.off()