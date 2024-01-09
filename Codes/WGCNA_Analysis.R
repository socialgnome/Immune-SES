#### This file is a script that is responsible for the analysing WGCNA results

#### This script also creates Supplementary Figure S2

#### Prerequisites: 1. WGCNA.R

#### Run Order = 5


## Setup env

library(here)
dev.off()
source(here("R/startup.R"))

## Load clustering results
load(str_c(here("Res/WGCNA/"), "WGCNA_Res.RData"))
rm(list= ls()[!(ls() %in% c('net'))])

## Module characteristics and numbers
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
moduleInfo = data.frame(Gene = names(moduleLabels), Module = as.numeric(moduleLabels), Colors = moduleColors)
ColInfo = moduleInfo[!duplicated(moduleInfo[-c(1)]),-c(1)]
MEs = net$MEs;
# MEs0 = moduleEigengenes(t(exprs), moduleColors)$eigengenes
# MEs1 = orderMEs(MEs0)
mod_nos = as.data.frame(table(moduleLabels))
colnames(mod_nos) = c("Cluster", "Total")


## Load Pheno Data and test for Eigen gene - SES relationship
filt = readRDS(str_c(here("data/"), "Filtered_ExpressionSet_Clustering_AllSubData.rds"))
pheno = pData(filt)

# Define treatment, control, mediation and disease vars
treatment = c("ses_sss_composite","sss_5","SEI_ff5",
              "edu_max","income_hh_ff5")

controls = c(
  "sex_interv", "re","age_w5"
  ,"pregnant_biow5", "FastHrs","Plate"
  ,"H5INFECT", "H5SUBCLN", "H5CRP8"
)

imp_vars = c("batch")

mediators = c(
  "stress_perceived",
  "w5bmi",
  "bills",
  "currentsmoke",
  "insurance_lack",
  "alcohol_use"
)

disease = c("diabetes", "heartatk", "H5ID6F", "H5ID6FM", "H5ID6C", "H5ID6CM",
            "H5ID6Q","H5ID6QM", "H5ID6A", "H5ID6AM")

subData = pheno %>% dplyr::select(AID,all_of(treatment), all_of(controls), all_of(imp_vars), all_of(mediators), all_of(disease))
selsubData = pheno %>% dplyr::select(AID,all_of(treatment), all_of(controls), all_of(imp_vars))

expr = t(as.matrix(MEs))
rownames(expr) = as.numeric(substring(rownames(expr),3))
i = 1

rhs = str_c(c(treatment[[i]],controls), collapse = " + ")
model_formula = str_c(" ~ ",rhs) %>% as.formula()
sub = selsubData %>% dplyr::select(treatment[[i]], all_of(controls))
sub = droplevels(sub)

if (class(sub[,c(1)])=="factor") {
  sub[,c(1)] = as.numeric(sub[,c(1)])
}

design  = model.matrix(model_formula, data = sub)

## DE
fit = lmFit(expr, design) 
fit = eBayes(fit, trend = T)

res = topTable(fit, coef = treatment[[i]], n= Inf) %>%
  rownames_to_column(var = "gene")
res$treatment = treatment[[i]]

mod_nos$Adj.P = res$adj.P.Val[match(mod_nos$Cluster,res$gene)]
mod_nos$LogFC = res$logFC[match(mod_nos$Cluster,res$gene)]

res = res[res$adj.P.Val<0.05,]
res = res[order(-res$logFC),]

sig = moduleInfo[moduleInfo$Module %in% res$gene,]
sig = sig[order(sig$Module),]

#write.table(sig, str_c(here("Res/WGCNA/"), "WGCNA_SES_Modules.txt"), row.names = F, col.names = T, sep = "\t")

## Read DE results
DE = list()
df = readRDS(str_c(here("Res/DE/"), "ses_sss_composite" , ".rds"))
DE[[1]] = df$gene[df$adj.P.Val<0.05]
DE[[2]] = df$gene[df$adj.P.Val<0.05 & df$logFC>0]
DE[[3]] = df$gene[df$adj.P.Val<0.05 & df$logFC<0]

name = c("All","Upregulated","Downregulated")
DE = setNames(DE, name)

## Merge results
moduleInfo$LFC = df$logFC[match(moduleInfo$Gene, df$gene)]
moduleInfo$P = df$adj.P.Val[match(moduleInfo$Gene, df$gene)]

sel_up = moduleInfo[moduleInfo$LFC>0 & moduleInfo$P<0.05,]
mod_nos_up = as.data.frame(table(sel_up$Module))
mod_nos$Up = mod_nos_up$Freq[match(mod_nos$Cluster, mod_nos_up$Var1)]

sel_down = moduleInfo[moduleInfo$LFC<0 & moduleInfo$P<0.05,]
mod_nos_down = as.data.frame(table(sel_down$Module))
mod_nos$Down = mod_nos_down$Freq[match(mod_nos$Cluster, mod_nos_down$Var1)]

mod_nos[is.na(mod_nos)] = 0

mod_nos$P_Up = 0
mod_nos$Odds_Up = 0
mod_nos$P_Down = 0
mod_nos$Odds_Down = 0


## Fisher test for modules
for (i in 1:length(mod_nos$Cluster)) {
  
  ## Fisher test for up
  n1 = mod_nos$Up[i]  
  n2 = mod_nos$Total[i]
  n3 = sum(mod_nos$Up)
  n4 = sum(mod_nos$Total)
  f = fisher.test(matrix(c(n1,n3-n1,n2-n1,n4-n3-n2+n1), nrow =2, ncol=2), alternative ="two.sided")
  mod_nos$P_Up[[i]] = as.numeric(f$p.value)
  mod_nos$Odds_Up[[i]] = as.numeric(f$estimate)
  
  ## Fisher test for down
  n1 = mod_nos$Down[i]  
  n2 = mod_nos$Total[i]
  n3 = sum(mod_nos$Down)
  n4 = sum(mod_nos$Total)
  f = fisher.test(matrix(c(n1,n3-n1,n2-n1,n4-n3-n2+n1), nrow =2, ncol=2), alternative ="two.sided")
  mod_nos$P_Down[[i]] = as.numeric(f$p.value)
  mod_nos$Odds_Down[[i]] = as.numeric(f$estimate)
  
}

## Results for significantly enriched clusters
# Up
mod_nos$Cluster[mod_nos$P_Up<0.05 & mod_nos$Odds_Up>2]
mod_nos$Cluster[mod_nos$P_Down<0.05 & mod_nos$Odds_Down>2]
rm(list= ls()[!(ls() %in% c('moduleInfo','mod_nos'))])
#save.image(str_c(here("Res/WGCNA/"), "WGCNA_Res_Enrichment.RData"))


## Plot
data = data.frame(Cluster = rep(mod_nos$Cluster,2), Adj.P = rep(mod_nos$Adj.P,2), 
                  Odds = c(mod_nos$Odds_Up, mod_nos$Odds_Down), P = c(mod_nos$P_Up, mod_nos$P_Down),
                  Reg = c(rep("Upregulated", length(mod_nos$Cluster)), rep("Downregulated", length(mod_nos$Cluster))))
data$Neg = -log10(data$Adj.P)
data$Neg[which(data$Neg>-log10(0.001))] = -log10(0.001)

data$neg2 = -log10(data$P)
data$neg2[which(data$neg2>-log10(0.001))] = -log10(0.001)

data$Odds[data$Odds>4] = 4
data$OddsCat = 0
data$OddsCat[data$Odds>1] = 1
data$OddsCat2 = "Underrepresented"
data$OddsCat2[data$Odds>1] = "Overrepresented"

tiff(str_c(here("Res/WGCNA/"), "WGCNA_Analysis.tiff"), units="px", width=(6*1100), height=(6*800), res=600)
p = ggplot(data, aes(x = Neg, y = Cluster, color = OddsCat2, fill = Odds, size = neg2)) +
  geom_point(alpha  = 0.75,shape = 21, stroke = 2) +
  theme_bw(base_size = 12) +
  
  #scale_color_distiller(palette = "RdBu", limits = c(0,1), breaks = c(0,1), direction = -1)+
  scale_color_manual(values = c(head(brewer.pal(11, "RdBu"), n=1),tail(brewer.pal(11, "RdBu"), n=1)), name = "Enrichment") +
  
  
  # scale_color_viridis(alpha = 0.9,
  #                     breaks = c(min(data$Odds),(min(data$Odds)+ max(data$Odds))/2, max(data$Odds)),
  #                     labels = c("0", "2", ">4"),
  #                     discrete = F, name = "Odds Ratio")+
  scale_fill_viridis(alpha = 0.9,
                     breaks = c(min(data$Odds),(min(data$Odds)+ max(data$Odds))/2, max(data$Odds)),
                     labels = c("0", "2", ">4"),
                     discrete = F, name = "Odds Ratio")+
  
  scale_size_continuous(range = c(4,10),
                        name = "Fisher's<br/>*p*-value", 
                        limits = c(-log10(1), -log10(0.001)),
                        breaks = c(-log10(1),-log10(0.05), -log10(0.01),-log10(0.001)), 
                        labels = c(expression(italic("p")~'>'~0.05), expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001))) +
  
  facet_grid(~Reg,scales="free",space="free")+
  
  geom_vline(xintercept = -log10(0.05),  linetype="dashed", color = "red", size= 1.5) +
  
  theme(strip.text = element_text(angle = 0, family = "Calibri", face = "bold", size = 15))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1.5), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1.5, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth=1,color='#476b6b')) +
  
  xlab("Cluster-SES Relationship<br/>-Log~10~ (Adj. *p*)") + ylab("WGCNA Cluster Label") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title.y = element_text(size = 15, face = "bold", family = "Calibri")) +
  theme(axis.title.x = element_markdown(size = 15, face = "bold", family = "Calibri")) +
  
  theme(axis.text.y = element_text(size=14, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=14,face = "bold" ,family = "Calibri")) +

  theme(legend.text=element_text(size=14, family = "Calibri")) +
  theme(legend.title = element_markdown(size = 14, face = "bold", family = "Calibri")) +
  #guides(color = "none")+
  guides(size = guide_legend(keyheight = 1.5)) +
  guides(color = guide_legend(override.aes = list(size = 6), keyheight = 1.5)) +
  scale_y_discrete(limits = rev)  

p
dev.off()

## Write genes
#Upregulated
sig = moduleInfo[moduleInfo$Module %in% c(7,13,17),]
sig = sig[sig$P>0.05,]
gene = as.character(sig$Gene)
write.table(gene,here("Res/Network/SES/Up/Genes/Wgcna.txt"), row.names = F, col.names = F, sep = "\t", quote = F)

#Downregulated
sig = moduleInfo[moduleInfo$Module %in% c(11),]
sig = sig[sig$P>0.05,]
gene = as.character(sig$Gene)
write.table(gene,here("Res/Network/SES/Down/Genes/Wgcna.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
