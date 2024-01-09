#### This file is a script that is responsible for performing the randomization analysis of the upstream regulators

#### This script also creates Supplementary Figure S10, S11

#### Prerequisites: 1. Randomization_Analysis.R

#### Run Order = 17


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
name = c("DE Genes", "Upstream")
grn_act = list()
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Up/Genes"),"/",name[[i]], ".txt"), header = F, sep = "\t")
  grn_act[[i]] = as.character(df$V1)
}

grn_act = setNames(grn_act, name)

## Extract Randomization results
library(metap)
iter = 1000
res = readRDS(str_c(here("Res/"), "/Randomization/Results_Up.rds"))
sumP = list()
P_up = list()
for (i in 1:length(grn_act)) {
  if (i>1) {
    test = lapply(res, `[[`, 1) 
    t1 = lapply(test, `[[`, 1)
    t2 = lapply(test, `[[`, 2)
    t3 = lapply(test, `[[`, 3)
    t4 = lapply(test, `[[`, 4)
    for (k in 1:iter) {
      if (!is.null(t3[[k]])) {
        t3[[k]] = setdiff(setdiff(t3[[k]],t1[[k]]),t4[[k]])
      }
      if (!is.null(t4[[k]])) {
        t4[[k]] = setdiff(t4[[k]],t1[[k]])
      }
    }
    test = table(c(unlist(t2),unlist(t3), unlist(t4)))
  } else {
    test = lapply(res, `[[`, 1) 
    test = lapply(test, `[[`, i)
    test = table(unlist(test))
  }
  P = c()
  for (j in 1:length(grn_act[[i]])) {
    t = grn_act[[i]][[j]]
    n = test[which(names(test)==t)]
    p = as.numeric((n+1)/(iter+1))
    P = c(P,p)
  }
  sumP[[i]] = data.frame(Condition = name[[i]], Regulation = "Upregulated", Value =  sumlog(P)$p)
  P_up[[i]] = data.frame(Condition = name[[i]], Regulation = "Upregulated", Value = P)
}

sumP_Up = sumP


## Extract true Downregulated GRN genes
name = c("DE Genes", "Upstream")
grn_act = list()
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Down/Genes"),"/",name[[i]], ".txt"), header = F, sep = "\t")
  grn_act[[i]] = as.character(df$V1)
}

grn_act = setNames(grn_act, name)

## Extract Randomization results
library(metap)
iter = 1000
res = readRDS(str_c(here("Res/"), "/Randomization/Results_Down.rds"))
sumP = list()
P_down = list()
for (i in 1:length(grn_act)) {
  if (i>1) {
    test = lapply(res, `[[`, 1)
    t1 = lapply(test, `[[`, 1)
    t2 = lapply(test, `[[`, 2)
    t3 = lapply(test, `[[`, 3)
    t4 = lapply(test, `[[`, 4)
    for (k in 1:iter) {
      if (!is.null(t3[[k]])) {
        t3[[k]] = setdiff(setdiff(t3[[k]],t1[[k]]),t4[[k]])
      }
      if (!is.null(t4[[k]])) {
        t4[[k]] = setdiff(t4[[k]],t1[[k]])
      }
    }
    test = table(c(unlist(t2),unlist(t3), unlist(t4)))
  } else {
    test = lapply(res, `[[`, 1) 
    test = lapply(test, `[[`, i)
    test = table(unlist(test))
  }
  P = c()
  for (j in 1:length(grn_act[[i]])) {
    t = grn_act[[i]][[j]]
    n = test[which(names(test)==t)]
    p = as.numeric((n+1)/(iter+1))
    P = c(P,p)
  }
  sumP[[i]] = data.frame(Condition = name[[i]], Regulation = "Downregulated", Value =  sumlog(P)$p)
  P_down[[i]] = data.frame(Condition = name[[i]], Regulation = "Downregulated", Value = P)
}

sumP_Down = sumP


## Plots
res = rbind(do.call("rbind",P_up), do.call("rbind",P_down))

res$Class = "Upregulated DE-genes"
res$Class[res$Condition=="DE Genes" & res$Regulation=="Downregulated"] = "Downregulated DE-genes"
res$Class[res$Condition=="Upstream" & res$Regulation=="Upregulated"] = "Upstream targets of\nupregulated DE-genes"
res$Class[res$Condition=="Upstream" & res$Regulation=="Downregulated"] = "Upstream targets of\ndownregulated DE-genes"

res$Class = factor(res$Class, levels = c("Upregulated DE-genes","Downregulated DE-genes","Upstream targets of\nupregulated DE-genes","Upstream targets of\ndownregulated DE-genes"))


tiff(str_c(here("Res/Randomization/"), "Density.tiff"), units="px", width=(6*800), height=(6*600), res=600)
p = ggplot(res, aes(x = Value, color = Class, fill = Class)) +
  geom_density(alpha  = 0.2,linewidth = 0.8) +
  geom_vline(xintercept = 0.05,  linetype="dashed", color = "red", size= 1.5) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#D73027", "#1A9850", "#FFD320", "#0072B2")) +
  scale_fill_manual(values = c("#D73027", "#1A9850", "#FFD320", "#0072B2")) +
  
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1.5), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1.5, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth=1,color='#476b6b')) +
  xlim(0,0.5) +
  xlab("*p* - value") + ylab("Frequency") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title.y = element_text(size = 15, face = "bold", family = "Calibri")) +
  theme(axis.title.x = element_markdown(size = 15, face = "bold", family = "Calibri")) +
  
  theme(axis.text.y = element_text(size=14, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=14,face = "bold" ,family = "Calibri")) +
  
  theme(legend.text=element_text(size=14, family = "Calibri")) +
  theme(legend.title = element_markdown(size = 14, face = "bold", family = "Calibri")) +
  theme(legend.position = c(.75, .75)) +
  guides(fill = guide_legend(override.aes = list(size = 6), keyheight = 2.5)) +
  guides(color = guide_legend(override.aes = list(size = 6), keyheight = 2.5)) 
#guides(color = "none")+
#guides(fill = "none")
p
dev.off()

res2 = rbind(do.call("rbind",sumP_Up), do.call("rbind",sumP_Down))
res2$Class = "Upregulated DE-genes"
res2$Class[res2$Condition=="DE Genes" & res2$Regulation=="Downregulated"] = "Downregulated DE-genes"
res2$Class[res2$Condition=="Upstream" & res2$Regulation=="Upregulated"] = "Upstream targets of\nupregulated DE-genes"
res2$Class[res2$Condition=="Upstream" & res2$Regulation=="Downregulated"] = "Upstream targets of\ndownregulated DE-genes"
res2$Class = factor(res2$Class, levels = c("Upregulated DE-genes","Downregulated DE-genes","Upstream targets of\nupregulated DE-genes","Upstream targets of\ndownregulated DE-genes"))

res2 = res2[order(res2$Class),]
res2$Total = c(length(which(res$Value[res$Class=="Upregulated DE-genes"]<1)),
               length(which(res$Value[res$Class=="Downregulated DE-genes"]<1)),
               length(which(res$Value[res$Class=="Upstream targets of\nupregulated DE-genes"]<1)),
               length(which(res$Value[res$Class=="Upstream targets of\ndownregulated DE-genes"]<1)))

res2$Sig = c(length(which(res$Value[res$Class=="Upregulated DE-genes"]<0.1)),
               length(which(res$Value[res$Class=="Downregulated DE-genes"]<0.1)),
               length(which(res$Value[res$Class=="Upstream targets of\nupregulated DE-genes"]<0.1)),
               length(which(res$Value[res$Class=="Upstream targets of\ndownregulated DE-genes"]<0.1)))
res2$Prop = res2$Sig/res2$Total

res2$Neg = -log10(res2$Value)
res2$Neg[which(res2$Neg>-log10(0.001))] = -log10(0.001)

tiff(str_c(here("Res/Randomization/"), "Fishers.tiff"), units="px", width=(6*900), height=(6*500), res=600)
p = ggplot(res2, aes(x = Neg, y = Class, color = Prop, fill = Prop, size = Total)) +
  geom_point(alpha  = 0.75,shape = 21) +
  theme_bw(base_size = 12) +
  scale_color_viridis(alpha = 0.9,
                      limits = c(0,1),
                     breaks = c(0,0.5,1),
                     labels = c("0", "0.5", "1"),
                     discrete = F, name = "Proportion of\nSignificant genes")+
  scale_fill_viridis(alpha = 0.9,
                     limits = c(0,1),
                     breaks = c(0,0.5,1),
                     labels = c("0", "0.5", "1"),
                     discrete = F, name = "Proportion of\nSignificant genes")+
  
  scale_size_continuous(range = c(4,10),
                        name = "Total number\nof Genes", 
                        limits = c(50,450),
                        breaks = c(50, 250, 450), 
                        labels = c("50", "250", "450")) +
  
  geom_vline(xintercept = -log10(0.05),  linetype="dashed", color = "red", size= 1.5) +
  
  theme(strip.text = element_text(angle = 0, family = "Calibri", face = "bold", size = 15))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1.5), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1.5, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth=1,color='#476b6b')) +
  xlim(0,3)+
  
  xlab("Fisher's Combined *p*-value<br/>-Log~10~ (*p*)") + ylab("Class") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title.y = element_text(size = 15, face = "bold", family = "Calibri")) +
  theme(axis.title.x = element_markdown(size = 15, face = "bold", family = "Calibri")) +
  
  theme(axis.text.y = element_text(size=14, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=14,face = "bold" ,family = "Calibri")) +
  
  theme(legend.text=element_text(size=14, family = "Calibri")) +
  theme(legend.title = element_text(size = 14, face = "bold", family = "Calibri")) +
  #guides(color = "none")+
  #guides(size = guide_legend(keyheight = 1.5)) +
  #guides(color = guide_legend(override.aes = list(size = 6), keyheight = 1.5)) +
  scale_y_discrete(limits = rev)  

p
dev.off()