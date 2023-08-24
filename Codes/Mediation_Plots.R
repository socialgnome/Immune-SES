#### This file is a script that is responsible for plotting mediation plots for SES

#### This script also creates Figure 5, Supplementary Figure S6, S7

#### Prerequisites: 1. Mediation.R

## Setup env
library(here)
dev.off()
source(here("R/startup.R"))


## Read Mediation results
res = readRDS(str_c(here("Res/"), "/Mediation/SES/Filtered_genes.rds"))
meds = readRDS(str_c(here("Res/"), "/Mediation/SES/Mediators.rds"))
genes = readRDS(str_c(here("Res/"), "/Mediation/SES/Genes.rds"))


## Extracting genes
grn_up = list()
#name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes", "All")
name = c("DE Genes", "Upstream")
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Up/Genes/"),name[[i]], ".txt"))
  df = as.character(df$V1)
  grn_up[[i]] = df
}

grn_down = list()
#name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes", "All")
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Down/Genes/"),name[[i]], ".txt"))
  df = as.character(df$V1)
  grn_down[[i]] = df
}


## Calculate mediational results
mediators = unique(meds)
final = c()
for (i in 1:length(grn_up)) {
  sel = which(genes %in% grn_up[[i]])
  for (j in 1:length(mediators)) {
    sel_m = which(meds %in% mediators[[j]])
    sel_f = intersect(sel, sel_m)
    av = c()
    p = c()
    if (length(sel_f)>0) {
      for (k in 1:length(sel_f)) {
        av[k] = res[[sel_f[[k]]]][1]
        p[k] = res[[sel_f[[k]]]][2]
        p[p==0] = 1e-16
      }
      dat = data.frame(Set = name[i], Mediator = mediators[[j]],  Av.Prop = median(av), P = combine.test(p, method = "fisher"))
      final = rbind(final,dat)
    }
  }
}
final$Reg = "Upregulated"
sav = final

final = c()
for (i in 1:length(grn_down)) {
  sel = which(genes %in% grn_down[[i]])
  for (j in 1:length(mediators)) {
    sel_m = which(meds %in% mediators[[j]])
    sel_f = intersect(sel, sel_m)
    av = c()
    p = c()
    if (length(sel_f)>0) {
      for (k in 1:length(sel_f)) {
        av[k] = res[[sel_f[[k]]]][1]
        p[k] = res[[sel_f[[k]]]][2]
        p[p==0] = 1e-16
      }
      dat = data.frame(Set = name[i], Mediator = mediators[[j]],  Av.Prop = median(av), P = combine.test(p, method = "fisher"))
      final = rbind(final,dat)
    }
  }
}
final$Reg = "Downregulated"
final = rbind(final,sav)
final = final[final$P<0.05,]
#name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes", "All")


## Prepare data for plotting
final = final %>%
  mutate(
    # Set  = case_when(Set == "DE Genes" ~ "Set A",
    #                  Set == "Directly affected TF's" ~ "Set B", 
    #                  Set == "First Neighbors of affected TF's" ~ "Set C",
    #                  Set == "TF's from DE genes" ~ "Set D",
    #                  Set == "All" ~ "All"),
    Set = case_when(Set == "DE Genes" ~ "DE Genes",
                    Set == "Upstream" ~ "Upstream\ntargets"),
    Mediator = case_when(Mediator == "stress_perceived" ~ "Stress",
                         Mediator == "w5bmi" ~ "BMI",
                         Mediator == "bills" ~ "Bills",
                         Mediator == "currentsmoke" ~ "Smoking",
                         Mediator == "insurance_lack" ~ "Lack of Insurance",
                         Mediator == "alcohol_use" ~ "Alcohol")
  )

#final$Set = factor(final$Set, levels = c("Set A", "Set B", "Set C", "Set D", "All"))
final$Mediator = factor(final$Mediator, levels = c("BMI", "Stress", "Smoking", "Alcohol", "Bills","Lack of Insurance"))
final$Reg = factor(final$Reg,levels = c("Upregulated", "Downregulated"))

final$Neg = -log10(final$P)
final$Neg[which(final$Neg>-log10(0.001))] = -log10(0.001)

final$Av.Prop  = round(final$Av.Prop*100, digits = 1)
final$Av.Prop = format(final$Av.Prop,nsmall = 1)

tiff(str_c(here("Res/Figures/"), "Network_SES_Mediation.tiff"), units="px", width=(6*750), height=(6*600), res=600)
p = ggplot(final[final$Reg=="Downregulated",], aes(x = Mediator,y  = Set, color = Reg, fill = Reg,  size = Neg)) +
  geom_point(alpha  = 0.90, shape="\u25D6", stroke = 1,position = position_nudge(x = -0.13)) +
  geom_point(data = final[final$Reg=="Upregulated",], shape="\u25D7",inherit.aes = T, position = position_nudge(x = +0.13), alpha = 0.90, stroke = 1,show.legend = F) +
  
  geom_text(data = final[final$Reg=="Downregulated",], aes(x = Mediator,y  = Set, label = Av.Prop), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3.5, inherit.aes = FALSE,position = position_nudge(x = -0.23)) +
  geom_text(data = final[final$Reg=="Upregulated",], aes(x = Mediator,y  = Set, label = Av.Prop), color = "black", family = "Calibri", fontface="bold",alpha=1, size=3.5, inherit.aes = FALSE,position = position_nudge(x = +0.18)) +
  
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("#1A9850","#D73027"), name = "Direction\nof change",
                     #guide = guide_legend(override.aes = list(shape = c("\u25D0","\u25D1"), size =7))) +
                     guide = guide_legend(override.aes = list(shape = c("\u25D6","\u25D7"), size =7))) +
  scale_fill_manual(values = c("#1A9850","#D73027"), name = "Direction\nof change",
                    #guide = guide_legend(override.aes = list(shape = c("\u25D0","\u25D1"), size =7))) +
                    guide = guide_legend(override.aes = list(shape = c("\u25D6","\u25D7"), size =7))) +
  scale_size_continuous(range = c(6,12),
                        name = "Aggregated\nadjusted\np-value", 
                        limits = c(-log10(0.05), -log10(0.001)),
                        breaks = c(-log10(0.05),-log10(0.01),-log10(0.001)), 
                        labels = c(expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001)),
                        guide = guide_legend(override.aes = list(shape = c(21), fill = "NA", colour = "grey30"))) + 
  xlab("Mediators") + ylab("") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70", linewidth = 1.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_markdown(size=12, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=12,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text = element_text(size=12, family = "Calibri")) +
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(legend.position=c(1.21,0.55)) +
  theme(plot.margin = unit(c(5, 120 , 5, 5), "pt")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(limits=rev)
p
dev.off()


### Plot results for gene layers
library(here)
dev.off()
source(here("R/startup.R"))


## Read Mediation results
res = readRDS(str_c(here("Res/"), "/Mediation/SES/Filtered_genes.rds"))
meds = readRDS(str_c(here("Res/"), "/Mediation/SES/Mediators.rds"))
genes = readRDS(str_c(here("Res/"), "/Mediation/SES/Genes.rds"))


## Read selected nodes
sel = read.table(here("Res/Network/SES/All/Cytoscape/Selected_Node.txt"), sep = "\t", header = T)

glist = split(sel$gene, sel$Class)
name = c("Layer 1","Layer 2", "Layer 3", "Layer 4")
glist = setNames(glist, name)

## Calculate mediational results
mediators = unique(meds)
final = c()
for (i in 1:length(glist)) {
  sel = which(genes %in% glist[[i]])
  for (j in 1:length(mediators)) {
    sel_m = which(meds %in% mediators[[j]])
    sel_f = intersect(sel, sel_m)
    av = c()
    p = c()
    if (length(sel_f)>0) {
      for (k in 1:length(sel_f)) {
        av[k] = res[[sel_f[[k]]]][1]
        p[k] = res[[sel_f[[k]]]][2]
        p[p==0] = 1e-16
      }
      dat = data.frame(Set = name[i], Mediator = mediators[[j]],  Av.Prop = median(av), P = combine.test(p, method = "fisher"))
      final = rbind(final,dat)
    }
  }
}

final = final[final$P<0.05,]
#name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes", "All")


## Prepare data for plotting
final = final %>%
  mutate(
    Mediator = case_when(Mediator == "stress_perceived" ~ "Stress",
                         Mediator == "w5bmi" ~ "BMI",
                         Mediator == "bills" ~ "Bills",
                         Mediator == "currentsmoke" ~ "Smoking",
                         Mediator == "insurance_lack" ~ "Lack of Insurance",
                         Mediator == "alcohol_use" ~ "Alcohol")
  )

final$Mediator = factor(final$Mediator, levels = c("BMI", "Stress", "Smoking", "Alcohol", "Bills","Lack of Insurance"))


final$Neg = -log10(final$P)
final$Neg[which(final$Neg>-log10(0.001))] = -log10(0.001)

final$Av.Prop  = round(final$Av.Prop*100, digits = 1)
#$Av.Prop[final$Av.Prop<0] = 0
breaks_col = c(round_any(min(final$Av.Prop),10), round_any(max(final$Av.Prop),10,f = ceiling))

tiff(str_c(here("Res/Figures/"), "Network_SES_Mediation_Onion.tiff"), units="px", width=(6*800), height=(6*600), res=600)
p = ggplot(final, aes(x = Mediator,y  = Set, color = Av.Prop, fill = Av.Prop,  size = Neg)) +
  geom_point(shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  scale_color_viridis(alpha = 0.9,
                      breaks = breaks_col, limits = breaks_col,
                      labels = c("-10", "30"),
                      discrete = F, name = "Median Percent\nProportion Mediated\n")+
  scale_fill_viridis(alpha = 0.9,
                     breaks = breaks_col,limits = breaks_col,
                     labels = c("-10", "30"),
                     discrete = F, name = "Median Percent\nProportion Mediated\n")+
  scale_size_continuous(range = c(6,12),
                        name = "Aggregated\nadjusted\np-value", 
                        limits = c(-log10(0.05), -log10(0.001)),
                        breaks = c(-log10(0.05),-log10(0.01),-log10(0.001)), 
                        labels = c(expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001)),
                        guide = guide_legend(override.aes = list(shape = c(21), fill = "NA", colour = "grey30"))) + 
  xlab("Mediators") + ylab("Gene Class") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70", linewidth = 1.5)) +
  theme(axis.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(axis.text.y = element_markdown(size=13, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=13,face = "bold" ,family = "Calibri", angle = 30,hjust = 1)) +
  theme(legend.text = element_text(size=12, family = "Calibri")) +
  theme(legend.text.align = 0) +
  theme(legend.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(legend.position="right") +
  scale_x_discrete(drop = F) 
p
dev.off()

###
library(reshape)
library(ggradar)
test = final
test = droplevels(test)

tiff(str_c(here("Res/Figures/"), "Network_SES_Mediation_Onion_Spyder_Alt.tiff"), units="px", width=(6*750), height=(6*720), res=600)

p = ggplot(test, 
           aes(x = Mediator, y = Av.Prop, fill = Set)) +
  geom_bar(position="dodge", stat="identity", alpha = 0.8) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c("#4379a7", "#F28E2B", "#e15759", "#59A14F"), name = "Gene Class") +
  #scale_color_manual(values = c("#4379a7", "#F28E2B", "#e15759", "#59A14F"), name = "Gene Class") +
  xlab("Mediator") + ylab("Median Percent Mediated Ratio") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  ylim(-10,30) +
  theme(axis.title = element_text(size = 15, face = "bold", family = "Calibri")) +

  theme(axis.text.y = element_text(size=14, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=14,face = "bold" ,family = "Calibri")) +
  
  theme(legend.text=element_text(size=13, family = "Calibri")) +
  theme(legend.title = element_text(size = 14, face = "bold", family = "Calibri")) 
  

p

dev.off()

finalT = cast(final[,c(1,2,3)], Set~Mediator)
finalT[is.na(finalT)] = 0

breaks_col = c(round_any(min(final$Av.Prop),10), round_any(max(final$Av.Prop),10,f = ceiling))

tiff(str_c(here("Res/Figures/"), "Network_SES_Mediation_Onion_Spyder.tiff"), units="px", width=(6*900), height=(6*800), res=600)

p = ggradar(finalT, 
            base.size = 13,
            font.radar = "Calibri",
            values.radar = c("-10%", "10%", "30%"),
            grid.min = breaks_col[[1]],
            grid.mid = mean(breaks_col),
            grid.max = breaks_col[[2]],
            gridline.min.linetype = "solid",
            gridline.mid.linetype = "longdash",
            gridline.max.linetype = "solid",
            gridline.min.colour = "#476b6b",
            gridline.mid.colour = "grey",
            gridline.max.colour = "#476b6b",
            grid.label.size = 6,
            gridline.label.offset = -4,
            axis.label.offset = 1.1,
            axis.label.size = 7,
            axis.line.colour = "grey",
            group.line.width =  1,
            group.point.size = 4,
            background.circle.colour = "#d6dce4",
            background.circle.transparency = 0.15,
            legend.title = "Gene Class"
            )
p = p + 
  scale_color_tableau(palette = "Classic 10 Medium") +
  theme(axis.title = element_text(size = 14, face = "bold", family = "Calibri")) +
  theme(legend.text = element_text(size=13, family = "Calibri")) +
  theme(legend.title = element_text(size = 13, face = "bold", family = "Calibri")) +
  theme(legend.position="bottom")
p

dev.off()