#### This file is a script that is responsible for legend creation for Figure 3,4

#### Prerequisites: 1. Network_SES_Genes_Creation.R

## Setup env
library(here)
dev.off()
source(here("R/startup.R"))

## Read nodes and edge data from pathway files
node = read.csv(here("Res/Cytoscape/SES/Immune/Pathways/Selected_nodes.csv"))
node = node[which(node$GeneType==""),]

node$Class = "Upregulated DE-genes"
node$Class[19:36] = "Downregulated DE-genes"
node$Class[37:54] = "Upstream targets of\nupregulated DE-genes"
node$Class[55:72] = "Upstream targets of\ndownregulated DE-genes"
node$Class = factor(node$Class, levels = c("Upregulated DE-genes","Downregulated DE-genes","Upstream targets of\nupregulated DE-genes","Upstream targets of\ndownregulated DE-genes"))

node$FDR = -log10(node$Term.PValue.Corrected.with.Benjamini.Hochberg)
node$FDR[which(node$FDR>-log10(0.0001))] = -log10(0.0001)


tiff(str_c(here("Res/Cytoscape/SES/Immune/Pathways/"), "Legend.tiff"), units="px", width=(6*1050), height=(6*850), res=600)
p = ggplot(node, aes(x = Class, y = UNIQUE_ID, color = Class, fill = Class, size = FDR)) +
  geom_point(alpha  = 0.9,shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  
  scale_color_manual(values = c("#D73027", "#1A9850", "#FFD320", "#0072B2"), name = "Gene Class") +
  
  scale_fill_manual(values = c("#D73027", "#1A9850", "#FFD320", "#0072B2"), name = "Gene Class") +
  
  scale_size_continuous(range = c(6,20),
                        name = "Adjusted *p*-value<br/>", 
                        limits = c(-log10(1), -log10(0.0001)),
                        breaks = c(-log10(1),0.7, 2.25,4),
                        labels = c(expression(italic("p")~'>'~0.05), expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001))) +
  
  
  theme(strip.text = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth=1,color='#476b6b')) +
  
  xlab("Module-Trait Relationship<br/>-Log~10~ (Adj. *p*)") + ylab("WGCNA Clusters") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title.y = element_text(size = 14, face = "bold", family = "Calibri")) +
  theme(axis.title.x = element_markdown(size = 14, face = "bold", family = "Calibri")) +
  
  theme(axis.text.y = element_text(size=13, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=13,face = "bold" ,family = "Calibri")) +
  
  theme(legend.text=element_text(size=17, family = "Calibri")) +
  theme(legend.title = element_markdown(size = 18, face = "bold", family = "Calibri")) +
  scale_y_discrete(limits = rev)  +
  guides(color = guide_legend(override.aes = list(size = 11),  keywidth = 1.5, keyheight = 1.5, default.unit = "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 11), keywidth = 1.5, keyheight = 1.5, default.unit = "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 7))) 

p
dev.off()

tiff(str_c(here("Res/Cytoscape/SES/Immune/Genes//"), "Legend.tiff"), units="px", width=(6*775), height=(6*600), res=600)
p = ggplot(node, aes(x = Class, y = UNIQUE_ID, color = Class, fill = Class, size = FDR)) +
  geom_point(alpha  = 0.9,shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  
  scale_color_manual(values = c("#D73027", "#1A9850", "#FFD320", "#0072B2"), name = "Gene Class") +
  
  scale_fill_manual(values = c("#D73027", "#1A9850", "#FFD320", "#0072B2"), name = "Gene Class") +
  
  scale_size_continuous(range = c(6,20),
                        name = "Adjusted *p*-value<br/>", 
                        limits = c(-log10(1), -log10(0.0001)),
                        breaks = c(-log10(1),0.7, 2.25,4),
                        labels = c(expression(italic("p")~'>'~0.05), expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001))) +
  
  
  theme(strip.text = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth=1,color='#476b6b')) +
  
  xlab("Module-Trait Relationship<br/>-Log~10~ (Adj. *p*)") + ylab("WGCNA Clusters") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title.y = element_text(size = 14, face = "bold", family = "Calibri")) +
  theme(axis.title.x = element_markdown(size = 14, face = "bold", family = "Calibri")) +
  
  theme(axis.text.y = element_text(size=13, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=13,face = "bold" ,family = "Calibri")) +
  
  theme(legend.text=element_text(size=17, family = "Calibri")) +
  theme(legend.title = element_markdown(size = 18, face = "bold", family = "Calibri", hjust = 0.5)) +
  scale_y_discrete(limits = rev)  +
  guides(color = guide_legend(nrow = 2,byrow = T, title.position = "top",override.aes = list(size = 11),  keywidth = 1.5, keyheight = 1.5, default.unit = "cm")) +
  guides(fill = guide_legend(nrow = 2,byrow = T, title.position = "top",override.aes = list(size = 11), keywidth = 1.5, keyheight = 1.5, default.unit = "cm")) +
  guides(size = FALSE) + 
  theme(legend.position = "bottom") 

p
dev.off()

node = node[37:72,]
node = droplevels(node)
tiff(str_c(here("Res/Cytoscape/SES/Immune/Pathways/"), "Legend_Rev.tiff"), units="px", width=(6*1050), height=(6*850), res=600)
p = ggplot(node, aes(x = Class, y = UNIQUE_ID, color = Class, fill = Class, size = FDR)) +
  geom_point(alpha  = 0.9,shape = 21, stroke = 1) +
  theme_bw(base_size = 12) +
  
  scale_color_manual(values = c( "#FFD320", "#0072B2"), name = "Gene Class") +
  
  scale_fill_manual(values = c("#FFD320", "#0072B2"), name = "Gene Class") +
  
  scale_size_continuous(range = c(6,20),
                        name = "Adjusted *p*-value<br/>", 
                        limits = c(-log10(1), -log10(0.0001)),
                        breaks = c(-log10(1),0.7, 2.25,4),
                        labels = c(expression(italic("p")~'>'~0.05), expression(italic("p")~'<'~0.05), expression(italic("p")~'<'~0.01),expression(italic("p")~'<'~0.001))) +
  
  
  theme(strip.text = element_text(angle = 0, family = "Calibri", face = "bold", size = 12))+
  theme(panel.spacing =unit(0.1, "lines"),
        panel.border = element_rect(color = "#476b6b", fill = NA, linewidth = 1), 
        strip.background = element_rect(color = "#476b6b", linewidth = 1, fill = "#d6dce4"))+
  theme(axis.ticks = element_line(linewidth=1,color='#476b6b')) +
  
  xlab("Module-Trait Relationship<br/>-Log~10~ (Adj. *p*)") + ylab("WGCNA Clusters") + ggtitle("") + theme(panel.background = element_rect(colour = "grey70")) +
  theme(axis.title.y = element_text(size = 14, face = "bold", family = "Calibri")) +
  theme(axis.title.x = element_markdown(size = 14, face = "bold", family = "Calibri")) +
  
  theme(axis.text.y = element_text(size=13, face = "bold", family = "Calibri")) +
  theme(axis.text.x = element_text(size=13,face = "bold" ,family = "Calibri")) +
  
  theme(legend.text=element_text(size=17, family = "Calibri")) +
  theme(legend.title = element_markdown(size = 18, face = "bold", family = "Calibri")) +
  scale_y_discrete(limits = rev)  +
  guides(color = guide_legend(override.aes = list(size = 11),  keywidth = 1.5, keyheight = 1.5, default.unit = "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 11), keywidth = 1.5, keyheight = 1.5, default.unit = "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 7))) 

p
dev.off()