#### This file is a script that is responsible for the analysing the functional 
#### roles of significant (To ses) cluster that is identified by WGCNA

#### This script also creates Figure 1

#### Prerequisites: 1. WGCNA_Analysis.R

#### Run Order = 7


## Setup env
library(here)
#dev.off()
source(here("R/startup.R"))

## Load clustering results
load(str_c(here("Res/WGCNA/"), "WGCNA_Res_Enrichment.RData"))


## Find enriched pathways
glists = split(moduleInfo$Gene, moduleInfo$Module)
glists = glists[mod_nos$Cluster[mod_nos$Adj.P<0.05]]
sig = mod_nos[mod_nos$Adj.P<0.05,]

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


relation = read.table("~/Projects/Helper files/Reactome/ReactomePathwaysRelation.txt", sep = "\t", stringsAsFactors = F)
names = read.table("~/Projects/Helper files/Reactome/ReactomePathways.txt", sep = "\t", stringsAsFactors = F, quote = "")
names = names[names$V3=="Homo sapiens",]
test = eres
test = test[test$Count>5,]
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
Final = Final[order(Final$Cluster),]


##### Chord diagram
library(RColorBrewer)
library(viridis)
library(circlize)
library(grid)
library(ComplexHeatmap)
library(gridBase)
library(extrafont)

## Cluster legend
rdims = data.frame(Label = unique(Final$Cluster), Name = c("(I)", "(II)", "(III)"), stringsAsFactors = F)
rcols = c('#4E79A7', '#F28E2B', '#B07AA1')
rname = c("I", "II", "III")

## Reactome legend
dims = data.frame(Label = LETTERS[seq( from = 1, to = length(unique(Final$Level1)))], Name = sort(unique(Final$Level1)), stringsAsFactors = F)
dims$Label[dims$Label=="I"] = " I "
ccols = viridis(length(dims$Label), alpha = 0.8)
mycolor <- c(rcols, ccols)


## Link colors
lcols = rep(ccols, each = length(rcols))
lcols = brewer.pal(10, "RdYlGn")
lcols = c(lcols[1], lcols[10], lcols[1])
lcols = rep(lcols, times = length(dims$Label))


## Link data - Original
df = data.frame(from = rep(rname, times = length(dims$Label)),
                to = rep(dims$Label, each = length(rname)),
                value = 0, value2 = 0,
                stringsAsFactors = FALSE)

for (i in 1:length(df$value)) {
  clus = which(Final$Cluster==rdims$Label[rname==df$from[i]])
  reac = which(Final$Level1==dims$Name[dims$Label==df$to[i]])
  sel = intersect(clus, reac)
  if (length(sel)>0) {
    valp = Final$p.adjust[sel]
    valc = Final$Cramer[sel]
    valg =  length(unique(unlist(strsplit(Final$geneID[sel],"/"))))
    valp = combine.test(valp, method = c("fisher"))
    valc = median(valc)
    
    #df$value[i] = valp
    df$value2[i] = valg
    df$value[i] = length(sel)
  }
}
#df$value[df$value!=0] = -log10(df$value[df$value!=0])

df = df[,c(1:3)]
## Draw 
circos.clear()
circos.par(start.degree = 90, clock.wise = FALSE)
circos.initializeWithIdeogram(plotType = NULL)
circlize_plot_dollar = function() {
  #circos.initialize(NULL)
  chordDiagram(df, annotationTrack = c("grid"),# "axis"),
               grid.col = mycolor,
               col = lcols,
               transparency = 0.3,
               diffHeight = mm_h(3),
               annotationTrackHeight = mm_h(c(3, 0.01)), #2 is default; 2.5 or 3 is great
               big.gap = 25, directional = 2, direction.type = c("arrows", "diffHeight"), link.arr.lwd = 2.2,
               link.arr.length = 0.25, link.arr.col = lcols, link.arr.width = 0.25
               ,link.arr.type = "triangle",
               link.target.prop = F, small.gap = 1.2
               , link.lwd = 1.2, link.border = "black"
               #, link.zindex = rank(df[[3]]))
               , link.sort = T, link.decreasing = T)
  
  #abline(v = 0, lty = 2, col = "#00000080")
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    colidx = get.cell.meta.data("sector.numeric.index")
    
    circos.axis(h = "bottom",direction = "inside", minor.ticks = 1, #h determines if axis is under or over grid
                major.at = seq(from = 0, to = xlim[2]+5, by = 5), 
                labels.cex = 0.9, labels.niceFacing = T)
    
    
    circos.text(mean(xlim), 
                ylim[1], 
                sector.name,
                facing = "inside", # facing determines if text is under or over axis
                cex = 1.6, 
                adj = c(0.5, -1.3), #1.3 adjust is the verticle space
                col = mycolor[colidx],
                #font = par(family = "Serif"),
                niceFacing = T)}, 
    bg.border = NA)
}
circlize_plot_dollar()
opar = par()

tiff(str_c(here("Res/Review/"), "WGCNA_Reactome_Circle_1.tiff"), units="px", width=(6*1200), height=(6*900), res=600)
plot.new()
par(family = "Cambria", font = 2)
circlize_plot_dollar()
dev.off()



## Draw 
circos.clear()
circos.par(start.degree = 90, clock.wise = FALSE)
circos.initializeWithIdeogram(plotType = NULL)
circlize_plot_dollar = function() {
  #circos.initialize(NULL)
  chordDiagram(df, annotationTrack = c("grid"),# "axis"),
               grid.col = mycolor,
               col = lcols,
               transparency = 0.3,
               diffHeight = mm_h(3),
               annotationTrackHeight = mm_h(c(3, 0.01)), #2 is default; 2.5 or 3 is great
               big.gap = 25, directional = 2, direction.type = c("arrows", "diffHeight"), link.arr.lwd = 2.2,
               link.arr.length = 0.25, link.arr.col = lcols, link.arr.width = 0.25
               ,link.arr.type = "none",
               link.target.prop = F, small.gap = 1.2
               , link.lwd = 1.2, link.border = "black"
               #, link.zindex = rank(df[[3]]))
               , link.sort = T, link.decreasing = T)
  
  #abline(v = 0, lty = 2, col = "#00000080")
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    colidx = get.cell.meta.data("sector.numeric.index")
    
    circos.axis(h = "bottom",direction = "inside", minor.ticks = 1, #h determines if axis is under or over grid
                major.at = seq(from = 0, to = xlim[2]+5, by = 5), 
                labels.cex = 0.9, labels.niceFacing = T)
    
    
    circos.text(mean(xlim), 
                ylim[1], 
                sector.name,
                facing = "inside", # facing determines if text is under or over axis
                cex = 1.6, 
                adj = c(0.5, -1.3), #1.3 adjust is the verticle space
                col = mycolor[colidx],
                #font = par(family = "Serif"),
                niceFacing = T)}, 
    bg.border = NA)
}
circlize_plot_dollar()
opar = par()

tiff(str_c(here("Res/Review/"), "WGCNA_Reactome_Circle_2.tiff"), units="px", width=(6*1200), height=(6*900), res=600)
plot.new()
par(family = "Cambria", font = 2)
circlize_plot_dollar()
dev.off()


lgd_list = Legend(
  labels = dims$Name,
  labels_gp = gpar(fontsize = 14, fontfamily = "Calibri"),
  title_gp = gpar(fontsize = 15, font = 2, fontfamily = "Calibri"),
  title = "Reactome Pathways",
  type = "points",
  grid_height = unit(6.5, "mm"), grid_width = unit(6.5, "mm"),
  pch = paste0(dims$Label,""),legend_gp = gpar(col = "white", font = 2, fontfamily = "Cambria"),
  background = ccols, title_gap = unit(4, "mm"))

lgd_groups = Legend(
  labels = paste0("Cluster ",as.character(unique(Final$Cluster))),
  labels_gp = gpar(fontsize = 14, fontfamily = "Calibri"),
  title_gp = gpar(fontsize = 15, font = 2, fontfamily = "Calibri"),
  title = "WGCNA Clusters",
  type = "points",
  grid_height = unit(6.5, "mm"), grid_width = unit(6.5, "mm"),
  pch = c("I", "II", "III"),legend_gp = gpar(col = "white", font = 2, fontfamily = "Cambria"),
  background = rcols, title_gap = unit(4, "mm"))

lgd_reg = Legend(
  labels = c("Downregulated", "Upregulated"),
  labels_gp = gpar(fontsize = 14, fontfamily = "Calibri"),
  title_gp = gpar(fontsize = 15, font = 2, fontfamily = "Calibri"),
  title = "Regulation",
  type = "lines",
  grid_height = unit(6.5, "mm"), grid_width = unit(6.5, "mm"),
  legend_gp = gpar(col = c(lcols[1], lcols[2]),lwd = 4, alpha = 0.75, lineend = "butt"),
  background = "white", title_gap = unit(4, "mm"))



par = opar
tiff(str_c(here("Res/Review/"), "WGCNA_Reactome_Circle2.tiff"), units="px", width=(6*1200), height=(6*900), res=600)
plot.new()
#par(family = "Calibri", font = 2)
grid.draw(lgd_list)
dev.off()

tiff(str_c(here("Res/Review/"), "WGCNA_Reactome_Circle3.tiff"), units="px", width=(6*1200), height=(6*900), res=600)
plot.new()
#par(family = "Calibri", font = 2)
grid.draw(lgd_groups)
dev.off()

tiff(str_c(here("Res/Review/"), "WGCNA_Reactome_Circle4.tiff"), units="px", width=(6*1200), height=(6*900), res=600)
plot.new()
#par(family = "Calibri", font = 2)
grid.draw(lgd_reg)
dev.off()

