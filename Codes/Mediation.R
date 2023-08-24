#### This file is a script that is responsible for performing the mediation analysis

#### Prerequisites: 1. Expression_Set_Creation.R

## Setup env
rm(list=ls(all=TRUE))
library(here)
library(mediation)
library(stringr)
library(purrr)
library(foreach)
library(parallel)
library(doParallel)
library(bigmemory)

## Load expression data
dat = readRDS(str_c(here("data/"), "/Filtered_ExpressionSet_Clustering_SelSubData.rds"))
exprs = Biobase::exprs(dat)

## Define treatment, control, mediation and disease vars
subData = Biobase::pData(dat)
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

subData = subData %>% dplyr::select(AID,all_of(treatment), all_of(controls), all_of(mediators))
non_missing = complete.cases(subData)
subData = subData[non_missing, ]
subData = droplevels(subData)
exprs = exprs[,non_missing]
all.equal(colnames(exprs), subData$AID)

subData$alcohol_use = factor(subData$alcohol_use, levels = c("Never", "Mild", "Moderate", "Severe"))
subData[,which(colnames(subData) %in% treatment | colnames(subData) %in% mediators)] <- sapply(subData[,which(colnames(subData) %in% treatment | colnames(subData) %in% mediators)], as.numeric)

## Mediation function and extraction
med_fun = function(treatment, controls, mediator, subjectData, expressionData) {
  ## Data
  data = data.frame(Y = expressionData, subjectData)
  ## fit.x
  mod_formula = as.formula(str_c(mediator," ~ ",str_c(c(treatment,controls), collapse = " + ")))
  fit_x = lm(mod_formula, data)
  ## fit.y
  mod_formula = as.formula(str_c("Y ~ ",str_c(c(treatment,controls, mediator), collapse = " + ")))
  fit_y = lm(mod_formula,data)
  ## Mediate
  m = mediate(fit_x, fit_y, treat = treatment, mediator = mediator)
  ## Extract results
  ACME_sd = m$d0.sims %>% sd
  ADE_sd = m$z0.sims %>% sd
  Total_sd = m$tau.sims %>% sd 
  y_sd = summary(m$model.y)$sigma
  p = m$d1.p
  med_prop = m$n1
  med_ACME = m$d1
  med_ADE = m$z1
  med_ACME_p = m$d1.p
  med_ADE_p = m$z1.p
  med_Total = m$tau.coef
  allres = c(med_prop,p,med_ACME, med_ADE,med_Total, med_ACME_p, med_ADE_p,ACME_sd, ADE_sd,Total_sd, y_sd)
  return(allres)
}

## Extract genes for mediation testing ##
name = c("DE Genes", "Directly affected TF's", "First Neighbors of affected TF's", "TF's from DE genes")
glist = list()
for (i in 1:length(name)) {
  df = read.table(str_c(here("Res/Network/SES/Up/Genes"),"/",name[[i]], ".txt"), sep = "\t")
  glist[[i]] = df$V1
  df = read.table(str_c(here("Res/Network/SES/Down/Genes"),"/",name[[i]], ".txt"), sep = "\t")
  glist[[i+4]] = df$V1
}

sel = unlist(glist)
exprs = exprs[rownames(exprs) %in% sel,]
genes = rownames(exprs)

treatment = treatment[[1]]
treat = rep(treatment, length(mediators))
treat = rep(treat, nrow(exprs))
meds = rep(mediators, each = length(treatment))
meds = rep(meds, nrow(exprs))
exprs = do.call(rbind, replicate((length(treatment)*length(mediators)), exprs, simplify=FALSE))
genes = rep(genes,length(treatment)*length(mediators))
exp = as.big.matrix(exprs)
exp1 = describe(exp)
rm(exprs)
strt<-Sys.time()
cl <- parallel::makeCluster(24)
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, library(stringr))
clusterEvalQ(cl, library(mediation))
clusterEvalQ(cl, library(bigmemory))

#clusterExport(cl=cl, c('exprs', 'subData','treat','controls','meds'))
fullout = foreach(i=1:length(treat) ,.packages = "foreach") %dopar% {  
  exp = attach.big.matrix(exp1)
  t1 = as.numeric(exp[i,])
  out = med_fun(treat[[i]], controls, meds[[i]], subData, t1)
  out
}
parallel::stopCluster(cl)
print(Sys.time()-strt)
fullout = setNames(fullout, genes)
saveRDS(fullout, str_c(here("Res/"), "/Mediation/SES/Filtered_genes.rds"))
saveRDS(meds, str_c(here("Res/"), "/Mediation/SES/Mediators.rds"))
saveRDS(genes, str_c(here("Res/"), "/Mediation/SES/Genes.rds"))

