#### This file is a script that is responsible for performing the DE analysis

#### Prerequisites: 1. BioMart Download of gene symbols (Biomart.R)

#### Run Order = 3

## Setup env
library(here)
dev.off()
source(here("R/startup.R"))


## Biomart (Downloaded 11.12.22)
sapiens_ensembl <- readRDS(str_c(here("data/"), "sapiens_ensembl.rds"))


## Read raw counts from all batches for all subjects <- Brandts data
raw = readRDS("~/Data/wave5/all.batches.expression.set.070121.Rds")
temp = pData(raw)
raw = exprs(raw)
genes = rownames(raw)


## Get pheno data for all subjects from the RDS file created using preprocess <- Wenjia's recoded phenotypes
df = readRDS("/home/share/preprocessing/preprocessed_two_batches/all.batches.expression.set.rawcount_waves_05.08.2021.rds")
pheno = pData(df)
rownames(pheno) = pheno$AID
pheno = pheno[rownames(temp),]
all.equal(rownames(temp), pheno$AID)
pheno$AnyFlag = temp$AnyFlag
pheno$batch = temp$batch
all.equal(colnames(raw), pheno$AID)


# Inflammation report from AddHealth
data_input = "~/Data/wave5/"
df = read.xport(str_c(data_input, "bcrp5.xpt"))
pheno$H5INFECT = as.factor(df$H5INFECT[match(pheno$AID, df$AID)])
pheno$H5SUBCLN = as.factor(df$H5SUBCLN[match(pheno$AID, df$AID)])
pheno$H5CRP8 = df$H5CRP8[match(pheno$AID, df$AID)]

# Alcohol use report from AddHealth
df = read_xpt(str_c(data_input, "WAVE5.xpt"))
pheno$alcohol_use = df$H5TO14[match(pheno$AID, df$AID)]

## Wave 5 Biocovariates for flagging of subjects
df = read_sas(str_c(data_input, "w5biocovars.sas7bdat"))
pheno$KITCOND = df$KITCOND[match(pheno$AID, df$AID)]
pheno$TUBECOND = df$TUBECOND[match(pheno$AID, df$AID)]
pheno$BloodQual = df$Q086[match(pheno$AID, df$AID)]
pheno$QualCode = df$QualCode[match(pheno$AID, df$AID)]
pheno$AlqQual = df$AlqQuality[match(pheno$AID, df$AID)]


## Convert gene names to HGNC
genedat = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters ="ensembl_gene_id", values = genes, mart = sapiens_ensembl)
rawgenedata = data.frame(Ensembl = genes, stringsAsFactors = F)
rawgenedata$HGNC = NA
for (i in 1:length(rawgenedata$Ensembl)) {
  match = which(genedat$ensembl_gene_id==rawgenedata$Ensembl[i])
  if  (!is_empty(match)) {
    rawgenedata$HGNC[i] = genedat$hgnc_symbol[match[1]]
  }
}
keepgenes = which(rawgenedata$HGNC!="" | rawgenedata$HGNC!=NA)
raw = raw[keepgenes,]
genes = genes[keepgenes]
rawgenedata = rawgenedata[keepgenes,]


## Count exploration and data filtering
keepgenes = filterByExpr(raw, min.count = 10, min.pro = 0.05)
raw = raw[keepgenes,]
genes = genes[keepgenes]
rawgenedata = rawgenedata[keepgenes,]
rownames(raw) = rawgenedata$HGNC


## Filter subjects based on QC metrics (AnyFlag for subjects)
sel = which(pheno$AnyFlag==0)
pheno = pheno[sel,]
raw = raw[,sel]


## Recode variables in Pheno data
pheno = pheno %>%
  mutate_at(
    .vars = vars(matches("^edu_p$|^edu_max$")),
    .funs = list(~ .x %>%
                   factor() %>%
                   fct_collapse(
                     "high or less" = "high",
                     "more than high" = c("votec","college","post")
                   ))
  ) %>% 
  mutate_at(vars(c("H5SUBCLN",
                   "H5CRP8",
                   "H5INFECT")), .funs = list(~ ifelse(.x %in% c(0,1,2,3), .x, NA))) %>%
  
  mutate(alcohol_use = case_when(alcohol_use %in% c(0,997) ~ "Never",
                                 alcohol_use %in% c(1,2) ~ "Mild",
                                 alcohol_use %in% c(3,4,5) ~ "Moderate",
                                 alcohol_use > 5 ~ "Severe"                                 
  ))

## Define treatment, control, mediation and disease vars
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


## Filter subjects and genes based on missing values
non_missing = complete.cases(selsubData)
subData = subData[non_missing, ]
subData = droplevels(subData)
selsubData = selsubData[non_missing, ]
selsubData = droplevels(selsubData)
pheno = pheno[non_missing, ]

keepgenes = filterByExpr(raw[,non_missing], min.count = 10, min.pro = 0.05)
dge = DGEList(counts = raw[keepgenes,non_missing])


## Normalization - TMM
dge = calcNormFactors(dge, method = "TMM")


## Linear modeling
for (i in 1:length(treatment)) {
  rhs = str_c(c(treatment[[i]],controls), collapse = " + ")
  model_formula = str_c(" ~ ",rhs) %>% as.formula()
  sub = selsubData %>% dplyr::select(treatment[[i]], all_of(controls))
  sub = droplevels(sub)
  
  if (class(sub[,c(1)])=="factor") {
    sub[,c(1)] = as.numeric(sub[,c(1)])
  }
  
  design  = model.matrix(model_formula, data = sub)
  
  ## voom based linear modeling
  v = voom(dge, design, plot = T)
  
  ## Batch correct using sva::ComBat
  v$E = ComBat(v$E, selsubData$batch)
  
  ## DE
  fit = lmFit(v$E, design) 
  fit = eBayes(fit, trend = T)
  
  res = topTable(fit, coef = treatment[[i]], n= Inf) %>%
    rownames_to_column(var = "gene")
  res$treatment = treatment[[i]]
  saveRDS(res, str_c(here("Res/DE/"), treatment[[i]],".rds"))
  
  res = res[res$adj.P.Val<0.05,]
  res = res[order(-res$logFC),]
  write.table(res, str_c(here("Res/DE/"), treatment[[i]],".txt"), col.names = T, row.names = F, sep = "\t", quote = F)
}
