Colon = readRDS("2019_Colon.rds")
BI = readRDS("2020_BI.rds")
CCRCC = readRDS("2020_CCRCC.rds")
HNSCC = readRDS("2020_HNSCC.rds")
LSCC = readRDS("2020_LSCC.rds")
LUAD = readRDS("2020_LUAD.rds")
OVA = readRDS("2020_OVA.rds")
UCEC = readRDS("2020_UCEC.rds")
GBM = readRDS("2021_GBM.rds")
PDAC = readRDS("2022_PDAC.rds")

#PDC下载的蛋白组数据是可以直接处理的，因为已经经过常规的预处理(除了缺失值插补)，而RNA数据是FPKM的，需要log(x + 1)后使用

# require("cluster")
# require("survival")
# require("randomForest")
# require("missForest") #需要R 4.0.5版本
# require("glmnet")
# require("Rcpp")
# require("foreach")
# require("itertools")
# require("iterators")
# require("Matrix")
# require("devtools")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("impute")
# require("impute")
# require("remotes")
# install_github("WangLab-MSSM/DreamAI/Code")
# BiocManager::install("GSVA")

####################################################################缺失值插补##############################################################################################

library(DreamAI)
total = list("2020_BI"=BI[[1]],"2019_Colon"=Colon[[1]],"2020_CCRCC"=CCRCC[[1]],"2020_HNSCC"=HNSCC[[1]],"2020_LSCC"=LSCC[[1]],"2020_LUAD"=LUAD[[1]],"2020_OVA"=OVA[[1]],"2020_UCEC"=UCEC[[1]],"2021_GBM"=GBM[[1]],"2021_PDAC"=PDAC[[1]])
for(i in names(total)){
  tmp = apply(total[[i]]$protein,1,function(x){
    if(sum(is.na(x)) < length(x)/2){
      return(x)
    }
  })
  tmp = do.call("rbind",tmp)
  total[[i]]$protein = DreamAI(data = tmp)$Ensemble
}

saveRDS(object = total,file = "CPTAC.rds")

total = readRDS("CPTAC.rds")

#################################################################RNA数值预处理########################################################################################

for(i in names(total)){
  total[[i]]$rna = log2(total[[i]]$rna + 1)
  
}
####################################################################免疫相关通路##############################################################################################
immune = list(
  "6-gene IFN signature" = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG"),
  "18-gene IFN signature" = c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB"),
  "Gene expression profile" = c("CXCR6", "TIGIT", "CD27", "CD274", "PDCD1LG2", "LAG3", "NKG7", "PSMB10", "CMKLR1", "CD8A", "IDO1", "CCL5", "CXCL9", "HLA-DQA1", "CD276", "HLA-DRB1", "STAT1", "HLA-E"),
  "Cytolytic activity" = c("GZMA", "PRF1"),
  "13 T-cell signature" = c("CD8A", "CCL2", "CCL3", "CCL4", "CXCL9", "CXCL10", "ICOS", "GZMK", "IRF1","HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB"),
  "Effective T cell score" = c("CD8A", "EOMES", "PRF1", "IFNG", "CD274"),
  "Immune checkpoint expression" = c("CD274", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "PDCD1LG2", "TIGIT"),
  "TLS" = c("CCL21", "CCL19", "CXCL13", "CXCL11", "CCL8", "CXCL10", "CXCL9", "CCL2", "CCL3", "CCL18", "CCL5")
)

for(n in names(total)){
  immune_pathway = vector()
  for( i in names(immune)){
    data = total[[n]]$rna
    tmp_genes = immune[[i]][immune[[i]] %in% rownames(data)]
    tmp_data = data[tmp_genes,]
    
    immune_pathway = rbind(immune_pathway,colMeans(tmp_data))
    
  } 
  rownames(immune_pathway) = names(immune)
  total[[n]]$immune_pathway = immune_pathway
}

#######################################################################免疫浸润###################################################################################
immune_infiltration_data = read.table(file = "immune_infiltration.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")
celltype_geneset = list()
for( i in immune_infiltration_data$`Cell type`){
  celltype_geneset[[i]] = immune_infiltration_data$Metagene[ immune_infiltration_data$`Cell type` == i]
}

library(GSVA)
for(i in names(total)){
  
  ssgsea_score = gsva(as.matrix(total[[i]]$rna), celltype_geneset, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE,parallel.sz=12)
  total[[i]]$immune_infiltration = ssgsea_score
  
}

####################################################################免疫相关通路(基于蛋白质组的)##############################################################################################
immune = list(
  "6-gene IFN signature" = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG"),
  "18-gene IFN signature" = c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB"),
  "Gene expression profile" = c("CXCR6", "TIGIT", "CD27", "CD274", "PDCD1LG2", "LAG3", "NKG7", "PSMB10", "CMKLR1", "CD8A", "IDO1", "CCL5", "CXCL9", "HLA-DQA1", "CD276", "HLA-DRB1", "STAT1", "HLA-E"),
  "Cytolytic activity" = c("GZMA", "PRF1"),
  "13 T-cell signature" = c("CD8A", "CCL2", "CCL3", "CCL4", "CXCL9", "CXCL10", "ICOS", "GZMK", "IRF1","HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB"),
  "Effective T cell score" = c("CD8A", "EOMES", "PRF1", "IFNG", "CD274"),
  "Immune checkpoint expression" = c("CD274", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "PDCD1LG2", "TIGIT"),
  "TLS" = c("CCL21", "CCL19", "CXCL13", "CXCL11", "CCL8", "CXCL10", "CXCL9", "CCL2", "CCL3", "CCL18", "CCL5")
)

for(n in names(total)){
  immune_pathway = vector()
  rn = vector()
  for( i in names(immune)){
    data = total[[n]]$protein
    tmp_genes = immune[[i]][immune[[i]] %in% rownames(data)]
    tmp_data = data[tmp_genes,]
    if(length(tmp_genes) >= 2){  #每个通路至少包含2个基因才纳入计算中
      immune_pathway = rbind(immune_pathway,colMeans(tmp_data))
      rn = c(rn,i)
    }
    
    
  } 
  rownames(immune_pathway) = rn
  total[[n]]$immune_pathway_protein = immune_pathway
}

#######################################################################免疫浸润(基于蛋白质组的)###################################################################################
immune_infiltration_data = read.table(file = "immune_infiltration.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")
celltype_geneset = list()
for( i in immune_infiltration_data$`Cell type`){
  celltype_geneset[[i]] = immune_infiltration_data$Metagene[ immune_infiltration_data$`Cell type` == i]
}

library(GSVA)
for(i in names(total)){
  
  ssgsea_score = gsva(as.matrix(total[[i]]$protein), celltype_geneset, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE,parallel.sz=12)
  total[[i]]$immune_infiltration_protein = ssgsea_score
  
}

saveRDS(object = total,file = "CPTAC.rds")
##########细节##################

CPTAC = readRDS("CPTAC.rds")
for(i in names(CPTAC)){
  CPTAC[[i]]$protein = as.data.frame(CPTAC[[i]]$protein)
}
tmp = CPTAC[[i]]$protein

saveRDS(object = CPTAC,file = "CPTAC.rds")
###################################################################整理临床信息#########################################################################

total = readRDS("CPTAC.rds")

####2020_BI####

clinical = total[["2020_BI"]]$clinical_detial
OS.time = clinical[,3]
OS.time[is.na(OS.time)] = clinical[,5][is.na(OS.time)]
OS = clinical$`Vital.Status.(at.time.of.last.contact)`
OS = ifelse(OS == "Living",0,ifelse(OS == "Deceased",1,NA))
survival_data = data.frame("PatientID" = clinical$Participant.ID,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2020_BI"]]$Survival = survival_data

####2019_Colon####
clinical = total[["2019_Colon"]]$clinical_detial
OS.time = clinical[,3]
OS.time[is.na(OS.time)] = clinical[,5][is.na(OS.time)]
OS = clinical$`Vital.Status.(at.time.of.last.contact)`
OS = ifelse(OS == "Living",0,ifelse(OS == "Deceased",1,NA))
survival_data = data.frame("PatientID" = clinical$Participant.ID,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2019_Colon"]]$Survival = survival_data

####2020_CCRCC###
clinical = total[["2020_CCRCC"]]$clinical_detial

tmp = colnames(clinical)[grepl(pattern = "follow_up",colnames(clinical))]

c(        tmp[grepl("vital_status",tmp)],
          tmp[grepl("initial.*contact",tmp)],
          tmp[grepl("initial.*death",tmp)]
  )

tmp_OS = clinical[,tmp[grepl("vital_status",tmp)]]
m0 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(is.na(x)) == length(x) })
m1 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") > 0})
m2 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") == 0})
tmp_OS$OS = NA
tmp_OS$OS[m1] = 1
tmp_OS$OS[m2] = 0
tmp_OS$OS[m0] = NA
OS = tmp_OS$OS

tmp_OS.time = clinical[,c(tmp[grepl("initial.*contact",tmp)],tmp[grepl("initial.*death",tmp)])]
m = apply(tmp_OS.time,1,function(x){sum(is.na(x)) == length(x)}) 
tmp_OS.time$OS.time = NA
tmp_OS.time$OS.time = apply(tmp_OS.time,1,function(x){max(x[!is.na(x)])})
tmp_OS.time$OS.time[m] = NA
OS.time = tmp_OS.time$OS.time
OS.time = as.numeric(OS.time) 
OS.time[OS.time < 0] = 0
survival_data = data.frame("PatientID" = clinical$case_id,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2020_CCRCC"]]$Survival = survival_data

#####2020_HNSCC#####
clinical = total[["2020_HNSCC"]]$clinical_detial

tmp = colnames(clinical)[grepl(pattern = "follow_up",colnames(clinical))]


c(        tmp[grepl("vital_status",tmp)],
          tmp[grepl("initial.*contact",tmp)],
          tmp[grepl("initial.*death",tmp)]
)


tmp_OS = clinical[,tmp[grepl("vital_status",tmp)]]
m0 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(is.na(x)) == length(x) })
m1 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") > 0})
m2 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") == 0})
tmp_OS$OS = NA
tmp_OS$OS[m1] = 1
tmp_OS$OS[m2] = 0
tmp_OS$OS[m0] = NA
OS = tmp_OS$OS

tmp_OS.time = clinical[,c(tmp[grepl("initial.*contact",tmp)],tmp[grepl("initial.*death",tmp)])]
m = apply(tmp_OS.time,1,function(x){sum(is.na(x)) == length(x)}) 
tmp_OS.time$OS.time = NA
tmp_OS.time$OS.time = apply(tmp_OS.time,1,function(x){max(x[!is.na(x)])})
tmp_OS.time$OS.time[m] = NA
OS.time = tmp_OS.time$OS.time
OS.time = as.numeric(OS.time)
OS.time[OS.time < 0] = 0
survival_data = data.frame("PatientID" = clinical$case_id,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2020_HNSCC"]]$Survival = survival_data

#####2020_LSCC#####
clinical = total[["2020_LSCC"]]$clinical_detial

tmp = colnames(clinical)[grepl(pattern = "follow_up",colnames(clinical))]

c(        tmp[grepl("vital_status",tmp)],
          tmp[grepl("initial.*contact",tmp)],
          tmp[grepl("initial.*death",tmp)]
)

tmp_OS = clinical[,tmp[grepl("vital_status",tmp)]]
m0 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(is.na(x)) == length(x) })
m1 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") > 0})
m2 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") == 0})
tmp_OS$OS = NA
tmp_OS$OS[m1] = 1
tmp_OS$OS[m2] = 0
tmp_OS$OS[m0] = NA
OS = tmp_OS$OS

tmp_OS.time = clinical[,c(tmp[grepl("initial.*contact",tmp)],tmp[grepl("initial.*death",tmp)])]
m = apply(tmp_OS.time,1,function(x){sum(is.na(x)) == length(x)}) 
tmp_OS.time$OS.time = NA
tmp_OS.time$OS.time = apply(tmp_OS.time,1,function(x){max(x[!is.na(x)])})
tmp_OS.time$OS.time[m] = NA
OS.time = tmp_OS.time$OS.time
OS.time = as.numeric(OS.time) #在里面有一个unknow，但是在强制转数值型变量时，会被变成NA
OS.time[OS.time < 0] = 0
survival_data = data.frame("PatientID" = clinical$case_id,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2020_LSCC"]]$Survival = survival_data

#####2020_LUAD#####
clinical = total[["2020_LUAD"]]$clinical_detial

tmp = colnames(clinical)[grepl(pattern = "follow_up",colnames(clinical))]

c(        tmp[grepl("vital_status",tmp)],
          tmp[grepl("initial.*contact",tmp)],
          tmp[grepl("initial.*death",tmp)]
)

tmp_OS = clinical[,tmp[grepl("vital_status",tmp)]]
m0 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(is.na(x)) == length(x) })
m1 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") > 0})
m2 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") == 0})
tmp_OS$OS = NA
tmp_OS$OS[m1] = 1
tmp_OS$OS[m2] = 0
tmp_OS$OS[m0] = NA
OS = tmp_OS$OS

tmp_OS.time = clinical[,c(tmp[grepl("initial.*contact",tmp)],tmp[grepl("initial.*death",tmp)])]
m = apply(tmp_OS.time,1,function(x){sum(is.na(x)) == length(x)}) 
tmp_OS.time$OS.time = NA
tmp_OS.time$OS.time = apply(tmp_OS.time,1,function(x){max(x[!is.na(x)])})
tmp_OS.time$OS.time[m] = NA
OS.time = tmp_OS.time$OS.time
OS.time = as.numeric(OS.time) 
OS.time[OS.time < 0] = 0
survival_data = data.frame("PatientID" = clinical$case_id,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2020_LUAD"]]$Survival = survival_data

#####2020_OVA#####
clinical = total[["2020_OVA"]]$clinical_detial

OS.time = clinical[,3]
OS.time[is.na(OS.time)] = clinical[,5][is.na(OS.time)]
OS = clinical$`Vital.Status.(at.time.of.last.contact)`
OS = ifelse(OS == "Living",0,ifelse(OS == "Deceased",1,NA))
survival_data = data.frame("PatientID" = clinical$Participant.ID,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2020_OVA"]]$Survival = survival_data

#####2020_UCEC#####
clinical = total[["2020_UCEC"]]$clinical_detial

tmp = colnames(clinical)[grepl(pattern = "follow_up",colnames(clinical))]

c(        tmp[grepl("vital_status",tmp)],
          tmp[grepl("initial.*contact",tmp)],
          tmp[grepl("initial.*death",tmp)]
)

tmp_OS = clinical[,tmp[grepl("vital_status",tmp)]]
m0 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(is.na(x)) == length(x) })
m1 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") > 0})
m2 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") == 0})
tmp_OS$OS = NA
tmp_OS$OS[m1] = 1
tmp_OS$OS[m2] = 0
tmp_OS$OS[m0] = NA
OS = tmp_OS$OS

tmp_OS.time = clinical[,c(tmp[grepl("initial.*contact",tmp)],tmp[grepl("initial.*death",tmp)])]
m = apply(tmp_OS.time,1,function(x){sum(is.na(x)) == length(x)}) 
tmp_OS.time$OS.time = NA
tmp_OS.time$OS.time = apply(tmp_OS.time,1,function(x){max(x[!is.na(x)])})
tmp_OS.time$OS.time[m] = NA
OS.time = tmp_OS.time$OS.time
OS.time = as.numeric(OS.time) 
OS.time[OS.time < 0] = 0
survival_data = data.frame("PatientID" = clinical$case_id,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2020_UCEC"]]$Survival = survival_data

#####2021_GBM#####
clinical = total[["2021_GBM"]]$clinical_detial

tmp = colnames(clinical)[grepl(pattern = "follow_up",colnames(clinical))]

c(        tmp[grepl("vital_status",tmp)],
          tmp[grepl("initial.*contact",tmp)],
          tmp[grepl("initial.*death",tmp)]
)

tmp_OS = clinical[,tmp[grepl("vital_status",tmp)]]
m0 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(is.na(x)) == length(x) })
m1 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") > 0})
m2 = apply(X = tmp_OS,MARGIN = 1,FUN = function(x){sum(x[!is.na(x)] == "Deceased") == 0})
tmp_OS$OS = NA
tmp_OS$OS[m1] = 1
tmp_OS$OS[m2] = 0
tmp_OS$OS[m0] = NA
OS = tmp_OS$OS

tmp_OS.time = clinical[,c(tmp[grepl("initial.*contact",tmp)],tmp[grepl("initial.*death",tmp)])]
m = apply(tmp_OS.time,1,function(x){sum(is.na(x)) == length(x)}) 
tmp_OS.time$OS.time = NA
tmp_OS.time$OS.time = apply(tmp_OS.time,1,function(x){max(x[!is.na(x)])})
tmp_OS.time$OS.time[m] = NA
OS.time = tmp_OS.time$OS.time
OS.time = as.numeric(OS.time) #在里面有一个unknow，但是在强制转数值型变量时，会被变成NA
OS.time[OS.time < 0] = 0
survival_data = data.frame("PatientID" = clinical$case_id,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2021_GBM"]]$Survival = survival_data

#####2021_PDAC#####
clinical = total[["2021_PDAC"]]$clinical_detial

tmp = colnames(clinical)

c(        tmp[grepl("vital_status",tmp)],
          tmp[grepl("initial.*contact",tmp)],
          tmp[grepl("initial.*death",tmp)]
)

OS = clinical[,tmp[grepl("vital_status",tmp)]]
OS = ifelse(OS == "Deceased",1,ifelse(OS == "Living",0,NA))

OS.time = clinical[,tmp[grepl("initial.*contact",tmp)]][,1]

survival_data = data.frame("PatientID" = clinical$case_id,"OS" = OS,"OS.time" = OS.time)
survival_data = survival_data[!apply(is.na(survival_data),1,any),]

total[["2021_PDAC"]]$Survival = survival_data

saveRDS(total,"CPTAC.rds")

#############################################NOTCH4在2021GBM中是明显的离群值###########################

total = readRDS("CPTAC.rds")
tmp = total$`2021_GBM`$protein

tmp = tmp[-which(rownames(tmp) == "NOTCH4"),]
total$`2021_GBM`$protein = tmp

saveRDS(total,"CPTAC.rds")
