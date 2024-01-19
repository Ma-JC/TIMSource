library(openxlsx)
library(maftools)
library(limma)

#######################################################################2021_GBM#################################################################
protein = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(1)2021_GBM/pdc/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(1)2021_GBM/pdc/S057_S048_CPTAC_GBM_Discovery_Cohort_Clinical_Data_Feb2021_r2.xlsx",sheet = 2)
cptac_meta = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(1)2021_GBM/pdc/S057_S048_CPTAC_GBM_Discovery_Cohort_TMT11_CaseID_SampleID_AliquotID_Map_Feb2021_r2.xlsx",startRow = 7)
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",header = T,quote = "",check.names = F)
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")

extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$`Case.ID.(Participant.ID)`,"AliquoID" = cptac_meta$Aliquot.ID, "Sample" = cptac_meta$Sample.type)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID)) & cptac_meta$Sample == "tumor",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]

################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)]
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna)

#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$case_id %in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$case_id %in% unique(submaf@data$TCGA),]


res = list("2021_GBM"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2021_GBM.rds")



#######################################################################2022_PDAC#################################################################
protein = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(2)2022_PDAC(embargo)/pdc/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.tmt11.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(2)2022_PDAC(embargo)/S061_CPTAC_PDA_Discovery_Cohort_Clinical_Data_r1_Feb2021.xlsx",sheet = 2)
cptac_meta = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(2)2022_PDAC(embargo)/S061_CPTAC_PDA_Discovery_Cohort_Specimens_r1_Feb2021.xlsx")
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",check.names = F,header = T,quote = "")
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")


extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$Case.ID,"AliquoID" = cptac_meta$Aliquot.ID, "Sample" = cptac_meta$`Tumor/Normal`,"withdrawn" = cptac_meta$Withdrawn)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID)) & cptac_meta$withdrawn == "no" & cptac_meta$Sample == "Tumor",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]

################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)]
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna)
#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$case_id %in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$case_id %in% unique(submaf@data$TCGA),]

res = list("2022_PDAC"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2022_PDAC.rds")



#######################################################################2020_BI#################################################################
protein = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(3)2020_BI/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(3)2020_BI/S060_S039_CPTAC_Breast_Confirmatory_Study_1year_Followup.xlsx",sheet = 1)
cptac_meta = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(3)2020_BI/S060_S039_Breast_Cancer_Confirmatory_Collection_Specimens_r1.xlsx")
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",header = T,quote = "",check.names = F)
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")


extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$`Participant.Protocol.Identifier.:.Collection.Protocol.Registration`,"AliquoID" = cptac_meta$Specimen.Label, "Sample" = cptac_meta$Sample.Type)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID)) & cptac_meta$Sample == "Tumor",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]

################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA(RNA 或 DNA)

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会根据患者统一为一个

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)]
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna)
#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$Participant.ID%in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$Participant.ID %in% unique(submaf@data$TCGA),]

res = list("2020_BI"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2020_BI.rds")



#######################################################################2020_OVA#################################################################
#该数据集的特殊性，由两个中心合作完成，且两个中心所测患者样本是有重复的！！
protein1 = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(4)2020_OVA/CPTAC2_Ovarian_Prospective_Collection_JHU_Proteome.tmt10.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
protein2 = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(4)2020_OVA/CPTAC2_Ovarian_Prospective_Collection_PNNL_Proteome.tmt10.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
protein1$genes = rownames(protein1)
protein2$genes = rownames(protein2)
protein = merge(x = protein1,y = protein2,by = "genes",all = TRUE)
rownames(protein) = protein$genes
protein$genes = NULL

clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(4)2020_OVA/CPTAC_Ovarian_Cancer_Confirmatory_Study_1Year_Followup_Clinical_Data_S038.xlsx",sheet = 1)
cptac_meta = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(4)2020_OVA/CPTAC_S038_OV_prospective_sample_tumor_normal_status_r1.xlsx")
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",header = T,quote = "",check.names = F)
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")


extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$Participant_ID,"AliquoID" = cptac_meta$Specimen_Lable, "Sample" = cptac_meta$`Tumor/Normal`)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID)) & cptac_meta$Sample == "Tumor",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]

################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA(RNA 或 DNA)

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会根据患者统一为一个

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)]
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna)
#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$Participant.ID%in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$Participant.ID %in% unique(submaf@data$TCGA),]

res = list("2020_OVA"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2020_OVA.rds")


#######################################################################2020_LSCC#################################################################
protein = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(5)2020_LSCC(embargo)/CPTAC3_Lung_Squamous_Cell_Carcinoma_Proteome.tmt11.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(5)2020_LSCC(embargo)/S063_S058_CPTAC_LSCC_Discovery_Cohort_Clinical_Data_r2_July2021.xlsx",sheet = 2,check.names = F)
cptac_meta = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(5)2020_LSCC(embargo)/S063_S058updated_BI_CPTAC3_LSCC_Tumor_Normal_Mapping_r2.xlsx",check.names = F)
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",header = T,quote = "",check.names = F)
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")


extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$Participant_ID,"AliquoID" = cptac_meta$Specimen_Label, "Sample" = cptac_meta$Tumor.or.Normal)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID)) & cptac_meta$Sample == "Tumor",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]

################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA(RNA 或 DNA)

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会根据患者统一为一个

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)] #subaliquot是唯一批次，而此处是唯一患者
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna) #这个数据集中一个测了超过两次的情况比较多
#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$case_id%in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$case_id %in% unique(submaf@data$TCGA),]

res = list("2020_LSCC"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2020_LSCC.rds")



#######################################################################2020_LUAD#################################################################
protein = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(6)2020_LUAD/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(6)2020_LUAD/S046_S056_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r2_July2020.xlsx",sheet = 2,check.names = F)
cptac_meta = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(6)2020_LUAD/S046_S056_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r2_July2020.xlsx",check.names = F)
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",header = T,quote = "",check.names = F)
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")


extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$`Participant.ID.(case_id)`,"AliquoID" = cptac_meta$`Aliquot.(Specimen.Label)`, "Sample" = cptac_meta$Type)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID) | is.na(cptac_meta$Sample)) & cptac_meta$Sample == "Tumor",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]

################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA(RNA 或 DNA)

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会根据患者统一为一个

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)] #subaliquot是唯一批次，而此处是唯一患者
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna) #这个数据集中一个测了超过两次的情况比较多
#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$case_id%in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$case_id %in% unique(submaf@data$TCGA),]

res = list("2020_LUAD"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2020_LUAD.rds")



#######################################################################2020_HNSCC#################################################################
#这个数据集没有现场的样本信息的文件(cptac_meta)，所以我是直接复制PDC官网的表格内容，所以会有点瑕疵~ ~ ~
protein = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(7)2021_HNSCC/CPTAC3_Head_and_Neck_Carcinoma_Proteome.tmt11.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(7)2021_HNSCC/S054_CPTAC_HNSCC_Discovery_Cohort_Clinical_Data_r1_May2020.xlsx",sheet = 2,check.names = F)
cptac_meta = read.csv("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(7)2021_HNSCC/2020_HNSCC_sample.csv",check.names = F)
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",header = T,quote = "",check.names = F)
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")


extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$CaseID,"AliquoID" = cptac_meta$AliquoID, "Sample" = cptac_meta$Sample)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID) | is.na(cptac_meta$Sample)) & cptac_meta$Sample == "Primary Tumor",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]
################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA(RNA 或 DNA)

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会根据患者统一为一个

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)] #subaliquot是唯一批次，而此处是唯一患者
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna) #这个数据集中一个测了超过两次的情况比较多
#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$case_id%in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$case_id %in% unique(submaf@data$TCGA),]

res = list("2020_HNSCC"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2020_HNSCC.rds")



#######################################################################2020_UCEC#################################################################
protein = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(8)2020_UCEC/CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Proteome.tmt10.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(8)2020_UCEC/S053_S043_CPTAC_UCEC_Discovery_Cohort_Clinical_Data_r2_Feb2020.xlsx",sheet = 2,check.names = F)
cptac_meta = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(8)2020_UCEC/S053_S043_CPTAC_UCEC_Discovery_Cohort_Study_Specimens_r2_Feb2020.xlsx",startRow = 6,check.names = F)
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",header = T,quote = "",check.names = F)
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")


extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$`ParticipantID.(Case_ID)`,"AliquoID" = cptac_meta$Aliquot.ID, "Sample" = cptac_meta$Group)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID) | is.na(cptac_meta$Sample)) & cptac_meta$Sample == "Tumor ",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]
################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA(RNA 或 DNA)

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会根据患者统一为一个

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)] #subaliquot是唯一批次，而此处是唯一患者
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna) #这个数据集中一个测了超过两次的情况比较多
#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$case_id%in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$case_id %in% unique(submaf@data$TCGA),]

res = list("2020_UCEC"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2020_UCEC.rds")



#######################################################################2020_CCRCC#################################################################
protein = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(9)2019_CCRCC/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(9)2019_CCRCC/S050_S044_CPTAC_ccRCC_Discovery_Cohort_Clinical_Data_r4_Sept2019.xlsx",sheet = 2,check.names = F)
cptac_meta = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(9)2019_CCRCC/S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.xlsx",check.names = F)
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",header = T,quote = "",check.names = F)
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")


extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$ParticipantID,"AliquoID" = cptac_meta$Aliquot.ID, "Sample" = cptac_meta$Group)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID) | is.na(cptac_meta$Sample)) & cptac_meta$Sample == "Tumor",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]
################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA(RNA 或 DNA)

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会根据患者统一为一个

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)] #subaliquot是唯一批次，而此处是唯一患者
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna) #这个数据集中一个测了超过两次的情况比较多
#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$case_id%in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$case_id %in% unique(submaf@data$TCGA),]

res = list("2020_CCRCC"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2020_CCRCC.rds")



#######################################################################2020_Colon#################################################################
protein = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(11)2019_Colon/CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F,quote = "")
clinical = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(11)2019_Colon/CPTAC_Colon_Cancer_Prospective_Collection_1Year_Followup_S037.xlsx",sheet = 1,check.names = F)
cptac_meta = read.xlsx("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/(11)2019_Colon/CPTAC_Colon_Cancer_Prospective_Collection_Biospecimens_S037.xlsx",check.names = F)
maf = read.maf("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/SNV_CPTAC.maf")
rna = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/result/mRNA_FPKM_CPTAC.csv",sep = ",",header = T,quote = "",check.names = F)
TCGA_meta = read.table("E:/1.MyDataBase/Program/jx/web/SNVIO_V1_complete_data/data/CPTAC/Genomics/aliquot.tsv",sep = "\t",header = T,check.names = F,stringsAsFactors = F,quote = "")


extract_cptac_data = function(data,used_pattern){
  if(used_pattern == "Unshared"){
    m = grepl(pattern = "(.*) Unshared Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }else if(used_pattern == "shared"){
    m = grepl(pattern = "(.*) Log Ratio",colnames(data))
    return(data[4:nrow(data),colnames(data)[m]])
  }
}

##########CPTAC的预处理#############
colnames(protein)

used_pattern = "Unshared"
protein = extract_cptac_data(data = protein,used_pattern = used_pattern)
colnames(protein) = gsub(pattern = "(.*) Unshared Log Ratio",replacement = "\\1",x = colnames(protein))
colnames(protein) = sapply(X = strsplit(colnames(protein),"\\."),FUN = function(x){x[1]})
protein = protein[,!duplicated(colnames(protein))]

cptac_meta = data.frame("CaseID" = cptac_meta$Participant.ID,"AliquoID" = cptac_meta$Specimen.Label, "Sample" = cptac_meta$Pathological.Status)
cptac_meta = cptac_meta[!(is.na(cptac_meta$AliquoID) | is.na(cptac_meta$CaseID) | is.na(cptac_meta$Sample)) & cptac_meta$Sample == "Malignant",]
cptac_meta = cptac_meta[!duplicated(cptac_meta),]
rownames(cptac_meta) = cptac_meta$AliquoID

sum(colnames(protein) %in% rownames(cptac_meta)) #有三个患者没有蛋白质组表达数据

share_name = intersect(colnames(protein),rownames(cptac_meta))
protein = protein[,share_name]
colnames(protein) = cptac_meta[share_name,"CaseID"]

protein_tumor = protein[,!duplicated(colnames(protein))]
################TCGA_maf的预处理################################
################SNVIO的分析都是以突变为大前提！！！#############
patientID = colnames(protein_tumor)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA(RNA 或 DNA)

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID] #此时的protein带有的PatientID将大于或等于MAF的patientID

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% maf@data$TCGA]
subaliquot = unique(tmp_meta[tmp_meta %in% maf@data$TCGA]) #不同批次的测序数据最后都会根据患者统一为一个

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

submaf = maf@data[maf@data$TCGA %in% subaliquot,]
sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
submaf$TCGA = sapply(submaf$TCGA,FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

write.table(submaf,file = "tmp.maf",quote = F,col.names = T,row.names = F,sep = "\t")
submaf = read.maf("tmp.maf")
plotmafSummary(submaf)

protein_tumor = protein_tumor[,unique(submaf@data$TCGA)] #subaliquot是唯一批次，而此处是唯一患者
#####################################
###########TCGA_mRNA预处理###########
patientID = unique(submaf@data$TCGA)
sum(patientID %in% TCGA_meta$case_submitter_id) #有多少CPTAC患者同时在TCGA

#提取与CPTAC患者对应的TCGA相应的aliquot_ID(注意：TCGA_meta包含了RNA和DNA等多个数据类型的ID关系，所以会出现一个patientID对应多个aliquot_ID的情况)
tmp_meta = TCGA_meta$aliquot_submitter_id[TCGA_meta$case_submitter_id %in% patientID]

#下面是检验对应情况，结果看到存在一小部分患者会有多次测序数据
tmp_meta[tmp_meta %in% colnames(rna)]
subaliquot = unique(tmp_meta[tmp_meta %in% colnames(rna)])

sub_meta = TCGA_meta[ TCGA_meta$aliquot_submitter_id %in% subaliquot,c("case_submitter_id","aliquot_submitter_id")]
sub_meta = sub_meta[!duplicated(sub_meta),]

rna_data = rna[,2:ncol(rna)]
rna_data = as.matrix(rna_data)
rownames(rna_data) = rna[,1]
rna = avereps(rna_data)
remove(rna_data)
rna = as.data.frame(rna)

subrna = rna[,subaliquot]
sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})
colnames(subrna) = sapply(colnames(subrna),FUN = function(x){unique(sub_meta$case_submitter_id[sub_meta$aliquot_submitter_id == x])})

subrna = t(subrna)
subrna = avereps(subrna)
subrna = t(subrna)
subrna = as.data.frame(subrna) #这个数据集中一个测了超过两次的情况比较多
#############clinical的预处理#############
sum(colnames(protein_tumor) %in% unique(submaf@data$TCGA)) #基因组和蛋白组的对应
sum(colnames(subrna) %in%  unique(submaf@data$TCGA)) #基因组和转录组的对应
sum(clinical$Participant.ID%in% unique(submaf@data$TCGA)) #与临床信息的对应关系(注意：会出现某些患者被撤销的情况)

clinical = clinical[clinical$Participant.ID %in% unique(submaf@data$TCGA),]

res = list("2019_Colon"=list("protein"=protein_tumor,"rna"=subrna,"maf"=submaf,"clinical_detial"=clinical))

saveRDS(res,"2019_Colon.rds")

