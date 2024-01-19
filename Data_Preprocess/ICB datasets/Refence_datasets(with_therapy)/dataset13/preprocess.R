library(openxlsx)
data1 = read.xlsx("origin_data/aan5951_TableS1.xlsx",sheet = 2,startRow = 2)
data2 = read.xlsx("origin_data/aan5951_TableS1.xlsx",sheet = 3,startRow = 2)
data3 = read.xlsx("origin_data/aan5951_TableS1.xlsx",sheet = 1,startRow = 2)
data4 = read.table("origin_data/ccrcc_dfci_2019/data_clinical_patient.txt",header = T,sep = "\t",row.names = 1)
data5 = read.table("origin_data/ccrcc_dfci_2019/data_clinical_sample.txt",header = T,sep = "\t",row.names = 2)

data4 = data4[data1$patient_id,]
data5 = data5[data1$patient_id,]

data3 = data3[data3$exclusion_reason == 0,]

sum(data1$patient_id %in% data3$patient_id)

length(unique(data2$Tumor_Sample_Barcode))

data = data.frame("OS_TIME" = data1$os_days,
                  "OS_STATUS" = data1$os_censor,
                  "PFS_TIME" = data1$pfs_days,
                  "PFS_STATUS" = data1$pfs_censor,
                  "RECIST" = data1$best_RECIST,
                  "RESPONSE" = data1$response_category,
                  "TMB" = data5$TMB_NONSYNONYMOUS,
                  "THERAPY" = "anti-PD1/PDL1",
                  "origin_therapy" = data4$DRUG
                  )

rownames(data) = data1$patient_id

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)

data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE = ifelse(data$RESPONSE %in% c("clinical benefit"),"response",ifelse(data$RESPONSE %in% c("intermediate benefit","no clinical benefit"),"nonresponse",NA))
data$OS_TIME = data$OS_TIME/30
data$PFS_TIME = data$PFS_TIME/30

data2$PatientID = do.call(rbind,strsplit(x = data2$Tumor_Sample_Barcode,"-"))[,1]
sum(unique(data2$PatientID) %in% rownames(data))


gene_name = unique(data2$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data2[data2$Hugo_Symbol %in% x,"PatientID"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

#该文章的使用了CNSR，此处cnsr = 0才是发生事件(死亡或进展)
m = data$OS_STATUS == 1
data$OS_STATUS[m] = 0
data$OS_STATUS[!m] = 1

m = data$PFS_STATUS == 1
data$PFS_STATUS[m] = 0
data$PFS_STATUS[!m] = 1

write.table(x = data,file = "result_data/dataset14.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

###############################################3

WES = data.frame(
  "ID" = data2$Tumor_Sample_Barcode,
  "Hugo_Symbol" = data2$Hugo_Symbol,
  "Chrmosome" = data2$Chromosome,
  "Start_Position" = data2$Start_position,
  "End_Position" = data2$End_position,
  "Variant_Classification" = data2$Variant_Classification,
  "Variant_Type" = data2$Variant_Type,
  "Reference_Allele" = data2$Reference_Allele,
  "Tumor_Seq_Allele1" = data2$Tumor_Seq_Allele1,
  "Tumor_Seq_Allele2" = data2$Tumor_Seq_Allele2)

write.table(x = WES,file = "result_data/dataset14_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
