data1 = read.table("origin_data//data_clinical_patient.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data2 = read.table("origin_data/data_clinical_sample.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data3 = read.table("origin_data/data_mutations_mskcc.txt",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F)

data4 = merge(data1,data2,by = "PATIENT_ID")
rownames(data4) = data4$SAMPLE_ID
data4 = data4[data4$IMMUNE_TREATMENT == 1,]

sum(data4$SAMPLE_ID %in% data3$Tumor_Sample_Barcode)



length(unique(data4$PATIENT_ID))
length(unique(data4$SAMPLE_ID))
length(unique(data3$Tumor_Sample_Barcode))


# data4 = data4[rownames(data4) %in% data3$Tumor_Sample_Barcode,]  #确保了临床和突变检查的患者是对应的(这个不一定对，因为panel测序会存在某些患者一个基因都不突变)


# 文章虽然提及了治疗药物的情况，但是数据中并没有明确是哪种免疫检查点抑制剂
data = data.frame("OS_TIME" = data4$OS_IO,
                  "OS_STATUS" = data4$OS_STATUS_IO,
                  "PFS_TIME" = data4$PFS_IO,
                  "PFS_STATUS" = data4$PFS_STATUS_IO,
                  "RECIST" = data4$IO_RESPONSE,
                  "TMB" = data4$TMB_NONSYNONYMOUS,
                  "THERAPY" = NA,
                  "origin" = NA
                  )
rownames(data) = rownames(data4) 

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$TMB)
# table(data$RESPONSE)
table(data$RECIST)
sum(is.na(as.matrix(data)))



sum(is.na(data$OS_TIME))
sum(is.na(data$OS_STATUS))
sum(is.na(data$PFS_TIME))
sum(is.na(data$PFS_STATUS))
sum(is.na(data$RECIST))
sum(is.na(data$RESPONSE))

data = data[!data$RECIST == "NE",] #特殊情况



data$OS_STATUS = ifelse(data$OS_STATUS == "DECEASED",1,ifelse(data$OS_STATUS == "LIVING",0,NA))
data$PFS_STATUS = ifelse(data$PFS_STATUS == "Progressed or Deceased",1,ifelse(data$PFS_STATUS == "Alive without progression",0,NA))
# data$RESPONSE = ifelse(data$RESPONSE == "DCB","response",ifelse(data$RESPONSE == "NDB","nonresponse",NA))
data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))

data3 = data3[data3$Tumor_Sample_Barcode %in% rownames(data),] # 不仅去除了其他治疗的患者突变数据，而且避免测了突变，但是没有临床数据

gene_name = unique(data3$Hugo_Symbol)

for(x in gene_name){
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
}

write.table(x = data,file = "result_data/dataset17.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

###############################################

WES = data.frame(
  "ID" = data3$Tumor_Sample_Barcode,
  "Hugo_Symbol" = data3$Hugo_Symbol,
  "Chrmosome" = data3$Chromosome,
  "Start_Position" = data3$Start_Position,
  "End_Position" = data3$End_Position,
  "Variant_Classification" = data3$Variant_Classification,
  "Variant_Type" = data3$Variant_Type,
  "Reference_Allele" = data3$Reference_Allele,
  "Tumor_Seq_Allele1" = data3$Tumor_Seq_Allele1,
  "Tumor_Seq_Allele2" = data3$Tumor_Seq_Allele2)

write.table(x = WES,file = "result_data/dataset17_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
