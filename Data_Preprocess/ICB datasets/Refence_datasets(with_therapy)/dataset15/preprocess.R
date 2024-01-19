data1 = read.table("origin_data/data_clinical_patient.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data2 = read.table("origin_data/data_clinical_sample.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data3 = read.table("origin_data/data_mutations_mskcc.txt",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F)

data4 = merge(data1,data2,by = "PATIENT_ID")
rownames(data4) = data4$PATIENT_ID


sum(data4$PATIENT_ID %in% data3$Tumor_Sample_Barcode)



length(unique(data4$PATIENT_ID))
length(unique(data3$Tumor_Sample_Barcode))


# data4 = data4[rownames(data4) %in% data3$Tumor_Sample_Barcode,]  #确保了临床和突变检查的患者是对应的

data = data.frame("PFS_TIME" = data4$PFS_MONTHS,
                  "PFS_STATUS" = data4$PFS_STATUS,
                  "RECIST" = data4$OVERALL_RESPONSE,
                  "RESPONSE" = data4$DURABLE_CLINICAL_BENEFIT,
                  "TMB" = data4$TMB_NONSYNONYMOUS,
                  "THERAPY" = "anti-PD1/PDL1",
                  "origin_therapy" = NA
                  )
rownames(data) = rownames(data4) 

# table(data$OS_TIME)
# table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$TMB)
table(data$RESPONSE)
table(data$RECIST)
sum(is.na(as.matrix(data)))



sum(is.na(data$OS_TIME))
sum(is.na(data$OS_STATUS))
sum(is.na(data$PFS_TIME))
sum(is.na(data$PFS_STATUS))
sum(is.na(data$RECIST))
sum(is.na(data$RESPONSE))


data$PFS_STATUS = ifelse(data$PFS_STATUS == "1:Event",1,ifelse(data$PFS_STATUS == "0:Censure",0,NA))
data$RESPONSE = ifelse(data$RESPONSE == "DCB","response",ifelse(data$RESPONSE == "NDB","nonresponse",NA))
data$RECIST = ifelse(data$RECIST %in% c("Partial Response"),"CR/PR",ifelse(data$RECIST %in% c("Progression of Disease","Stable Response"),"PD/SD",NA))

gene_name = unique(data3$Hugo_Symbol)

for(x in gene_name){
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
}

write.table(x = data,file = "result_data/dataset16.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

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

write.table(x = WES,file = "result_data/dataset16_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
