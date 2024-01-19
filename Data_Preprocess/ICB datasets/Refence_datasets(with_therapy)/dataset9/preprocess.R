data1 = read.table(file = "origin_data/data_clinical_patient.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
data2 = read.table(file = "origin_data/data_clinical_sample.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
data3 = read.table(file = "origin_data/data_mutations_mskcc.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")

sum(data1$PATIENT_ID %in% data2$PATIENT_ID)
sum(data1$PATIENT_ID %in% unique(data3$Tumor_Sample_Barcode))



rownames(data1) = data1$PATIENT_ID
rownames(data2) = data2$PATIENT_ID
nn = intersect(data1$PATIENT_ID,data2$PATIENT_ID)
data1 = data1[nn,]
data2 = data2[nn,]
data = data.frame("OS_TIME" = data1$OS_MONTHS,
                  "OS_STATUS" = data1$OS_STATUS,
                  "RESPONSE" = data1$DURABLE_CLINICAL_BENEFIT,
                  "TMB" = data2$TMB_NONSYNONYMOUS,
                  "THERAPY" = "anti-CTLA4",
                  "origin_therapy" = data1$TREATMENT
                  )

data$RESPONSE = ifelse(data$RESPONSE == "LB","response",ifelse(data$RESPONSE == "NB","nonresponse",NA))
data$OS_STATUS = ifelse(data$OS_STATUS == "0:LIVING",0,ifelse(data$OS_STATUS == "1:DECEASED",1,NA))
rownames(data) = data1$PATIENT_ID

gene_name = unique(data3$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
}


#这里的OS时间(cbioportal下载的)与附录PDF整理出来的不一样，其他的大致相同(其实突变列表也不尽相同)
write.table(x = data,file = "result_data/dataset10.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


###############################################3

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

write.table(x = WES,file = "result_data/dataset10_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
