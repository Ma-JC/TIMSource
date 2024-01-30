data1 = read.table("origin_data/data_clinical_patient.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data2 = read.table("origin_data/data_clinical_sample.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data3 = read.table("origin_data/data_mutations_mskcc.txt",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F)
data4 = read.table("origin_data/nsclc_pd1_msk_2018_clinical_data .tsv",row.names = 3,sep = "\t",header = T,quote = "",check.names = F,stringsAsFactors = F)

sum(rownames(data1) %in% rownames(data2))
sum(rownames(data4) %in% data2$SAMPLE_ID)
sum(data2$SAMPLE_ID %in% data3$Tumor_Sample_Barcode)
sum(rownames(data4) %in% data3$Tumor_Sample_Barcode)
length(unique(data2$SAMPLE_ID))
length(unique(data3$Tumor_Sample_Barcode))





rownames(data1) = data2[rownames(data1),1]

# data1 = data1[rownames(data1) %in% data3$Tumor_Sample_Barcode,] #确保了临床和突变检查的患者是对应的(排除，这是panel测序，但是此处的结果是一样的)
# data4 = data4[rownames(data4) %in% data3$Tumor_Sample_Barcode,] 

data = data.frame("PFS_TIME" = data4$`Progress Free Survival (Months)`,
                  "PFS_STATUS" = data4$`Progression Free Status`,
                  "RESPONSE" = data4$`Durable Clinical Benefit`,
                  "TMB" = data4$`Mutation Rate`, # 也不是真实意义上的TMB
                  "THERAPY" = data4$`Treatment Type`,
                  "origin_therapy" = data4$`Treatment Type`
                  )
rownames(data) = rownames(data4) 


table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RESPONSE)
table(data$TMB)
table(data$THERAPY)
sum(is.na(as.matrix(data)))


data$PFS_STATUS= ifelse(data$PFS_STATUS == "Progressed",1,0)
data$RESPONSE = sapply(data$RESPONSE,FUN = function(x){switch(x,"YES"="response","NO"="nonresponse","NE" = NA)})
data$THERAPY = ifelse(data$THERAPY == "Combination","anti-PD1/PDL1 + anti-CTLA4",
                      ifelse(data$THERAPY == "Monotherapy","anti-PD1/PDL1",NA))

gene_name = unique(data3$Hugo_Symbol)

for(x in gene_name){
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
}




write.table(x = data,file = "./result_data//dataset5.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


#########################################

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

write.table(x = WES,file = "result_data/dataset5_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
