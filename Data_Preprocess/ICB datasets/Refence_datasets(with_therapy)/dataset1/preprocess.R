data1 = read.table("./origin_data/data_clinical_patient.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data2 = read.table("./origin_data/data_clinical_sample.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data3 = read.table("./origin_data/data_mutations_mskcc.txt",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F)


sum(rownames(data1) %in% rownames(data2))
sum(data2$SAMPLE_ID %in% data3$Tumor_Sample_Barcode)
length(unique(data2$SAMPLE_ID))
length(unique(data3$Tumor_Sample_Barcode))

pie(table(data1$DRUG_TYPE))

rownames(data1) = data2[rownames(data1),1]

# data1 = data1[rownames(data1) %in% data3$Tumor_Sample_Barcode,] #确保了临床和突变检查的患者是对应的(排除这一步，这里是panel测序)

data = data.frame("OS_TIME" = data1$OS_MONTHS,
                  "OS_STATUS" = data1$OS_STATUS,
                  "TMB" = data1$TMB_SCORE,
                  "THERAPY" = data1$DRUG_TYPE,
                  "origin_therapy" = data1$DRUG_TYPE)

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$TMB)
table(data$THERAPY)
sum(is.na(as.matrix(data)))

data$OS_STATUS = ifelse(data$OS_STATUS == "DECEASED",1,ifelse(data$OS_STATUS == "LIVING",0,NA))
data$THERAPY = ifelse(data$THERAPY == "CTLA4","anti-CTLA4",
                      ifelse(data$THERAPY == "PD-1/PDL-1","anti-PD1/PDL1",
                             ifelse(data$THERAPY == "Combo","anti-PD1/PDL1 + anti-CTLA4",NA)))
rownames(data) = rownames(data1) 
gene_name = unique(data3$Hugo_Symbol)

for(x in gene_name){
  
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

write.table(x = data,file = "./result_data/dataset1.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


########################

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

write.table(x = WES,file = "./result_data/dataset1_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)


##############################################

data = read.table(file = "./result_data/dataset1.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
WES = read.table(file = "./result_data/dataset1_mutation.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data2 = read.table("./origin_data//data_clinical_sample.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")

item_id = c("OS_TIME","OS_STATUS","PFS_TIME","PFS_STATUS","TMB","RECIST","RESPONSE","THERAPY","origin_therapy")

for(i in unique(data2$CANCER_TYPE)[1:10]){
  
  tmp = data[data2$SAMPLE_ID[ data2$CANCER_TYPE == i],]
  tmp = tmp[,setdiff(colnames(tmp),setdiff(colnames(tmp)[colSums(tmp == "Mutation") == 0],item_id))]
  write.table(x = tmp,
              file = paste("./result_data/dataset1_",paste(strsplit(i," ")[[1]],collapse = "_"),".txt",sep = ""),sep = "\t",quote = FALSE,row.names = T,col.names = T)
  
  write.table(x = WES[ WES$ID %in% data2$SAMPLE_ID[ data2$CANCER_TYPE == i],],
              file = paste("./result_data/dataset1_mutation_",paste(strsplit(i," ")[[1]],collapse = "_"),".txt",sep = ""),sep = "\t",quote = FALSE,row.names = F,col.names = T)
}



