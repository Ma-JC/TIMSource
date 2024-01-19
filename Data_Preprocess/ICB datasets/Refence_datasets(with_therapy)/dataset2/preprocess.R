############ 该文章数据可从cbioportal和文章附录获得，各信息能对应，但是TMB和治疗类型只在cbioportal中 ##############

data1 = read.csv("./origin_data/Allen TableS2_Revised.csv",row.names = 1)
data2 = read.csv("./origin_data/Allen TableS1.Mutation_list_all_patients.csv")
data3 = read.table("origin_data/skcm_dfci_2015/data_clinical_patient.txt",sep = "\t",header = T,row.names = 1)
data4 = read.table("origin_data/skcm_dfci_2015/data_clinical_sample.txt",sep = "\t",header = T,row.names = 1)

sum(table(c(rownames(data1),rownames(data3),rownames(data4))) == 3)
sum(table(c(rownames(data1),rownames(data3),rownames(data4))) != 3)

sum(rownames(data1) %in% data2$patient)
sum(unique(data2$patient) %in% rownames(data1)) #确保了临床和突变检查的患者是对应的
sum(rownames(data1) %in% rownames(data3))
sum(rownames(data1) %in% rownames(data4))

data3 = data3[rownames(data1),]
data4 = data4[rownames(data1),]

data = data.frame("OS_TIME" = data1$overall_survival,
                  "OS_STATUS" = data1$dead,
                  "PFS_TIME" = data1$progression_free,
                  "PFS_STATUS" = data1$progression,
                  "RECIST" = data1$RECIST,
                  "RESPONSE" = data1$group,
                  "TMB" = data4$TMB_NONSYNONYMOUS,
                  "THERAPY" = data3$TREATMENT,
                  "origin_therapy" = data3$TREATMENT
                  )
rownames(data) = rownames(data1) 

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$THERAPY)

data$OS_TIME = data$OS_TIME/30
data$PFS_TIME = data$PFS_TIME/30
data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE[data$RESPONSE == "long-survival"] = "nonresponse"
data$THERAPY = "anti-CTLA4"



gene_name = unique(data2$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data2[data2$Hugo_Symbol %in% x,"patient"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}
write.table(x = data,file = "./result_data/dataset2.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


#########################################

WES = data.frame(
  "ID" = data2$patient,
  "Hugo_Symbol" = data2$Hugo_Symbol,
  "Chrmosome" = data2$Chromosome,
  "Start_Position" = data2$Start_position,
  "End_Position" = data2$End_position,
  "Variant_Classification" = data2$Variant_Classification,
  "Variant_Type" = data2$Variant_Type,
  "Reference_Allele" = data2$Reference_Allele,
  "Tumor_Seq_Allele1" = data2$Tumor_Seq_Allele1,
  "Tumor_Seq_Allele2" = data2$Tumor_Seq_Allele2)

write.table(x = WES,file = "./result_data/dataset2_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
