data3 = read.table("origin_data/data_mutations_mskcc.txt",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F)
data4 = read.table(file = "origin_data/egc_msk_2017_clinical_data (2).tsv",sep = "\t",header = T,quote = "",check.names = F,stringsAsFactors = F)
data4 = data4[!is.na(data4$`IO Response`),] #只留接受免疫治疗的患者

table(data4$Treatment,data4$`IO Regimen Monotherapy Vs Combination`)
table(data4$`IO PDL1 Tested`,data4$`IO Regimen Monotherapy Vs Combination`)
table(data4$`Details on First Line & Comments`,data4$`IO Regimen Monotherapy Vs Combination`)


sum(data4$`Sample ID` %in% data3$Tumor_Sample_Barcode)
data3 = data3[data3$Tumor_Sample_Barcode %in% data4$`Sample ID`,] # 取了患者子集以后记得突变数据也要取子集！！！

data = data.frame("OS_TIME" = data4[,"IO OS Months From Start of IO"],
                  "OS_STATUS" = data4[,"Patient's Vital Status"],
                  "PFS_TIME" = data4[,"IO PFS from Start of IO to POD"] ,
                  "PFS_STATUS" = data4[,"IO Progression Free Status"],
                  "RECIST" = data4[,"IO Response"],
                  "TMB" = data4[,"Mutation Rate"], ## 这也是不真正意义上的的TMB
                  "THERAPY" = data4$`IO Regimen Monotherapy Vs Combination`, # 没有说明具体的免疫检查点抑制剂
                  "origin_therapy" = data4$`IO Regimen Monotherapy Vs Combination`
                  
                  )
rownames(data) = data4$`Sample ID`

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$TMB)
table(data$THERAPY)
sum(is.na(as.matrix(data)))



sum(is.na(data$OS_TIME))
sum(is.na(data$OS_STATUS))
sum(is.na(data$PFS_TIME))
sum(is.na(data$PFS_STATUS))
sum(is.na(data$RECIST))
sum(is.na(data$RESPONSE))

data$OS_STATUS = ifelse(data$OS_STATUS == "DOD",1,ifelse(data$OS_STATUS == "AWD",0,NA))
data$PFS_STATUS = ifelse(data$PFS_STATUS == "Progressed",1,0)
data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("POD","SD"),"PD/SD",NA))
data$THERAPY = ifelse(data$THERAPY == "Monotherapy","anti-PD1/PDL1",
                      ifelse(data$THERAPY == "Combination","anti-PD1/PDL1 + anti-CTLA4",NA))

gene_name = unique(data3$Hugo_Symbol)

for(x in gene_name){
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
}

write.table(x = data,file = "result_data/dataset4.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

###########################################

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

write.table(x = WES,file = "result_data/dataset4_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
