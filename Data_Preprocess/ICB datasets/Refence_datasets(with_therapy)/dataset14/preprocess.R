data1 = read.table("origin_data/data_clinical_patient.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data2 = read.table("origin_data/data_clinical_sample.txt",header = T,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data3 = read.table("origin_data/data_mutations_mskcc.txt",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F)
data5 = openxlsx::read.xlsx("origin_data/41591_2019_349_MOESM1_ESM.xlsx",sheet = 1,startRow = 2)

data4 = merge(data1,data2,by = "PATIENT_ID")
rownames(data4) = data4$PATIENT_ID


# 这篇文章只有一部分患者测了WES,文章附件有
sn = data4$PATIENT_ID[which(gsub(pattern = "_",replacement = " ",do.call(rbind,strsplit(rownames(data4),"2019_"))[,2])
      %in% data5$`Patient.#`[data5$`Sequenced?` == "Yes"])]

sum(data4$PATIENT_ID %in% data3$Tumor_Sample_Barcode)# 结果表明，文章提到测了WES的患者有17个
length(unique(data3$Tumor_Sample_Barcode)) # 但是WES数据中检测到有突变的患者只有15个
sum(sn %in% data3$Tumor_Sample_Barcode)  # 而这15个患者都是文章提到测序的患者，可能有两个患者不存在任何突变？？(GBM突变频率低)


data4 = data4[rownames(data4) %in% data3$Tumor_Sample_Barcode,]  #这一步本不应该执行，但是那两个不存在任何突变的患者不清楚是否正常，所以还是去除吧

data = data.frame("OS_TIME" = data4$OS_FROM_PD1I_MONTHS,
                  "OS_STATUS" = data4$OS_FROM_PD1I_STATUS,
                  "PFS_TIME" = data4$PFS_MONTHS,
                  "PFS_STATUS" = data4$PFS_STATUS,
                  "RESPONSE" = data4$RESPONSE,
                  "TMB" = data4$TMB_NONSYNONYMOUS,
                  "THERAPY" = "anti-PD1/PDL1",
                  "origin_therapy" = data4$PD1_INHIBITOR_DRUG
                  )
rownames(data) = rownames(data4) 

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$TMB)
table(data$RESPONSE)
sum(is.na(as.matrix(data)))



sum(is.na(data$OS_TIME))
sum(is.na(data$OS_STATUS))
sum(is.na(data$PFS_TIME))
sum(is.na(data$PFS_STATUS))
sum(is.na(data$RECIST))
sum(is.na(data$RESPONSE))

data$OS_STATUS = ifelse(data$OS_STATUS == "1:DECEASED",1,ifelse(data$OS_STATUS == "0:LIVING",0,NA))
data$PFS_STATUS = ifelse(data$PFS_STATUS == "1:Yes",1,ifelse(data$PFS_STATUS == "0:No",0,NA))
data$RESPONSE = ifelse(data$RESPONSE == "Yes","response",ifelse(data$RESPONSE == "No","nonresponse",NA))


gene_name = unique(data3$Hugo_Symbol)

for(x in gene_name){
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
}

apply(data[,9:742
           ],1,function(x){sum(x == "Mutation")})

write.table(x = data,file = "result_data/dataset15.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

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

write.table(x = WES,file = "result_data/dataset15_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
