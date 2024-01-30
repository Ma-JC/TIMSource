library(openxlsx)
data1 = read.xlsx("origin_data/MSK-GI_JP_PUCH_clinical_info_with_GIPS.xlsx",sheet = 3)
data2 = read.table("origin_data/data_mutations_PUCH.txt",header = T,sep = "\t",quote = "",check.names = F,stringsAsFactors = F)

sum(data1$Sample_ID %in% data2$Tumor_Sample_Barcode)
sum(unique(data2$Tumor_Sample_Barcode) %in% data1$Sample_ID)

data = data.frame("PatientID" = data1$Sample_ID,
                  "OS_TIME" = data1$OS_time,
                  "OS_STATUS" = data1$OS_status,
                  "PFS_TIME" = data1$PFS_time,
                  "PFS_STATUS" = data1$PFS_status,
                  "RESPONSE" = data1$DCB,
                  "TMB" = data1$TMB_Score,
                  "THERAPY" = data1$Immunotherapy_regimen,
                  "origin_therapy" = data1$Immunotherapy_regimen
                  )
rownames(data) = data$PatientID
data = data[,-1]

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RESPONSE)
table(data$TMB)
pie(table(data$THERAPY))
sum(is.na(data))

data$RESPONSE = ifelse(data$RESPONSE == "DCB","response",ifelse(data$RESPONSE == "NDB","nonresponse",NA))
data$THERAPY = ifelse(data$THERAPY %in% c("PD1","PDL1"),"anti-PD1/PDL1",
                      ifelse(data$THERAPY %in% c("PD1+CTLA4","PDL1+CTLA4"),"anti-PD1/PDL1 + anti-CTLA4",NA))

gene_name = unique(data2$Gene.refGene)
for(x in gene_name){
  tmp = unique(data2[data2$Gene.refGene %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}


write.table(x = data,file = "result_data/dataset9.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

######################################
WES = data.frame(
  "ID" = data2$Tumor_Sample_Barcode,
  "Hugo_Symbol" = data2$Gene.refGene,
  "Chrmosome" = data2$Chr,
  "Start_Position" = data2$Start,
  "End_Position" = data2$End,
  "Variant_Classification" = data2$type_specific,
  "Variant_Type" = data2$type,
  "Reference" = data2$ref,
  "Alteration" = data2$alt)

write.table(x = WES,file = "result_data/dataset9_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
