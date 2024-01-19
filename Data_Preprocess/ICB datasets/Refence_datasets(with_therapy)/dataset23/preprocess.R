############ 该文章数据可从cbioportal和文章附录获得，各信息能对应，但是TMB和治疗类型只在cbioportal中 ##############

data1 = openxlsx::read.xlsx("./origin_data/clinical.xlsx")
data2 = read.table("./origin_data/mutation.csv",sep = ",",header = T)

unique(data2$ID)
sum(data1$ID %in% data2$ID)

sum(unique(data2$ID) %in% data1$ID)

length(unique(data2$ID))

data = data.frame("OS_TIME" = data1$OS_time_20210201,
                  "OS_STATUS" = data1$OS_status_20210201,
                  "PFS_TIME" = data1$PFS_time_20210201,
                  "PFS_STATUS" = data1$PFS_status,
                  "RECIST" = data1$`OR_RECIST-20210201`,
                  "RESPONSE" = data1$DCB_20210201,
                  "TMB" = data1$TMB,
                  "THERAPY" = data1$Regimen_of_immune_therapy,
                  "origin_therapy" = data1$Regimen_of_immune_therapy
                  )
rownames(data) = data1$ID

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$THERAPY)

data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE = ifelse(data$RESPONSE %in% 1,"response",ifelse(data$RESPONSE %in% 0,"nonresponse",NA))
data$THERAPY = ifelse(data$THERAPY %in% c("Combined therapy"),"anti-PD1/PDL1 + anti-CTLA4",
                      ifelse(data$THERAPY %in% c("Monotherapy"),"anti-PD1/PDL1",NA)
                      )

data2 = data2[ data2$ID %in% data1$ID,]

gene_name = unique(data2$Hugo_Symbol)
x = gene_name[1]
for(x in gene_name){
  tmp = unique(data2[data2$Hugo_Symbol %in% x,"ID"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}
write.table(x = data,file = "./result_data/dataset24.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


#########################################

WES = data.frame(
  "ID" = data2$ID,
  "Hugo_Symbol" = data2$Hugo_Symbol,
  "Chrmosome" = data2$Chrmosome,
  "Start_Position" = data2$Start_Position,
  "End_Position" = data2$End_Position,
  # "Variant_Classification" = data2$Alternation.type2,
  # "Variant_Type" = data2$Variant_Type,
  "Reference" = data2$Reference,
  "Alteration" = data2$Tumor_Seq_Allele1
  # "Tumor_Seq_Allele2" = data2$Tumor_Seq_Allele2
  )

write.table(x = WES,file = "./result_data/dataset24_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)


