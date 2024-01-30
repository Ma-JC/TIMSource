############ 该文章数据可从cbioportal和文章附录获得，各信息能对应，但是TMB和治疗类型只在cbioportal中 ##############

data1 = openxlsx::read.xlsx("./origin_data/clinical.xlsx")
data2 = openxlsx::read.xlsx("./origin_data/mutation.xlsx")


sum(data1$`3DID` %in% data2$ID)
sum(unique(data2$ID) %in% data1$`3DID`)
length(unique(data2$ID))

unique(data2$ID[!data2$ID %in% data1$`3DID`]) #有一例患者有测序数据，都是没有临床资料
data2 = data2[ data2$ID %in% data1$`3DID`,]
length(unique(data2$ID))

data1$`CTLA-4`[ is.na(data1$`CTLA-4`) ] = ""
data = data.frame("OS_TIME" = data1$OS_time_20210421,
                  "OS_STATUS" = data1$OS_status,
                  "PFS_TIME" = data1$PFS_time_20210421,
                  "PFS_STATUS" = data1$PFS_status,
                  "RECIST" = data1$Best_ORR,
                  # "RESPONSE" = data1$group,
                  # "TMB" = data4$TMB_NONSYNONYMOUS,
                  "THERAPY" = data1$ICI_strategy,
                  "origin_therapy" = paste(data1$`PD1/-L1`,data1$`CTLA-4`)
                  )
rownames(data) = paste("P",data1$`3DID`,sep = "_")
data2$ID = paste("P",data2$ID,sep = "_")

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$THERAPY)

data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$THERAPY = ifelse(data$THERAPY == "Combined therapy","anti-PD1/PDL1 + anti-CTLA4",
                      ifelse(data$THERAPY == "Monotherapy","anti-PD1/PDL1",NA))



gene_name = unique(data2$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data2[data2$Hugo_Symbol %in% x,"ID"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}
write.table(x = data,file = "./result_data/dataset22.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


#########################################

WES = data.frame(
  "ID" = data2$ID,
  "Hugo_Symbol" = data2$Hugo_Symbol,
  "Chrmosome" = data2$Chrmosome,
  "Start_Position" = data2$Start_Position,
  "End_Position" = data2$End_Position,
  "Variant_Classification" = data2$Variant_Classification,
  # "Variant_Type" = data2$Variant_Type,
  "Reference" = data2$Reference,
  "Alteration" = data2$Alteration
  # "Tumor_Seq_Allele2" = data2$Tumor_Seq_Allele2
  )

WES$Variant_Classification[ is.na(WES$Variant_Classification) ] = "Other"
write.table(x = WES,file = "./result_data/dataset22_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
