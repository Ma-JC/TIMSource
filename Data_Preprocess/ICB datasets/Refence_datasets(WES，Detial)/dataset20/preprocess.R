data1 = read.table(file = "origin_data/data_clinical_patient.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
data2 = read.table(file = "origin_data/data_clinical_sample.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
data3 = read.table(file = "origin_data/data_mutations_mskcc.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
# data4 = read.table(file = "origin_data/data_RNA_Seq_expression_tpm.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")

sum(data1$PATIENT_ID %in% data2$PATIENT_ID)
sum(data1$PATIENT_ID %in% unique(data3$Tumor_Sample_Barcode))
sum(data2$SAMPLE_ID %in% unique(data3$Tumor_Sample_Barcode))


rownames(data2) = data2$PATIENT_ID
rownames(data1) = data2[data1$PATIENT_ID,"SAMPLE_ID"]
rownames(data2) = data2$SAMPLE_ID

nn = intersect(rownames(data1),rownames(data2))
data1 = data1[nn,]
data2 = data2[nn,]
data = data.frame("OS_TIME" = data1$OS_MONTHS,
                  "OS_STATUS" = data1$OS_STATUS,
                  "PFS_TIME" = data1$PFS_MONTHS,
                  "PFS_STATUS" = data1$PFS_STATUS,
                  "RECIST" = data1$BR,
                  "TMB" = data2$TMB_NONSYNONYMOUS,
                  "THERAPY" = "anti-PD1/PDL1",
                  "origin_therapy" = data1$IO_THERAPY,
                  "BIOPSY_SITE" = data2$BIOPSY_SITE,
                  "IO_THERAPY" = data1$IO_THERAPY,
                  "PRIOR_MAPK_TX" = data1$PRIOR_MAPK_TX,
                  "PRIOR_CTLA4" = data1$PRIOR_CTLA4,
                  "POST_CTLA4" = data1$POST_CTLA4,
                  "POST_MAPK_TX" = data1$POST_MAPK_TX,
                  "POAR_COMBINED_CTLA_PD1" = data1$POAR_COMBINED_CTLA_PD1,
                  "Gender" = data1$SEX,
                  "M" = data1$M_STAGE
                  )

data$RECIST = ifelse(data$RECIST %in% c("Complete Response","Mixed Response","Partial Response"),"CR/PR",ifelse(data$RECIST %in% c("Progressive Disease","Stable Disease"),"PD/SD",NA))
data$OS_STATUS = ifelse(data$OS_STATUS == "0:LIVING",0,ifelse(data$OS_STATUS == "1:DECEASED",1,NA))
data$PFS_STATUS = ifelse(data$PFS_STATUS == "0:CENSORED",0,ifelse(data$PFS_STATUS == "1:PROGRESSION",1,NA))
rownames(data) = data2$SAMPLE_ID

gene_name = unique(data3$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
}

#这里的OS时间(cbioportal下载的)与附录PDF整理出来的不一样，其他的大致相同(其实突变列表也不尽相同)
write.table(x = data,file = "result_data/dataset21.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("G:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])

newdata = data[,1:20]
newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.8),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")


cox = coxph(Surv(PFS_TIME,PFS_STATUS)~Complement+BTMB+Gender+M,data = newdata,)

ggforest(cox)
