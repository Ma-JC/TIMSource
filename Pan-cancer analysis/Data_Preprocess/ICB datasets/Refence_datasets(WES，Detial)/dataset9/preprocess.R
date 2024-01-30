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
                  "origin_therapy" = data1$TREATMENT,
                  "Age" = data1$AGE,
                  "Gender" = data1$SEX,
                  "M" = data1$M_STAGE
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


################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("G:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])

newdata = data[,1:9]
newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.8),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")


cox = coxph(Surv(OS_TIME,OS_STATUS)~Complement+BAge+Gender+BTMB+M,data = newdata,)

ggforest(cox)
