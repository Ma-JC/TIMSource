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
                  "origin_therapy" = data1$Immunotherapy_regimen,
                  "Age" = data1$Age,
                  "Gender" = data1$Sex
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

################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("G:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])

newdata = data[,1:10]
newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.9),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")


cox = coxph(Surv(OS_TIME,OS_STATUS)~Complement+Age+Gender+BTMB+THERAPY,data = newdata,)

ggforest(cox)
