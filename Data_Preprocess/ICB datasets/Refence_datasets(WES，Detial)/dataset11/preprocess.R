library(openxlsx)
data1 = read.xlsx(xlsxFile = "origin_data/41591_2020_839_MOESM2_ESM.xlsx",sheet = 1,startRow = 2)
data2 = read.xlsx(xlsxFile = "origin_data/41591_2020_839_MOESM2_ESM.xlsx",sheet = 2,startRow = 2)
data3 = read.xlsx(xlsxFile = "origin_data/41591_2020_839_MOESM2_ESM.xlsx",sheet = 5,startRow = 2)

pie(table(data1$Arm))
data1 = data1[data1$Arm == "NIVOLUMAB",] #有一部分患者是anti-mTOR治疗的

#跟dataset13不同，这里同样是CNSR(censor),但是这里的cnsr = 1 表示的是事件(即死亡或进展)，也因此提醒我以后对于类似的信息(如生存状态)，如果只给了数字来代表其分类，应该确认不同数字对应的分类(文章也可能未明确提出，此时应该通过复现文章的图来确认)
data = data.frame("OS_TIME" = data1$OS,
                  "OS_STATUS" = data1$OS_CNSR,
                  "PFS_TIME" = data1$PFS,
                  "PFS_STATUS" = data1$PFS_CNSR,
                  "RECIST" = data1$ORR,
                  "RESPONSE" = data1$Benefit,
                  "TMB" = data1$TMB_Counts,
                  "THERAPY" = "anti-PD1/PDL1",
                  "origin_therapy" = data1$Arm,
                  "Age" = data1$Age,
                  "Gender" = data1$Sex,
                  "Site" = data1$Tumor_Sample_Primary_or_Metastasis
                  )

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$TMB)

data$RECIST = ifelse(data$RECIST %in% c("CR","PR","CRPR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE = ifelse(data$RESPONSE == "CB","response",ifelse(data$RESPONSE %in% c("ICB","NCB"),"nonresponse",NA))
data$Gender = ifelse(data$Gender %in% c("Female","F","FEMALE"),"Female",ifelse(data$Gender %in% c("Male","M","MALE"),"Male",NA))

data = data[!is.na(data1$MAF_Tumor_ID),]
rownames(data) = data1[!is.na(data1$MAF_Tumor_ID),"MAF_Tumor_ID"]

sum(unique(data2$Tumor_Sample_Barcode) %in% rownames(data))

data2 = data2[data2$Tumor_Sample_Barcode %in% rownames(data),] # 一部分患者过滤了，突变数据也需要过滤

gene_name = unique(data2$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data2[data2$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

write.table(x = data,file = "result_data/dataset12.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("G:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])

newdata = data[,1:12]
newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.8),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")

cox = coxph(Surv(PFS_TIME,PFS_STATUS)~Complement+BAge+Gender+BTMB,data = newdata,)

ggforest(cox)

unique(data$Gender)
