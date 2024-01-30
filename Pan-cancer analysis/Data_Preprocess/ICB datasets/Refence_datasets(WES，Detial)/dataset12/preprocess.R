library(openxlsx)
data1 = read.xlsx(xlsxFile = "origin_data/41591_2020_1044_MOESM3_ESM.xlsx",sheet = 11,startRow = 2)
data2 = read.xlsx(xlsxFile = "origin_data/41591_2020_1044_MOESM3_ESM.xlsx",sheet = 18,startRow = 2)
data3 = read.xlsx(xlsxFile = "origin_data/41591_2020_1044_MOESM3_ESM.xlsx",sheet = 21,startRow = 2)

data3 = data3[!duplicated(data3$HUGO),]
rownames(data3) = data3$HUGO
data3 = data3[,-1]

pie(table(data1$TRT01P))
data1 = data1[data1$TRT01P == "Avelumab+Axitinib",]
data2 = data2[data2$ID %in% data1$ID,]
data2 = data2[!is.na(data2$nonsyn_var_MB),] # 去除那些没有测WES的患者
data1 = data1[data1$ID %in% data2$ID,]
data3 = data3[,data1$ID] # 突变数据也对齐患者了

share_id = intersect(data1$ID,data2$ID)
rownames(data1) = data1$ID
rownames(data2) = data2$ID
data1 = data1[share_id,]
data2 = data2[share_id,]

data = data.frame("PFS_TIME" = data1$PFS_P,
                  "PFS_STATUS" = data1$PFS_P_CNSR,
                  "TMB" = data2$nonsyn_var_MB,
                  "THERAPY" = "Anti-PD1/PDL1+Axitinib",
                  "origin_therapy" = data1$TRT01P,
                  "Age" = data1$AGE,
                  "Gender" = data1$SEX
                  )
rownames(data) = data1$ID

table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$TMB)

###这里需要注意的是cnsr是censor,此处cnsr = 0才表示事件(进展)，否则表示censor(即在研究时间内未发生事件，后续是否发生就未知了)
m = data$PFS_STATUS == 1
data$PFS_STATUS[m] = 0
data$PFS_STATUS[!m] = 1

library(reshape2)
data3$Gene = rownames(data3)
data3 = melt(data3)
data3 = data3[data3$value == 1,]
data3$variable = as.character(data3$variable)

gene_name = unique(data3$Gene)
for(x in gene_name){
  tmp = unique(data3[data3$Gene %in% x,"variable"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

write.table(x = data,file = "result_data/dataset13.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("G:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])

newdata = data[,1:7]
newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.8),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")


cox = coxph(Surv(PFS_TIME,PFS_STATUS)~Complement+BAge+Gender+BTMB,data = newdata)

ggforest(cox)
