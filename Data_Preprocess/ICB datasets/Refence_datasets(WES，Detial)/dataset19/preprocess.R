library(openxlsx)

tmp1 = read.xlsx(xlsxFile = "origin_data/mmc2 (1).xlsx",sheet = 1,startRow = 2,check.names = F)
tmp2 = read.csv(file = "origin_data/bms038_clinical_data.csv",header = T,quote = "",stringsAsFactors = F,check.names = F)

rownames(tmp1) = paste(tmp1$Patient,"_Pre",sep = "")
rownames(tmp2) = paste(tmp2$PatientID,"_Pre",sep = "")



share_id = intersect(rownames(tmp1),rownames(tmp2))

tmp1 = tmp1[share_id,]
tmp2 = tmp2[share_id,]

data1 = cbind(tmp1,tmp2)


data = data.frame("OS_TIME" = data1$OS,
                  "OS_STATUS" = data1$OS_SOR,
                  "PFS_TIME" = data1$PFS,
                  "PFS_STATUS" = data1$PFS_SOR,
                  "RECIST" = data1$BOR,
                  "TMB" = data1$Mutation.Load,
                  "THERAPY" = data1$Cohort,
                  "origin_therapy" = data1$Cohort
                  )



table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
pie(table(data$THERAPY))
table(data$TMB)


data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$OS_STATUS = abs(data$OS_STATUS - 1)  #该数据集将事件定义为0，所以这里将0,1反过来
data$PFS_STATUS = abs(data$PFS_STATUS - 1)
data$OS_TIME = data$OS_TIME/30
data$PFS_TIME = data$PFS_TIME/30
data$THERAPY = ifelse(data$THERAPY %in% c("NIV3-NAIVE"),"anti-PD1/anti-PDL1",NA)

rownames(data) = rownames(data1)


wes_sample = read.csv(file = "origin_data/genomic_data_per_case.csv",header = T,quote = "",stringsAsFactors = F,check.names = F)
wes_sample$id = paste(wes_sample$Patient,"_Pre",sep = "")

WES_id = wes_sample$id[wes_sample$`Pre-treatment Exome` == 1] #73个患者不都测了WES或RNA，所以要删除一些患者的临床信息

data = data[WES_id,]


data3 = read.csv(file = "origin_data/pre_therapy_nonsynonmous_mutations.csv",header = T,quote = "",stringsAsFactors = F,check.names = F)
data3$Patient = paste(data3$Patient,"_Pre",sep = "")

length(unique(data3$Patient))

length(intersect(unique(data3$Patient),rownames(data)))



TMB = read.xlsx("origin_data/mmc1.xlsx",sheet = 1,startRow = 4,check.names = F)
rownames(TMB) = paste(TMB$Patient,"_Pre",sep = "")

data$TMB = (data$TMB*10e5)/(TMB[rownames(data),]$ON_TARGET_BASES/TMB[rownames(data),]$MEAN_TARGET_COVERAGE)

gene_name = unique(data3$`Hugo Symbol`)
for(x in gene_name){
  tmp = unique(data3[data3$`Hugo Symbol` %in% x,"Patient"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

write.table(x = data,file = "result_data/dataset20.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("F:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])


data2 = data[data$origin_therapy == "NIV3-NAIVE",]
# data2 = data
newdata = data2[,1:8]
newdata$Complement = ifelse(rowSums(data2[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.35),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

# newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")


cox = coxph(Surv(OS_TIME,OS_STATUS)~Complement,data = newdata,)

ggforest(cox)
