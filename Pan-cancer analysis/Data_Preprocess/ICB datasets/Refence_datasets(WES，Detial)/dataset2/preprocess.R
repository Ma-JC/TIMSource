############ 该文章数据可从cbioportal和文章附录获得，各信息能对应，但是TMB和治疗类型只在cbioportal中 ##############

data1 = read.csv("./origin_data/Allen TableS2_Revised.csv",row.names = 1)
data2 = read.csv("./origin_data/Allen TableS1.Mutation_list_all_patients.csv")
data3 = read.table("origin_data/skcm_dfci_2015/data_clinical_patient.txt",sep = "\t",header = T,row.names = 1)
data4 = read.table("origin_data/skcm_dfci_2015/data_clinical_sample.txt",sep = "\t",header = T,row.names = 1)

sum(table(c(rownames(data1),rownames(data3),rownames(data4))) == 3)
sum(table(c(rownames(data1),rownames(data3),rownames(data4))) != 3)

sum(rownames(data1) %in% data2$patient)
sum(unique(data2$patient) %in% rownames(data1)) #确保了临床和突变检查的患者是对应的
sum(rownames(data1) %in% rownames(data3))
sum(rownames(data1) %in% rownames(data4))

data3 = data3[rownames(data1),]
data4 = data4[rownames(data1),]

data = data.frame("OS_TIME" = data1$overall_survival,
                  "OS_STATUS" = data1$dead,
                  "PFS_TIME" = data1$progression_free,
                  "PFS_STATUS" = data1$progression,
                  "RECIST" = data1$RECIST,
                  "RESPONSE" = data1$group,
                  "TMB" = data4$TMB_NONSYNONYMOUS,
                  "THERAPY" = data3$TREATMENT,
                  "origin_therapy" = data3$TREATMENT,
                  "Age" = data1$age_start,
                  "Gender" = data1$gender,
                  "Primary" = data1$primary,
                  "Stage" = data1$stage,
                  "M" = data1$M
                  
                  )
rownames(data) = rownames(data1) 

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$THERAPY)

data$OS_TIME = data$OS_TIME/30
data$PFS_TIME = data$PFS_TIME/30
data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE[data$RESPONSE == "long-survival"] = "nonresponse"
data$THERAPY = "anti-CTLA4"



gene_name = unique(data2$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data2[data2$Hugo_Symbol %in% x,"patient"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}
write.table(x = data,file = "./result_data/dataset2.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("G:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])

newdata = data[,1:14]
newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.75),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")


cox = coxph(Surv(OS_TIME,OS_STATUS)~Complement+BAge+Gender+BTMB,data = newdata,)

ggforest(cox)



