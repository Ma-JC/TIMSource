############ 该文章数据可从cbioportal和文章附录获得，各信息能对应，但是TMB和治疗类型只在cbioportal中 ##############

data1 = openxlsx::read.xlsx("./origin_data/clinical.xlsx")
data2 = read.table("./origin_data/mutation.txt",sep = "\t",header = T)
TMB = read.csv("./origin_data/dataset23_TMB.txt",sep = "\t",row.names = 1)

sum(data1$OM_ID %in% data2$OM_ID)

sum(unique(data2$OM_ID) %in% data1$OM_ID)

length(unique(data2$OM_ID)) #测序的患者不止有35个，多出来的患者需要剔除

share_names = intersect(data1$OM_ID,data2$OM_ID)

setdiff(data1$OM_ID,share_names)

data1 = data1[ data1$OM_ID %in% share_names,]
data2 = data2[ data2$OM_ID %in% share_names,]

data = data.frame("OS_TIME" = data1$OS_TIME,
                  "OS_STATUS" = data1$OS_STATUS,
                  "PFS_TIME" = data1$PFS_TIME,
                  "PFS_STATUS" = data1$PFS_STATUS,
                  "RECIST" = data1$RECIST,
                  # "RESPONSE" = data1$group,
                  # "TMB" = data4$TMB_NONSYNONYMOUS,
                  "THERAPY" = "anti-PD1/PDL1",
                  "origin_therapy" = "anti-PD1/PDL1",
                  "Age" = data1$age,
                  "Gender" = data1$gender
                  )
rownames(data) = data1$OM_ID

data$TMB = TMB[rownames(data),1]

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$THERAPY)

data$RECIST = ifelse(data$RECIST %in% c("CR","PR","PR "),"CR/PR",ifelse(data$RECIST %in% c("PD","PD ","SD","SD "),"PD/SD",NA))

gene_name = unique(data2$Gene)
for(x in gene_name){
  tmp = unique(data2[data2$Gene %in% x,"OM_ID"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}
write.table(x = data,file = "./result_data/dataset23.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("G:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])

newdata = data[,1:11]
newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.75),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")


cox = coxph(Surv(PFS_TIME,PFS_STATUS)~Complement+BAge+Gender+BTMB,data = newdata,)

ggforest(cox)

