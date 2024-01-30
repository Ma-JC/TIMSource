library(openxlsx)
data1 = read.xlsx(xlsxFile = "origin_data/43018_2019_8_MOESM2_ESM (1).xlsx",sheet = 1,startRow = 2)
data3 = read.xlsx(xlsxFile = "origin_data/43018_2019_8_MOESM2_ESM (1).xlsx",sheet = 3,startRow = 2)


data = data.frame("OS_TIME" = data1$OS,
                  "OS_STATUS" = data1$`OS.censor.(0=censored,.1=DOD)`,
                  "PFS_TIME" = data1$PFS,
                  "PFS_STATUS" = data1$`PFS.censor.(0=censored,.1=progressed)`,
                  "RESPONSE" = data1$Clinical.Benefit,
                  "THERAPY" = data1$Treatment,
                  "origin_therapy" = data1$Treatment,
                  "Age" = data1$Age.at.ICB.initiation,
                  "Gender" = data1$Gender,
                  "Smoke" = data1$Smoking.Status
                  )

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RESPONSE)
pie(table(data$THERAPY))
# table(data$RESPONSE)
# table(data$TMB)

# data$RECIST = ifelse(data$RECIST %in% c("Complete Response","Partial Response"),"CR/PR",ifelse(data$RECIST %in% c("Disease Progression","Stable Disease"),"PD/SD",NA))
data$RESPONSE = ifelse(data$RESPONSE == "DCB","response",ifelse(data$RESPONSE %in% c("NDB"),"nonresponse",NA))
data$THERAPY = ifelse(data$THERAPY %in% c("Anti-PD1 (nivolumab)","Anti-PD1 (pembrolizumab)"),"anti-PD1/PDL1",
                      ifelse(data$THERAPY %in% c("Dual ICB (anti-PD1+anti-CTLA4)"),"anti-PD1/PDL1 + anti-CTLA4",NA) )

rownames(data) = data1$Patient.ID

sum(unique(data3$Patient.ID) %in% rownames(data))

# data = data[rownames(data) %in% data3$Patient.ID,] 

gene_name = unique(data3$Gene)
for(x in gene_name){
  tmp = unique(data3[data3$Gene %in% x,"Patient.ID"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

data = cbind(data[,1:10],"TMB" = rowSums(data == "Mutation",na.rm = T),data[,11:8379])

TMB = read.xlsx("origin_data/43018_2019_8_MOESM2_ESM (1).xlsx",sheet = 2,startRow = 2)
data$TMB = (data$TMB*10e5)/(TMB$Sequenced.Bases.Mapped.to.Genome * TMB$Percent.Mapped.to.Target.Regions / TMB$`Average.High.Quality.Total.Coverage.(Tumor)`)

write.table(x = data,file = "result_data/dataset19.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("F:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])

newdata = data[,1:11]
newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.75),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")


cox = coxph(Surv(OS_TIME,OS_STATUS)~Complement+Age+BTMB,data = newdata,)

ggforest(cox)
