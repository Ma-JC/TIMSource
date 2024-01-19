library(openxlsx)
data1 = read.xlsx("origin_data/Hugo mmc1 (1).xlsx",sheet = 2,startRow = 3)
data2 = read.xlsx("origin_data/Hugo mmc1 (1).xlsx",sheet = 5,startRow = 3)
data3 = read.table("origin_data/mel_ucla_2016/data_clinical_patient.txt",sep = "\t",header = T)
data4 = read.table("origin_data/mel_ucla_2016/data_clinical_sample.txt",sep = "\t",header = T)

share_pn = intersect(data1$Patient.ID,data3$PATIENT_ID)
data1 = data1[data1$Patient.ID %in% share_pn,]
rownames(data1) = data1$Patient.ID
rownames(data3) = data3$PATIENT_ID
rownames(data4) = data4$PATIENT_ID
data1 = data1[share_pn,]
data3 = data3[share_pn,]
data4 = data4[share_pn,]

sum(data1$Patient.ID %in% data2$Sample)
sum(unique(data2$Sample) %in% data1$Patient.ID)

# data1 = data1[data1$Patient.ID %in% data2$Sample,] #确保了临床和突变检查的患者是对应的(排除，这是panel测序)
# data2 = data2[data2$Sample %in% data1$Patient.ID,]


data = data.frame("PatientID" = data1$Patient.ID,
                  "OS_TIME" = data1$Overall.Survival,
                  "OS_STATUS" = data1$Vital.Status,
                  "RECIST" = data1$irRECIST,
                  "TMB" = data4$TMB_NONSYNONYMOUS,
                  "THERAPY" = data3$TREATMENT,
                  "origin_therapy" = data3$TREATMENT,
                  "Age" = data1$Age,
                  "Gender" = data1$Gender,
                  "M" = data1$Disease.Status,
                  "MAPK" = data1$Previous.MAPKi
                  )
# data = data[!rowSums(is.na(data)) == 4,] # 前面已经删除了
rownames(data) = data$PatientID
data = data[,-1]

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$RECIST)

sum(is.na(as.matrix(data)))

data$OS_STATUS = ifelse(data$OS_STATUS == "Dead",1,0)
data$RECIST = ifelse(data$RECIST %in% c("Complete Response","Partial Response"),"CR/PR",ifelse(data$RECIST %in% c("Progressive Disease"),"PD/SD",NA))
data$OS_TIME = data$OS_TIME/30
data$THERAPY = "anti-PD1/PDL1"


gene_name = unique(data2$Gene)
for(x in gene_name){
  tmp = unique(data2[data2$Gene %in% x,"Sample"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}


write.table(x = data,file = "result_data/dataset8.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

################# 多因素 #####################

library("survival")
library("survminer")

pathway_list = readRDS("G:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])


data = data[!is.na(data$OS_TIME),]
newdata = data[,1:10]



newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.9),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")


cox = coxph(Surv(OS_TIME,OS_STATUS)~Complement+BAge+Gender+BTMB+MAPK+M,data = newdata,)

ggforest(cox)
