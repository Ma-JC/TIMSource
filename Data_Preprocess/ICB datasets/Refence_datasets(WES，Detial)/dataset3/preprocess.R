data1 = read.table("origin_data/data_clinical_patient.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data2 = read.table("origin_data/data_clinical_sample.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data3 = read.table("origin_data/data_mutations_mskcc.txt",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F)
data4 = read.table("origin_data/mixed_allen_2018_clinical_data (1).tsv",sep = "\t",header = T,row.names = 3)
TMB = read.table("origin_data/41588_2018_200_MOESM3_ESM.csv",skip = 28,sep = ",",header = T,row.names = 1)

sum(rownames(data1) %in% rownames(data2))
sum(data2$SAMPLE_ID %in% data3$Tumor_Sample_Barcode)
length(unique(data2$SAMPLE_ID))
length(unique(data3$Tumor_Sample_Barcode))



rownames(data1) = data2[rownames(data1),1]

# data1 = data1[rownames(data1) %in% data3$Tumor_Sample_Barcode,]  #确保了临床和突变检查的患者是对应的
data4 = data4[rownames(data1),]

rownames(TMB)[ rownames(TMB) == "BLCA-DFCI_11-104_011"] = "DFCI_11-104_011"
a = na.omit(unlist(sapply(rownames(TMB),function(x){rownames(data1)[grepl(x,rownames(data1))][1]})))
TMB = TMB[names(a),]
rownames(TMB) = a
TMB = TMB[ rownames(data1),]

data = data.frame("OS_TIME" = data1$OS_MONTHS,
                  "OS_STATUS" = data1$OS_STATUS,
                  "PFS_TIME" = data1$PFS_MONTHS,
                  "PFS_STATUS" = data1$PFS_STATUS,
                  "RECIST" = data1$RECIST,
                  "RESPONSE" = data1$ROH_RESPONSE,
                  "TMB" = (data4$Mutation.Count*10e5)/TMB$somatic_mutation_covered_bases_capture, # 该数据没有提供TMB，但是提供了mutation count，跟自己统计的不一样，总体偏小
                  "THERAPY" = data1$DRUG_TYPE,
                  "origin_therapy" = data1$DRUG_TYPE,
                  "Age" = data4$Age.Start.IO,
                  "Gender" = data4$Sex,
                  "Cancer_Type" = data4$Cancer.Type
                  )
rownames(data) = rownames(data1) 

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$THERAPY)
sum(is.na(as.matrix(data)))



sum(is.na(data$OS_TIME))
sum(is.na(data$OS_STATUS))
sum(is.na(data$PFS_TIME))
sum(is.na(data$PFS_STATUS))
sum(is.na(data$RECIST))
sum(is.na(data$RESPONSE))

data$OS_STATUS = ifelse(data$OS_STATUS == "DECEASED",1,0)
data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE = sapply(data$RESPONSE,FUN = function(x){switch(x,"clinical benefit"="response","no clinical benefit"="nonresponse")})
data$THERAPY = ifelse(data$THERAPY == "anti-CTLA-4","anti-CTLA4",
                      ifelse(data$THERAPY == "anti-PD-1/anti-PD-L1","anti-PD1/PDL1",
                             ifelse(data$THERAPY == "anti-CTLA-4 + anti-PD-1/PD-L1","anti-PD1/PDL1 + anti-CTLA4",NA)))

gene_name = unique(data3$Hugo_Symbol)

for(x in gene_name){
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
}



data$TMB
as.numeric(rowSums(data[,13:17936] == "Mutation")) # 自己统计的会偏大

write.table(x = data,file = "./result_data/dataset3.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

################################################

data = read.table(file = "result_data/dataset3.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")
data2 = read.table("origin_data/data_clinical_sample.txt",header = T,row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F,quote = "")

item_id = c("OS_TIME","OS_STATUS","PFS_TIME","PFS_STATUS","TMB","RECIST","RESPONSE","THERAPY","origin_therapy","Age","Gender","Cancer_Type")


for(i in unique(data2$CANCER_TYPE)[2:4]){
  
  tmp = data[data2$SAMPLE_ID[ data2$CANCER_TYPE == i],]
  tmp = tmp[,setdiff(colnames(tmp),setdiff(colnames(tmp)[colSums(tmp == "Mutation") == 0],item_id))]
  write.table(x = tmp,
              file = paste("result_data/dataset3_",paste(strsplit(i," ")[[1]],collapse = "_"),".txt",sep = ""),sep = "\t",quote = FALSE,row.names = T,col.names = T)
  
}



##### 多因素 #####
library("survival")
library("survminer")

pathway_list = readRDS("F:/web/Refence_datasets/Refence_datasets(COX)/all_pathway.rds")

genes = intersect(gene_name,pathway_list[["GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]])

newdata = data[,1:12]


newdata$Complement = ifelse(rowSums(data[,genes] == "Mutation") >0, "Mutation","Wildtype")
newdata$Complement = factor(newdata$Complement,levels = c("Wildtype","Mutation"))

newdata$BTMB = ifelse(newdata$TMB >= quantile(newdata$TMB,probs = 0.75),"TMB high","TMB low") # 二分化TMB
newdata$BTMB = factor(newdata$BTMB,levels = c("TMB low","TMB high"))

newdata$BAge = ifelse(newdata$Age >= 65,"≥65","<65")

newdata = newdata[!newdata$Cancer_Type %in% c("Anal Cancer","Soft Tissue Sarcoma","Small Cell Lung Cancer"),]
newdata$Cancer_Type = factor(newdata$Cancer_Type,levels = c("Melanoma","Non-Small Cell Lung Cancer","Head and Neck Cancer","Bladder Cancer"))
# newdata = newdata[newdata$Cancer_Type %in% c("Non-Small Cell Lung Cancer"),]

cox = coxph(Surv(OS_TIME,OS_STATUS)~Complement+Gender+BTMB+THERAPY,data = newdata) #这里的多因素只能选择瘤种和治疗的其中之一，因为这两者是高度相关的，

ggforest(cox,data = newdata)

table(newdata$Cancer_Type)
