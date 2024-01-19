library(openxlsx)
data1 = read.xlsx(xlsxFile = "origin_data/Hellmann nsclc_mskcc_2018_clinical_data.xlsx")
data2 = read.xlsx(xlsxFile = "origin_data/mmc3 (1).xlsx",sheet = 2,startRow = 2)
data4 = read.table("origin_data/nsclc_mskcc_2018/data_clinical_sample.txt",sep = "\t",header = T,row.names = 1)

data4 = data4[data1$Sample.ID,]
plot(data1$Nonsynonymous.Mutation.Burden,data4$TMB_NONSYNONYMOUS)

data = data.frame("PFS_TIME" = data1$`Progress.Free.Survival.(Months)`,
                  "PFS_STATUS" = data1$Progression.Free.Status,
                  "RECIST" = data1$Best.Overall.Response,
                  "RESPONSE" = data1$Durable.Clinical.Benefit,
                  "TMB" = data4$TMB_NONSYNONYMOUS, # 这里跟SNVIO的用的突变负荷不一样
                  "THERAPY" = "anti-PD1/PDL1 + anti-CTLA4",
                  "origin_therapy" = NA,
                  "Age" = data1$`Age.(yrs)`,
                  "Gender" = data1$Sex,
                  "Smoke" = data1$Smoking.Status
                  )


rownames(data) = data1$Sample.ID


table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$TMB)
table(data$THERAPY)

data$PFS_STATUS = ifelse(data$PFS_STATUS == "Event",1,0)
data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE = ifelse(data$RESPONSE == "Durable Clinical Benefit","response",ifelse(data$RESPONSE == "No Durable Benefit","nonresponse",NA))

data2$Patient.ID = paste("nsclc_mskcc_2018s",data2$Patient.ID,sep = "")
sum(rownames(data) %in% data2$Patient.ID) #确保了临床和突变检查的患者是对应的
length(unique(data2$Patient.ID))

data3 = vector()
for(i in 1:nrow(data2)){
  tmp = data2[i,]
  tmp_genes = unlist(strsplit(tmp[,"gene_name"],";"))
  for(j in tmp_genes){
    tmp[,"gene_name"] = j
    data3 = rbind(data3,tmp)
  }
  }


gene_name = unique(data3$gene_name)
gene_name = gene_name[!is.na(gene_name)]
for(x in gene_name){
  tmp = unique(data3[data3$gene_name %in% x,"Patient.ID"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}


write.table(x = data,file = "result_data/dataset7.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)


###################################
#之前只是考虑了网站的构建，所以没有用data3，而实际上要用于分析数据，就应该是dataset3
data3$Variant_Classification = "Other"
vars = c("missense_snv","nonsynonymous_snv","nonsynonymous_indel","frameshift")
data3$Variant_Classification[data3$nonsynonymous_snv] = "nonsynonymous_snv"
data3$Variant_Classification[data3$nonsynonymous_indel] = "nonsynonymous_indel"
data3$Variant_Classification[data3$frameshift] = "frameshift"
data3$Variant_Classification[data3$missense_snv] = "missense_snv"

WES = data.frame(
  "ID" = data3$Patient.ID,
  "Hugo_Symbol" = data3$gene_name,
  "Chrmosome" = data3$chr,
  "Position" = data3$start,
  "Variant_Classification" = data3$Variant_Classification,
  "Reference" = data3$ref,
  "Alteration" = data3$alt
)

write.table(x = WES,file = "result_data/dataset7_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)

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


cox = coxph(Surv(PFS_TIME,PFS_STATUS)~Complement+BAge+Gender+BTMB+Smoke,data = newdata,)

ggforest(cox)

