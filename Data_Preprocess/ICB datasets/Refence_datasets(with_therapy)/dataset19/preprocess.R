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


gene_name = unique(data3$`Hugo Symbol`)
for(x in gene_name){
  tmp = unique(data3[data3$`Hugo Symbol` %in% x,"Patient"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

write.table(x = data,file = "result_data/dataset20.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

#########################################################

WES = data.frame(
  "ID" = data3$Patient,
  "Hugo_Symbol" = data3$`Hugo Symbol`,
  "Chrmosome" = data3$Chromosome,
  "Start_Position" = data3$Start,
  "End_Position" = data3$End,
  "Variant_Classification" = data3$`Variant Classification`,
  "Alter" = data3$HGVS_c)

write.table(x = WES,file = "result_data/dataset20_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
