library(openxlsx)
data1 = read.xlsx(xlsxFile = "origin_data/1-s2.0-S266637912030183X-mmc2(1).xlsx",sheet = 1,startRow = 2)
data3 = read.xlsx(xlsxFile = "origin_data/1-s2.0-S266637912030183X-mmc2(1).xlsx",sheet = 2,startRow = 2)
data4 = read.xlsx(xlsxFile = "origin_data/1-s2.0-S266637912030183X-mmc2(1).xlsx",sheet = 15,startRow = 2)

data1 = data1[data1$Patient.ID %in% data4$Patient.ID,]

data1$Patient.ID = paste("P",data1$Patient.ID,sep = "")
data3$Patient.ID = paste("P",data3$Patient.ID,sep = "")
#跟dataset13不同，这里同样是CNSR(censor),但是这里的cnsr = 1 表示的是事件(即死亡或进展)，也因此提醒我以后对于类似的信息(如生存状态)，如果只给了数字来代表其分类，应该确认不同数字对应的分类(文章也可能未明确提出，此时应该通过复现文章的图来确认)
data = data.frame("OS_TIME" = data1$OS,
                  "OS_STATUS" = data1$`OS.censor.(1=censored,.0=DOD)`,
                  "PFS_TIME" = data1$`PFS.(months)`,
                  "PFS_STATUS" = data1$`PFS.censor.(1=censored,.0=progressed)`,
                  "RECIST" = data1$BOR,
                  "THERAPY" = data1$Treatment.Group,
                  "origin_therapy" = data1$Treatment.Group
                  )

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
pie(table(data$THERAPY))
# table(data$TMB)

data$RECIST = ifelse(data$RECIST %in% c("Complete Response","Partial Response"),"CR/PR",ifelse(data$RECIST %in% c("Disease Progression","Stable Disease"),"PD/SD",NA))
data$THERAPY = ifelse(data$THERAPY %in% c("NIV1+IPI3 P2","NIV1+IPI3 P3","NIV1+IPI3 P4"),"anti-PD1/PDL1 + anti-CTLA4",
                      ifelse(data$THERAPY %in% c("NIV3-Q2W P3","NIV3-Q2W P4"),"anti-PD1/PDL1",
                             ifelse(data$THERAPY %in% "IPI3-Q3W P3","anti-CTLA1",NA)))

# data$RESPONSE = ifelse(data$RESPONSE == "CB","response",ifelse(data$RESPONSE %in% c("ICB","NCB"),"nonresponse",NA))

data$OS_STATUS = abs(data$OS_STATUS - 1)  #该数据集将事件定义为0，所以这里将0,1反过来
data$PFS_STATUS = abs(data$PFS_STATUS - 1)


rownames(data) = data1$Patient.ID

sum(unique(data3$Patient.ID) %in% rownames(data))

# data = data[rownames(data) %in% data3$Patient.ID,] #此处不是panel测序,而且根据文章来看就是46人的WES测序

gene_name = unique(data3$Gene)
for(x in gene_name){
  tmp = unique(data3[data3$Gene %in% x,"Patient.ID"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

write.table(x = data,file = "result_data/dataset18.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

######################################

tmp_word = strsplit(data3$`Genomic.Change.(Hg19)`,"_|-")

for(i in 1:length(tmp_word)){
  if(length(tmp_word[[i]]) < 5){
    tmp_word[[i]] = c(tmp_word[[i]][1:3],"",tmp_word[[i]][4])
  }
}

tmp_word = do.call(rbind,tmp_word)
tmp_word = as.data.frame(tmp_word)

WES = data.frame(
  "ID" = data3$Patient.ID,
  "Hugo_Symbol" = data3$Gene,
  "Chrmosome" = tmp_word$V1,
  "Start_Position" = tmp_word$V2,
  "End_Position" = tmp_word$V3,
  "Variant_Classification" = data3$Consequence,
  "Variant_Type" = data3$Mutation.Type,
  "Reference" = tmp_word$V4,
  "Alter" = tmp_word$V5)

write.table(x = WES,file = "result_data/dataset18_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
