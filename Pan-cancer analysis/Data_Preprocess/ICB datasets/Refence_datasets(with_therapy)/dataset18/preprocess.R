library(openxlsx)
data1 = read.xlsx(xlsxFile = "origin_data/43018_2019_8_MOESM2_ESM (1).xlsx",sheet = 1,startRow = 2)
data3 = read.xlsx(xlsxFile = "origin_data/43018_2019_8_MOESM2_ESM (1).xlsx",sheet = 3,startRow = 2)


data = data.frame("OS_TIME" = data1$OS,
                  "OS_STATUS" = data1$`OS.censor.(0=censored,.1=DOD)`,
                  "PFS_TIME" = data1$PFS,"PFS_STATUS" = data1$`PFS.censor.(0=censored,.1=progressed)`,
                  "RESPONSE" = data1$Clinical.Benefit,
                  "THERAPY" = data1$Treatment,
                  "origin_therapy" = data1$Treatment
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

write.table(x = data,file = "result_data/dataset19.txt",sep = "\t",quote = FALSE,row.names = T,col.names = T)

######################################

tmp_word = strsplit(data3$`Nucleotide.Position.(hg19,.Genomic)`,"_|-")

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

write.table(x = WES,file = "result_data/dataset19_mutation.txt",quote = F,sep = "\t",col.names = T,row.names = F)
