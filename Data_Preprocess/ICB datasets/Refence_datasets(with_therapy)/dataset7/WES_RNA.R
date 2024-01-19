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


raw_RNA = read.xlsx("origin_data/GSE78220_PatientFPKM.xlsx",check.names = F)

rownames(raw_RNA) = raw_RNA$Gene
raw_RNA = raw_RNA[,-1]
colnames(raw_RNA) = do.call(rbind,strsplit(colnames(raw_RNA),"\\."))[,1]
raw_RNA$Pt27 = rowMeans(raw_RNA[,c("Pt27A","Pt27B")]) #必须清楚没有匹配到的样本是什么情况，不要随意就丢弃而不顾原因
raw_RNA$Pt27A = NULL
raw_RNA$Pt27B = NULL

sum(data1$Patient.ID %in% data2$Sample)
sum(unique(data2$Sample) %in% data1$Patient.ID)

# data1 = data1[data1$Patient.ID %in% data2$Sample,] #确保了临床和突变检查的患者是对应的
# data2 = data2[data2$Sample %in% data1$Patient.ID,]

sum(colnames(raw_RNA) %in% data1$Patient.ID)
sum(data1$Patient.ID %in% colnames(raw_RNA))
sum(unique(data2$Sample) %in% colnames(raw_RNA)) #27个患者三种信息全匹配

data1 = data1[data1$Patient.ID %in% colnames(raw_RNA),] 
data2 = data2[data2$Sample %in% colnames(raw_RNA),]

data3 = data3[rownames(data1),]
data4 = data4[rownames(data1),]

data = data.frame("PatientID" = data1$Patient.ID,
                  "OS_TIME" = data1$Overall.Survival,
                  "OS_STATUS" = data1$Vital.Status,
                  "RECIST" = data1$irRECIST,
                  "TMB" = data4$TMB_NONSYNONYMOUS,
                  "THERAPY" = data3$TREATMENT,
                  "origin_therapy" = data3$TREATMENT
)
# data = data[!rowSums(is.na(data)) == 4,]
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

share_id = intersect(colnames(raw_RNA),rownames(data)) #已预先对齐了~ ~
RNA = raw_RNA[,share_id]
data = data[share_id,]

RNA = log2(RNA + 1)
#########################
immune = list(
  "6-gene IFN signature" = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG"),
  "18-gene IFN signature" = c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB"),
  "Gene expression profile" = c("CXCR6", "TIGIT", "CD27", "CD274", "PDCD1LG2", "LAG3", "NKG7", "PSMB10", "CMKLR1", "CD8A", "IDO1", "CCL5", "CXCL9", "HLA-DQA1", "CD276", "HLA-DRB1", "STAT1", "HLA-E"),
  "Cytolytic activity" = c("GZMA", "PRF1"),
  "13 T-cell signature" = c("CD8A", "CCL2", "CCL3", "CCL4", "CXCL9", "CXCL10", "ICOS", "GZMK", "IRF1","HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB"),
  "Effective T cell score" = c("CD8A", "EOMES", "PRF1", "IFNG", "CD274"),
  "Immune checkpoint expression" = c("CD274", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "PDCD1LG2", "TIGIT"),
  "TLS" = c("CCL21", "CCL19", "CXCL13", "CXCL11", "CCL8", "CXCL10", "CXCL9", "CCL2", "CCL3", "CCL18", "CCL5")
)


immune_pathway = vector()
for( i in names(immune)){
  tmp_genes = immune[[i]][immune[[i]] %in% rownames(RNA)]
  tmp_data = RNA[tmp_genes,]
  
  immune_pathway = rbind(immune_pathway,colMeans(tmp_data))
  
} 
rownames(immune_pathway) = names(immune)


immune_infiltration_data = read.table(file = "origin_data/immune_infiltration.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")
celltype_geneset = list()
for( i in immune_infiltration_data$`Cell type`){
  celltype_geneset[[i]] = immune_infiltration_data$Metagene[ immune_infiltration_data$`Cell type` == i]
}

library(GSVA)
ssgsea_score = gsva(as.matrix(RNA), celltype_geneset, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE,parallel.sz=12)


datasets_rna_wes = list("dataset8"=list("WES" = data,"RNA" = RNA,"immune_infiltration" = ssgsea_score,"immune_pathway" = immune_pathway))

saveRDS(datasets_rna_wes,"result_data/datasets8_rna_wes.rds")
