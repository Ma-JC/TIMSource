library(openxlsx)
data1 = read.xlsx("origin_data/aan5951_TableS1.xlsx",sheet = 2,startRow = 2)
data2 = read.xlsx("origin_data/aan5951_TableS1.xlsx",sheet = 3,startRow = 2)
data3 = read.xlsx("origin_data/aan5951_TableS1.xlsx",sheet = 1,startRow = 2)
data4 = read.xlsx("origin_data/aan5951_TableS8.xlsx",startRow = 2)

data5 = read.table("origin_data/ccrcc_dfci_2019/data_clinical_patient.txt",header = T,sep = "\t",row.names = 1)
data6 = read.table("origin_data/ccrcc_dfci_2019/data_clinical_sample.txt",header = T,sep = "\t",row.names = 2)

library(stringi)
tmp_ensembl = stri_sub(data4$gene_id,1,15)
data4$gene_id = tmp_ensembl

gtf_data = rtracklayer::import('origin_data/Homo_sapiens.GRCh37.87.gtf')
gtf_data = as.data.frame(gtf_data)
gtf_data2 = gtf_data[,c("type","gene_id","gene_name")]
gtf_data2 = gtf_data2[!duplicated( gtf_data2 ),]
gtf_data2 = gtf_data2[gtf_data2$type == "gene",]
rownames(gtf_data2) = gtf_data2$gene_id
share_gene = intersect(tmp_ensembl,gtf_data2$gene_id)

data4 = data4[data4$gene_id %in% share_gene,]
rownames(data4) = data4$gene_id
data4 = data4[,-1]
data4 = data4[share_gene,]
gtf_data2 = gtf_data2[share_gene,]

data4 = as.matrix(data4)
rownames(data4) = gtf_data2$gene_name
data4 = limma::avereps(data4)
data4 = as.data.frame(data4)

data3 = data3[data3$exclusion_reason == 0,]

sum(data1$patient_id %in% data3$patient_id)

data5 = data5[data1$patient_id,]
data6 = data6[data1$patient_id,]


data = data.frame("OS_TIME" = data1$os_days,
                  "OS_STATUS" = data1$os_censor,
                  "PFS_TIME" = data1$pfs_days,
                  "PFS_STATUS" = data1$pfs_censor,
                  "RECIST" = data1$best_RECIST,
                  "RESPONSE" = data1$response_category,
                  "TMB" = data6$TMB_NONSYNONYMOUS,
                  "THERAPY" = "anti-PD1/PDL1",
                  "origin_therapy" = data5$DRUG
)

rownames(data) = data1$patient_id

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)

data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE = ifelse(data$RESPONSE %in% c("clinical benefit"),"response",ifelse(data$RESPONSE %in% c("intermediate benefit","no clinical benefit"),"nonresponse",NA))
data$OS_TIME = data$OS_TIME/30
data$PFS_TIME = data$PFS_TIME/30

data2$PatientID = do.call(rbind,strsplit(x = data2$Tumor_Sample_Barcode,"-"))[,1]
sum(unique(data2$PatientID) %in% rownames(data))


gene_name = unique(data2$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data2[data2$Hugo_Symbol %in% x,"PatientID"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

colnames(data4)
rownames(data)
###转录组数据中超过一半患者是来自与验证集的，而该文章并未提供验证集的突变数据(只提供发现集的)，所以最后匹配WES和RNA-seq的数据少之又少

data4 = data4[,grepl('RCC',colnames(data4))]
colnames(data4) = do.call(rbind,strsplit(colnames(data4)[grepl('RCC',colnames(data4))],"_T"))[,1]

data = data[colnames(data4),] #匹配了~ ~

RNA = data4

RNA = log2(RNA + 1)
#################################################

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


datasets_rna_wes = list("dataset14"=list("WES" = data,"RNA" = RNA,"immune_infiltration" = ssgsea_score,"immune_pathway" = immune_pathway))

saveRDS(datasets_rna_wes,"result_data/datasets14_rna_wes.rds")
