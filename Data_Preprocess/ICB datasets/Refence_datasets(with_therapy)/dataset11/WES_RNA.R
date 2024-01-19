library(openxlsx)
data1 = read.xlsx(xlsxFile = "origin_data/41591_2020_839_MOESM2_ESM.xlsx",sheet = 1,startRow = 2)
data2 = read.xlsx(xlsxFile = "origin_data/41591_2020_839_MOESM2_ESM.xlsx",sheet = 2,startRow = 2)
data3 = read.xlsx(xlsxFile = "origin_data/41591_2020_839_MOESM2_ESM.xlsx",sheet = 5,startRow = 2) #此处的标准化RNA-seq数据并非FPKM，而是CPM的标准化

data1 = data1[data1$Arm == "NIVOLUMAB",]

#跟dataset13不同，这里同样是CNSR(censor),但是这里的cnsr = 1 表示的是事件(即死亡或进展)，也因此提醒我以后对于类似的信息(如生存状态)，如果只给了数字来代表其分类，应该确认不同数字对应的分类(文章也可能未明确提出，此时应该通过复现文章的图来确认)
data = data.frame("OS_TIME" = data1$OS,
                  "OS_STATUS" = data1$OS_CNSR,
                  "PFS_TIME" = data1$PFS,
                  "PFS_STATUS" = data1$PFS_CNSR,
                  "RECIST" = data1$ORR,
                  "RESPONSE" = data1$Benefit,
                  "TMB" = data1$TMB_Counts,
                  "THERAPY" = "anti-PD1/PDL1",
                  "origin_therapy" = data1$Arm
)

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$TMB)

data$RECIST = ifelse(data$RECIST %in% c("CR","PR","CRPR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE = ifelse(data$RESPONSE == "CB","response",ifelse(data$RESPONSE %in% c("ICB","NCB"),"nonresponse",NA))

data = data[!is.na(data1$MAF_Tumor_ID),]
rownames(data) = data1[!is.na(data1$MAF_Tumor_ID),"MAF_Tumor_ID"]

sum(unique(data2$Tumor_Sample_Barcode) %in% rownames(data))

data2 = data2[data2$Tumor_Sample_Barcode %in% rownames(data),] #确保了~

gene_name = unique(data2$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data2[data2$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

tmp_RNA = data3[,2:312]
tmp_RNA = as.matrix(tmp_RNA)
rownames(tmp_RNA) = data3$gene_name
library(limma)
tmp_RNA = avereps(tmp_RNA)
tmp_RNA = as.data.frame(tmp_RNA)
sum(duplicated(colnames(tmp_RNA))) #样本没有重复

tmp = data1[data1$MAF_Tumor_ID %in% rownames(data) & !is.na(data1$RNA_ID),] #通过meta信息来保留既有WES又有RNA的样本数据

WES = data[tmp$MAF_Tumor_ID,] #对齐
RNA = tmp_RNA[,tmp$RNA_ID]
colnames(RNA) = tmp$MAF_Tumor_ID #将RNA的ID换成对应的DNA的ID

###############################################

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

datasets_rna_wes = list("dataset12"=list("WES" = WES,"RNA" = RNA,"immune_infiltration" = ssgsea_score,"immune_pathway" = immune_pathway))

saveRDS(datasets_rna_wes,"result_data/datasets12_rna_wes.rds")
