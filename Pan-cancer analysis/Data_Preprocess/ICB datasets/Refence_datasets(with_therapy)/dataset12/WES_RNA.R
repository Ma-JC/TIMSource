library(openxlsx)
data1 = read.xlsx(xlsxFile = "origin_data/41591_2020_1044_MOESM3_ESM.xlsx",sheet = 11,startRow = 2)
data2 = read.xlsx(xlsxFile = "origin_data/41591_2020_1044_MOESM3_ESM.xlsx",sheet = 18,startRow = 2)
data3 = read.xlsx(xlsxFile = "origin_data/41591_2020_1044_MOESM3_ESM.xlsx",sheet = 21,startRow = 2)
data4 = read.xlsx(xlsxFile = "origin_data/41591_2020_1044_MOESM3_ESM.xlsx",sheet = 13,startRow = 2)

data3 = data3[!duplicated(data3$HUGO),]
rownames(data3) = data3$HUGO
data3 = data3[,-1]

tmp_RNA = data4[,2:727]
tmp_RNA = as.matrix(tmp_RNA)
rownames(tmp_RNA) = data4$HUGO
library(limma)
tmp_RNA = avereps(tmp_RNA)
tmp_RNA = as.data.frame(tmp_RNA)

#以下是为了筛选具有免疫治疗和WES数据的患者
data1 = data1[data1$TRT01P == "Avelumab+Axitinib",] 
data2 = data2[data2$ID %in% data1$ID,]
data2 = data2[!is.na(data2$nonsyn_var_MB),]
data1 = data1[data1$ID %in% data2$ID,]
data3 = data3[,data1$ID]

share_id = intersect(colnames(tmp_RNA),colnames(data3)) #对齐
RNA = tmp_RNA[,share_id]
WES = data3[,share_id]
tmp = apply(X = WES,MARGIN = 1,FUN = function(x){
  x = x + 1
  c("Wildtype","Mutation")[x]
})
rownames(tmp) = colnames(WES)
WES = as.data.frame(tmp)


############################################

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

datasets_rna_wes = list("dataset13"=list("WES" = WES,"RNA" = RNA,"immune_infiltration" = ssgsea_score,"immune_pathway" = immune_pathway))

saveRDS(datasets_rna_wes,"result_data/datasets13_rna_wes.rds")
