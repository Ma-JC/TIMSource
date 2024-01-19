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
# table(data$RESPONSE)
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

WES_id = wes_sample$id[wes_sample$`Pre-treatment Exome` == 1]

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



data4 = read.csv(file = "origin_data/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv",header = T,row.names = 1,stringsAsFactors = F,check.names = F)

library(org.Hs.eg.db)
convert_gene = clusterProfiler::bitr(rownames(data4),fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db)

data4 = data4[convert_gene$ENTREZID,]
rownames(data4) = convert_gene$SYMBOL

colnames(data4) = paste(do.call(rbind,strsplit(colnames(data4),"_"))[,1],do.call(rbind,strsplit(colnames(data4),"_"))[,2],sep = "_")

WES_RNA_id = wes_sample$id[wes_sample$`Pre-treatment Exome` == 1 & wes_sample$`Pre-treatment RNA-Seq` == 1]

WES_RNA_id2 = intersect(rownames(data),colnames(data4))

share_id = intersect(WES_RNA_id,WES_RNA_id2)


WES = data[share_id,] #对齐
RNA = data4[,share_id]

RNA = log2(RNA + 1)
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

datasets_rna_wes = list("dataset20"=list("WES" = WES,"RNA" = RNA,"immune_infiltration" = ssgsea_score,"immune_pathway" = immune_pathway))

saveRDS(datasets_rna_wes,"result_data/datasets20_rna_wes.rds")
