data1 = read.table(file = "origin_data/data_clinical_patient.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
data2 = read.table(file = "origin_data/data_clinical_sample.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
data3 = read.table(file = "origin_data/data_mutations_mskcc.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
data4 = read.table(file = "origin_data/data_RNA_Seq_expression_median.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")

sum(data1$PATIENT_ID %in% data2$PATIENT_ID)
sum(data1$PATIENT_ID %in% unique(data3$Tumor_Sample_Barcode))

library(org.Hs.eg.db)
convert_gene = clusterProfiler::bitr(data4$Entrez_Gene_Id,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db)
data4 = data4[data4$Entrez_Gene_Id %in% convert_gene$ENTREZID,]
rownames(convert_gene) = convert_gene$ENTREZID
data4$Entrez_Gene_Id = convert_gene[as.character(data4$Entrez_Gene_Id),"SYMBOL"]

tmp_data4 = data4[,2:22]
tmp_data4 = as.matrix(tmp_data4)
rownames(tmp_data4) = data4$Entrez_Gene_Id
tmp_data4 = limma::avereps(tmp_data4)
data4 = as.data.frame(tmp_data4)


rownames(data1) = data1$PATIENT_ID
rownames(data2) = data2$PATIENT_ID
nn = intersect(data1$PATIENT_ID,data2$PATIENT_ID)
data1 = data1[nn,]
data2 = data2[nn,]
data = data.frame("OS_TIME" = data1$OS_MONTHS,
                  "OS_STATUS" = data1$OS_STATUS,
                  "RESPONSE" = data1$DURABLE_CLINICAL_BENEFIT,
                  "TMB" = data2$TMB_NONSYNONYMOUS,
                  "THERAPY" = "anti-CTLA4",
                  "origin_therapy" = data1$TREATMENT
)

data$RESPONSE = ifelse(data$RESPONSE == "LB","response",ifelse(data$RESPONSE == "NB","nonresponse",NA))
data$OS_STATUS = ifelse(data$OS_STATUS == "0:LIVING",0,ifelse(data$OS_STATUS == "1:DECEASED",1,NA))
rownames(data) = data1$PATIENT_ID

gene_name = unique(data3$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data3[data3$Hugo_Symbol %in% x,"Tumor_Sample_Barcode"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
}


share_id = intersect(rownames(data),colnames(data4))
data = data[share_id,]
data4 = data4[,share_id]

RNA = data4

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


datasets_rna_wes = list("dataset10"=list("WES" = data,"RNA" = RNA,"immune_infiltration" = ssgsea_score,"immune_pathway" = immune_pathway))

saveRDS(datasets_rna_wes,"result_data/datasets10_rna_wes.rds")
