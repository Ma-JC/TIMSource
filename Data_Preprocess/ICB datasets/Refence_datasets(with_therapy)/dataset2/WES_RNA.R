data1 = read.csv("./origin_data/Allen TableS2_Revised.csv",row.names = 1)
data2 = read.csv("./origin_data/Allen TableS1.Mutation_list_all_patients.csv")
data4 = read.table("./origin_data/TPM_RSEM_VAScience2015.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
data5 = read.table("origin_data/skcm_dfci_2015/data_clinical_patient.txt",sep = "\t",header = T,row.names = 1)
data6 = read.table("origin_data/skcm_dfci_2015/data_clinical_sample.txt",sep = "\t",header = T,row.names = 1)

library(stringi)
tmp_ensembl = stri_sub(data4$NAME,1,15)
data4$NAME = tmp_ensembl

gtf_data = rtracklayer::import('/web/Refence_datasets/dataset14/Homo_sapiens.GRCh37.87.gtf')
gtf_data = as.data.frame(gtf_data)
gtf_data2 = gtf_data[,c("type","gene_id","gene_name")]
gtf_data2 = gtf_data2[!duplicated( gtf_data2 ),]
gtf_data2 = gtf_data2[gtf_data2$type == "gene",]
rownames(gtf_data2) = gtf_data2$gene_id
share_gene = intersect(tmp_ensembl,gtf_data2$gene_id)

data4 = data4[data4$NAME %in% share_gene,]
rownames(data4) = data4$NAME
data4 = data4[,-1]
data4 = data4[share_gene,]
gtf_data2 = gtf_data2[share_gene,]

data4 = as.matrix(data4)
rownames(data4) = gtf_data2$gene_name
data4 = limma::avereps(data4)
data4 = as.data.frame(data4)

colnames(data4) = do.call(rbind,strsplit(do.call(rbind,strsplit(x = colnames(data4),split = "-"))[,2],"_"))[,2]

colnames(data4)[!colnames(data4) %in% rownames(data1)]

#############
sum(rownames(data1) %in% data2$patient)
sum(unique(data2$patient) %in% rownames(data1)) #确保了临床和突变检查的患者是对应的

data = data.frame("OS_TIME" = data1$overall_survival,
                  "OS_STATUS" = data1$dead,
                  "PFS_TIME" = data1$progression_free,
                  "PFS_STATUS" = data1$progression,
                  "RECIST" = data1$RECIST,
                  "RESPONSE" = data1$group,
                  "TMB" = data6$TMB_NONSYNONYMOUS,
                  "THERAPY" = data5$TREATMENT,
                  "origin_therapy" = data5$TREATMENT
)
rownames(data) = rownames(data1) 

table(data$OS_TIME)
table(data$OS_STATUS)
table(data$PFS_TIME)
table(data$PFS_STATUS)
table(data$RECIST)
table(data$RESPONSE)
table(data$THERAPY)

data$OS_TIME = data$OS_TIME/30
data$PFS_TIME = data$PFS_TIME/30
data$RECIST = ifelse(data$RECIST %in% c("CR","PR"),"CR/PR",ifelse(data$RECIST %in% c("PD","SD"),"PD/SD",NA))
data$RESPONSE[data$RESPONSE == "long-survival"] = "nonresponse"
data$THERAPY = "anti-CTLA4"



gene_name = unique(data2$Hugo_Symbol)
for(x in gene_name){
  tmp = unique(data2[data2$Hugo_Symbol %in% x,"patient"])
  data[,x] = NA
  data[tmp,x] = "Mutation"
  data[is.na(data[,x]),x] = "Wildtype"
  
  
  
}

share_id = intersect(rownames(data),colnames(data4)) #已预先对齐了~ ~
RNA = data4[,share_id]
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


immune_infiltration_data = read.table(file = "./origin_data/immune_infiltration.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")
celltype_geneset = list()
for( i in immune_infiltration_data$`Cell type`){
  celltype_geneset[[i]] = immune_infiltration_data$Metagene[ immune_infiltration_data$`Cell type` == i]
}

library(GSVA)
ssgsea_score = gsva(as.matrix(RNA), celltype_geneset, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE,parallel.sz=12)


datasets_rna_wes = list("dataset2"=list("WES" = data,"RNA" = RNA,"immune_infiltration" = ssgsea_score,"immune_pathway" = immune_pathway))

saveRDS(datasets_rna_wes,"./result_data/datasets2_rna_wes.rds")

