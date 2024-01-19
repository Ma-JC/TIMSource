library(GSVA)

###### RNA-seq Data preprocess #########
data = read.table("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
rownames(data) = data[,1]
data = data[,-1]

genes_name = do.call(rbind,strsplit(rownames(data),"\\|"))[,1] 
m = genes_name != "?"
data = data[m,]

data = as.matrix(data)
rownames(data) = genes_name[m]
data = limma::avereps(data)
data = as.data.frame(data)

data = t(data)
rownames(data) = gsub(pattern = "(TCGA-..-....)-.*",replacement = "\\1",x = rownames(data))
data[ is.na(data) ] = 0
data = limma::avereps(data)
data = as.data.frame(t(data))

data = log2(data + 1)

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

immune_infiltration_data = read.table(file = "immune_infiltration.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")


immune_pathway = vector()
for( i in names(immune)){
  
  tmp_genes = immune[[i]][immune[[i]] %in% rownames(data)]
  tmp_data = data[tmp_genes,]
  
  immune_pathway = rbind(immune_pathway,colMeans(tmp_data))
  
}
rownames(immune_pathway) = names(immune)

celltype_geneset = list()
for( i in immune_infiltration_data$`Cell type`){
  celltype_geneset[[i]] = immune_infiltration_data$Metagene[ immune_infiltration_data$`Cell type` == i]
}

ssgsea_score = gsva(as.matrix(data), celltype_geneset, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE,parallel.sz=16)

######## WES Data preprocess ############
maf1 = read.maf("mc3.v0.2.8.PUBLIC.maf",isTCGA = F,useAll = T,removeDuplicatedVariants = T)
maf2 = read.maf("mc3.v0.2.8.PUBLIC.maf",isTCGA = T,useAll = T,removeDuplicatedVariants = T)

maf_patient_id = unique(maf2@data$Tumor_Sample_Barcode)
rna_patient_id = unique(gsub(pattern = "(TCGA-..-....)-.*",replacement = "\\1",x = colnames(data)))

share_patient_id = intersect(maf_patient_id,rna_patient_id)

clinical = read.table("clinical_PANCAN_patient_with_followup.tsv",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F)
clinical_endpoint = read.table("clinical.txt",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F)

total_maf_matrix = maf2@data

creat_tumor_maf_rna = function(tumor_type,clinical,total_maf_matrix,total_mRNA_matrix,total_immune_pathway,total_immune_infiltration){
  
  patient_id = clinical$bcr_patient_barcode[ clinical$type == tumor_type]
  
  maf_patient_id = as.character(total_maf_matrix$Tumor_Sample_Barcode)
  maf_tumor_barcode =  unique(maf_patient_id[maf_patient_id %in% patient_id])
  rna_tumor_bardoce = colnames(total_mRNA_matrix)
  
  sub_clinical = clinical[ clinical$bcr_patient_barcode %in% maf_tumor_barcode,]
  sub_maf = total_maf_matrix[maf_patient_id %in% patient_id,]
  sub_rna = total_mRNA_matrix[,maf_tumor_barcode[ maf_tumor_barcode %in% rna_tumor_bardoce]]
  sub_immune_pathway = total_immune_pathway[,maf_tumor_barcode[ maf_tumor_barcode %in% rna_tumor_bardoce]]
  sub_immune_infiltration = total_immune_infiltration[,maf_tumor_barcode[ maf_tumor_barcode %in% rna_tumor_bardoce]]
  
  write.table(sub_maf,file = "tmp.maf",sep = "\t",quote = F,col.names = T,row.names = F)
  sub_maf = read.maf("tmp.maf")
  return(list("clinical" = sub_clinical,"maf" = sub_maf,"rna" = sub_rna,"immune_pathway" = sub_immune_pathway,"immune_infiltration" = sub_immune_infiltration))
}

result_list = list()
for( i in unique(clinical_endpoint$type)){
  result_list[[i]] = creat_tumor_maf_rna(tumor_type = i,clinical = clinical_endpoint,total_maf_matrix = total_maf_matrix,total_mRNA_matrix = data,total_immune_pathway = immune_pathway,total_immune_infiltration = ssgsea_score)
}

real_clinical = read.table("real_clinical.txt",sep = "\t",quote = "",header = T,check.names = F,stringsAsFactors = F,row.names = 2)
for( i in unique(clinical_endpoint$type)){
  result_list[[i]][["clinical_detial"]] = real_clinical[result_list[[i]][["clinical"]]$bcr_patient_barcode,]
}


saveRDS(result_list,"panacanlt_TCGA_log2.rds")


