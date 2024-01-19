library(IMvigor210CoreBiologies)
library(limma)

###################RNA#######################
data(cds)

RNA = counts(cds)
clinical1 = pData(cds)
fmeta = fData(cds)
clinical1$ANONPT_ID = paste("P_",clinical1$ANONPT_ID,sep = "")

sum(rownames(RNA) == fmeta$entrez_id)
rownames(RNA) = fmeta$symbol
m = !is.na(rownames(RNA))
RNA = RNA[m,]
fmeta = fmeta[m,]
RNA = as.matrix(RNA)
total = colSums(RNA)
len = fmeta$length
tmp = RNA/len
tmp = t(tmp)
tmp = tmp/total
tmp = t(tmp)
tmp = 10^9 * tmp
FPKM = tmp
remove(tmp)
FPKM = avereps(FPKM)

FPKM = as.data.frame(FPKM)
colnames(FPKM) == rownames(clinical1)
colnames(FPKM) = clinical1$ANONPT_ID

sum(duplicated(colnames(FPKM))) #有一个重复

FPKM = t(FPKM)
FPKM = avereps(FPKM) #去重复
FPKM = as.data.frame(t(FPKM))

################################################

data(fmone)

clinical2 = pData(fmone)
Mutation = any_mutation(fmone)
clinical2$ANONPT_ID = paste("P_",clinical2$ANONPT_ID,sep = "")
colnames(Mutation) = clinical2$ANONPT_ID
tmp = apply(Mutation,1,function(x){
  m = as.integer(x)
  m = m +1
  c("Wildtype","Mutation")[m]
  })

rownames(tmp) = colnames(Mutation)


share_id = intersect(rownames(tmp),colnames(FPKM))

WES = tmp[share_id,] #已预先对齐了~ ~
WES = as.data.frame(WES)
FPKM = FPKM[,share_id] 

FPKM = log2(FPKM + 1)
################################################

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
  tmp_genes = immune[[i]][immune[[i]] %in% rownames(FPKM)]
  tmp_data = FPKM[tmp_genes,]
  
  immune_pathway = rbind(immune_pathway,colMeans(tmp_data))
  
} 
rownames(immune_pathway) = names(immune)


immune_infiltration_data = read.table(file = "origin_data/immune_infiltration.txt",header = T,check.names = F,stringsAsFactors = F,quote = "",sep = "\t")
celltype_geneset = list()
for( i in immune_infiltration_data$`Cell type`){
  celltype_geneset[[i]] = immune_infiltration_data$Metagene[ immune_infiltration_data$`Cell type` == i]
}

library(GSVA)
ssgsea_score = gsva(as.matrix(FPKM), celltype_geneset, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE,parallel.sz=12)

datasets_rna_wes = list("dataset11"=list("WES" = WES,"RNA" = FPKM,"immune_infiltration" = ssgsea_score,"immune_pathway" = immune_pathway))

saveRDS(datasets_rna_wes,"result_data/datasets11_rna_wes.rds")
