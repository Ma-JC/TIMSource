datasets = readRDS("data/Refence_dataset/Refence_datasets_with_therapy.rds")
datasets_mu = readRDS("data/Refence_dataset/Refence_datasets_mutation_with_therapy.rds")
datasets_rna_wes = readRDS("data/Refence_dataset/Refence_datasets_rna_wes_with_therapy.rds")
gene_fre = readRDS("data/Refence_dataset/Refence_datasets_gene_maf_with_therapy.rds")
datasets_overview = read.csv("data/Refence_dataset/datasets_overview.csv",quote = "",row.names = NULL)
tooltip_text = read.csv("data/tooltip.csv",row.names = 1)
###### 森林图所需数据 #########
dataset_meta = read.csv("data/Refence_dataset/dataset_name.csv",header = F,row.names = NULL)
rownames(datasets_overview) = dataset_meta$V1

dataset_name2 = dataset_meta$V2
names(dataset_name2) = dataset_meta$V1

dataset_name = dataset_meta$V1
names(dataset_name) = dataset_meta$V2

# datasets_overview$OS = FALSE
# datasets_overview$PFS = FALSE
# datasets_overview$RECIST = FALSE
# datasets_overview$RESPONSE = FALSE
# datasets_overview$TMB = FALSE
datasets_overview$RNA = FALSE
# for(i in rownames(datasets_overview)){
#   if( "OS_TIME" %in% colnames(datasets[[i]]) ){datasets_overview[i,"OS"] = TRUE}
#   if( "PFS_TIME" %in% colnames(datasets[[i]]) ){datasets_overview[i,"PFS"] = TRUE}
#   if( "RECIST" %in% colnames(datasets[[i]]) ){datasets_overview[i,"RECIST"] = TRUE}
#   if( "RESPONSE" %in% colnames(datasets[[i]]) ){datasets_overview[i,"RESPONSE"] = TRUE}
#   if( "TMB" %in% colnames(datasets[[i]]) ){datasets_overview[i,"TMB"] = TRUE}
# }
datasets_overview[c("dataset2","dataset6","dataset8","dataset10","dataset11","dataset12","dataset13","dataset14","dataset20"),"RNA"] = TRUE



total_genes = unlist(gene_fre)
total_genes = table(total_genes)
total_genes = total_genes[ total_genes>= 3 ] # 作用不大，对于只在个别数据集存在，并且突变频率极低的基因，肯定报错
total_genes = names(total_genes)[order(total_genes,decreasing = T)]



TCGA = readRDS("data/TCGA/panacanlt_TCGA_log2.rds")
CPTAC = readRDS("data/CPTAC/CPTAC.rds")





pathway_list = readRDS("data/all_pathway.rds")

ref_total_OS_single = readRDS("data/Explore/ref_OS_total_survival.rds")
ref_total_PFS_single = readRDS("data/Explore/ref_PFS_total_survival.rds")
ref_total_RECIST_single = readRDS("data/Explore/ref_total_RECIST.rds")
ref_total_RESPONSE_single = readRDS("data/Explore/ref_total_RESPONSE.rds")
ref_total_immune_infiltrating_single = readRDS("data/Explore/ref_total_immune_infiltration.rds")
ref_total_immune_pathway_single = readRDS("data/Explore/ref_total_immune_pathway.rds")

ref_total_OS_pm = readRDS("data/Explore/ref_OS_total_survival_pm.rds")
ref_total_PFS_pm = readRDS("data/Explore/ref_PFS_total_survival_pm.rds")
ref_total_RECIST_pm = readRDS("data/Explore/ref_total_RECIST_pm.rds")
ref_total_RESPONSE_pm = readRDS("data/Explore/ref_total_RESPONSE_pm.rds")
ref_total_immune_infiltrating_pm = readRDS("data/Explore/ref_total_immune_infiltration_pm.rds")
ref_total_immune_pathway_pm = readRDS("data/Explore/ref_total_immune_pathway_pm.rds")

# TCGA_overview = read.csv("data/TCGA/TCGA_overview.csv",quote = "")
# TCGA_total_immune_infiltration_single = readRDS("data/Explore/TCGA_total_immune_infiltration.rds")
# TCGA_total_immune_pathway_single = readRDS("data/Explore/TCGA_total_immune_pathway.rds")
# TCGA_total_immune_infiltration_pm = readRDS("data/Explore/TCGA_total_immune_infiltration_pm.rds")
# TCGA_total_immune_pathway_pm = readRDS("data/Explore/TCGA_total_immune_pathway_pm.rds")
# 
# CPTAC_overview = read.csv("data/CPTAC/CPTAC_overview.csv",quote = "")
# CPTAC_total_immune_infiltration_rna_single = readRDS("data/Explore/CPTAC_total_immune_infiltration_rna.rds")
# CPTAC_total_immune_pathway_rna_single = readRDS("data/Explore/CPTAC_total_immune_pathway_rna.rds")
# CPTAC_total_immune_infiltration_rna_pm = readRDS("data/Explore/CPTAC_total_immune_infiltration_rna_pm.rds")
# CPTAC_total_immune_pathway_rna_pm = readRDS("data/Explore/CPTAC_total_immune_pathway_rna_pm.rds")
# 
# CPTAC_total_immune_infiltration_protein_single = readRDS("data/Explore/CPTAC_total_immune_infiltration_protein.rds")
# CPTAC_total_immune_pathway_protein_single = readRDS("data/Explore/CPTAC_total_immune_pathway_protein.rds")
# CPTAC_total_immune_infiltration_protein_pm = readRDS("data/Explore/CPTAC_total_immune_infiltration_protein_pm.rds")
# CPTAC_total_immune_pathway_protein_pm = readRDS("data/Explore/CPTAC_total_immune_pathway_protein_pm.rds")

###################Explore data clean & preprocess ######################################

for(i in names(ref_total_OS_single)){
  ref_total_OS_single[[i]] = round(ref_total_OS_single[[i]],3)
}

for(i in names(ref_total_PFS_single)){
  ref_total_PFS_single[[i]] = round(ref_total_PFS_single[[i]],3)
}

for(i in names(ref_total_RECIST_single)){
  tmp = as.data.frame(ref_total_RECIST_single[[i]])
  ref_total_RECIST_single[[i]] = cbind(tmp[,1:2],apply(tmp[3:8],2,function(x){round(as.numeric(x),3)}))
}

for(i in names(ref_total_RECIST_single)){


  a = ref_total_RECIST_single[[i]]$`CR & PR(Mutation/total)`
  b = ref_total_RECIST_single[[i]]$`SD & PD(Mutation/total)`

  ref_total_RECIST_single[[i]]$`CR & PR(Mutation/total)` = b
  ref_total_RECIST_single[[i]]$`SD & PD(Mutation/total)` = a
  # ref_total_RECIST_single[[i]]$OR = 1/ref_total_RECIST_single[[i]]$OR
  # a = 1/ref_total_RECIST_single[[i]]$`Lower(95%)`
  # b = 1/ref_total_RECIST_single[[i]]$`Upper(95%)`
  # ref_total_RECIST_single[[i]]$`Lower(95%)` = b
  # ref_total_RECIST_single[[i]]$`Upper(95%)` = a
}

for(i in names(ref_total_RESPONSE_single)){
  tmp = as.data.frame(ref_total_RESPONSE_single[[i]])
  ref_total_RESPONSE_single[[i]] = cbind(tmp[,1:2],apply(tmp[3:8],2,function(x){round(as.numeric(x),3)}))
}

for(i in names(ref_total_RESPONSE_single)){
  
  
  ref_total_RESPONSE_single[[i]]$OR = 1/ref_total_RESPONSE_single[[i]]$OR
  a = 1/ref_total_RESPONSE_single[[i]]$`Lower(95%)`
  b = 1/ref_total_RESPONSE_single[[i]]$`Upper(95%)`
  ref_total_RESPONSE_single[[i]]$`Lower(95%)` = b
  ref_total_RESPONSE_single[[i]]$`Upper(95%)` = a
}


for(i in names(ref_total_immune_infiltrating_single)){
  for(j in names(ref_total_immune_infiltrating_single[[i]])){
    tmp = ref_total_immune_infiltrating_single[[i]][[j]]
    tmp[ tmp == "None" ] = NA
    tmp = as.data.frame(tmp)
    tmp = as.data.frame(apply(tmp,2,function(x){round(as.numeric(x),3)}))
    rownames(tmp) = rownames(ref_total_immune_infiltrating_single[[i]][[j]])
    ref_total_immune_infiltrating_single[[i]][[j]] = tmp
  }
}

for(i in names(ref_total_immune_pathway_single)){
  for(j in names(ref_total_immune_pathway_single[[i]])){
    tmp = ref_total_immune_pathway_single[[i]][[j]]
    tmp[ tmp == "None" ] = NA
    tmp = as.data.frame(tmp)
    tmp = as.data.frame(apply(tmp,2,function(x){round(as.numeric(x),3)}))
    rownames(tmp) = rownames(ref_total_immune_pathway_single[[i]][[j]])
    ref_total_immune_pathway_single[[i]][[j]] = tmp
  }
}

for(i in names(ref_total_OS_pm)){
  tmp = as.data.frame(ref_total_OS_pm[[i]])
  ref_total_OS_pm[[i]] = cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:8],2,function(x){round(as.numeric(x),3)}))
}

for(i in names(ref_total_PFS_pm)){
  tmp = as.data.frame(ref_total_PFS_pm[[i]])
  ref_total_PFS_pm[[i]] = cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:8],2,function(x){round(as.numeric(x),3)}))
}

for(i in names(ref_total_RECIST_pm)){
  tmp = as.data.frame(ref_total_RECIST_pm[[i]])
  ref_total_RECIST_pm[[i]] = cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3:5],apply(tmp[6:11],2,function(x){round(as.numeric(x),3)}))
}

for(i in names(ref_total_RECIST_pm)){


  a = ref_total_RECIST_pm[[i]]$`CR & PR(Mutation/total)`
  b = ref_total_RECIST_pm[[i]]$`SD & PD(Mutation/total)`

  ref_total_RECIST_pm[[i]]$`CR & PR(Mutation/total)` = b
  ref_total_RECIST_pm[[i]]$`SD & PD(Mutation/total)` = a
  # ref_total_RECIST_pm[[i]]$OR = 1/ref_total_RECIST_pm[[i]]$OR
  # a = 1/ref_total_RECIST_pm[[i]]$`Lower(95%)`
  # b = 1/ref_total_RECIST_pm[[i]]$`Upper(95%)`
  # ref_total_RECIST_pm[[i]]$`Lower(95%)` = b
  # ref_total_RECIST_pm[[i]]$`Upper(95%)` = a
}

for(i in names(ref_total_RESPONSE_pm)){
  tmp = as.data.frame(ref_total_RESPONSE_pm[[i]])
  ref_total_RESPONSE_pm[[i]] = cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3:5],apply(tmp[6:11],2,function(x){round(as.numeric(x),3)}))
}

for(i in names(ref_total_RESPONSE_pm)){
  
  
  ref_total_RESPONSE_pm[[i]]$OR = 1/ref_total_RESPONSE_pm[[i]]$OR
  a = 1/ref_total_RESPONSE_pm[[i]]$`Lower(95%)`
  b = 1/ref_total_RESPONSE_pm[[i]]$`Upper(95%)`
  ref_total_RESPONSE_pm[[i]]$`Lower(95%)` = b
  ref_total_RESPONSE_pm[[i]]$`Upper(95%)` = a
}

for(i in names(ref_total_immune_infiltrating_pm)){
  for(j in names(ref_total_immune_infiltrating_pm[[i]])){
    tmp = ref_total_immune_infiltrating_pm[[i]][[j]]
    tmp[ tmp == "None" ] = NA
    tmp = as.data.frame(tmp)
    tmp = as.data.frame(cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:6],2,function(x){round(as.numeric(x),3)})))
    rownames(tmp) = rownames(ref_total_immune_infiltrating_pm[[i]][[j]])
    ref_total_immune_infiltrating_pm[[i]][[j]] = tmp
  }
}

for(i in names(ref_total_immune_pathway_pm)){
  for(j in names(ref_total_immune_pathway_pm[[i]])){
    tmp = ref_total_immune_pathway_pm[[i]][[j]]
    tmp[ tmp == "None" ] = NA
    tmp = as.data.frame(tmp)
    tmp = as.data.frame(cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:6],2,function(x){round(as.numeric(x),3)})))
    rownames(tmp) = rownames(ref_total_immune_pathway_pm[[i]][[j]])
    ref_total_immune_pathway_pm[[i]][[j]] = tmp
  }
}

# for(i in names(TCGA_total_immune_infiltration_single)){
#   for(j in names(TCGA_total_immune_infiltration_single[[i]])){
#     tmp = TCGA_total_immune_infiltration_single[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(apply(tmp,2,function(x){round(as.numeric(x),3)}))
#     rownames(tmp) = rownames(TCGA_total_immune_infiltration_single[[i]][[j]])
#     TCGA_total_immune_infiltration_single[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(TCGA_total_immune_pathway_single)){
#   for(j in names(TCGA_total_immune_pathway_single[[i]])){
#     tmp = TCGA_total_immune_pathway_single[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(apply(tmp,2,function(x){round(as.numeric(x),3)}))
#     rownames(tmp) = rownames(TCGA_total_immune_pathway_single[[i]][[j]])
#     TCGA_total_immune_pathway_single[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(TCGA_total_immune_infiltration_pm)){
#   for(j in names(TCGA_total_immune_infiltration_pm[[i]])){
#     tmp = TCGA_total_immune_infiltration_pm[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:6],2,function(x){round(as.numeric(x),3)})))
#     rownames(tmp) = rownames(TCGA_total_immune_infiltration_pm[[i]][[j]])
#     TCGA_total_immune_infiltration_pm[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(TCGA_total_immune_pathway_pm)){
#   for(j in names(TCGA_total_immune_pathway_pm[[i]])){
#     tmp = TCGA_total_immune_pathway_pm[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:6],2,function(x){round(as.numeric(x),3)})))
#     rownames(tmp) = rownames(TCGA_total_immune_pathway_pm[[i]][[j]])
#     TCGA_total_immune_pathway_pm[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(CPTAC_total_immune_infiltration_rna_single)){
#   for(j in names(CPTAC_total_immune_infiltration_rna_single[[i]])){
#     tmp = CPTAC_total_immune_infiltration_rna_single[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(apply(tmp,2,function(x){round(as.numeric(x),3)}))
#     rownames(tmp) = rownames(CPTAC_total_immune_infiltration_rna_single[[i]][[j]])
#     CPTAC_total_immune_infiltration_rna_single[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(CPTAC_total_immune_pathway_rna_single)){
#   for(j in names(CPTAC_total_immune_pathway_rna_single[[i]])){
#     tmp = CPTAC_total_immune_pathway_rna_single[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(apply(tmp,2,function(x){round(as.numeric(x),3)}))
#     rownames(tmp) = rownames(CPTAC_total_immune_pathway_rna_single[[i]][[j]])
#     CPTAC_total_immune_pathway_rna_single[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(CPTAC_total_immune_infiltration_protein_single)){
#   for(j in names(CPTAC_total_immune_infiltration_protein_single[[i]])){
#     tmp = CPTAC_total_immune_infiltration_protein_single[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(apply(tmp,2,function(x){round(as.numeric(x),3)}))
#     rownames(tmp) = rownames(CPTAC_total_immune_infiltration_protein_single[[i]][[j]])
#     CPTAC_total_immune_infiltration_protein_single[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(CPTAC_total_immune_pathway_protein_single)){
#   for(j in names(CPTAC_total_immune_pathway_protein_single[[i]])){
#     tmp = CPTAC_total_immune_pathway_protein_single[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(apply(tmp,2,function(x){round(as.numeric(x),3)}))
#     rownames(tmp) = rownames(CPTAC_total_immune_pathway_protein_single[[i]][[j]])
#     CPTAC_total_immune_pathway_protein_single[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(CPTAC_total_immune_infiltration_rna_pm)){
#   for(j in names(CPTAC_total_immune_infiltration_rna_pm[[i]])){
#     tmp = CPTAC_total_immune_infiltration_rna_pm[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:6],2,function(x){round(as.numeric(x),3)})))
#     rownames(tmp) = rownames(CPTAC_total_immune_infiltration_rna_pm[[i]][[j]])
#     CPTAC_total_immune_infiltration_rna_pm[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(CPTAC_total_immune_pathway_rna_pm)){
#   for(j in names(CPTAC_total_immune_pathway_rna_pm[[i]])){
#     tmp = CPTAC_total_immune_pathway_rna_pm[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:6],2,function(x){round(as.numeric(x),3)})))
#     rownames(tmp) = rownames(CPTAC_total_immune_pathway_rna_pm[[i]][[j]])
#     CPTAC_total_immune_pathway_rna_pm[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(CPTAC_total_immune_infiltration_protein_pm)){
#   for(j in names(CPTAC_total_immune_infiltration_protein_pm[[i]])){
#     tmp = CPTAC_total_immune_infiltration_protein_pm[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:6],2,function(x){round(as.numeric(x),3)})))
#     rownames(tmp) = rownames(CPTAC_total_immune_infiltration_protein_pm[[i]][[j]])
#     CPTAC_total_immune_infiltration_protein_pm[[i]][[j]] = tmp
#   }
# }
# 
# for(i in names(CPTAC_total_immune_pathway_protein_pm)){
#   for(j in names(CPTAC_total_immune_pathway_protein_pm[[i]])){
#     tmp = CPTAC_total_immune_pathway_protein_pm[[i]][[j]]
#     tmp[ tmp == "None" ] = NA
#     tmp = as.data.frame(tmp)
#     tmp = as.data.frame(cbind(apply(tmp[1:2],2,function(x){round(as.numeric(x),3)}),tmp[3],apply(tmp[4:6],2,function(x){round(as.numeric(x),3)})))
#     rownames(tmp) = rownames(CPTAC_total_immune_pathway_protein_pm[[i]][[j]])
#     CPTAC_total_immune_pathway_protein_pm[[i]][[j]] = tmp
#   }
# }
##########################################################################################
library("msigdbr")
KEGG = msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:KEGG") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
REACTOME = msigdbr(species="Homo sapiens",category="C2",subcategory = "CP:REACTOME") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
HALLMARK = msigdbr(species="Homo sapiens",category="H") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
GO_BP = msigdbr(species="Homo sapiens",category="C5",subcategory = "GO:BP") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
GO_CC = msigdbr(species="Homo sapiens",category="C5",subcategory = "GO:CC") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
GO_MF = msigdbr(species="Homo sapiens",category="C5",subcategory = "GO:MF") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)

pathway_database = list("GO_BP"=GO_BP,"GO_CC"=GO_CC,"GO_MF"=GO_MF,"KEGG"=KEGG,"REACTOME"=REACTOME,"HALLMARK"=HALLMARK)


CPTAC_name = c("2021_GBM","2021_PDAC","2020_BI","2020_OVA","2020_LSCC","2020_LUAD","2020_HNSCC","2020_UCEC","2020_CCRCC","2019_Colon")
TCGA_name = names(TCGA)
ref_name = names(datasets)

#TCGA mutation type as same as CPTAC
TCGA_mutation_type = c("All","Missense_Mutation","Frame_Shift_Del","Splice_Site","Nonsense_Mutation","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonstop_Mutation","Translation_Start_Site")

Mcolor = c("red","blue")
names(Mcolor) = c("Mutation","Wildtype")
inf_names = factor(
  c('Activated B cell',
    'Activated CD4 T cell',
    'Activated CD8 T cell',
    'Central memory CD4 T cell',
    'Central memory CD8 T cell',
    'Effector memeory CD4 T cell',
    'Effector memeory CD8 T cell',
    'Gamma delta T cell',
    'Immature  B cell',
    'Memory B cell',
    'Regulatory T cell',
    'T follicular helper cell',
    'Type 1 T helper cell',
    'Type 17 T helper cell',
    'Type 2 T helper cell',
    'Activated dendritic cell',
    'CD56bright natural killer cell',
    'CD56dim natural killer cell',
    'Eosinophil',
    'Immature dendritic cell',
    'Macrophage','Mast cell',
    'MDSC','Monocyte',
    'Natural killer cell',
    'Natural killer T cell',
    'Neutrophil',
    'Plasmacytoid dendritic cell')
)

sig_names = factor(
  c(
    '6-gene IFN signature',
    '18-gene IFN signature',
    'Gene expression profile',
    'Cytolytic activity',
    '13 T-cell signature',
    'Effective T cell score',
    'Immune checkpoint expression',
    'TLS'
  )
)

# for(i in names(datasets_mu)){
#   for(j in colnames(datasets_mu[[i]])){
#     if(!j %in% c("Start_Position","End_Position","Position")){
#       datasets_mu[[i]][,j] = as.factor(datasets_mu[[i]][,j])
#     }
#   }
# }

