scdata <- readRDS('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/scdata_all_anno_lv1.rds')
meta_data <- read.csv('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/meta_data_macmono-nkt.csv',row.names=1)
scdata@meta.data <- meta_data
colnames(scdata@meta.data)
meta_data <- scdata@meta.data %>% 
  select(-c('anno_res_data2','anno_res_data2_data2','anno_res_data2_data2_data2'))
scdata@meta.data <- meta_data
scdata$anno_res_lv2 %>% unique()
celltype <- scdata$anno_res_lv2 %>% unique() %>% .[!is.na(.)]
#### Step2.1.1.2 GSEA analysis --------------------------------------------------------
#> Iteratively calculate the treatment results between two groups across different cell types.
#> 
celltype_list <-
  scdata$anno_res_lv2 %>% 
  unique() %>% .[!is.na(.)] %>% 
  as.character()


## differential gene calculation ----------------------------------------------------------------
for(celltype in celltype_list){
  print(celltype)
  res <- fun_diffgene(scdata_use = scdata,celltype)
}

fun_diffgene <- function(scdata_use,celltype){
  if(celltype != 'allcell'){
    scdata_sub <- subset(
      x = scdata_use,
      subset = anno_res_lv2 %in% c(celltype)
    )
    file_name = str_replace_all(
      string = celltype,
      pattern = ' ',
      replacement = '-'
    ) %>% str_replace_all(pattern = '/',replacement = '_')
  }else{
    scdata_sub <- scdata_use
    file_name = str_replace_all(
      string = celltype,
      pattern = ' ',
      replacement = '-'
    ) %>% str_replace_all(pattern = '/',replacement = '_')
  }
  table(scdata_sub$group_maf) %>% print()
  # table(scdata_sub$batch) %>% print()
  diff_gene_sub <- FindMarkers(
    object = scdata_sub,
    ident.1 = 'complement-MUT',
    ident.2 = 'complement-WT',
    min.pct = 0,
    group.by = 'group_maf',
    only.pos = FALSE
  ) %>% 
    mutate(
      group = ifelse(p_val_adj < 0.05,yes  = 'credible',no = NA)
    ) %>% 
    arrange(group,desc(avg_log2FC))
  file_path <- 
    paste(
      '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/cd45_diffgene/diffgene_',file_name,'.csv',
      sep = ''
    )
  write.csv(diff_gene_sub,file_path)
}




## 进行gsea ------------------------------------------------------------------
library(msigdbr)

# anno_res_lv1 ------------------------------------------------------------
celltype_list <- scdata@meta.data %>% 
  pull(anno_res_lv1) %>% 
  unique() %>% .[!is.na(.)] %>% 
  as.character()
celltype_list
# scdata_use <- subset(scdata,cell = colnames(scdata)[!is.na(scdata$anno_res_lv2)])
for(celltype in celltype_list[8:10]){
  # celltype <- celltype_list[23]
  print(celltype)
  res <- fun_gsea(scdata_use = scdata,celltype)
}

# anno_res_lv2 ------------------------------------------------------------
celltype_list <- scdata@meta.data %>% 
  pull(anno_res_lv2) %>% 
  unique() %>% .[!is.na(.)] %>% 
  as.character()
celltype_list
# scdata_use <- subset(scdata,cell = colnames(scdata)[!is.na(scdata$anno_res_lv2)])
for(celltype in celltype_list){
  # celltype <- celltype_list[7]
  print(celltype)
  res <- fun_gsea(
    scdata_use = scdata,
    celltype = celltype,
    batch = 'group_maf',
    group_control = 'complement-MUT',
    cellanno = 'anno_res_lv1',
    species = 'Homo sapiens',
    # savefile_path = 
    savefile = TRUE
  )
}

res_gsea <- fun_gsea(
  scdata_use = scdata,
  celltype = celltype,
  batch = 'group_maf',
  group_control = 'complement-MUT',
  cellanno = 'anno_res_lv2',
  species = 'Homo sapiens',
  savefile = TRUE
)
fun_gsea <- function(
    scdata_use,
    celltype,batch = 'group_maf',
    group_control = 'complement-MUT', 
    cellanno = 'anno_res_lv2',
    species = "Homo sapiens",
    savefile = TRUE){
  library(clusterProfiler)
  library(purrr)
  library(msigdbr)
  meta_data <- scdata_use@meta.data %>% 
    dplyr::select(all_of(c(batch, cellanno))) %>% 
    rename(batch = !!sym(batch), cellanno = !!sym(cellanno))
  scdata_use@meta.data <- meta_data
  if(celltype != 'allcell'){
    scdata_sub <- subset(
      x = scdata_use,
      subset = cellanno %in% c(celltype)
    )
    file_name = str_replace_all(
      string = celltype,
      pattern = ' ',
      replacement = '-'
    ) %>% str_replace_all(pattern = '/',replacement = '_')
  }else{
    scdata_sub <- scdata_use
    file_name = str_replace_all(
      string = celltype,
      pattern = ' ',
      replacement = '-'
    ) %>% str_replace_all(pattern = '/',replacement = '_')
  }
  table(scdata_sub$cellanno) %>% print()
  table(scdata_sub$batch) %>% print()
  library(presto)
  exp.genes <- wilcoxauc(scdata_sub, 'batch')
  # dplyr::count(exp.genes, group) %>% print()
  cluster0.genes<- exp.genes %>%
    filter(group == group_control) %>%
    arrange(desc(logFC),desc(auc)) %>%
    dplyr::select(feature,logFC,auc)
  
  ranks<- cluster0.genes$logFC
  names(ranks) <- cluster0.genes$feature
  geneList <- ranks
  
  gene_list <- msigdbr(species = species)
  
  ##### Hallmark -------------------------------------------------------------
  # gene_list$gs_cat %>% unique()
  cat('Hallmark...')
  gene_list_sub <- gene_list %>% 
    dplyr::filter((gs_cat == 'H'))%>% 
    dplyr::select(gs_name,gene_symbol) %>% 
    rename_all(~c('term','gene'))
  egmt_hallmark <- clusterProfiler::GSEA(geneList, 
                                         TERM2GENE=gene_list_sub, 
                                         minGSSize = 1,
                                         pvalueCutoff = 0.5,
                                         verbose=FALSE)
  if(nrow(egmt_hallmark@result)>0){
    gsea_results_hallmark <- egmt_hallmark@result %>% 
      arrange(desc(NES)) %>%
      mutate(Description = str_replace(Description,'HALLMARK_','') %>% 
               str_split("_") %>% 
               map_chr(~ str_to_title(.x) %>% paste(collapse = " "))) %>% 
      mutate(Description = factor(
        Description,levels = Description %>% unique() %>% rev())
      ) %>% 
      mutate(group = ifelse(test = p.adjust<0.05,yes = 'Credible','NS')) %>% 
      arrange(desc(NES)) 
    
    if(savefile){
      file_path <- 
        paste(
          '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/gsea/gsea_Hallmark_',file_name,'.csv',
          sep = ''
        )
      write.csv(gsea_results_hallmark,file_path)
    }
  }else{
    gsea_results_hallmark <- NULL
  }
  
  
  ##### M5 Gobp-------------------------------------------------------------
  cat('M5 Gobp...')
  gene_list_sub <- gene_list %>% 
    dplyr::filter((gs_subcat == 'GO:BP'))%>% 
    dplyr::select(gs_name,gene_symbol) %>% 
    rename_all(~c('term','gene'))
  tmp <- gene_list_sub %>% 
    group_by(term) %>% 
    summarise(counts = n()) %>% 
    arrange(counts) %>% 
    filter(counts>10) %>% 
    pull(term)
  egmt_gobp <- GSEA(geneList, 
                    TERM2GENE = gene_list_sub[gene_list_sub$term %in% tmp,], 
                    minGSSize = 1,
                    pvalueCutoff = 0.5,
                    verbose=FALSE)
  if(nrow(egmt_gobp@result)>0){
    gsea_results_gobp <- egmt_gobp@result %>% 
      arrange(desc(NES)) %>%
      mutate(group = ifelse(test = p.adjust<0.05,yes = 'Credible','NS')) %>% 
      mutate(Description = str_replace(Description,'GOBP_','') %>% 
               str_split("_") %>% 
               map_chr(~ str_to_title(.x) %>% paste(collapse = " "))) %>% 
      mutate(Description = factor(
        Description,levels = Description %>% unique() %>% rev())
      )
    if(savefile){
      file_path <- 
        paste(
          '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/gsea/gsea_gobp_',file_name,'.csv',
          sep = ''
        )
      write.csv(gsea_results_gobp,file_path)
    }
  }else{
    gsea_results_gobp <- NULL
  }
  
  ##### KEGG --------------------------------------------------------------------
  cat('M5 KEGG...')
  gene_list_sub <- gene_list %>% 
    dplyr::filter((gs_subcat == 'CP:KEGG'))%>% 
    dplyr::select(gs_name,gene_symbol) %>% 
    rename_all(~c('term','gene'))
  egmt_kegg <- GSEA(geneList = geneList, 
                    TERM2GENE=gene_list_sub, 
                    minGSSize = 1,
                    pvalueCutoff = 0.5,
                    verbose=FALSE)
  if(nrow(egmt_kegg@result)>0){
    gsea_results_kegg <- egmt_kegg@result %>% 
      arrange(desc(NES)) %>%
      mutate(group = ifelse(test = pvalue<0.05,yes = 'credible','NS')) %>% 
      mutate(Description = str_replace(Description,'KEGG_','') %>% 
               str_split("_") %>% 
               map_chr(~ str_to_title(.x) %>% paste(collapse = " "))) %>% 
      mutate(Description = factor(Description,levels = Description %>% unique() %>% rev()))
    if(savefile){
      file_path <- 
        paste(
          '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/gsea/gsea_kegg_',file_name,'.csv',
          sep = ''
        )
      write.csv(gsea_results_kegg,file_path)
    }
  }else{
    gsea_results_kegg <- NULL
  }
  
  ##### CP:REACTOME --------------------------------------------------------------------
  cat('M5 CP:REACTOME...')
  gene_list_sub <- gene_list %>% 
    dplyr::filter((gs_subcat == 'CP:REACTOME'))%>% 
    dplyr::select(gs_name,gene_symbol) %>% 
    rename_all(~c('term','gene'))
  egmt_reactome <- GSEA(geneList = geneList, 
                        TERM2GENE=gene_list_sub, 
                        minGSSize = 1,
                        pvalueCutoff = 0.5,
                        verbose=FALSE)
  if(nrow(egmt_reactome@result)>0){
    gsea_results_reactome <- egmt_reactome@result %>% 
      arrange(desc(NES)) %>%
      mutate(group = ifelse(test = pvalue<0.05,yes = 'credible','NS')) %>% 
      mutate(Description = str_replace(Description,'REACTOME_','') %>% 
               str_split("_") %>% 
               map_chr(~ str_to_title(.x) %>% paste(collapse = " "))) %>% 
      mutate(Description = factor(Description,levels = Description %>% unique() %>% rev()))
    if(savefile){
      file_path <- 
        paste(
          '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/gsea/gsea_reactome_',file_name,'.csv',
          sep = ''
        )
      write.csv(gsea_results_reactome,file_path)
    }
  }else{
    gsea_results_reactome <- NULL
  }
  
  ##### CP:WIKIPATHWAYS --------------------------------------------------------------------
  cat('CP:WIKIPATHWAYS...')
  gene_list_sub <- gene_list %>% 
    dplyr::filter((gs_subcat == 'CP:WIKIPATHWAYS'))%>% 
    dplyr::select(gs_name,gene_symbol) %>% 
    rename_all(~c('term','gene'))
  egmt_wikipathways <- GSEA(geneList = geneList, 
                            TERM2GENE=gene_list_sub, 
                            minGSSize = 10,
                            pvalueCutoff = 0.99,
                            verbose=FALSE)
  if(nrow(egmt_wikipathways@result)>0){
    gsea_results_wikipathways <- egmt_wikipathways@result %>% 
      arrange(desc(NES)) %>%
      mutate(group = ifelse(test = pvalue<0.05,yes = 'credible','NS')) %>% 
      mutate(Description = str_replace(Description,'KEGG_','') %>% 
               str_split("_") %>% 
               map_chr(~ str_to_title(.x) %>% paste(collapse = " "))) %>% 
      mutate(Description = factor(Description,levels = Description %>% unique() %>% rev()))
    if(savefile){
      file_path <- 
        paste(
          '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/gsea/gsea_wikipathways_',file_name,'.csv',
          sep = ''
        )
      write.csv(gsea_results_wikipathways,file_path)
    }
  }else{
    gsea_results_wikipathways <- NULL
  }
  return(list(
    'gsea_results_hallmark' = gsea_results_hallmark,
    'gsea_results_gobp' = gsea_results_gobp,
    'gsea_results_kegg' = gsea_results_kegg,
    'gsea_results_reactome' = gsea_results_reactome,
    'gsea_results_wikipathways' = gsea_results_wikipathways,
    'egmt_hallmark' = egmt_hallmark,
    'egmt_gobp' = egmt_gobp,
    'egmt_kegg' = egmt_kegg,
    'egmt_reactome' = egmt_reactome,
    'egmt_wikipathways' = egmt_wikipathways
  ))
}

# Organize the results of the GSEA enrichment analysis -------------------------------------------------------------
getwd()
dir.create('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/gsea_anno_res_lv1/')
file_list <- dir('./gsea')
celltype_list <- scdata@meta.data %>% 
  pull(anno_res_lv1) %>% 
  unique() %>% .[!is.na(.)] %>% 
  as.character()
celltype_list
data_gsea <- lapply(file_list, function(file_name){
  # file_name <- file_list[1]
  data <- read.csv(file = file.path('./gsea/',file_name),row.names = 1)
})
names(data_gsea) <- file_list
data_res <- lapply(celltype_list, function(celltype){
  celltype <- celltype %>%
    str_replace_all(
      pattern = ' ',
      replacement = '-'
    ) %>% str_replace_all(pattern = '/',replacement = '_')
  data_use_name <- file_list[grepl(celltype,file_list)]
  data_res <- data_gsea[data_use_name] %>% 
    do.call(rbind,.) %>% 
    ungroup() %>% 
    mutate(
      group_p = ifelse(pvalue<0.05,yes = 'Credible','NS'),
      group_p.adj = ifelse(p.adjust<0.05,yes = 'Credible','NS')
    ) %>% 
    select(-group) %>% 
    arrange(group_p,group_p.adj,desc(NES))
  file_save <- paste('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/gsea_anno_res_lv1/',celltype,'.csv',sep = '')
  write.csv(data_res,file_save)
})

# Plotting of specific GSEA pathways. -------------------------------------------------------------
pathway_select <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_KRAS_SIGNALING_DN",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)
library(presto)
scdata_use <- subset(scdata,subset = anno_res_lv1 == 'Epithelial/Malignant cell')
exp.genes <- wilcoxauc(scdata_use, 'group_maf')
dplyr::count(exp.genes, group) %>% print()
cluster0.genes<- exp.genes %>%
  filter(group == 'complement-MUT') %>%
  arrange(desc(logFC),desc(auc)) %>%
  dplyr::select(feature,logFC,auc)

ranks<- cluster0.genes$logFC
names(ranks) <- cluster0.genes$feature
geneList <- ranks
gene_list <- msigdbr(species = 'Homo sapiens')
gene_list_sub <- gene_list %>% 
  dplyr::filter((gs_name %in% pathway_select))%>% 
  dplyr::select(gs_name,gene_symbol) %>% 
  rename_all(~c('term','gene'))
gene_list_sub$term %>% unique()
egmt_hallmark <- GSEA(geneList, 
                      TERM2GENE=gene_list_sub, 
                      minGSSize = 5,
                      pvalueCutoff = 1,
                      verbose=FALSE)
enrichplot::gseaplot2(
  egmt_hallmark,
  pathway_select,
  title = "",
  base_size = 11,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = TRUE
)
ggsave(
  filename = '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs_0814/gsea_malignant_hallmarker.png',
  width = 16,height = 8,dpi = 300
)
data_plot <- egmt_hallmark@result %>% 
  arrange(desc(NES)) %>%
  mutate(Description = str_replace(Description,'HALLMARK_','') %>% 
           str_split("_") %>% 
           map_chr(~ str_to_title(.x) %>% paste(collapse = " "))) %>% 
  mutate(Description = factor(
    Description,levels = Description %>% unique() %>% rev())
  ) %>% 
  mutate(group = ifelse(test = p.adjust<0.05,yes = 'Credible','NS')) %>% 
  arrange(desc(NES)) %>% 
  mutate(
    x_loc = case_when(
      NES>0 ~ -0.1,
      NES<0 ~ 0.1
    )
  )
colnames(data_plot)
ggplot(
  data = data_plot,
  mapping = aes(x = NES,y = Description,fill = -log10(p.adjust))
) +
  geom_bar(stat = 'identity') +
  annotate(
    "segment",
    x = 0,
    y = 0,
    xend = 0,
    yend = 6.6,
    color = 'black',
    linewidth = 0.3
  ) +
  geom_text(
    data = data_plot[data_plot$NES>0,],
    aes(x = x_loc,y = Description,label = Description),
    hjust = 1,size = 6
  ) +
  geom_text(
    data = data_plot[data_plot$NES<0,],
    aes(x = x_loc,y = Description,label = Description),
    hjust = 0,size = 6
  ) +
  annotate(
    "segment",
    x = -0.1,
    y = 6.8,
    xend = -abs(data_plot$NES) %>% max(),
    yend = 6.8,
    arrow = arrow(type = "closed", length = unit(0.1, "inches")),
    color = 'black',
    size = 0.3
  ) +
  annotate(
    "text",
    x = -0.7,
    y = 7.2,
    label = 'Enriched in \nComplement MUT',
    size = 6
  ) +
  annotate(
    "segment",
    x = 0.1,
    y = 6.8,
    xend = abs(data_plot$NES) %>% max(),
    yend = 6.8,
    arrow = arrow(type = "closed", length = unit(0.1, "inches")),
    color = 'black',
    size = 0.3
  ) +
  annotate(
    "text",
    x = 0.8,
    y = 7.2,
    label = 'Enriched in \nComplement WT',
    size = 6
  ) +
  scale_fill_gradient(low = '#B2D8EE',high = '#3B6895',breaks = seq(2, 2.5, by = 0.2)) +
  coord_cartesian(ylim = c(1, 7),xlim = c(-1.8,1.8)) +
  # scale_y_continuous(limits = c(0,8)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(t = 10),
    panel.border = element_blank(),
    axis.line.x = element_line(color = 'black',linewidth = 0.3),
    axis.ticks.x = element_line(linewidth = 0.3),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 18,hjust = 0.3),
    axis.title.y = element_blank()
  )
ggsave(
  filename = '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs_0814/barplot_malignant_hallmarker.png',
  width = 9,height = 7,bg = 'white'
)
ggsave(
  filename = '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs_0814/barplot_malignant_hallmarker.pdf',
  width = 9,height = 7,bg = 'white'
)
# Generate a classic GSEA plot for specific GSEA_GO pathways ---------------------------------------------------

pathway_select <- c(
  "GOBP_Antigen_Receptor_Mediated_Signaling_Pathway",
  "GOBP_Lymphocyte_Mediated_Immunity",
  "GOBP_Immune_Response_Regulating_Signaling_Pathway",
  "GOBP_Activation_of_Immune_Response"
) %>% str_to_upper()
library(presto)
scdata_use <- subset(scdata,subset = anno_res_lv1 == 'Epithelial/Malignant cell')
exp.genes <- wilcoxauc(scdata_use, 'group_maf')
dplyr::count(exp.genes, group) %>% print()
cluster0.genes<- exp.genes %>%
  filter(group == 'complement-MUT') %>%
  arrange(desc(logFC),desc(auc)) %>%
  dplyr::select(feature,logFC,auc)
ranks<- cluster0.genes$logFC
names(ranks) <- cluster0.genes$feature
geneList <- ranks

gene_list <- msigdbr(species = 'Homo sapiens')
gene_list_sub <- gene_list %>% 
  dplyr::filter((gs_name %in% pathway_select))%>% 
  dplyr::select(gs_name,gene_symbol) %>% 
  rename_all(~c('term','gene'))
gene_list_sub$term %>% unique()
egmt_hallmark <- GSEA(geneList, 
                      TERM2GENE=gene_list_sub, 
                      minGSSize = 1,
                      pvalueCutoff = 1,
                      verbose=FALSE)
enrichplot::gseaplot2(
  egmt_hallmark,
  pathway_select,
  title = "",
  base_size = 11,
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:3,
  pvalue_table = TRUE
)
ggsave(
  filename = '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs_0814/gsea_malignant_gobp.png',
  width = 16,height = 6,dpi = 300
)
ggsave(
  filename = '/mnt/data3/jiaoxi/新辅助治疗队列/result_figs_0814/gsea_malignant_gobp.pdf',
  width = 16,height = 6#,dpi = 300
)
# write.csv(  tmp,'E:/zhuona/res/gsea_PI3K_AKR.csv')
# gseaplot2(Go_gseresult, 1:3, title = "Specific GO Biological Process in T2D group", pvalue_table = FALSE)  #1:3：这表示在图上显示前3条富集结果，也可以根据自己分析需要指定输出某一条结果；Go_gseresult：GO富集分析结果；title：加上标题；pvalue_table：是否在图上显示P值列表。
tmp <- egmt_hallmark@result
# endline -----------------------------------------------------------------


