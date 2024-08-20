rm(list = ls())
library(AUCell)
library(ggalluvial)
library(cowplot)
library(patchwork)
library(ggpubr)
library(cowplot)
library(tibble)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(biomaRt)
library(tidyverse)
library(GSVA)
library(msigdbr)
library(ggpubr)
library(entropy)
library(immunarch)

# Figs7 A dotplot of CD8+ T cell-----------------------------------------------------------------
scRNA_T_sub <- readRDS("D:/jiaoxi/scRNA_T_cd8_再注释.rds")
scRNA_T_sub <- subset(scRNA_T_sub,subset = anno_lv3 != 'low_qc cell' )
scRNA_T_sub$anno_lv3 %>% unique()
cell_type <- scRNA_T_sub$anno_lv3 %>% unique() %>% .[c(1,2,4,5,8)]
cell_type
scRNA_T_sub <- subset(scRNA_T_sub,subset = anno_lv3 %in% cell_type)
scRNA_T_sub$batch <- 
  scRNA_T_sub$batch %>% 
  factor(levels = c('MOCK','C5aRa','PD1','C5aRA_PD1'))
meta_data <- scRNA_T_sub@meta.data
meta_data <- meta_data %>%
  mutate(
    batch = case_when(
      batch == 'MOCK' ~ 'Mock',
      batch == 'C5aRa' ~ 'C5aRA',
      batch == 'PD1' ~ 'PD1',
      batch == 'C5aRA_PD1' ~ 'C5aRA_PD1'
    )
  )
meta_data$batch %>% unique()
scRNA_T_sub@meta.data <- meta_data

Idents(scRNA_T_sub) <- scRNA_T_sub$anno_lv3
scRNA_T_sub$anno_lv3 %>% table()
gene_list <- c(
  'Havcr2','Tigit','Pdcd1','Ctla4',
  'Gzmb','Gzmk','Icos','Cd28','Nkg7','Ifng',
  'Il10ra','Ccl4','Ccl3','Ccl5','Stat1','Irf7',
  'Bhlhe40','Runx2','Isg15','Ly6c1','Ly6c2'
)
cell_level <- c('Mock','C5aRA','PD1','C5aRA_PD1')
data.usage <- DotPlot(
  object = scRNA_T_sub,
  features = gene_list,
  split.by = 'batch',
  group.by = 'anno_lv3',
  cols = c("lightgrey", "blue",'Set1','Set2')
)$data
data.usage <-  data.usage %>%
  mutate(
    celltype = id %>% 
      str_remove('_C5aRA_PD1') %>% 
      str_remove('_Mock') %>% 
      str_remove('_PD1') %>% 
      str_remove('_C5aRA')
  ) %>% 
  mutate(
    batch = case_when(
      grepl("_C5aRA_PD1", id) ~ 'C5aRA_PD1',
      grepl("_Mock", id) ~ 'Mock',
      grepl("_PD1", id) ~ 'PD1',
      grepl("_C5aRA", id) ~ 'C5aRA'
    ) %>% factor(.,levels = cell_level)
  ) %>% 
  mutate(celltype = factor(celltype,levels = c(
    "Cytolytic_T cell","Gzmb_T cell","Dysfunction_T cell",
    "Proliferative CD8+ T","CD8+ T memory"
  )))
data.usage <- data.usage %>% 
  rownames_to_column('gene')
ggplot(data.usage,aes(x=batch,y =  features.plot))+
  geom_point(
    mapping = aes(size = pct.exp, color = avg.exp.scaled)
  ) +
  scale_size("% detected", range = c(0,10)) + 
  scale_color_gradientn(
    colours = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    guide = guide_colorbar(
      ticks.colour = "black",
      frame.colour = "black"
    ),
    name = "Average\nexpression"
  ) +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("Markers") + 
  theme_bw() +
  facet_grid(~celltype, scales="free",space = "free")+
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=16, angle=-45, hjust=0,vjust = 0.5, color="black",face="bold"),
    axis.text.y = element_text(size=16, color="black",face="bold"),
    axis.ticks.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.2),
    panel.spacing=unit(0, "mm"), 
    panel.border = element_rect(fill = NA,linetype = 'solid',linewidth = 0.5),
    panel.grid.major.y = element_line(),
    strip.text.x = element_text(size=18, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),
    strip.background = element_rect(colour="grey10", fill="grey40",size = 1)
  ) +
  geom_hline(yintercept = 4.5,size = 0.5) +
  geom_hline(yintercept = 10.5,size = 0.5)

# Figs7 B umap of T cells--------------------------------------------------
scRNA_11$anno_res %>% unique()
scRNA_11_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory"
))
cell_label_correct <- c(
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory"
)
names(cell_label_correct) <- c(
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory"
)
scRNA_11_tmp$anno_tmp <- cell_label_correct[scRNA_11_tmp$anno_res] %>% as.character()
scRNA_11_tmp$anno_tmp %>% unique()
scRNA_11_tmp <-
  FindVariableFeatures(object = scRNA_11_tmp, nfeatures = 2000)
scRNA_11_tmp <- RunPCA(scRNA_11_tmp,
                       features = VariableFeatures(object = scRNA_11_tmp))
scRNA_11_tmp <- RunUMAP(scRNA_11_tmp, dims = 1:16, label = T)
DimPlot(
  object = scRNA_11_tmp,
  group.by = c('anno_tmp'),
  label = TRUE
) + NoLegend()
cell_leves <- cell_label_correct %>% as.character()
label_legend <- paste(
  1:length(cell_leves),
  cell_leves,sep = ': '
)
p <- DimPlot(
  object = scRNA_11_tmp,
  group.by = c('anno_tmp'),
  label = TRUE,
  pt.size = 0.8,
  label.size = 6,
  label.box = TRUE,
  repel = 30,
  alpha = 0.8
)
data_plot <- p$data %>% 
  mutate(anno_tmp = factor(anno_tmp,levels = cell_leves %>% unique())) %>% 
  mutate(label_inner = as.numeric(anno_tmp) %>% as.factor()) %>% 
  mutate(label_legend = paste(label_inner,anno_tmp,sep = ': '))

cell_label_loc <- data_plot %>%
  group_by(anno_tmp) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )
colnames(data_plot)
data_plot$umap_1 %>% min()
data_plot$umap_2 %>% min()
my36colors <-c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) %>% sample(36,replace = FALSE)
ggplot() +
  geom_point(
    data = data_plot,
    aes(x = umap_1,y = umap_2,color = label_inner),
    size = 0.8) +
  geom_label_repel(data = cell_label_loc,
                   aes(
                     x = umap_1,
                     y = umap_2,
                     label = anno_tmp),
                   label.size = 0,
                   size = 8,
                   alpha = 0.8,
                   fill = 'grey90',
                   segment.color = NA,
                   show.legend = FALSE)+
  #x
  geom_segment(aes(
    x = (data_plot$umap_1 %>% min())-2,
    y = (data_plot$umap_2 %>% min())-2,
    xend = (data_plot$umap_1 %>% min())+
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2,
    yend = ((data_plot$umap_2 %>% min())-2)
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.5) +
  #Y
  geom_segment(aes(
    x = (data_plot$umap_1 %>% min())-2,
    y = (data_plot$umap_2 %>% min())-2,
    xend = ((data_plot$umap_1 %>% min())-2),
    yend = (data_plot$umap_2 %>% min()) +
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.5) +
  scale_color_manual(
    values = my36colors
  ) +
  labs(x = 'Umap 1',y= 'Umap 2')+
  guides(color = guide_legend(override.aes = list(size = 5),ncol = 1))+
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = 'none',
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_text(hjust = 0.1,vjust = 7),
    axis.title.y = element_text(vjust = -7,hjust = 0.1)
  )
# Figs7 B Ro/e T cells--------------------------------------------------
scRNA_11$anno_res %>% unique()
scRNA_11_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  "Cancer cell","Fibroblast",
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory","NK",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "S100A8-Neu","Ccl4-Neu","ISG-Neu",
  "Mast_Mcpt2","B cell","cDC1_Xcr1","tDC"
))
cell_label_correct <- c(
  "Cancer cell","Fibroblast",
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory","NK",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "S100A8-Neu","Ccl4-Neu","ISG-Neu",
  "Mast_Mcpt2","B cell","cDC1_Xcr1","tDC"
)
names(cell_label_correct) <- c(
  "Cancer cell","Fibroblast",
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory","NK",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "S100A8-Neu","Ccl4-Neu","ISG-Neu",
  "Mast_Mcpt2","B cell","cDC1_Xcr1","tDC"
)
scRNA_11_tmp$anno_tmp <- cell_label_correct[scRNA_11_tmp$anno_res] %>% as.character()
scRNA_11_tmp$anno_tmp %>% unique()
meta_data_sub <- scRNA_11_tmp@meta.data
meta_data_sub$anno_res  <- meta_data_sub$anno_tmp
celltype_select <- c(
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory"
)

chisq <- chisq.test(table(meta_data_sub$anno_res,meta_data_sub$batch))
data_plot <- chisq$observed/chisq$expected %>% 
  as.matrix() %>% as.data.frame()


data_ggplot <- data_plot %>%
  rownames_to_column('celltype') %>%
  reshape2::melt(id.var = 'celltype',variable.name = 'batch',value.name = 'Ro/e') %>% 
  mutate(
    facet_group = ifelse(batch %in% c('Mock','C5aRA'),'C5aRA vs Mock','C5aRA_PD1 vs PD1')
  ) %>% 
  filter(celltype %in% celltype_select) %>% 
  mutate(
    batch = factor(batch,levels = c('Mock','C5aRA','PD1','C5aRA_PD1')),
    celltype = factor(celltype,levels = celltype_select)
    ) %>% 
  group_by(celltype) %>% 
  mutate(roe_scale = scale(`Ro/e`,center = FALSE)[,1]) 
data_ggplot$roe_scale <- data_ggplot$`Ro/e`
limit_min <- data_ggplot$roe_scale %>% min()
limit_max <- data_ggplot$roe_scale %>% max()

ggplot(data = data_ggplot,
       mapping = aes(x = batch,y = celltype)) +
  geom_tile(aes(fill = roe_scale)) +
  geom_text(aes(label = round(roe_scale,digits = 2))) +
  guides(
    fill = guide_colorbar(title = 'Ro/e',title.vjust = 1)
  ) +
  scale_fill_gradient(low = 'grey90',high =  '#3131F2',limit = c(limit_min,limit_max)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 30,hjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.text = element_text(size = 16),
    legend.position = 'right',
    legend.justification = c(0,1),
    legend.title = element_text(size = 20),
    strip.text = element_blank(),
    strip.background = element_blank()
  )

# Figs7 C GSEA of T cells within AKR-----------------------------------------------------------------
scRNA_T$anno_res_1 %>% unique()
cell_type <- scRNA_T$anno_res_1 %>% unique() %>% .[c(1,3)]
cell_type

## C5aRA_PD1 ---------------------------------------------------------------
scRNA_T_sub <- subset(scRNA_T,subset = anno_res_1%in%cell_type)
Idents(scRNA_T_sub) <- scRNA_T_sub$batch
scRNA_T_sub$batch %>% table()
scRNA_T_tmp <-
  subset(scRNA_T_sub, 
         subset = batch %in% c('C5aRA_PD1', 'PD1')
  )

marker_T_sub <- FindMarkers(
  object = scRNA_T_tmp,
  ident.1 = 'C5aRA_PD1',
  ident.2 = 'PD1',
  logfc.threshold = 0,
  min.pct = 0.25,
  only.pos = FALSE)
cluster0.genes<- marker_T_sub %>%
  arrange(desc(avg_log2FC)) %>%
  rownames_to_column() %>% 
  rename(feature = names(.)[1]) %>% 
  dplyr::select(feature,avg_log2FC) %>% 
  rename(logFC = avg_log2FC)
ranks<- cluster0.genes$logFC
names(ranks) <- cluster0.genes$feature
head(ranks)
geneList <- ranks
gmtfile ='D:/jiaoxi/mh.all.v2023.1.Mm.symbols.gmt'
geneset <- read.gmt(gmtfile)
egmt <- GSEA(geneList, 
             TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
gsea_results_df <- egmt@result %>% 
  arrange(desc(NES)) %>%
  mutate(group = ifelse(test = pvalue<0.05,yes = 'Credible','NS')) %>% 
  mutate(Description = str_replace(Description,'HALLMARK_','') %>% 
           str_split("_") %>% 
           map_chr(~ str_to_title(.x) %>% paste(collapse = " "))) %>% 
  mutate(Description = factor(Description,levels = Description %>% unique() %>% rev()))
gsea_results_df$group %>% table()
gsea_results_df_sub_des_list <- 
  gsea_results_df$Description[gsea_results_df$group == 'Credible']
gsea_results_df_sub <- gsea_results_df %>% 
  filter(group == 'Credible') %>% 
  arrange(desc(NES)) %>% 
  filter(Description %in% gsea_results_df_sub_des_list) %>% 
  mutate(Description = factor(gsea_results_df_sub_des_list))
gsea_results_df_sub$Description
gsea_results_df_sub_des_list <-
  gsea_results_df_sub$Description[c(1,3,4,8,10)]
gsea_results_df_sub <- gsea_results_df %>% 
  filter(group == 'Credible') %>% 
  arrange(desc(NES)) %>% 
  filter(Description %in% gsea_results_df_sub_des_list) %>% 
  mutate(Description = factor(gsea_results_df_sub_des_list))
gsea_plot_self(
  gsea_results_df_sub = gsea_results_df_sub,
  title_plot = "MisgDB_H terms of CD8+ T(C5aRA_PD1)",
  title_vjust = -10,
  title_hjust = -0.1,
  color_name = 'hallmark'
)

## C5aRA -------------------------------------------------------------------


scRNA_T_tmp <-
  subset(scRNA_T_sub, 
         subset = batch %in% c('C5aRA', 'Mock')
  )
marker_T_sub <- FindMarkers(
  object = scRNA_T_tmp,
  ident.1 = 'C5aRA',
  ident.2 = 'Mock',
  logfc.threshold = 0,
  min.pct = 0.25,
  only.pos = FALSE)
cluster0.genes<- marker_T_sub %>%
  arrange(desc(avg_log2FC)) %>%
  rownames_to_column() %>% 
  rename(feature = names(.)[1]) %>% 
  dplyr::select(feature,avg_log2FC) %>% 
  rename(logFC = avg_log2FC)
ranks<- cluster0.genes$logFC
names(ranks) <- cluster0.genes$feature
head(ranks)
geneList <- ranks
gmtfile ='D:/jiaoxi/mh.all.v2023.1.Mm.symbols.gmt'
geneset <- read.gmt(gmtfile)
egmt <- GSEA(geneList, 
             TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
gsea_results_df <- egmt@result %>% 
  arrange(desc(NES)) %>%
  mutate(group = ifelse(test = pvalue<0.05,yes = 'Credible','NS')) %>% 
  mutate(Description = str_replace(Description,'HALLMARK_','') %>% 
           str_split("_") %>% 
           map_chr(~ str_to_title(.x) %>% paste(collapse = " "))) %>% 
  mutate(Description = factor(Description,levels = Description %>% unique() %>% rev()))
gsea_results_df$group %>% table()
gsea_results_df_sub_des_list <- 
  gsea_results_df$Description[gsea_results_df$group == 'Credible']
gsea_results_df_sub <- gsea_results_df %>% 
  filter(group == 'Credible') %>% 
  arrange(desc(NES)) %>% 
  filter(Description %in% gsea_results_df_sub_des_list) %>% 
  mutate(Description = factor(gsea_results_df_sub_des_list))
gsea_results_df_sub$Description
gsea_results_df_sub_des_list <-
  gsea_results_df_sub$Description[c(1,2,3,4,7)]
gsea_results_df_sub <- gsea_results_df %>% 
  filter(group == 'Credible') %>% 
  arrange(desc(NES)) %>% 
  filter(Description %in% gsea_results_df_sub_des_list) %>% 
  mutate(Description = factor(gsea_results_df_sub_des_list))
gsea_plot_self(
  gsea_results_df_sub = gsea_results_df_sub,
  title_plot = "MisgDB_H terms of CD8+ T(C5aRA)",
  title_vjust = -10,
  title_hjust = -0.1,
  color_name = 'hallmark'
)
# Figs7 E TCR clone score-----------------------------------------------------------------
library(entropy)
library(immunarch)
library(ggplot2)
library(Seurat)
library(dplyr)
library(stringr)
library(ggpubr)
file_path <- 'D:/jiaoxi_12/TCR/immunarch_PD1'
immdata_10x <- repLoad(file_path)
scRNA_T <- readRDS('D:/jiaoxi_12/scRNA_T_TCR.rds')
meta_data_T <- scRNA_T@meta.data
colnames(meta_data_T)
tmp <- meta_data_T %>% 
  mutate(
    anno_res_1 = ifelse(
      test = anno_res_1 %in% c("Proliferative CD4+ T",'Proliferative CD8+ T','Treg_Proliferating'),
      yes = 'Proliferative T',
      no = anno_res_1
    )
  ) %>% 
  filter(!is.na(t_cdr3s_aa)) %>% 
  group_by(batch,anno_res_1) %>% 
  summarise(cellnum = n()) %>% 
  filter(cellnum>5) %>% 
  mutate(group = paste(batch,anno_res_1,sep = '_'))
tmp
cell_type <- meta_data_T$anno_res_1 %>% unique()
cell_type
cell_type_select <- cell_type#[-c(4,8)]
cell_type_select
data_tcr <- meta_data_T %>%
  filter(anno_res_1 %in% cell_type_select) %>%
  mutate(
    anno_res_1 = ifelse(
      test = anno_res_1 %in% c("Proliferative CD4+ T",'Proliferative CD8+ T','Treg_Proliferating'),
      yes = 'Proliferative T',
      no = anno_res_1
    )
  ) %>%
  filter(!is.na(t_cdr3s_aa)) %>%
  group_by(batch, anno_res_1, t_cdr3s_aa) %>%
  summarise(Clones = n()) %>%
  group_by(batch, anno_res_1) %>%
  mutate(p = Clones / sum(Clones),
         N = sum(Clones)) %>%
  summarise(clonality = 1 + sum(p * log(p)) / N[1]) %>%
  ungroup() %>%
  mutate(
    batch = factor(batch, levels = c("Mock","C5aRA","PD1","C5aRA_PD1")),
    group = paste(batch,anno_res_1,sep = '_')
  ) %>%
  filter(group %in% tmp$group)


data_tcr$batch %>% unique()
data_tcr$anno_res_1 %>% unique()
data_tcr_sub <- data_tcr %>% 
  filter(batch %in% c("Mock","C5aRA")) %>% 
  group_by(anno_res_1) %>% 
  filter(n() == 2) %>% 
  filter(!(anno_res_1 %in% c('Dysfunction_T cell')))
data_tcr_sub$anno_res_1 %>% unique()
p1 <- ggplot(data_tcr_sub, aes(x = batch, y = clonality)) +
  geom_boxplot(
    aes(fill = batch), 
    show.legend = F,
    width = 0.6) +  
  scale_fill_manual(values = c('#00468B', '#ED0000')) + 
  geom_point(size = 3,fill='#374E54',shape=21) +
  geom_line(
    aes(group = anno_res_1), 
    color = 'gray',
    lwd = 0.5) +  
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(size=16, angle=-45, hjust=0, color="black",face="bold"),
    axis.text.y = element_text(size=25, color="black",face="bold"),
    axis.title = element_text(size=25,colour = 'black',vjust = -0.8,hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.y.right = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.2)
  ) +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,
    label.x.npc = 0.5,
    label = 'p.format',
    fontface = 'bold'
  ) + 
  labs(title = 'C5aRa  vs Mock')
p1

# Figs7 E Monocle3 of T cell within AKR -----------------------------------
scRNA_t_CD8_rep <- readRDS('D:/jiaoxi_12/scRNA_t_CD8_rep.rds')
data <- GetAssayData(scRNA_t_CD8_rep, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNA_t_CD8_rep@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds,num_dim = 70)
cds <- reduce_dimension(
  cds,
  preprocess_method = 'PCA',
  max_components = 2,
  umap.n_neighbors = 5)
plot_cells(
  cds = cds,
  color_cells_by = 'predicted.id',
  label_groups_by_cluster = TRUE,
  label_cell_groups = TRUE,
  group_cells_by = 'cluster',
  group_label_size = 5,
  cell_size = 1,
  reduction_method = 'UMAP'
)
colnames(colData(cds))
cds <- cluster_cells(
  cds = cds,
  resolution=1e-5,
  k = 10
)
plot_cells(
  cds = cds,
  group_label_size = 10,
  cell_size = 1
)

cds <- learn_graph(
  cds = cds,
  use_partition = TRUE,
  learn_graph_control = list(
    minimal_branch_len = 6
  )
)
my36colors <-c('#14517c', "#96C37D", 
               "#D8383A","#BD956A", '#585658')
plot_cells(
  cds,
  color_cells_by = 'anno_pseudo',
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 6,
  cell_size = 1.4,
  trajectory_graph_segment_size = 1.4
) +
  scale_color_manual(values = my36colors) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 25),
    axis.text = element_blank(),
    axis.title = element_text(size = 20)
  )
cds <- order_cells(cds = cds)
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 10,
  cell_size = 1,
  reduction_method = 'UMAP',
  label_groups_by_cluster = FALSE
)

# endline -----------------------------------------------------------------
