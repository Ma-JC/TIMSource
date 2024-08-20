rm(list = ls())
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
# Fig3 J(left) umap of AKR ----------------------------------------------------
scRNA_12 <- readRDS('D:/jiaoxi_12/scRNA_12_fig_use.rds')
scRNA_12_sub <- scRNA_12
meta_data_sub <- scRNA_12_sub@meta.data
colnames(meta_data_sub)
meta_data_sub$anno_res_1 %>% unique()
scRNA_12_sub$anno_res_low <- scRNA_12_sub$anno_res_1
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Mac_C1qa","Mac_Mki67/C1qa","Mac_Vegfa")] = 
  'Mac/Mono Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Treg", "CD4+ T","Proliferative CD4+ T",
                              "Proliferative Treg")] = 
  'T Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Cytolytic/Gzmb_T cell","Proliferative CD8+ T",
                              "CD8+ T memory","Dysfunction_T cell")] = 
  'T Cells'

scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("MonoMac_Ccr2","Mono_Isg15")] =  
  'Mac/Mono Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("cDC","pDC","tDC")] =  'Dendritic Cells (DCs)'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Naive T")] = 
  'T Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("NK","Proliferative NK")] = 
  'Natural Killer (NK) Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Neutrophil")] = 
  'Neutrophil Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Mast")] = 
  'Mast Cells'
scRNA_12_sub$anno_res_low %>% unique()
scRNA_12_sub <- subset(scRNA_12_sub,subset = anno_res_low %in% c(
  "T Cells","Natural Killer (NK) Cells","Mac/Mono Cells",
  "Dendritic Cells (DCs)","Neutrophil Cells","Mast Cells"
))
meta_data <- scRNA_12_sub@meta.data %>% 
  select(-c("umap_1","umap_2"))
scRNA_12_sub@meta.data <- meta_data
scRNA_12_sub$anno_res_low %>% unique()
scRNA_12_sub <-
  FindVariableFeatures(object = scRNA_12_sub, nfeatures = 3000)
scRNA_12_sub <- RunPCA(scRNA_12_sub,
                       features = VariableFeatures(object = scRNA_12_sub))
scRNA_12_sub <- RunUMAP(scRNA_12_sub, dims = 1:15, label = T)
cell_leves <- c(
  "T Cells","Natural Killer (NK) Cells","Mac/Mono Cells",
  "Dendritic Cells (DCs)","Neutrophil Cells","Mast Cells"
)
cell_label_loc <- data_plot %>%
  group_by(anno_res_low) %>%
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
                     label = anno_res_low),
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
    values = my7colors
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
# Fig3 J(right) Ro/e of AKR ----------------------------------------------------
scRNA_12 <- readRDS('D:/jiaoxi_12/scRNA_12_fig_use.rds')
scRNA_12_sub <- scRNA_12
meta_data_sub <- scRNA_12_sub@meta.data
colnames(meta_data_sub)
meta_data_sub$anno_res_1 %>% unique()
scRNA_12_sub$anno_res_low <- scRNA_12_sub$anno_res_1
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Mac_C1qa","Mac_Mki67/C1qa","Mac_Vegfa")] = 
  'Mac/Mono Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Treg", "CD4+ T","Proliferative CD4+ T",
                              "Proliferative Treg")] = 
  'T Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Cytolytic/Gzmb_T cell","Proliferative CD8+ T",
                              "CD8+ T memory","Dysfunction_T cell")] = 
  'T Cells'

scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("MonoMac_Ccr2","Mono_Isg15")] =  
  'Mac/Mono Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("cDC","pDC","tDC")] =  'Dendritic Cells (DCs)'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Naive T")] = 
  'T Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("NK","Proliferative NK")] = 
  'Natural Killer (NK) Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Neutrophil")] = 
  'Neutrophil Cells'
scRNA_12_sub$anno_res_low[scRNA_12_sub$anno_res_1 %in% 
                            c("Mast")] = 
  'Mast Cells'
scRNA_12_sub$anno_res_low %>% unique()
scRNA_12_sub <- subset(scRNA_12_sub,subset = anno_res_low %in% c(
  "T Cells","Natural Killer (NK) Cells","Mac/Mono Cells",
  "Neutrophil Cells","Dendritic Cells (DCs)","Mast Cells"
))
meta_data <- scRNA_12_sub@meta.data %>% 
  dplyr::select(-c("umap_1","umap_2"))
scRNA_12_sub@meta.data <- meta_data
scRNA_12_sub$anno_res_low %>% unique()
cell_levels <- c(
  "T Cells","Natural Killer (NK) Cells","Mac/Mono Cells",
  "Neutrophil Cells","Dendritic Cells (DCs)","Mast Cells"
)
meta_data_sub <- scRNA_12_sub@meta.data
meta_data_sub$anno_res_low %>% unique()
meta_data_sub$batch %>% unique()
meta_data_sub$anno_res_low %>% unique()
meta_data_sub$anno_res <- meta_data_sub$anno_res_low
tmp <- table(meta_data_sub$anno_res,meta_data_sub$batch) %>%
  as.data.frame() %>% 
  pivot_wider(id_cols = 'Var1',names_from = 'Var2',values_from = 'Freq') %>% 
  as.data.frame() %>% 
  column_to_rownames('Var1')
celltype_select <- c(
  "T Cells","Natural Killer (NK) Cells","Mac/Mono Cells",
  "Neutrophil Cells","Dendritic Cells (DCs)","Mast Cells"
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
# Fig4 D(right) TCR clone proportion------------------------------------------------------------------
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
cell_type_select <- cell_type
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
  filter(batch %in% c("PD1","C5aRA_PD1")) %>% 
  group_by(anno_res_1) %>% 
  filter(n() == 2) 
p2 <- ggplot(data_tcr_sub, aes(x = batch, y = clonality)) +
  geom_boxplot(
    aes(fill = batch), 
    show.legend = F,
    width = 0.6) +  
  scale_fill_manual(values =  c('#42B540','#0099B4')) + 
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
    axis.line = element_line(colour = 'grey30',size = 0.2), 
  ) +
  stat_compare_means(
    method = "wilcox.test",
    paired = TRUE,
    label.x.npc = 0.5, 
    label = 'p.format',
    fontface = 'bold'
  ) +
  labs(title = 'C5aRa_PD1  vs PD1')
p2
# Fig4 F(right)TCR clone score ----------------------------------------------------


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
cell_type_select <- cell_type
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
  filter(batch %in% c("PD1","C5aRA_PD1")) %>% 
  group_by(anno_res_1) %>% 
  filter(n() == 2) 
p2 <- ggplot(data_tcr_sub, aes(x = batch, y = clonality)) +
  geom_boxplot(
    aes(fill = batch), 
    show.legend = F,
    width = 0.6) +  
  scale_fill_manual(values =  c('#42B540','#0099B4')) + 
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
    axis.line = element_line(colour = 'grey30',size = 0.2), 
  ) +
  stat_compare_means(
    method = "wilcox.test",
    paired = TRUE,
    label.x.npc = 0.5, 
    label = 'p.format',
    fontface = 'bold'
  ) +  
  labs(title = 'C5aRa_PD1  vs PD1')
p2