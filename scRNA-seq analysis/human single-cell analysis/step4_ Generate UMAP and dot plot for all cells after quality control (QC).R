library(Seurat)
library(stringr)
library(dplyr)
library(tibble)
library(harmony)
library(tibble)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
setwd('/mnt/data3/jiaoxi/新辅助治疗队列/')
scdata <- readRDS('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/scdata_all_anno_lv1.rds')
# scdata_res <- readRDS('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/scdata_all_anno_lv1.rds')
# identical(scdata_res@meta.data %>% colnames(),scdata_use@meta.data %>% colnames())
meta_data <- read.csv('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/meta_data_macmono-nkt.csv',row.names=1)
meta_data <- meta_data %>% 
  select(-c('anno_res_data2','anno_res_data2_data2','anno_res_data2_data2_data2'))
# scdata@meta.data <- meta_data
# colnames(scdata_res@meta.data)
identical(colnames(scdata),rownames(meta_data))
scdata@meta.data <- meta_data
scdata$anno_res_lv2 %>% unique()
celltype <- scdata$anno_res_lv2 %>% unique() %>% .[!is.na(.)]
meta_data <- scdata@meta.data %>% 
  mutate(
    anno_res_lv1_modif = case_when(
      anno_res_lv2 %in% c(
        'Cd8_T_Proliferative','Cd8_T_FOS','Cd8_T_GZMK','Cd8_T_METRNL',
        'Cd8_T_HOPX','Cd8_T_CXCL13'
      ) ~ 'CD8+ T',
      anno_res_lv2 %in% c(
        'Cd4_Tmem_IL7R','Cd4_Tmem_FOS','Cd4_FOXP3','Cd4_T_CXCL13',
        'Cd4_KLRB1','Cd4_T_Proliferative'
      ) ~ 'CD4+ T',
      anno_res_lv2 %in% c('NK_NCAM1','NK_FCGR3A') ~ 'NK',
      anno_res_lv2 %in% c('pDC-LILRA4') ~ 'pDC-LILRA4',
      anno_res_lv2 %in% c('cDC') ~ 'cDC',
      anno_res_lv2 %in% c(
        'Mac_RPL/RPS','Mac_MSR1','Macro_prolif',
        'Mac_RGS1','Macro_prolif','Mac_FLT1','Mac_RGS1','Mac_IL7R'
      ) ~ 'Macro',
      anno_res_lv2 %in% c('Mono_FCN1','Mono_FCGR3B') ~ 'Mono',
      NA~NA,
      TRUE ~ anno_res_lv1
    )
  )
meta_data$anno_res_lv1_modif %>% unique()
meta_data$anno_res_lv2[meta_data$anno_res_lv1_modif == 'Macrophage/Monocyte'] %>% 
  unique()
scdata@meta.data <- meta_data
scdata_use <- subset(scdata,subset = anno_res_lv1_modif != c('NK/T cell'))
scdata_use <- subset(scdata_use,subset = anno_res_lv1_modif != c('Macrophage/Monocyte'))
## Umap --------------------------------------------------------------------
p <- DimPlot(
  object = scdata_use,
  pt.size = 1,
  label.size = 6,
  group.by = c('anno_res_lv1_modif'),
  label = FALSE
) + 
  labs(title = 'Cell Annotation of all cells')
p

# getwd()
# dir.create('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs_0814')
## anno_res_lv1_modif ------------------------------------------------------------


my36colors <-c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) %>% sample(36,replace = FALSE)

meta_data_sub <- scdata_use@meta.data
meta_data_sub[,c('umap_1','umap_2')] <- 
  scdata_use@reductions$umap@cell.embeddings
cell_leves <- c(
  'Epithelial/Malignant cell','Fibroblast','Vascular Smooth Muscle cell',
  'Endothelial cell',
  'CD8+ T','CD4+ T',"NK",'Macro','Mono',"Neutrophil cell",
  'B cell','Plasma cell','Mast cell','pDC-LILRA4','cDC'
)

label_legend <- paste(
  1:length(cell_leves),
  cell_leves,sep = ': '
)
data_plot <- p$data %>% 
  mutate(anno_res_lv1_modif = factor(anno_res_lv1_modif,levels = cell_leves)) %>% 
  mutate(label_inner = as.numeric(anno_res_lv1_modif) %>% as.factor()) %>% 
  mutate(label_legend = paste(label_inner,anno_res_lv1_modif,sep = ': '))
# #筛选出来部分离群点
# filter(umap_1 > -10)
cell_label_loc <- data_plot %>%
  group_by(label_inner) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )
colnames(data_plot)
data_plot$umap_1 %>% min()
data_plot$umap_2 %>% min()
res_p <- ggplot() +
  geom_point(
    data = data_plot,
    aes(x = umap_1,y = umap_2,color = label_inner),
    size = 0.8) +
  geom_label_repel(data = cell_label_loc,
                   aes(
                     x = umap_1,
                     y = umap_2,
                     label = label_inner),
                   fill = NA,
                   label.size = 0,size = 8,
                   segment.color = NA,
                   show.legend = FALSE)+
  scale_color_manual(
    values = my36colors,
    breaks = 1:length(cell_leves),
    label = label_legend
  ) +
  labs(x = 'Umap 1',y= 'Umap 2')+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.border = element_blank(),
    axis.line =  element_line(arrow = arrow(length = unit(0.3, "cm"))),
    # axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 16)
    # axis.title.x = element_text(hjust = 0.1,vjust = 7),
    # axis.title.y = element_text(vjust = -7,hjust = 0.1)
  )
res_p
ggsave('./result_figs_0814/allcell_umap.png',width = 9,height = 6)
ggsave('./result_figs_0814/allcell_umap.pdf',width = 9,height = 6)

#### gene list ---------------------------------------------------------------
cell_leves <- c(
  'Epithelial/Malignant cell','Fibroblast','Vascular Smooth Muscle cell',
  'Endothelial cell',
  'CD4+ T','CD8+ T',"NK",'Macro','Mono',"Neutrophil cell",
  'B cell','Plasma cell','Mast cell','pDC-LILRA4','cDC'
)
scdata_use$anno_res_lv1_modif <- factor(
  scdata_use$anno_res_lv1_modif,
  levels = cell_leves
)
gene_select <- list(
  'Epithelial/Malignant cell' = c('FXYD3','S100A14','PERP','SFN'),
  'Fibroblast' = c('MMP2','RARRES2','LUM','DCN'),
  'Vascular Smooth Muscle cell' = c('CALD1','MYL9','RGS5','IGFBP7'),
  'Endothelial cell' = c('PECAM1','VWF'),#"SPARCL1","A2M",
  'CD4+ T' = c('CD3D','CD3E','CD4'),
  'CD8+ T' = c('CD8A'),
  'NK' = c('KIR2DL4','KLRC1','FGFBP2','FCGR3A'),
  'Macrophage' = c('LYZ','CD68'),
  'Monocyte' = c('CD14','AIF1'),
  'Neutrophil cell' = c(,'CXCR2','CSF3R'),#'S100A8','S100A9'
  'B cell' = c('CD37','MS4A1','BANK1','CD79A'),
  'Plasma cell' = c('MZB1','JCHAIN','DERL3'),
  'Mast cell' = c("TPSB2",'TPSAB1'),
  'pDC-LILRA4' = c('BST2','LILRA4'),
  'cDC' = c('XCR1','CLEC9A')
)
gene_select %>% unlist()
#### dotplot ====
DotPlot(object = scdata_use,
        features = gene_select %>% unlist() %>% unique(),
        group.by = 'anno_res_lv1_modif') +
  scale_color_gradient(
    name = 'Average Expression',
    low = 'grey90',
    high = '#3131F2',
    breaks = c(-1,0,1),
    labels = c(-1,0,1)) +
  scale_size_continuous(name = 'Percent Expressed') +
  guides(
    color = guide_colourbar(nrow = 1,title.position = 'top'),
    size = guide_legend(ncol = 1,title.position = 'top')
  ) +
  labs(x = '', y= '') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0),
    axis.text = element_text(size = 18),
    legend.direction = 'vertical',
    axis.title.x = element_blank()
  )
ggsave('./result_figs_0814/allcell_dotplot.png',width = 22,height = 6)
ggsave('./result_figs_0814/allcell_dotplot.pdf',width = 22,height = 6)
### boxplot ----------------------------------------------------------------
library(ggpubr)
`%noin%` <- Negate(`%in%`)
meta_data <- scdata_use@meta.data 
colnames(meta_data)
data_plot <- meta_data %>%
  dplyr::rename(Tumor_Sample_Barcode = orig.ident) %>% 
  dplyr::select(c("Tumor_Sample_Barcode","group_maf","anno_res_lv2")) %>%
  filter(
    anno_res_lv2 %noin% c(
      'Epithelial/Malignant cell','Fibroblast','Vascular Smooth Muscle cell',
      'Endothelial cell'
    )
  ) %>% 
  group_by(group_maf,Tumor_Sample_Barcode,anno_res_lv2) %>%
  summarise(
    cell_num = n(),.groups = 'drop'
  ) %>%
  group_by(group_maf,Tumor_Sample_Barcode) %>%
  mutate(
    percent = (cell_num/sum(cell_num)) * 100
  ) %>% 
  filter(
    anno_res_lv2 %in% c(
      'Cd8_T_FOS', 'Cd8_T_GZMK', 'Cd8_T_CXCL13', 'Cd8_T_METRNL', 'Cd8_T_HOPX','NK_NCAM1','NK_FCGR3A',
      "cDC","pDC-LILRA4"
    )
  ) %>% 
  mutate(
    anno_res_lv2 = factor(anno_res_lv2,levels = c(
      'Cd8_T_FOS', 'Cd8_T_GZMK', 'Cd8_T_CXCL13', 'Cd8_T_METRNL', 'Cd8_T_HOPX','NK_NCAM1','NK_FCGR3A',
      "cDC","pDC-LILRA4"
    )),
    group_maf = factor(group_maf,levels = c("complement-WT","complement-MUT"))
  )
colnames(data_plot)
colour=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
  "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
  "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
  "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887"
)

colnames(data_plot)
data_plot$group_maf %>% unique()
ggboxplot(
  data_plot, x="group_maf", y="percent",
  color ="group_maf",
  width = 0.6,
  palette = colour[c(2,1)],add = "jitter",
  xlab = F,  bxp.errorbar=T,
  bxp.errorbar.width=0.5, #误差条大小
  size=0.5, outlier.shape=NA,legend = "right") +
  # ggrepel::geom_text_repel(aes(label = Tumor_Sample_Barcode)) +
  guides(color = guide_legend(title = 'Group of maf'))+
  # scale_y_continuous(limits = c(-0.5,0.7)) +
  stat_compare_means(
    aes(group = group_maf),
    label = "p.format",
    comparisons = list(c(
      "complement-MUT","complement-WT"
    )),
    method = "t.test",
  ) +
  facet_wrap(~anno_res_lv2,nrow = 1,scales = 'free',switch = 'y') +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10,face = 'italic',vjust = 0.7),
    strip.placement = 'outside'
  )
ggsave('./result_figs_0814/allcell_boxplot.png',width = 16,height = 4.5)
ggsave('./result_figs_0814/allcell_boxplot.pdf',width = 16,height = 4.5)

### malignant cells ----------------------------------------------------------------
library(ggpubr)
`%noin%` <- Negate(`%in%`)
meta_data <- scdata_use@meta.data 
colnames(meta_data)
data_plot <- meta_data %>%
  dplyr::rename(Tumor_Sample_Barcode = orig.ident) %>% 
  dplyr::select(c("Tumor_Sample_Barcode","group_maf","anno_res_lv2")) %>%
  group_by(group_maf,Tumor_Sample_Barcode,anno_res_lv2) %>%
  summarise(
    cell_num = n(),.groups = 'drop'
  ) %>%
  group_by(group_maf,Tumor_Sample_Barcode) %>%
  mutate(
    percent = (cell_num/sum(cell_num)) * 100
  ) %>% 
  filter(
    anno_res_lv2 %in% c(
      # 'Epithelial/Malignant cell',
      'Fibroblast'
      # 'Vascular Smooth Muscle cell',
      # 'Endothelial cell'
    )
  ) %>% 
  mutate(
    group_maf = factor(group_maf,levels = c("complement-WT","complement-MUT"))
  )
colnames(data_plot)
colour=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
  "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
  "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
  "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887"
)

colnames(data_plot)
data_plot$group_maf %>% unique()
ggboxplot(
  data_plot, x="group_maf", y="percent",
  color ="group_maf",
  width = 0.6,
  palette = colour[c(2,1)],add = "jitter",
  xlab = F,  bxp.errorbar=T,
  bxp.errorbar.width=0.5, #误差条大小
  size=0.5, outlier.shape=NA,legend = "right") +
  # ggrepel::geom_text_repel(aes(label = Tumor_Sample_Barcode)) +
  guides(color = guide_legend(title = 'Group of maf'))+
  # scale_y_continuous(limits = c(-0.5,0.7)) +
  stat_compare_means(
    aes(group = group_maf),
    label = "p.format",
    comparisons = list(c(
      "complement-MUT","complement-WT"
    )),
    method = "t.test",
  ) +
  facet_wrap(~anno_res_lv2,nrow = 1,scales = 'free',switch = 'y') +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10,face = 'italic',vjust = 0.7),
    strip.placement = 'outside'
  )
ggsave('./result_figs_0814/caf_boxplot.png',width = 4,height = 4.5)
ggsave('./result_figs_0814/caf_boxplot.pdf',width = 4,height = 4.5)
### boxplot allcells ----------------------------------------------------------------
library(ggpubr)
`%noin%` <- Negate(`%in%`)
meta_data <- scdata_use@meta.data 
colnames(meta_data)
data_plot <- meta_data %>%
  dplyr::rename(Tumor_Sample_Barcode = orig.ident) %>% 
  dplyr::select(c("Tumor_Sample_Barcode","group_maf","anno_res_lv2")) %>%
  filter(
    anno_res_lv2 %noin% c(
      'Epithelial/Malignant cell','Fibroblast','Vascular Smooth Muscle cell',
      'Endothelial cell'
    )
  ) %>% 
  group_by(group_maf,Tumor_Sample_Barcode,anno_res_lv2) %>%
  summarise(
    cell_num = n(),.groups = 'drop'
  ) %>%
  group_by(group_maf,Tumor_Sample_Barcode) %>%
  mutate(
    percent = (cell_num/sum(cell_num)) * 100
  ) %>% 
  mutate(
    anno_res_lv2 = factor(anno_res_lv2,levels = c(
      "Cd4_FOXP3","Cd4_T_CXCL13","Cd4_Tmem_IL7R","Cd4_KLRB1","Cd4_Tmem_FOS",
      "Cd4_T_Proliferative",'Cd8_T_Proliferative',"Cd8_T_FOS","Cd8_T_GZMK","Cd8_T_CXCL13","Cd8_T_METRNL","Cd8_T_HOPX",
      "NK_NCAM1","NK_FCGR3A",
      'Mac_MSR1','Mac_RGS1',
      'Mac_IL7R','Mac_FLT1','Macro_prolif',
      'Mono_FCN1','Mono_FCGR3B','Neutrophil cell',
      'B cell','Plasma cell','Mast cell','pDC-LILRA4','cDC'
    )),
    group_maf = factor(group_maf,levels = c("complement-WT","complement-MUT"))
  )
colnames(data_plot)
colour=c(
  "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
  "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
  "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
  "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887"
)

colnames(data_plot)
data_plot$group_maf %>% unique()
ggboxplot(
  data_plot, x="group_maf", y="percent",
  color ="group_maf",
  width = 0.6,
  palette = colour[c(2,1)],add = "jitter",
  xlab = F,  bxp.errorbar=T,
  bxp.errorbar.width=0.5, #误差条大小
  size=0.5, outlier.shape=NA,legend = "right") +
  # ggrepel::geom_text_repel(aes(label = Tumor_Sample_Barcode)) +
  guides(color = guide_legend(title = 'Group of maf'))+
  # scale_y_continuous(limits = c(-0.5,0.7)) +
  stat_compare_means(
    aes(group = group_maf),
    label = "p.format",
    comparisons = list(c(
      "complement-MUT","complement-WT"
    )),
    method = "t.test",
  ) +
  facet_wrap(~anno_res_lv2,nrow = 2,scales = 'free',switch = 'y') +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10,face = 'italic',vjust = 0.7),
    strip.placement = 'outside'
  )
ggsave('./result_figs_0814/allcell_boxplot_all.png',width = 19,height = 6.5)
ggsave('./result_figs_0814/allcell_boxplot_all.pdf',width = 19,height = 6.5)



## anno_res_lv2 ------------------------------------------------------------
scdata_use$anno_res_lv2 %>% unique()
scdata_use <- subset(scdata_use,anno_res_lv2 != 'Mac_RPL/RPS')
my36colors <-c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) %>% sample(36,replace = FALSE)

meta_data_sub <- scdata_use@meta.data
meta_data_sub[,c('umap_1','umap_2')] <- 
  scdata_use@reductions$umap@cell.embeddings
cell_leves <- c(
  'Epithelial/Malignant cell','Fibroblast','Vascular Smooth Muscle cell',
  'Endothelial cell',
  "Cd4_FOXP3","Cd4_T_CXCL13","Cd4_Tmem_IL7R","Cd4_KLRB1","Cd4_Tmem_FOS",
  "Cd4_T_Proliferative",'Cd8_T_Proliferative',"Cd8_T_FOS","Cd8_T_GZMK","Cd8_T_CXCL13","Cd8_T_METRNL","Cd8_T_HOPX",
  "NK_NCAM1","NK_FCGR3A",
  'Mac_MSR1','Mac_RGS1',
  # 'Mac_RPL/RPS',
  'Mac_IL7R','Mac_FLT1','Macro_prolif',
  'Mono_FCN1','Mono_FCGR3B','Neutrophil cell',
  'B cell','Plasma cell','Mast cell','pDC-LILRA4','cDC'
)

label_legend <- paste(
  1:length(cell_leves),
  cell_leves,sep = ': '
)
scdata_tmp <- subset(scdata_use,group_maf == "complement-WT")

p <- DimPlot(
  object = scdata_use,
  pt.size = 1,
  label.size = 6,
  group.by = c('anno_res_lv2'),
  label = FALSE
) + 
  labs(title = 'Cell Annotation of all cells')
data_plot <- p$data %>% 
  mutate(anno_res_lv2 = factor(anno_res_lv2,levels = cell_leves)) %>% 
  mutate(label_inner = as.numeric(anno_res_lv2) %>% as.factor()) %>% 
  mutate(label_legend = paste(label_inner,anno_res_lv2,sep = ': '))
# filter(umap_1 > -10)
cell_label_loc <- data_plot %>%
  group_by(label_inner) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )
colnames(data_plot)
data_plot$umap_1 %>% min()
data_plot$umap_2 %>% min()
res_p <- ggplot() +
  geom_point(
    data = data_plot,
    aes(x = umap_1,y = umap_2,color = label_inner),
    size = 0.8) +
  geom_label_repel(data = cell_label_loc,
                   aes(
                     x = umap_1,
                     y = umap_2,
                     label = label_inner),
                   fill = NA,
                   label.size = NA,size = 8,max.overlaps = 1000,
                   segment.color = NA,
                   show.legend = FALSE)+
  scale_color_manual(
    values = my36colors,
    breaks = 1:length(cell_leves),
    label = label_legend
  ) +
  labs(x = 'Umap 1',y= 'Umap 2',title = 'UMAP Plot of the Complement Wild Type (WT) Group')+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.border = element_blank(),
    axis.line =  element_line(arrow = arrow(length = unit(0.3, "cm"))),
    # axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18)
    # axis.title.x = element_text(hjust = 0.1,vjust = 7),
    # axis.title.y = element_text(vjust = -7,hjust = 0.1)
  )
res_p
ggsave('./result_figs_0814/allcell_umap_anno_lv2.png',width = 11,height = 6)
ggsave('./result_figs_0814/allcell_umap_anno_lv2.pdf',width = 11,height = 6)
# ggsave('./result_figs_0814/allcell_umap_anno_lv2_complement-WT.png',width = 11,height = 6)
# ggsave('./result_figs_0814/allcell_umap_anno_lv2_complement-WT.pdf',width = 11,height = 6)
#### gene list ---------------------------------------------------------------
cell_leves <- c(
  'Epithelial/Malignant cell','Fibroblast','Vascular Smooth Muscle cell',
  'Endothelial cell',
  "Cd4_FOXP3","Cd4_T_CXCL13","Cd4_Tmem_IL7R","Cd4_KLRB1","Cd4_Tmem_FOS",
  "Cd4_T_Proliferative",'Cd8_T_Proliferative',"Cd8_T_FOS","Cd8_T_GZMK","Cd8_T_CXCL13","Cd8_T_METRNL","Cd8_T_HOPX",
  "NK_NCAM1","NK_FCGR3A",
  'Mac_MSR1','Mac_RGS1',
  # 'Mac_RPL/RPS',
  'Mac_IL7R','Mac_FLT1','Macro_prolif',
  'Mono_FCN1','Mono_FCGR3B','Neutrophil cell',
  'B cell','Plasma cell','Mast cell','pDC-LILRA4','cDC'
)
scdata_use$anno_res_lv2 <- factor(
  scdata_use$anno_res_lv2,
  levels = cell_leves
)
table(is.na(scdata_use$anno_res_lv2))
gene_select <- list(
  'Epithelial/Malignant cell' = c('FXYD3','S100A14','PERP','SFN'),
  'Fibroblast' = c('MMP2','RARRES2','LUM','DCN'),
  'Vascular Smooth Muscle cell' = c('CALD1','MYL9','RGS5','IGFBP7'),
  'Endothelial cell' = c('PECAM1','VWF'),#"SPARCL1","A2M",
  "Cd4_FOXP3" = c('FOXP3','CTLA4','ICOS','BATF'),
  "Cd4_T_CXCL13" = c('PDCD1','NR3C1','CXCL13'),
  "Cd4_Tmem_IL7R" = c('IL7R','CCR7','ANXA1'),
  "Cd4_KLRB1" = c('KLRB1','JAML','FKBP11','SYTL2'),
  "Cd4_Tmem_FOS" = c('FOS','FOSB','JUN','CD69'),
  "Cd4_T_Proliferative" = c('STMN1','MCM7','CLSPN','HELLS'),
  "Cd8_T_Proliferative" = c('PCNA','MCM7','MCM5','CLSPN'),
  "Cd8_T_FOS" = c('FOSB','FOS','UBE2S'),
  "Cd8_T_GZMK" = c('GZMK','GZMM'),
  "Cd8_T_CXCL13" = c('CXCL13','GNLY','SOX4'),
  "Cd8_T_METRNL" = c('METRNL','VCAM1','HAVCR2','LAG3'),
  "Cd8_T_HOPX" = c('HOPX',"CAPG","GPR65"),
  "NK_NCAM1" = c('FCER1G','KIR2DL4','KLRC1'),
  "NK_FCGR3A" = c('GZMH','KLF2','FGFBP2','FCGR3A'),
  'Mac_MSR1' = c('MSR1','C2','STAB1','GPR34','CD84'),
  'Mac_RGS1' = c('PTGER4','GPR183','PRDM1','RGS1'),
  # 'Mac_RPL/RPS' = c("RPS18","RPS11","RPL30","RPS23"),
  'Mac_IL7R' = c('ANPEP','MMP9','MMP12','IL7R'),
  'Mac_FLT1' = c('TNIP3','INHBA','FLT1','CCL20'),
  'Macro_prolif' = c('PCLAF','TYMS','PCNA','MKI67'),
  'Mono_FCN1' = c('FCN1','VCAN','S100A12','LYZ'),
  'Mono_FCGR3B' = c('FCGR3B','CXCR2','CXCR1','ADGRG3'),
  'Neutrophil cell' = c('CXCR2','CSF3R'),#'S100A8','S100A9'
  'B cell' = c('CD37','MS4A1','BANK1','CD79A'),
  'Plasma cell' = c('MZB1','JCHAIN','DERL3'),
  'Mast cell' = c("TPSB2",'TPSAB1'),
  'pDC-LILRA4' = c('SIGLEC10','BST2','LILRA4'),
  'cDC' = c('XCR1','CLEC9A')
)
gene_select %>% unlist()
### dotplot ====
# tmp <- FetchData(scdata_use,vars = c(
#   gene_select %>% unlist() %>% unique(),
#   'anno_res_lv2'
# ),slot = 'data')
# cell_select <- colnames(scdata_use)[is.na(scdata_use$anno_res_lv2)]
# scdata_tmp <- subset(scdata_use,cells = cell_select)
# scdata_use <- subset(scdata_use,!is.na())
DotPlot(object = scdata_use,
        features = gene_select %>% unlist() %>% unique(),
        group.by = 'anno_res_lv2',assay = 'RNA',) +
  scale_color_gradient(
    name = 'Average Expression',
    low = 'grey90',
    high = '#3131F2',
    breaks = c(-1,0,1),
    labels = c(-1,0,1)) +
  scale_size_continuous(name = 'Percent Expressed') +
  guides(
    color = guide_colourbar(nrow = 1,title.position = 'top'),
    size = guide_legend(ncol = 1,title.position = 'top')
  ) +
  labs(x = '', y= '') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0),
    axis.text = element_text(size = 18),
    legend.direction = 'vertical',
    axis.title.x = element_blank()
  )
ggsave('./result_figs_0814/allcell_dotplot_anno_lv2.png',width = 40,height = 10,bg = 'white')
ggsave('./result_figs_0814/allcell_dotplot_anno_lv2.pdf',width = 40,height = 10,bg = 'white')

# endline -----------------------------------------------------------------


