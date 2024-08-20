library(tibble)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
`%noin%` <- Negate(`%in%`)
data_meta <- read.csv('E:/jiaoxi_数据/新辅助治疗队列/vcf/data_meta.csv')
scdata <- readRDS('E:/jiaoxi_数据/新辅助治疗队列/seurat_merge.rds')
Idents(scdata) <- scdata$sample_group
scdata$sample_group %>% table()
DimPlot(scdata,group.by = c('total_medium_v2','sample_group'),label = T)

# 筛选样本 --------------------------------------------------------------------

meta_data <- scdata@meta.data %>% 
  filter(
    sample_type == 'Pre_treatment' , 
    sample_group == 'Tumor'
  ) %>% 
  rownames_to_column('cell_id') %>% 
  rename(Tumor_Sample_Barcode = 'orig.ident') %>% 
  left_join(data_meta,by = 'Tumor_Sample_Barcode') %>% 
  filter(!is.na(group_maf)) %>% 
  column_to_rownames('cell_id') %>% 
  mutate(
    group_maf = ifelse(
      test = group_maf == 'mutation',
      yes = 'complement-MUT',
      no = 'complement-WT'
    )
  ) %>% 
  arrange(desc(group_maf)) %>% 
  mutate(
    anno_use = case_when(
      total_v2 %in% c(
        "pB-JCHAIN","Bmem-CD24","pB-SSR4","Bn-FCER2",
        "pB-STMN1"
      ) ~ 'B cells',
      total_v2 %in% c(
        # "CD4+Treg-FOXP3",
        "CD4+Tem-FOS",
        "CD4+Tm-IL7R", "CD4+Th17-IL17A","CD4+Tem-GZMK",
        "CD4+Th1-IFNG", "CD4+Temra-GNLY","CD4+Tfh-CXCL13",
        "CD4+ISG-IFIT3","CD4+HS-HSPA1A"
      ) ~ 'CD4+ T cells',
      total_v2 %in% c(
        "cDC-CD1A","cDC-CLEC9A", "cDC-LAMP3"
      ) ~ 'cDC',
      total_v2 %in% c("NK-XCL1","NK-FGFBP2") ~ 'NK',
      total_v2 %in% c("Neu-CCL4","Neu-ISG15","Neu-CXCR2",
                      "Neu-S100A8") ~ 'Neu',
      total_v2 %in% c(
        "CD8+ISG-IFIT3","CD8+Tex-IFNG", "CD8+Tem-RPL13", 
        "CD8+Tex-GNLY", "CD8+Tex-CXCL13","CD8+Tem-FOS", 
        "CD8+HS-HSPA1A","CD8+Tm-IL7R"
      ) ~ 'CD8+ T',
      total_v2 %in% c(
        "Macro-SPP1","Macro-MMP9"
        # "Macro-IL1B","Macro-APOE",
        # "Macro-EGR1", "Macro-IFIT2"
      ) ~ 'Macro_M2',
      total_v2 %in% c(
        "Macro-EGR1", "Macro-IFIT2"
      ) ~ 'Macro_M1',
      total_v2 %in% c("Mono-CD16","Mono-CD14") ~ 'Mono',
      TRUE ~ total_v2
    )
  )
scdata_use <- subset(scdata,cells = rownames(meta_data))
meta_data_use <- scdata_use@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(
    meta_data %>% 
      select("group_maf",'anno_use')%>% 
      rownames_to_column('cell_id'),
    by = 'cell_id'
  ) %>% 
  column_to_rownames('cell_id')
scdata_use@meta.data <- meta_data_use
DimPlot(scdata_use,group.by = c('total_medium_v2','sample_group'),label = T)
scdata_use$total_medium_v2 %>% table()

# batch integration ---------------------------------------------------------------

scdata_use@meta.data %>% colnames()
Idents(scdata_use) <- scdata_use$orig.ident

scdata_use$total_v2 %>% table()

all.genes <- rownames(scdata_use)
scdata_use <- ScaleData(
  object = scdata_use,features = all.genes
)
scdata_use <-
  FindVariableFeatures(
    object = scdata_use,
    nfeatures = 1500
  )
scdata_use <- RunPCA(
  scdata_use,
  features = VariableFeatures(object = scdata_use)
)
scdata_use <-
  RunHarmony(
    scdata_use,
    group.by.vars = "orig.ident",
    # assay.use = "SCT",
    max_iter = 3#max_iter
  )
ElbowPlot(scdata_use, ndims = 50)
ElbowPlot(scdata_use,reduction = "harmony", ndims = 50)
pc.num = 1:37
scdata_use <- RunUMAP(
  scdata_use,
  reduction = "harmony", 
  #n.neighbors = 10,
  dims = pc.num)
scdata_use <- FindNeighbors(
  scdata_use, 
  reduction = "harmony", 
  dims = pc.num)
scdata_use <- FindClusters(scdata_use,resolution = 0.5)

DimPlot(
  object = scdata_use,
  pt.size = 1,
  label.size = 9,
  group.by = 'RNA_snn_res.0.5',#'anno_use',#c('total_v2'),
  label = TRUE
) + NoLegend()
scdata_use$RNA_snn_res.0.5 %>% table()
FeaturePlot(
  object = scdata_use,
  # group.by = 'RNA_snn_res.0.5',
  # pt.size = 0,
  features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt')
)
DimPlot(
  object = scdata_use,
  pt.size = 2,
  label.size = 7,
  group.by = 'total_medium_v2',#'total_medium_v2',#'anno_use',#c('total_v2'),
  label = TRUE
) + NoLegend()



# cell annotation lv1 ------------------------------------------------------------------

meta_data_use <- scdata_use@meta.data %>% 
  mutate(
    anno_res_lv1 = RNA_snn_res.0.5 %>% as.character()
  )
meta_data_use$RNA_snn_res.0.5 %>% unique()
meta_data_use$RNA_snn_res.0.5 %>% table()
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(0,2,5,15,16),'anno_res_lv1'] = 'NK/T cell'
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(1,4,11,17),'anno_res_lv1'] = 'Epithelial/Malignant cell'
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(3),'anno_res_lv1'] = 'Neutrophil cell'
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(8),'anno_res_lv1'] = 'Endothelial cell'
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(9),'anno_res_lv1'] = 'Fibroblast'
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(12),'anno_res_lv1'] = 'Vascular Smooth Muscle cell'
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(14),'anno_res_lv1'] = 'Mast cell'
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(10),'anno_res_lv1'] = 'B cell'
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(7,18),'anno_res_lv1'] = 'Plasma cell'
meta_data_use[meta_data_use$RNA_snn_res.0.5 %in% c(6,13),'anno_res_lv1'] = 'Macrophage/Monocyte'
meta_data_use$anno_res_lv1 <- factor(
  meta_data_use$anno_res_lv1,levels = c(
    'Epithelial/Malignant cell','Fibroblast','Vascular Smooth Muscle cell','Endothelial cell',
    'NK/T cell','B cell','Plasma cell','Macrophage/Monocyte',"Neutrophil cell",'Mast cell'
  )
)
scdata_use@meta.data <- meta_data_use
scdata_use$anno_res_lv1 %>% table()
scdata_use$anno_res_lv1 %>% unique()

### gene list ---------------------------------------------------------------


gene_select <- list(
  'Epithelial/Malignant cell' = c('FXYD3','S100A14','PERP','SFN'),
  'Fibroblast' = c('MMP2','RARRES2','LUM','DCN'),
  'Vascular Smooth Muscle cell' = c('CALD1','MYL9','IGFBP7','RGS5'),
  'Endothelial cell' = c('PECAM1','VWF'),#"SPARCL1","A2M",
  'NK/T cell' = c('CD3D','CD3E','CXCR4','CD52'),
  'B cell' = c('CD37','MS4A1','BANK1','CD79A'),
  'Plasma cell' = c('MZB1','JCHAIN','DERL3'),
  'Macrophage/Monocyte' = c('LYZ','CD68','CD14','AIF1'),
  'Neutrophil cell' = c('S100A8','S100A9','CXCR2','CSF3R'),
  'Mast cell' = c("TPSB2",'TPSAB1')
)
gene_select %>% unlist()
### dotplot ====
scdata_use <- scdata
DotPlot(object = scdata_use,
        features = gene_select %>% unlist() %>% unique(),
        group.by = 'anno_res_lv1') +
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
DimPlot(
  object = scdata_use,
  pt.size = 1,
  label.size = 8,
  group.by = c('anno_res_lv1'),
  label = TRUE
) + NoLegend() +
  labs(title = 'Cell Annotation of all cells')
## boxplot ----------------------------------------------------------------
library(ggpubr)
meta_data <- scdata_use@meta.data %>% 
  filter(anno_res_lv1 %in% c(
    'Macrophage/Monocyte','NK/T cell','B cell','Plasma cell','Neutrophil cell'
  ))
colnames(meta_data)
data_plot <- meta_data %>%
  rename(Tumor_Sample_Barcode = orig.ident) %>% 
  dplyr::select(c("Tumor_Sample_Barcode","group_maf","anno_res_lv1")) %>%
  group_by(group_maf,Tumor_Sample_Barcode,anno_res_lv1) %>%
  summarise(
    cell_num = n(),.groups = 'drop'
  ) %>%
  group_by(group_maf,Tumor_Sample_Barcode) %>%
  mutate(
    percent = (cell_num/sum(cell_num)) * 100
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
  palette = colour,add = "jitter",
  #palette =c("#E7B800", "#00AFBB"),#分组着色
  xlab = F,  bxp.errorbar=T,
  bxp.errorbar.width=0.5, #误差条大小
  size=1, outlier.shape=NA,legend = "right") +
  # ggrepel::geom_text_repel(aes(label = Tumor_Sample_Barcode)) +
  guides(fill = guide_legend(title = 'Group of maf'))+
  # scale_y_continuous(limits = c(-0.5,0.7)) +
  stat_compare_means(
    aes(group = group_maf),
    label = "p.format",
    comparisons = list(c(
      "complement-MUT","complement-WT"
    )),
    method = "wilcox.test",
  ) +
  facet_wrap(~anno_res_lv1,nrow = 1,scales = 'free') +
  theme(
    # axis.text.x = element_text(
    #   size = 20,angle = 90,
    #   hjust = 1,vjust = 0.5
    # ),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )

# datasave --------------------------------------------------------------------
write.csv(scdata_use@meta.data,'E:/jiaoxi_数据/新辅助治疗队列/result_figs/meta_data_all_anno_lv1.csv')
saveRDS(scdata_use,'E:/jiaoxi_数据/新辅助治疗队列/result_figs/scdata_all_anno_lv1.rds')

# NT&T cellanno -----------------------------------------------------------

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
`%noin%` <- Negate(`%in%`)
scdata <- readRDS('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/scdata_all_anno_lv1.rds')

scdata$anno_res_lv1 %>% unique()
scdata_nkt <- subset(
  scdata,subset = anno_res_lv1 %in% c(
    "NK/T cell"
  )
)

scdata_nkt$total_v2 %>% table()
all.genes <- rownames(scdata_nkt)
scdata_nkt <- ScaleData(
  object = scdata_nkt,features = all.genes
)
scdata_nkt <-
  FindVariableFeatures(
    object = scdata_nkt,
    nfeatures = 1500
  )
scdata_nkt <- RunPCA(
  scdata_nkt,
  features = VariableFeatures(object = scdata_nkt)
)
scdata_nkt <-
  RunHarmony(
    scdata_nkt,
    group.by.vars = "orig.ident",
    # assay.use = "SCT",
    max_iter = 2#max_iter
  )

ElbowPlot(scdata_nkt, ndims = 50)
ElbowPlot(scdata_nkt,reduction = "harmony", ndims = 50)
pc.num = 1:32
scdata_nkt <- RunUMAP(
  scdata_nkt,
  reduction = "harmony", 
  #n.neighbors = 10,
  dims = pc.num)
scdata_nkt <- FindNeighbors(
  scdata_nkt, 
  reduction = "harmony", 
  dims = pc.num)
scdata_nkt <- FindClusters(scdata_nkt,resolution = 0.5)
scdata_nkt <- FindClusters(scdata_nkt,resolution = 1)
scdata_nkt <- FindClusters(scdata_nkt,resolution = 1.5)
scdata_nkt <- FindClusters(scdata_nkt,resolution = 2)
scdata_nkt$RNA_snn_res.1 %>% table()
# scdata_nkt <- FindClusters(scdata_nkt,resolution = 3)
# scdata_nkt <- RunUMAP(scdata_nkt, dims = 1:20, label = T)
DimPlot(
  object = scdata_nkt,
  pt.size = 1,
  label.size = 6,
  group.by = c('RNA_snn_res.2'),#'anno_use',#c('total_v2'),
  label = TRUE
) + NoLegend()

DimPlot(
  object = scdata_nkt,
  pt.size = 1,
  label.size = 5,
  group.by = 'anno_use',#'anno_use',#c('total_v2'),
  label = TRUE
) + NoLegend()


meta_data_nkt <- scdata_nkt@meta.data %>% 
  mutate(
    anno_res = NA
  )
meta_data_nkt$RNA_snn_res.2 %>% unique()
meta_data_nkt$RNA_snn_res.2 %>% table()
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(5),'anno_res'] = 'NK_FCGR3A'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(6),'anno_res'] = 'NK_NCAM1'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(0,3,4,11,23),'anno_res'] = 'Cd4_FOXP3'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(2,15),'anno_res'] = 'Cd4_Tmem_IL7R'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(12),'anno_res'] = 'Cd4_Tmem_FOS'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(14),'anno_res'] = 'Cd4_KLRB1'
# meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(15),'anno_res'] = 'Cd4_暂定'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(16,18),'anno_res'] = 'Cd4_T_CXCL13'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(1,17),'anno_res'] = 'Cd8_T_GZMK'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(7),'anno_res'] = 'Cd8_T_FOS'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(8,10),'anno_res'] = 'Cd8_T_HOPX'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(9),'anno_res'] = 'Cd8_T_CXCL13'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(13),'anno_res'] = 'Cd8_T_METRNL'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(20,21),'anno_res'] = 'Cd8_T_Proliferative'
meta_data_nkt[meta_data_nkt$RNA_snn_res.2 %in% c(22),'anno_res'] = 'Cd4_T_Proliferative'
meta_data_nkt$anno_res %>% table()
meta_data_nkt$anno_res %>% unique()
meta_data_nkt$anno_res <- factor(
  meta_data_nkt$anno_res,levels = c(
    "Cd4_FOXP3","Cd4_T_CXCL13","Cd4_Tmem_IL7R","Cd4_KLRB1","Cd4_Tmem_FOS",
    "Cd4_T_Proliferative",'Cd8_T_Proliferative',"Cd8_T_FOS","Cd8_T_GZMK","Cd8_T_CXCL13","Cd8_T_METRNL","Cd8_T_HOPX",
    "NK_NCAM1","NK_FCGR3A"
  )
)
scdata_nkt@meta.data <- meta_data_nkt
meta_data_nkt$total_v2 %>% table()
scdata_nkt$anno_res %>% table()
DimPlot(
  scdata_nkt,
  reduction = "umap",
  group.by = c('anno_res'),
  label = TRUE,
  pt.size = 1,label.size = 5
) +
  scale_color_manual(values = my36colors) +
  NoLegend() +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )
## anno_res ------------------------------------------------------------


my36colors <-c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) %>% sample(36,replace = FALSE)

meta_data_sub <- scdata_nkt@meta.data
meta_data_sub[,c('umap_1','umap_2')] <- 
  scdata_nkt@reductions$umap@cell.embeddings
cell_leves <- c(
  "Cd4_FOXP3","Cd4_T_CXCL13","Cd4_Tmem_IL7R","Cd4_KLRB1","Cd4_Tmem_FOS",
  "Cd4_T_Proliferative",'Cd8_T_Proliferative',"Cd8_T_FOS","Cd8_T_GZMK","Cd8_T_CXCL13","Cd8_T_METRNL","Cd8_T_HOPX",
  "NK_NCAM1","NK_FCGR3A"
)

label_legend <- paste(
  1:length(cell_leves),
  cell_leves,sep = ': '
)
p <- DimPlot(
  object = scdata_nkt,
  pt.size = 1,
  label.size = 6,
  group.by = c('anno_res'),
  label = FALSE
) 
data_plot <- p$data %>% 
  mutate(anno_res = factor(anno_res,levels = cell_leves)) %>% 
  mutate(label_inner = as.numeric(anno_res) %>% as.factor()) %>% 
  mutate(label_legend = paste(label_inner,anno_res,sep = ': '))

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
                   fill = NA,max.overlaps = 10000,
                   label.size = NA,size = 8,
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
ggsave('./result_figs_0814/nkt_umap.png',width = 9,height = 6)
ggsave('./result_figs_0814/nkt_umap.pdf',width = 9,height = 6)


# Transfer the annotations to the whole dataset ----------------------------------------------------------------
meta_data_nkt <- scdata_nkt@meta.data
meta_data <- scdata@meta.data
colnames(meta_data)
meta_data$anno_res_lv2 %>% unique()
meta_data$anno_res_lv2 %>% table()

meta_data_all <- meta_data %>% 
  rownames_to_column('cellid') %>% 
  mutate(anno_res_lv2 = anno_res_lv2 %>% as.character()) %>% 
  left_join(
    meta_data_nkt %>% rownames_to_column('cellid') %>% 
      dplyr::select(c('cellid','anno_res')) %>% 
      mutate(anno_res_lv2 = anno_res %>% as.character()), 
    by = "cellid", 
    suffix = c("", "_data2")) %>%
  mutate(anno_res_lv2 = ifelse(!is.na(anno_res_lv2_data2), anno_res_lv2_data2, anno_res_lv2)) %>%
  dplyr::select(-anno_res_lv2_data2) %>% 
  column_to_rownames('cellid') %>% 
  mutate(
    anno_res_lv2 = ifelse(
      test = anno_use %in% c('cDC','pDC-LILRA4'),
      yes = anno_use,
      no = anno_res_lv2
    )
  ) %>% 
  mutate(
    anno_res_lv2 = ifelse(
      test = anno_res_lv2 %in% c('CD8+ T','CD4+ T','NK/T cell','Macrophage/Monocyte'),
      yes = NA,
      no = anno_res_lv2
    )
  )
meta_data_all$anno_res_lv2 %>% table()

scdata@meta.data <- meta_data_all
## calculate Ro/e ------------------------------------------------------------------
meta_data <- scdata@meta.data %>% 
  filter(anno_res_lv2 %noin% c(
    "Epithelial/Malignant cell", "Fibroblast",
    "Endothelial cell","Vascular Smooth Muscle cell",
    'CD8+ T','CD4+ T','Mast cell'
  ))

meta_data$anno_res_lv2 %>% unique()
chisq <- chisq.test(
  table(
    meta_data$anno_res_lv2 %>% as.character(),
    meta_data$group_maf
  ))
data_plot <- chisq$observed/chisq$expected %>% 
  as.matrix() %>% as.data.frame()
pheatmap::pheatmap(
  mat = data_plot,
  scale = 'none',
  fontsize = 20,
  display_numbers = TRUE
)
## Box plot, statistical testing.----------------------------------------------------------------
library(ggpubr)
meta_data <- scdata@meta.data %>% 
  filter(
    anno_res_lv2 %noin% c(
      "Epithelial/Malignant cell", "Fibroblast",
      "Endothelial cell","Vascular Smooth Muscle cell",
      'CD8+ T','CD4+ T','Mast cell'
    )
  )
colnames(meta_data)
meta_data$anno_res_lv2 %>% unique()

data_plot <- meta_data %>%
  rename(Tumor_Sample_Barcode = orig.ident) %>% 
  mutate(anno_res = anno_res_lv2) %>% 
  dplyr::select(c("Tumor_Sample_Barcode","group_maf","anno_res")) %>%
  group_by(group_maf,Tumor_Sample_Barcode,anno_res) %>%
  summarise(
    cell_num = n(),.groups = 'drop'
  ) %>%
  group_by(group_maf,Tumor_Sample_Barcode) %>%
  mutate(
    percent = (cell_num/sum(cell_num)) * 100
  ) %>% 
  filter(
    anno_res %noin% c(
      'Plasma cell','B cell','Macrophage/Monocyte',NA
    )
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
  palette = colour,add = "jitter",
  #palette =c("#E7B800", "#00AFBB"),#分组着色
  xlab = F,  bxp.errorbar=T,
  bxp.errorbar.width=0.5, #误差条大小
  size=1, outlier.shape=NA,legend = "right") +
  # ggrepel::geom_nktext_repel(aes(label = Tumor_Sample_Barcode)) +
  guides(fill = guide_legend(title = 'Group of maf'))+
  # scale_y_continuous(limits = c(-0.5,0.7)) +
  stat_compare_means(
    aes(group = group_maf),
    label = "p.format",
    comparisons = list(c(
      "complement-MUT","complement-WT"
    )),
    method = "wilcox.test",
  ) +
  facet_wrap(~anno_res,nrow = 2,scales = 'free') +
  theme(
    # axis.text.x = element_nktext(
    #   size = 20,angle = 90,
    #   hjust = 1,vjust = 0.5
    # ),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )
# save data --------------------------------------------------------------------
# setwd('/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/')
saveRDS(scdata_nkt,'./scdata_nkt.rds')
# scdata_nkt <- readRDS('./scdata_nkt.rds')
write.csv(scdata@meta.data,'/mnt/data3/jiaoxi/新辅助治疗队列/result_figs/meta_data_macmono-nkt.csv')
# endline -----------------------------------------------------------------
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
`%noin%` <- Negate(`%in%`)
scdata <- readRDS('E:/jiaoxi_数据/新辅助治疗队列/result_figs/scdata_all_anno_lv1.rds')
scdata_mac <- subset(
  scdata,subset = anno_use %in% c(
    "Macro-APOE","Mono","Macro_M1",
    "Macro_M2","Macro-IL1B"
  )
)

scdata_mac$total_v2 %>% table()
all.genes <- rownames(scdata_mac)
scdata_mac <- ScaleData(
  object = scdata_mac,features = all.genes
)
scdata_mac <-
  FindVariableFeatures(
    object = scdata_mac,
    nfeatures = 1500
  )
scdata_mac <- RunPCA(
  scdata_mac,
  features = VariableFeatures(object = scdata_mac)
)
scdata_mac <-
  RunHarmony(
    scdata_mac,
    group.by.vars = "orig.ident",
    # assay.use = "SCT",
    max_iter = 3#max_iter
  )

ElbowPlot(scdata_mac, ndims = 50)
ElbowPlot(scdata_mac,reduction = "harmony", ndims = 50)
pc.num = 1:30
scdata_mac <- RunUMAP(
  scdata_mac,
  reduction = "harmony", 
  #n.neighbors = 10,
  dims = pc.num)
scdata_mac <- FindNeighbors(
  scdata_mac, 
  reduction = "harmony", 
  dims = pc.num)
scdata_mac <- FindClusters(scdata_mac,resolution = 0.5)
scdata_mac <- FindClusters(scdata_mac,resolution = 1)
scdata_mac$RNA_snn_res.1 %>% table()
DimPlot(
  object = scdata_mac,
  pt.size = 1,
  label.size = 9,
  group.by = 'RNA_snn_res.1',#'anno_use',#c('total_v2'),
  label = TRUE
) + NoLegend()

DimPlot(
  object = scdata_mac,
  pt.size = 2,
  label.size = 9,
  group.by = 'anno_use',#'anno_use',#c('total_v2'),
  label = TRUE
) + NoLegend()
## calculate different expression gene ----------------------------------------------------------------
scdata_mac$RNA_snn_res.1 %>% unique()
scdata_mac$RNA_snn_res.1 %>% table()
Idents(scdata_mac) <- scdata_mac$RNA_snn_res.1
diff_gene_mac <- FindAllMarkers(
  object = scdata_mac,
  min.pct = 0.1,
  only.pos = FALSE
)
colnames(diff_gene_mac)
marker_list_mac <- diff_gene_mac %>% 
  filter(p_val<0.05 & avg_log2FC>=1) %>% 
  group_by(cluster) %>% 
  arrange(desc(avg_log2FC)) %>% 
  dplyr::slice(1:20) %>% 
  dplyr::select(cluster,gene)
marker_mac <- lapply(
  marker_list_mac$cluster %>% unique(), 
  function(cluster){
    marker_list_mac$gene[marker_list_mac$cluster == cluster]
  })
names(marker_mac) <- marker_list_mac$cluster %>% unique()
marker_mac[[3]]

## cell anno ----------------------------------------------------------------
meta_data_mac <- scdata_mac@meta.data %>% 
  mutate(
    anno_res = anno_use %>% as.character()
  )
meta_data_mac$RNA_snn_res.1 %>% unique()
meta_data_mac$RNA_snn_res.1 %>% table()
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(0),'anno_res'] = 'Mac_MSR1'
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(1,10),'anno_res'] = 'Mac_RGS1'
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(2),'anno_res'] = 'Mac_RPL/RPS'
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(3,7),'anno_res'] = 'Mono_FCN1'
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(4,8),'anno_res'] = 'Mac_IL7R'
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(5,6),'anno_res'] = 'Mac_FLT1'
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(9),'anno_res'] = 'Macro_prolif'
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(11),'anno_res'] = 'CD4+ T'
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(12),'anno_res'] = 'Mono_FCGR3B'
meta_data_mac[meta_data_mac$RNA_snn_res.1 %in% c(13),'anno_res'] = 'CD8+ T'
meta_data_mac$anno_res <- factor(
  meta_data_mac$anno_res,levels = c(
    'Mac_MSR1','Mac_RGS1','Mac_RPL/RPS','Mac_IL7R','Mac_FLT1','Macro_prolif',
    'Mono_FCN1','Mono_FCGR3B','CD4+ T','CD8+ T'
  )
)
scdata_mac@meta.data <- meta_data_mac
scdata_mac$anno_res %>% table()
scdata_mac$anno_res %>% unique()
scdata_mac <- subset(scdata_mac,subset = anno_res %noin% c(
  'CD4+ T','CD8+ T'
))


my36colors <-c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) %>% sample(36,replace = FALSE)
scdata_mac <- subset(scdata_mac,subset = anno_res %noin% c(
  'CD4+ T','CD8+ T'
))
meta_data_sub <- scdata_mac@meta.data
meta_data_sub[,c('umap_1','umap_2')] <- 
  scdata_mac@reductions$umap@cell.embeddings
cell_leves <- c(
  'Mac_MSR1','Mac_RGS1','Mac_RPL/RPS','Mac_IL7R','Mac_FLT1','Macro_prolif',
  'Mono_FCN1','Mono_FCGR3B'
)

label_legend <- paste(
  1:length(cell_leves),
  cell_leves,sep = ': '
)
p <- DimPlot(
  object = scdata_mac,
  pt.size = 1,
  label.size = 6,
  group.by = c('anno_res'),
  label = FALSE
) 
data_plot <- p$data %>% 
  mutate(anno_res = factor(anno_res,levels = cell_leves)) %>% 
  mutate(label_inner = as.numeric(anno_res) %>% as.factor()) %>% 
  mutate(label_legend = paste(label_inner,anno_res,sep = ': '))

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
                   fill = NA,max.overlaps = 10000,
                   label.size = NA,size = 8,
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
ggsave('./result_figs_0814/macro_umap.png',width = 9,height = 6)
ggsave('./result_figs_0814/macro_umap.pdf',width = 9,height = 6)


# dotplot -----------------------------------------------------------------




gene_select <- list(
  'Mac_MSR1' = c('MSR1','C2','STAB1','GPR34','CD84'),
  'Mac_RGS1' = c('PTGER4','GPR183','PRDM1','RGS1'),
  'Mac_RPL/RPS' = c("RPS18","RPS11","RPL30","RPS23"),
  'Mac_IL7R' = c('ANPEP','MMP9','MMP12','IL7R'),
  'Mac_FLT1' = c('TNIP3','INHBA','FLT1','CCL20'),
  'Macro_prolif' = c('PCLAF','TYMS','PCNA','MKI67'),
  'Mono_FCN1' = c('FCN1','VCAN','S100A12','LYZ'),
  'Mono_FCGR3B' = c('FCGR3B','CXCR2','CXCR1','ADGRG3')
  # 'CD4+ T' = c('CD3D','CD3E','CD4'),
  # 'CD8+ T' = c('CD8A')
)
gene_select %>% unlist()
DotPlot(object = scdata_mac,
        features = gene_select %>% unlist() %>% unique(),
        group.by = 'anno_res') +
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
DimPlot(
  object = scdata_mac,
  pt.size = 1,
  label.size = 8,
  group.by = c('anno_res'),
  label = TRUE
) + NoLegend() +
  labs(title = 'Cell Annotation of Macro/Mono')
## transfer the annotations to the whole data ----------------------------------------------------------------
meta_data <- scdata@meta.data
colnames(meta_data)
meta_data$anno_res_lv1 %>% unique()

meta_data_all <- meta_data %>% 
  rownames_to_column('cellid') %>% 
  mutate(anno_res_lv2 = anno_res_lv1 %>% as.character()) %>% 
  left_join(
    meta_data_mac %>% rownames_to_column('cellid') %>% 
      dplyr::select(c('cellid','anno_res')) %>% 
      mutate(anno_res_lv2 = anno_res %>% as.character()), 
    by = "cellid", 
    suffix = c("", "_data2")) %>%
  mutate(anno_res_lv2 = ifelse(!is.na(anno_res_lv2_data2), anno_res_lv2_data2, anno_res_lv2)) %>%
  dplyr::select(-anno_res_lv2_data2) %>% 
  column_to_rownames('cellid')
meta_data_all$anno_res_lv2 %>% table()
scdata@meta.data <- meta_data_all

## calculate Ro/e ------------------------------------------------------------------
meta_data <- scdata@meta.data %>% 
  filter(anno_res_lv2 %noin% c(
    "Endothelial/Malignant cell", "Fibroblast",
    "Endothelial cell",
    'CD8+ T','CD4+ T','Mast cell'
  ))

meta_data$anno_res_lv2 %>% unique()
chisq <- chisq.test(
  table(
    meta_data$anno_res_lv2 %>% as.character(),
    meta_data$group_maf
  ))
data_plot <- chisq$observed/chisq$expected %>% 
  as.matrix() %>% as.data.frame()
pheatmap::pheatmap(
  mat = data_plot,
  scale = 'none',
  fontsize = 20,
  display_numbers = TRUE
)
## boxplot ----------------------------------------------------------------
library(ggpubr)
meta_data <- scdata@meta.data %>% 
  filter(
    anno_res_lv2 %noin% c(
      "Epithelial/Malignant cell", "Fibroblast",
      "Endothelial cell",'Vascular Smooth Muscle cell',
      'CD8+ T','CD4+ T'
    )
  )
colnames(meta_data)
meta_data$anno_res_lv2 %>% unique()

data_plot <- meta_data %>%
  rename(Tumor_Sample_Barcode = orig.ident) %>% 
  mutate(anno_res = anno_res_lv2) %>% 
  dplyr::select(c("Tumor_Sample_Barcode","group_maf","anno_res")) %>%
  group_by(group_maf,Tumor_Sample_Barcode,anno_res) %>%
  summarise(
    cell_num = n(),.groups = 'drop'
  ) %>%
  group_by(group_maf,Tumor_Sample_Barcode) %>%
  mutate(
    percent = (cell_num/sum(cell_num)) * 100
  ) %>% 
  filter(
    anno_res %noin% c(
      'Plasma cell','B cell','Macrophage/Monocyte'
    )
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
  palette = colour,add = "jitter",
  #palette =c("#E7B800", "#00AFBB"),#分组着色
  xlab = F,  bxp.errorbar=T,
  bxp.errorbar.width=0.5, #误差条大小
  size=1, outlier.shape=NA,legend = "right") +
  # ggrepel::geom_text_repel(aes(label = Tumor_Sample_Barcode)) +
  guides(fill = guide_legend(title = 'Group of maf'))+
  # scale_y_continuous(limits = c(-0.5,0.7)) +
  stat_compare_means(
    aes(group = group_maf),
    label = "p.format",
    comparisons = list(c(
      "complement-MUT","complement-WT"
    )),
    method = "wilcox.test",
  ) +
  facet_wrap(~anno_res,nrow = 2,scales = 'free') +
  theme(
    # axis.text.x = element_text(
    #   size = 20,angle = 90,
    #   hjust = 1,vjust = 0.5
    # ),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )
# save data --------------------------------------------------------------------
saveRDS(scdata,'../scdata_anno_myeloid.rds')
# scdata_mac <- readRDS('./scdata_myeloid.rds')
# endline -----------------------------------------------------------------






# endline -----------------------------------------------------------------


