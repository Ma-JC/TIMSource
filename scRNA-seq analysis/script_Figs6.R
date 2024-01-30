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
# Figs6 A  umap of CT26 ----------------------------------------------------------------
scRNA_11 <- readRDS('D:/jiaoxi/scRNA_11_res.rds')
cell_leves <- c(
  "Cancer cell","Fibroblast",
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory","NK",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "Mast_Mcpt2","B cell","cDC1_Xcr1","tDC"
)
scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% cell_leves)
scRNA_11_sub$anno_res %>% unique()
label_legend <- paste(
  1:length(cell_leves),
  cell_leves,sep = ': '
)
p <- DimPlot(
  object = scRNA_11_sub,
  group.by = c('anno_res'),
  label = TRUE,
  pt.size = 0.8,
  label.size = 6,
  label.box = TRUE,
  repel = 30,
  alpha = 0.8
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
my36colors <-c(
  "#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) 
ggplot() +
  geom_point(
    data = data_plot,
    aes(x = umap_1,y = umap_2,color = label_inner),
    size = 0.5) +
  geom_label_repel(data = cell_label_loc,
                   aes(
                     x = umap_1,
                     y = umap_2,
                     label = label_inner),
                   size = 5,
                   alpha = 0.8,
                   fill = NA,
                   label.padding = unit(0, "lines"),
                   label.r = unit(0,'lines'),
                   label.size = NA,
                   segment.color = NA,
                   show.legend = FALSE)+
  geom_segment(aes(
    x = data_plot$umap_1 %>% min() * 1.4, 
    y = data_plot$umap_2 %>% min() * 1.1 , 
    xend = (data_plot$umap_1 %>% min())* 1.4+
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2, 
    yend = (data_plot$umap_2 %>% min())* 1.1
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.5) +
  geom_segment(aes(
    x = data_plot$umap_1 %>% min()* 1.4, 
    y = data_plot$umap_2 %>% min() * 1.1, 
    xend = (data_plot$umap_1 %>% min())* 1.4, 
    yend = (data_plot$umap_2 %>% min()) + 
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2 
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.5) +
  scale_color_manual(
    values = my36colors,
    label = label_legend
  ) +
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 2)))+
  labs(x = 'Umap 1',y= 'Umap 2')+
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin =  margin(1, 1, 1,1, "lines"),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_text(size = 8,hjust = 0.1,vjust = 7),
    axis.title.y = element_text(size = 8,vjust = -6,hjust = 0.08)
  )

# Figs6 C dotplot of CT26 --------------------------------------------------

scRNA_11$anno_res %>% unique()
cell_type_select <- c(
  "Cancer cell","Fibroblast",
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory","NK",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "Mast_Mcpt2",
  "cDC1_Xcr1","tDC","B cell"
)
scRNA_11_tmp <- subset(
  scRNA_11,
  subset = anno_res %in% cell_type_select
)
scRNA_11_tmp$anno_res <- 
  scRNA_11_tmp$anno_res %>% factor(levels = cell_type_select)
scRNA_11_tmp$anno_res %>% unique()
gene_select <- list(
  'Cancer cell' = c('Esm1','Anln'),
  'Fibroblast' = c('Dcn','Aspn','Mmp3'),
  "Dysfunction_T cell" = c('Ctla2a','Ctsw','Pdcd4'),
  "Cytolytic_T cell" = c('Havcr2','Ccl4','Ccl5'),
  'Gzmb_T cell' = c('Gzmb','Gzme','Gzmk'),
  'Proliferative CD8+ T' = c('Top2a','Cdk1'),
  "CD8+ T memory" = c('Ccr7','Xcl1'),
  'Naive T' = c('Tcf7','Emb','Lef1'),
  'Treg' = c('Foxp3','Il2ra'),
  "Proliferative Treg" = c('Kif11','Aurkb'),
  "CD4+ T memory" = c('Cd4','Rora'),
  "NK" = c('Ncr1','Klrd1','Prf1'),
  "Mono_Arhgap26" = c('Arhgap26','Napsa'),
  "Mono_Ly6i" = c('Ly6i','Gpr141','Gbp5'),
  "Mono_Isg" = c('Csf1r','Isg15'),
  "Mac/Mono_Arg1" = c('Arg1','Mmp12'),
  "Mac/Mono_Cx3cr1" = c('Cx3cr1','Pmepa1'),
  "Mac/Mono_Ccl8" = c('Ccl8', 'Serpina3g'),
  "Mac/Mono_Mki67" = c('Mki67','Ccna2'),
  "NPBNs" = c('S100a8','S100a9'),
  "Exhausted TAN" = c('Ccl3','Ccl4','Cxcl2'),
  "interferon-stimulated NAN" = c('Gbp5','Ifit2'),
  'Mast_Mcpt2' = c('Mcpt2','Kit'),
  'cDC1_Xcr1' = c('Xcr1','Clec9a'),
  'tDC' = c('Flt3','Fscn1'),
  'B cell' = c('Cd79a','Cd19','Ms4a1')
)
gene_select %>% unlist()
data.usage <- DotPlot(
  object = scRNA_11_tmp,
  features = gene_select %>% unlist() %>% unique(),#gene_select %>% unique(),
  group.by = 'anno_res'
)$data
data.anno <- data.frame(
  features.plot = unique(data.usage$features.plot),
  label = lapply(gene_select %>% names(), function(x) {
    rep(x, length(gene_select[[x]]))
  }) %>% unlist()
)

df.plot <- plyr::join(data.usage,data.anno)
df.plot$label <- factor(
  df.plot$label,
  levels = cell_type_select
)
df.plot$id[1:3]
y_group <- c(
  'Cancer cell','Fibroblast',
  "CD8+ T","CD8+ T","CD8+ T",
  "CD8+ T","CD8+ T",
  "Naive T","CD4+ T","CD4+ T","CD4+ T","NK",
  "Mac","Mac","Mac","MonoMac",
  "Mono",
  'Isg15+ Neutrophil','Ccl4+ Neutrophil',
  "Mast",
  "DC","DC","B cell"
)
names(y_group) <- cell_type_select
df.plot$y_group <- y_group[df.plot$label]
p <- ggplot(df.plot,aes(x=features.plot,y =  as.numeric(id),size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("% detected", range = c(0,5)) + 
  scale_color_gradientn(colours = viridis::magma(30)[1:20],
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nexpression") +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("Markers") + theme_bw() +
  scale_y_continuous(breaks = 1:length(levels(df.plot$id)),labels = levels(df.plot$id),sec.axis = dup_axis())+ 
  facet_grid(~label, scales="free",space = "free")+theme_classic() +
  theme(
    axis.text.x = element_text(size=16, angle=-45, hjust=0, color="black",face="bold"),
    axis.text.y = element_text(size=16, color="black",face="bold"),
    axis.title.x = element_text(size=25,colour = 'black',vjust = -0.8,hjust = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x.top = element_line(linetype = 'solid'),
    axis.line = element_line(colour = 'grey30',size = 0.2),
    panel.spacing=unit(0, "mm"),
    panel.border = element_rect(fill = NA,linetype = 'dashed',linewidth = 0.4),
    panel.grid.major.y = element_line(),
    strip.text.x = element_blank(),
    strip.background = element_blank()
  )
p
# Figs6 D Ro/e of CT26 --------------------------------------------------
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
# Figs6 G ssGSEA of cancer cell -----------------------------------------------------------------
geneSet <- readRDS('D:/jiaoxi/ssgsea/gene_tran/geneSet.rds')
names(geneSet)
name_select <- c(
  "Neutrophil:ECM remodeling",
  "CancerCell:Tumor proliferation rate",
  "CancerCell:Checkpoint molecules",
  "self:Angiogenesis"
)
geneSet_sub <- geneSet[name_select]
names(geneSet_sub)
scRNA_caf <- subset(scRNA_11,subset = anno_res %in% c(
  "Cancer cell","Fibroblast"
))
cells_rankings <- AUCell_buildRankings(
  scRNA_caf@assays$RNA@data)
cells_AUC <- AUCell_calcAUC(
  geneSet_sub,
  cells_rankings, 
  aucMaxRank=nrow(cells_rankings)*0.1)
res <- data.frame(row.names = c('group','AUCS',"geneSet")) %>% t()
for(geneSet in names(geneSet_sub)){
  aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ]) %>% 
    as.data.frame() %>% 
    mutate(group = scRNA_caf$batch) %>% 
    rename_all(~c('AUCS','group')) %>% 
    mutate(geneSet = geneSet %>% str_replace_all('-',' ') %>% str_to_lower()) %>% 
    as.data.frame() %>% 
    dplyr::select('group','AUCS',"geneSet")
  res <- rbind(res,aucs)
}
colnames(res)
group_label <- c("Mock","C5aRA","PD1","C5aRA_PD1")
names(group_label) <- c("MOCK","C5aRa","PD1","C5aRA_PD1")
data_plot <- res %>% 
  filter(geneSet %in% c(
    "neutrophil:ecm remodeling","cancercell:tumor proliferation rate"
  )) %>% 
  mutate(group = group_label[group]) %>% 
  mutate(
    group = factor(group,levels = c(
      "Mock","C5aRA","PD1","C5aRA_PD1"
    ))
  ) %>% 
  mutate(
    geneSet = ifelse(test = geneSet == 'neutrophil:ecm remodeling',
                     yes = 'Ecm Remodeling',
                     no = 'Tumor Proliferation Rate')
  )
filtered_data_plot <- data_plot %>%
  group_by(geneSet,group) %>%
  mutate(Q1 = quantile(AUCS, .25),
         Q3 = quantile(AUCS, .75)) %>%
  ungroup() %>%
  mutate(IQR = Q3 - Q1) %>%
  filter(AUCS >= (Q1 - 1.5 * IQR) & AUCS <= (Q3 + 1.5 * IQR))
ggplot(filtered_data_plot, aes(x = group, y = AUCS,fill = group)) +
  geom_violin(trim=TRUE,color="white",scale = 'width') +
  geom_boxplot(width=0.2,position=position_dodge(0.9),show.legend = FALSE,size = 0.5,outlier.size = 0)+ 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x,width = 30,whitespace_only = FALSE)) +
  scale_fill_manual(values = c('lightslategrey','#0099B4', '#42B540', '#ED0000', '#00468B') %>% rev()) +
  stat_compare_means(aes(group =  as.factor(group),x = as.factor(group)),
                     label = paste0("p = ", after_stat("p.format")),
                     method = "t.test",
                     hide.ns = FALSE,
                     comparisons = 
                       list(
                         c('Mock','C5aRA'),
                         c('PD1','C5aRA_PD1')
                       ),
                     step.increase = 0
  ) +
  facet_wrap(~geneSet,scales = 'free_y') +
  labs(y = 'AUCell Score') +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16,color="black",face="bold",angle =45,hjust = 1,vjust = 1),
    axis.text.y = element_text(size=16, color="black",face="bold"),
    axis.title.y = element_text(size=16,colour = 'black',face="bold",vjust = 0,hjust = 0.5),#坐标轴标题
    axis.text.y.right = element_blank(),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 23),
    plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 24),
    plot.margin = margin(2,2,2,2),
    axis.line = element_line(colour = 'grey30',size = 0.2), 
    strip.text = element_text(size  = 18)
  )
# Figs6 B umap of AKR -----------------------------------------------------------------
scRNA_12 <- readRDS('D:/jiaoxi_12/scRNA_12_fig_use.rds')
scRNA_12_tmp <- scRNA_12
scRNA_12_tmp$anno_res %>% unique()
cell_type_select <- c(
  'Endothelial Cells',"Perithelial cell","Epithelial cell",
  'Cancer cell/Fibroblast',
  "Dysfunction_T cell","Cytolytic/Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg",
  'CD4+ T',"Proliferative CD4+ T",
  "NK",'Proliferative NK',
  "Mono_Arhgap26","Mac/Mono_Col1a1","Mac/Mono_Cx3cr1","Mac/Mono_Mki67",
  "Mono_Isg",
  "S100A8/ISG-Neu","Ccl4-Neu",
  "Mast","cDC","tDC",'pDC'
)
scRNA_12_tmp <- subset(
  scRNA_12,
  subset = anno_res %in% cell_type_select
)
scRNA_12_tmp$anno_res <- 
  scRNA_12_tmp$anno_res %>% factor(levels = cell_type_select)
scRNA_12_tmp$anno_res %>% unique()
cell_leves <- c(
  'Endothelial Cells',"Perithelial cell","Epithelial cell",
  'Cancer cell/Fibroblast',
  "Dysfunction_T cell","Cytolytic/Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg",
  'CD4+ T',"Proliferative CD4+ T",
  "NK",'Proliferative NK',
  "Mono_Arhgap26","Mac/Mono_Col1a1","Mac/Mono_Cx3cr1","Mac/Mono_Mki67",
  "Mono_Isg",
  "S100A8/ISG-Neu","Ccl4-Neu",
  "Mast","cDC","tDC",'pDC'
)
label_legend <- paste(
  1:length(cell_leves),
  cell_leves,sep = ': '
)
p <- DimPlot(
  object = scRNA_12_tmp,
  group.by = c('anno_res'),
  label = TRUE,
  pt.size = 0.8,
  label.size = 6,
  label.box = TRUE,
  repel = 30,
  alpha = 0.8
)
p
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
my36colors <-c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#ad494a","#8ca252","#e7ba52"
) %>% sample(26,replace = FALSE)
ggplot() +
  geom_point(
    data = data_plot,
    aes(x = umap_1,y = umap_2,color = label_inner),
    size = 0.5) +
  geom_label_repel(data = cell_label_loc,
                   aes(
                     x = umap_1,
                     y = umap_2,
                     label = label_inner),
                   size = 5,
                   alpha = 0.8,
                   fill = NA,
                   label.padding = unit(0, "lines"),
                   label.r = unit(0,'lines'),
                   label.size = NA,
                   segment.color = NA,
                   show.legend = FALSE)+
  geom_segment(aes(
    x = data_plot$umap_1 %>% min() * 1.4, 
    y = data_plot$umap_2 %>% min() * 1.1 , 
    xend = (data_plot$umap_1 %>% min())* 1.4+
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2, 
    yend = (data_plot$umap_2 %>% min())* 1.1
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.5) +
  geom_segment(aes(
    x = data_plot$umap_1 %>% min()* 1.4, 
    y = data_plot$umap_2 %>% min() * 1.1, 
    xend = (data_plot$umap_1 %>% min())* 1.4, 
    yend = (data_plot$umap_2 %>% min()) + 
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2 
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.5) +
  scale_color_manual(
    values = my36colors,
    label = label_legend
  ) +
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 2)))+
  labs(x = 'Umap 1',y= 'Umap 2')+
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin =  margin(1, 1, 1,1, "lines"),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_text(size = 8,hjust = 0.1,vjust = 7),
    axis.title.y = element_text(size = 8,vjust = -6,hjust = 0.08)
  )

# Figs6 E dotplot of AKR --------------------------------------------------

gene_select <- list(
  'Endothelial Cells' = c('Pecam1','Vwf'),
  "Perithelial cell" = c('Pdgfrb','Rgs5','Cspg4'),
  'Epithelial cell' = c('Epcam','Cdh1'),
  'Cancer cell/Fibroblast' = c('Dcn','Mmp3'),
  "Dysfunction_T cell" = c('Ctla2a','Ctsw','Pdcd4'),
  "Cytolytic/Gzmb_T cell" = c('Gzmk','Gzmb','Bcl2','Ccl5','Hcst',
                              'Havcr2','Tigit','Lag3','Ctla4'),
  'Proliferative CD8+ T' = c('Mcm2','Top2a','Cdk1'),
  "CD8+ T memory" = c('Ccr7','Xcl1'),
  'Naive T' = c('Tcf7','Emb','Lef1'),
  'Treg' = c('Foxp3','Il2ra'),
  "Proliferative Treg" = c('Kif11','Aurkb'),
  "Proliferative CD4+ T" = c('Cd4','Rora'),
  "NK" = c('Ncr1','Klrd1','Prf1'),
  "Mono_Arhgap26" = c('Arhgap26','Napsa'),
  "Mac/Mono_Col1a1" = c('Col1a1','Lamb1'),
  "Mono_Isg" = c('Csf1r','Isg15'),
  "Mac/Mono_Cx3cr1" = c('Cx3cr1','Pmepa1'),
  "Mac/Mono_Mki67" = c('Mki67','Ccna2'),
  "S100A8/ISG-Neu" = c('S100a8','S100a9','Gbp5','Ifit2'),
  "Ccl4-Neu" = c('Ccl3','Ccl4','Cxcl2'),
  'Mast' = c('Mcpt2','Kit'),
  'cDC' = c('Xcr1','Clec9a'),
  'tDC' = c('Flt3','Fscn1'),
  'pDC' = c('Siglech','Bst2')
)
gene_select %>% unlist()
tmp <- gene_select %>% unlist() %>% as.character()
tmp[duplicated(tmp)]
cell_leves <- c(
  'Endothelial Cells',"Perithelial cell","Epithelial cell",
  'Cancer cell/Fibroblast',
  "Dysfunction_T cell","Cytolytic/Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg",
  'CD4+ T',"Proliferative CD4+ T",
  "NK",'Proliferative NK',
  "Mono_Arhgap26","Mac/Mono_Col1a1","Mac/Mono_Cx3cr1","Mac/Mono_Mki67",
  "Mono_Isg",
  "S100A8/ISG-Neu","Ccl4-Neu",
  "Mast","cDC","tDC",'pDC'
)
scRNA_12_tmp <- subset(scRNA_12,subset = anno_res %in% cell_leves)
scRNA_12_tmp$anno_res <- factor(
  scRNA_12_tmp$anno_res,
  levels = cell_leves
)
DotPlot(object = scRNA_12_tmp,
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
    legend.direction = 'vertical',
    axis.title.x = element_blank()
  )

# Figs6 F Ro/e of AKR --------------------------------------------------
scRNA_12$anno_res %>% unique()
scRNA_12_tmp <- subset(scRNA_12,subset = anno_res %in% c(
  'Endothelial Cells',"Perithelial cell","Epithelial cell",
  'Cancer cell/Fibroblast',
  "Dysfunction_T cell","Cytolytic/Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg",
  'CD4+ T',"Proliferative CD4+ T",
  "NK",'Proliferative NK',
  "Mono_Arhgap26","Mac/Mono_Col1a1","Mac/Mono_Cx3cr1","Mac/Mono_Mki67",
  "Mono_Isg",
  "S100A8/ISG-Neu","Ccl4-Neu",
  "Mast","cDC","tDC",'pDC'
))
cell_label_correct <- c(
  'Endothelial Cells',"Perithelial cell","Epithelial cell",
  'Cancer cell/Fibroblast',
  "Dysfunction_T cell","Cytolytic/Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg",
  'CD4+ T',"Proliferative CD4+ T",
  "NK",'Proliferative NK',
  "Mono_Arhgap26","Mac/Mono_Col1a1","Mac/Mono_Cx3cr1","Mac/Mono_Mki67",
  "Mono_Isg",
  # "Neutrophil",
  "S100A8/ISG-Neu","Ccl4-Neu",
  "Mast","cDC","tDC",'pDC'
)
names(cell_label_correct) <- c(
  'Endothelial Cells',"Perithelial cell","Epithelial cell",
  'Cancer cell/Fibroblast',
  "Dysfunction_T cell","Cytolytic/Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg",
  'CD4+ T',"Proliferative CD4+ T", 
  "NK",'Proliferative NK',
  "Mono_Arhgap26","Mac/Mono_Col1a1","Mac/Mono_Cx3cr1","Mac/Mono_Mki67",
  "Mono_Isg",
  "S100A8/ISG-Neu","Ccl4-Neu",
  "Mast","cDC","tDC",'pDC'
)
scRNA_12_tmp$anno_tmp <- cell_label_correct[scRNA_12_tmp$anno_res] %>% as.character()
scRNA_12_tmp$anno_tmp %>% unique()
meta_data_sub <- scRNA_12_tmp@meta.data
meta_data_sub$anno_res  <- meta_data_sub$anno_tmp
celltype_select <- c(
  'Endothelial Cells',"Perithelial cell","Epithelial cell",
  'Cancer cell/Fibroblast',
  "Dysfunction_T cell","Cytolytic/Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg",
  'CD4+ T',"Proliferative CD4+ T", 
  "NK",'Proliferative NK',
  "Mono_Arhgap26","Mac/Mono_Col1a1","Mac/Mono_Cx3cr1","Mac/Mono_Mki67",
  "Mono_Isg",
  "S100A8/ISG-Neu","Ccl4-Neu",
  "Mast","cDC","tDC",'pDC'
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
    axis.text = element_text(size = 20,color = 'black'),
    axis.text.x = element_text(angle = 30,hjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.text = element_text(size = 16),
    legend.position = 'right',
    legend.justification = c(0,1),
    legend.title = element_text(size = 20),
    strip.text = element_blank(),
    strip.background = element_blank()
  )
# Figs6 H ssGSEA of cancer cell -----------------------------------------------------------------
geneSet <- readRDS('D:/jiaoxi/ssgsea/gene_tran/geneSet.rds')
names(geneSet)
name_select <- c(
 "Neutrophil:ECM remodeling",
  "CancerCell:Tumor proliferation rate",
  "CancerCell:Checkpoint molecules",
  "self:Angiogenesis")
geneSet_sub <- geneSet[name_select]
names(geneSet_sub)

scRNA_caf <- subset(scRNA_12,subset = anno_res_1 %in% c(
  "Cancer cell/Fibroblast"
))
scRNA_caf$anno_res <- scRNA_caf$anno_res_1
cells_rankings <- AUCell_buildRankings(
  scRNA_caf@assays$RNA@data)
cells_AUC <- AUCell_calcAUC(
  geneSet_sub,
  cells_rankings, 
  aucMaxRank=nrow(cells_rankings)*0.1)

res <- data.frame(row.names = c('group','AUCS',"geneSet")) %>% t()
for(geneSet in names(geneSet_sub)){
  aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ]) %>% 
    as.data.frame() %>% 
    mutate(group = scRNA_caf$batch) %>% 
    rename_all(~c('AUCS','group')) %>% 
    mutate(geneSet = geneSet %>% str_replace_all('-',' ') %>% str_to_lower()) %>% 
    as.data.frame() %>% 
    dplyr::select('group','AUCS',"geneSet")
  res <- rbind(res,aucs)
}
colnames(res)
data_plot <- res %>% 
  filter(geneSet %in% c(
    "neutrophil:ecm remodeling","cancercell:tumor proliferation rate"
  )) %>% 
  mutate(
    group = factor(group,levels = c(
      "Mock","C5aRA","PD1","C5aRA_PD1"
    ))
  ) %>% 
  mutate(
    geneSet = ifelse(test = geneSet == 'neutrophil:ecm remodeling',
                     yes = 'Ecm Remodeling',
                     no = 'Tumor Proliferation Rate')
  )
filtered_data_plot <- data_plot %>%
  group_by(geneSet,group) %>%
  mutate(Q1 = quantile(AUCS, .25),
         Q3 = quantile(AUCS, .75)) %>%
  ungroup() %>%
  mutate(IQR = Q3 - Q1) %>%
  filter(AUCS >= (Q1 - 1.5 * IQR) & AUCS <= (Q3 + 1.5 * IQR))
ggplot(filtered_data_plot, aes(x = group, y = AUCS,fill = group)) +
  geom_violin(trim=TRUE,color="white",scale = 'area') +
  geom_boxplot(width=0.2,position=position_dodge(0.9),show.legend = FALSE,size = 0.5,outlier.size = 0)+ 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x,width = 30,whitespace_only = FALSE)) +
  scale_fill_manual(values = c('lightslategrey','#0099B4', '#42B540', '#ED0000', '#00468B') %>% rev()) +
  stat_compare_means(aes(group =  as.factor(group),x = as.factor(group)),
                     label = paste0("p = ", after_stat("p.format")),
                     method = "t.test",
                     hide.ns = FALSE,
                     comparisons = 
                       list(
                         c('Mock','C5aRA'),
                         c('PD1','C5aRA_PD1')
                       ),
                     step.increase = 0
  ) +
  facet_wrap(~geneSet,scales = 'free_y') +
  labs(y = 'AUCell Score') +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16,color="black",face="bold",angle =45,hjust = 1,vjust = 1),
    axis.text.y = element_text(size=16, color="black",face="bold"),
    axis.title.y = element_text(size=16,colour = 'black',face="bold",vjust = 0,hjust = 0.5),
    axis.text.y.right = element_blank(),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 23),
    plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 24),
    plot.margin = margin(2,2,2,2),
    axis.line = element_line(colour = 'grey30',size = 0.2), 
    strip.text = element_text(size  = 18)
  )
# endline -----------------------------------------------------------------


