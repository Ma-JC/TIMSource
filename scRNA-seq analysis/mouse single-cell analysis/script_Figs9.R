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
# Figs9 A -----------------------------------------------------------------
scRNA_12 <- readRDS('D:/jiaoxi_12/scRNA_12_anno.rds')
scRNA_12$anno_res_1 %>% table()
geneSet <- readRDS('D:/jiaoxi_12/ssgsea/gene_tran/geneSet.rds')
scRNA_12$anno_res_1 %>% unique()
scRNA_use <- subset(
  scRNA_12,
  subset = anno_res_1 %in% c('Neutrophil') 
)
scRNA_use$anno_res_1 %>% unique()
scRNA_use$batch %>% table()
gene_exp <- scRNA_use@assays$RNA@counts
ssgsea <-
  gsva(
    gene_exp,
    geneSet,
    method = 'ssgsea',
    kcdf = 'Gaussian',
    abs.ranking = TRUE
  )
a <- ssgsea %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(group = scRNA_use$batch) %>% 
  rownames_to_column("sample")
names(a)
scRNA_use_ggsea <- gather(a,key = ssgsea, value = Expression, -c(group,sample)) %>% 
  mutate(
    group = factor(group,levels = c("Mock","C5aRA","PD1","C5aRA_PD1"))
  ) %>% 
  separate(
    ssgsea, 
    c("facet", "ssgsea_name"),
    sep = ":",
    remove = FALSE) %>% 
  filter(ssgsea %in% c(
    'self:18-gene IFN signature',
    'Neutrophil:ECM remodeling',
    "Neutrophil:Interferon signaling",
    "CancerDiscovery:Interleukin signaling pathway",
    "CancerDiscovery:Inflammation",
    "Neutrophil:Neutrophil cytotoxicity",
    "Neutrophil:Neutrophil",
    "self:Cytolytic activity",
    "Neutrophil:Angiogenesis",
    "CancerDiscovery:Phagocytosis"
  )) %>%
  mutate(ssgsea = ifelse(
    ssgsea == 'Neutrophil:Neutrophil',
    'Neutrophil:Neutrophil degranulation',
    ssgsea)
  ) %>% 
  mutate(
    ssgsea = factor(ssgsea)
  )
scRNA_use_ggsea$ssgsea %>% unique()



ssgsea_i <- 'Neutrophil:Interferon signaling'
ssgsea_i <- 'Neutrophil:Neutrophil degranulation'
ssgsea_violon_self(data = scRNA_use_ggsea,ssgsea_i = ssgsea_i)
scRNA_use_ggsea %>% colnames()
head(scRNA_use_ggsea)
data_plot <- scRNA_use_ggsea %>% 
  mutate(group =factor(group,levels = c('Mock', 'C5aRA','PD1','C5aRA_PD1'))) %>% 
  filter(group %in% c('PD1','C5aRA_PD1')) %>% 
  filter(ssgsea %in% c(
    'Neutrophil:Interferon signaling',
    'Neutrophil:Neutrophil degranulation'
  ))
ggplot(data_plot, aes(x = group, y = Expression,fill = group)) +
  geom_violin(trim=FALSE,color="white",scale = 'width') +
  geom_boxplot(width=0.2,position=position_dodge(0.9),show.legend = FALSE,size = 0.5)+ #绘制箱线图
  scale_x_discrete(labels = function(x) stringr::str_wrap(x,width = 30,whitespace_only = FALSE)) +
  scale_fill_manual(values = c('#0099B4', '#42B540') %>% rev()) +
  stat_compare_means(aes(group =  group),
                     label = "p.format",
                     method = "wilcox.test",
                     hide.ns = FALSE,
                     size = 7,
                     comparisons = list(c('PD1','C5aRA_PD1'))
  ) + 
  facet_wrap(~ssgsea_name,scale = 'free_y',ncol = 2) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NULL,colour = 'black'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=16, color="black",face="bold"),
    axis.title.y = element_text(size=22,colour = 'black',face="bold",vjust = 0,hjust = 0.5),#坐标轴标题
    axis.text.y.right = element_blank(),
    legend.position = 'top',
    legend.justification = c(0.5, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 24),
    plot.margin = margin(2,2,2,2),
    axis.line = element_line(colour = 'grey30',size = 0.2), 
    strip.text = element_text(
      size = 22,
      face = "bold",
      color = "#FFFFFF",
      vjust = 0.5,hjust = 0.5
    ),
    strip.background = element_rect(
      colour = 'black',
      fill = '#2878B5',
      size = 1
    )
  )
# Figs9 B -----------------------------------------------------------------
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

create_rotation_matrix <- function(angle_degrees) {
  angle_radians <- angle_degrees * (pi / 180)
  rotation_matrix <- matrix(c(cos(angle_radians), sin(angle_radians),
                              -sin(angle_radians), cos(angle_radians)),
                            nrow = 2, byrow = TRUE)
  return(rotation_matrix)
}

netVisual_circle_self <-
  function (net,
            color.use = NULL,
            title.name = NULL,
            sources.use = NULL,
            targets.use = NULL,
            idents.use = NULL,
            remove.isolate = FALSE,
            top = 1,
            weight.scale = FALSE,
            vertex.weight = 20,
            vertex.weight.max = NULL,
            vertex.size.max = NULL,
            vertex.label.cex = 1,
            vertex.label.color = "black",
            edge.weight.max = NULL,
            edge.width.max = 8,
            alpha.edge = 0.6,
            label.edge = FALSE,
            edge.label.color = "black",
            edge.label.cex = 0.8,
            edge.curved = 0.2,
            shape = "circle",
            layout = in_circle(),
            margin = 0.2,
            vertex.size = NULL,
            arrow.width = 1,
            arrow.size = 0.2,
            text.x = 0,
            text.y = 1.5,
            vertex.size.min = NULL,
            cell_type_order = NULL,
            angle_coords = NULL,
            label_dist = 2){
    if (!is.null(vertex.size)) {
      warning("'vertex.size' is deprecated. Use `vertex.weight`")
    }
    if (is.null(vertex.size.max)) {
      if (length(unique(vertex.weight)) == 1) {
        vertex.size.max <- 5
      }else {
        vertex.size.max <- 15
      }
    }
    options(warn = -1)
    thresh <- stats::quantile(net, probs = 1 - top)
    net[net < thresh] <- 0
    if ((!is.null(sources.use)) | (!is.null(targets.use)) |
        (!is.null(idents.use))) {
      if (is.null(rownames(net))) {
        stop("The input weighted matrix should have rownames!")
      }
      cells.level <- rownames(net)
      df.net <- reshape2::melt(net, value.name = "value")
      colnames(df.net)[1:2] <- c("source", "target")
      if (!is.null(sources.use)) {
        if (is.numeric(sources.use)) {
          sources.use <- cells.level[sources.use]
        }
        df.net <- subset(df.net, source %in% sources.use)
      }
      if (!is.null(targets.use)) {
        if (is.numeric(targets.use)) {
          targets.use <- cells.level[targets.use]
        }
        df.net <- subset(df.net, target %in% targets.use)
      }
      if (!is.null(idents.use)) {
        if (is.numeric(idents.use)) {
          idents.use <- cells.level[idents.use]
        }
        df.net <- filter(df.net, (source %in% idents.use) |
                           (target %in% idents.use))
      }
      df.net$source <- factor(df.net$source, levels = cells.level)
      df.net$target <- factor(df.net$target, levels = cells.level)
      df.net$value[is.na(df.net$value)] <- 0
      net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                            df.net[["target"]]), sum)
    }
    net[is.na(net)] <- 0
    if (remove.isolate) {
      idx1 <- which(Matrix::rowSums(net) == 0)
      idx2 <- which(Matrix::colSums(net) == 0)
      idx <- intersect(idx1, idx2)
      if (length(idx)>0) {
        net <- net[-idx,]
        net <- net[,-idx]
      }
    }
    if(length(cell_type_order) == 0){
      cell_type_order <- rownames(net)
    }else{
      cell_type_order <- cell_type_order[cell_type_order %in% colnames(net)]
    }
    net <- net[cell_type_order,cell_type_order]
    g <- graph_from_adjacency_matrix(net, mode = "directed",
                                     weighted = T)
    edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
    coords <- layout_(g, layout)
    if(!is.null(angle_coords)){
      rotation_matrix <- create_rotation_matrix(angle_coords)
      coords_rotated <- coords %*% rotation_matrix
      if (nrow(coords_rotated) != 1) {
        coords_scale = scale(coords_rotated)
      } else {
        coords_scale <- coords_rotated
      }
    }else{
      if (nrow(coords) != 1) {
        coords_scale = scale(coords)
      }else {
        coords_scale <- coords
      }
    }
    
    vertex_angles <- atan2(coords_scale[,2], coords_scale[,1])
    vertex_angles <- (vertex_angles + 2*pi) %% (2*pi)
    
    label.locs <- vertex_angles
    
    if (is.null(color.use)) {
      color.use = scPalette(length(igraph::V(g)))
    }
    if (is.null(vertex.weight.max)) {
      vertex.weight.max <- max(vertex.weight)
    }
    vertex.weight <-
      (vertex.weight[cell_type_order] %>% as.numeric()) / vertex.weight.max * vertex.size.max +
      vertex.size.max
    igraph::V(g)$size <- vertex.weight
    igraph::V(g)$color <- color.use[igraph::V(g)]
    igraph::V(g)$frame.color <- color.use[igraph::V(g)]
    igraph::V(g)$label.color <- vertex.label.color
    igraph::V(g)$label.cex <- vertex.label.cex
    if (label.edge) {
      igraph::E(g)$label <- igraph::E(g)$weight
      igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
    }
    if (is.null(edge.weight.max)) {
      edge.weight.max <- max(igraph::E(g)$weight)
    }
    if (weight.scale == TRUE) {
      igraph::E(g)$width <- 0.3 + igraph::E(g)$weight / edge.weight.max *
        edge.width.max
    }else {
      igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
    }
    igraph::E(g)$arrow.width <- arrow.width
    igraph::E(g)$arrow.size <- arrow.size
    igraph::E(g)$label.color <- edge.label.color
    igraph::E(g)$label.cex <- edge.label.cex
    igraph::E(g)$color <-
      grDevices::adjustcolor(igraph::V(g)$color[edge.start[,
                                                           1]], alpha.edge)
    label.dist <- vertex.weight / max(vertex.weight) + label_dist
    plot(
      g,
      edge.curved = edge.curved,
      vertex.shape = shape,
      layout = coords_scale,
      margin = margin,
      vertex.label.dist = label.dist,
      vertex.label.degree = label.locs,
      vertex.label.family = "Helvetica",
      edge.label.family = "Helvetica"
    )
    if (!is.null(title.name)) {
      text(text.x, text.y, title.name, cex = 2)
    }
    gg <- recordPlot()
    return(gg)
  }

netVisual_diffInteraction_self <- function(object, comparison = c(1, 2), measure = c("count", 
                                                                                     "weight", "count.merged", "weight.merged"), color.use = NULL, 
                                           color.edge = c("#b2182b", "#2166ac"), title.name = NULL, 
                                           sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
                                           top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
                                           vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black", 
                                           edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
                                           label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
                                           edge.curved = 0.2, shape = "circle", layout = in_circle(), 
                                           margin = 0.2, arrow.width = 1, arrow.size = 0.2,
                                           cells.level = cells.level,angle_coords = NULL) 
{
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level = cells.level
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top, 
                                 na.rm = T)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if(!is.null(angle_coords)){
    rotation_matrix <- create_rotation_matrix(angle_coords)
    coords_rotated <- coords %*% rotation_matrix
    if (nrow(coords_rotated) != 1) {
      coords_scale = scale(coords_rotated)
    } else {
      coords_scale <- coords_rotated
    }
  }else{
    if (nrow(coords) != 1) {
      coords_scale = scale(coords)
    }else {
      coords_scale <- coords
    }
  }
  
  vertex_angles <- atan2(coords_scale[,2], coords_scale[,1])
  vertex_angles <- (vertex_angles + 2*pi) %% (2*pi)
  label.locs <- vertex_angles
  
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, 
                       -atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g), 
                                                                        1]), pi - atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g), 
                                                                                                                                  1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  igraph::E(g)$loop.angle <- 0
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
                                                                1])] <- loop.angle[edge.start[which(edge.start[, 
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)),
                               direction = -1, start = 90)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}

scRNA_11 <- readRDS('D:/jiaoxi/scRNA_11_res.rds')
scRNA_use <- subset(scRNA_11,subset = anno_res %in% c(
  "Cancer cell","Fibroblast",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "Dysfunction_T cell",
  "Cytolytic_T cell",
  "Naive T","cDC1_Xcr1",
  "Proliferative CD8+ T","CD8+ T memory",
  "Treg",
  "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "Mast_Mcpt2","B cell",
  "Mono_Isg15",
  "Proliferative Treg","CD4+ T memory",
  "tDC","Gzmb_T cell"
))
rm(scRNA_11)
scRNA_use$res_use <- scRNA_use$anno_res
scRNA_use$res_use %>% unique()
scRNA_use$batch %>% unique()
group.cellType <- c(
  "Cancer cell","Fibroblast",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "Dysfunction_T cell",
  "Cytolytic_T cell",
  "Naive T","cDC1_Xcr1",
  "Proliferative CD8+ T","CD8+ T memory",
  "Treg",
  "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "Mast_Mcpt2","B cell",
  "Mono_Isg15",
  "Proliferative Treg","CD4+ T memory",
  "tDC","Gzmb_T cell"
) 
names(group.cellType) <- group.cellType
scRNA_use$res_low <- group.cellType[scRNA_use$res_use] %>% as.character()
table(scRNA_use$res_low == scRNA_use$anno_res)
scRNA_use$res_low[scRNA_use$res_low == 'NPBNs'] = 
  'S100A8-Neu'
scRNA_use$res_low[scRNA_use$res_low == 'Exhausted TAN'] = 
  'Ccl4-Neu'
scRNA_use$res_low[scRNA_use$res_low == 'interferon-stimulated NAN'] = 
  'ISG-Neu'
celltype <- 
  scRNA_use$res_low %>% unique() %>% .[!is.na(.)]
celltype
my36colors <-c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) %>% .[1:length(celltype)]
names(my36colors) <- celltype

scRNA_mock <- subset(scRNA_use,subset = batch ==  'MOCK')
Idents(scRNA_mock) <- scRNA_mock$res_low
scRNA_C5aRa <- subset(scRNA_use,subset = batch ==  'C5aRa')
Idents(scRNA_C5aRa) <- scRNA_C5aRa$res_low
scRNA_C5aRA_PD1 <- subset(scRNA_use,subset = batch ==  'C5aRA_PD1')
Idents(scRNA_C5aRA_PD1) <- scRNA_C5aRA_PD1$res_low
scRNA_PD1 <- subset(scRNA_use,subset = batch ==  'PD1')
Idents(scRNA_PD1) <- scRNA_PD1$res_low

scRNA_mock$anno_res %>% table()
data.input <- GetAssayData(scRNA_mock, assay = "RNA", slot = "data") 
labels <- Idents(scRNA_mock)
identity <- data.frame(group = labels, row.names = names(labels)) 
cellchat <-
  createCellChat(object = data.input,
                 meta = scRNA_mock@meta.data,
                 group.by = "res_low")
cellchat <- setIdent(cellchat, ident.use = "res_low")  
levels(cellchat@idents)  
groupSize <- as.numeric(table(cellchat@idents))  

CellChatDB <- CellChatDB.mouse 
colnames(CellChatDB$interaction) 
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(
  CellChatDB, 
  search = c("Cell-Cell Contact",'Secreted Signaling')
)
dplyr::glimpse(CellChatDB$interaction) 
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)   
cellchat@LR$LRsig %>% head()

cellchat <- computeCommunProb(cellchat)  
cellchat <- filterCommunication(cellchat, min.cells = 5)
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(
  cellchat,
  sources.use = c(1, 2),
  targets.use = c(4, 5))
cellchat <- computeCommunProbPathway(cellchat,thresh = 0.05)
cellchat <- aggregateNet(cellchat,thresh = 0.05) 
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(
  cellchat@net$weight, 
  vertex.weight = groupSize, 
  weight.scale = T, 
  label.edge= F, 
  title.name = "Interaction weights/strength")
cellchat_mock <- cellchat
df.net_mock <- subsetCommunication(cellchat_mock,thresh = 0.05)

data.input <- GetAssayData(scRNA_C5aRa, assay = "RNA", slot = "data")  
scRNA_C5aRa$res_use %>% table()
labels <- Idents(scRNA_C5aRa)
identity <- data.frame(group = labels, row.names = names(labels))  
cellchat <- createCellChat(object = data.input, meta = scRNA_C5aRa@meta.data,group.by = "res_low")
cellchat <- setIdent(cellchat, ident.use = "res_low") 
levels(cellchat@idents)  
groupSize <- as.numeric(table(cellchat@idents))  
CellChatDB <- CellChatDB.mouse 
colnames(CellChatDB$interaction) 
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact",'Secreted Signaling')) 
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)  
cellchat@LR$LRsig %>% head()


cellchat <- computeCommunProb(cellchat)  
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat,thresh = 0.05)

cellchat <- aggregateNet(cellchat,thresh = 0.05) 
groupSize <- as.numeric(table(cellchat@idents))
cellchat_C5aRa <- cellchat
df.net_C5aRa <- subsetCommunication(cellchat_C5aRa,thresh = 0.05)


data.input <-
  GetAssayData(scRNA_PD1, assay = "RNA", slot = "data")  
labels <- Idents(scRNA_PD1)
identity <-
  data.frame(group = labels, row.names = names(labels))  
cellchat <-
  createCellChat(object = data.input,
                 meta = scRNA_PD1@meta.data,
                 group.by = "res_low")
cellchat <-
  setIdent(cellchat, ident.use = "res_low")  
levels(cellchat@idents) 
groupSize <-
  as.numeric(table(cellchat@idents))  
CellChatDB <- CellChatDB.mouse
colnames(CellChatDB$interaction)  
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
CellChatDB.use <-
  subsetDB(CellChatDB, search = c("Cell-Cell Contact",'Secreted Signaling'))
cellchat@DB <- CellChatDB.use

cellchat <-
  subsetData(cellchat) 
cellchat <-
  identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <-
  projectData(cellchat, PPI.mouse)   
cellchat@LR$LRsig %>% head()


cellchat <- computeCommunProb(cellchat)  
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat,thresh = 0.05)
cellchat <- aggregateNet(cellchat,thresh = 0.05)
groupSize <- as.numeric(table(cellchat@idents))
cellchat_PD1 <- cellchat
df.net_PD1 <- subsetCommunication(cellchat_PD1,thresh = 0.05)

data.input <-
  GetAssayData(
    scRNA_C5aRA_PD1, 
    assay = "RNA", 
    slot = "data")  
labels <- Idents(scRNA_C5aRA_PD1)
identity <-
  data.frame(
    group = labels, 
    row.names = names(labels))  
cellchat <- createCellChat(
  object = data.input, 
  meta = scRNA_C5aRA_PD1@meta.data,
  group.by = "res_low")
cellchat <- setIdent(cellchat, ident.use = "res_low")  
levels(cellchat@idents)  
groupSize <- as.numeric(table(cellchat@idents)) 

CellChatDB <- CellChatDB.mouse 
colnames(CellChatDB$interaction)  
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact",'Secreted Signaling')) 
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)   
cellchat@LR$LRsig %>% head()


cellchat <- computeCommunProb(cellchat)  
cellchat <- filterCommunication(cellchat, min.cells = 5)

cellchat <- computeCommunProbPathway(cellchat,thresh = 1)
cellchat <- aggregateNet(cellchat,thresh = 0.05) 
cellchat_C5aRa_PD1 <- cellchat
df.net_C5aRa_PD1 <- subsetCommunication(cellchat_C5aRa_PD1,thresh = 0.05)

group.new <- levels(cellchat_C5aRa_PD1@idents)
unique(scRNA_use$res_use)
cellchat_C5aRa_PD1_mult <- liftCellChat(object = cellchat_C5aRa_PD1,group.new = group.new)
cellchat_PD1_mult <- liftCellChat(object = cellchat_PD1,group.new = group.new)
cellchat_C5aRa_mult <- liftCellChat(object = cellchat_C5aRa,group.new = group.new)
cellchat_mock_mult <- liftCellChat(object = cellchat_mock,group.new = group.new)
object.list <- list(
  PD1 = cellchat_PD1_mult,
  C5aRa_PD1 = cellchat_C5aRa_PD1_mult,
  mock = cellchat_mock_mult,
  C5aRa = cellchat_C5aRa_mult
)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
{
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_diffInteraction_self(
    object = cellchat,
    comparison = c(3,4),
    sources.use = c(
      "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
      "Proliferative CD8+ T","CD8+ T memory",
      'S100A8-Neu','Ccl4-Neu','ISG-Neu'
    ),
    targets.use = c(
      "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
      "Proliferative CD8+ T","CD8+ T memory",
      'S100A8-Neu','Ccl4-Neu','ISG-Neu'
    ),
    cells.level = c(
      "Cancer cell","Fibroblast",
      "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
      "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
      "Mac/Mono_Mki67",
      "S100A8-Neu","Ccl4-Neu","ISG-Neu",
      "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
      "Proliferative CD8+ T","CD8+ T memory",
      "CD4+ T memory","Treg","Naive T","Proliferative Treg",
      "cDC1_Xcr1","tDC","B cell","Mast_Mcpt2"
    ),
    remove.isolate = TRUE,
    title.name = 'C5aRA vs Mock',
    weight.scale = T)
  netVisual_diffInteraction_self(
    object = cellchat, 
    comparison = c(1,2),
    sources.use = c(
      "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
      "Proliferative CD8+ T","CD8+ T memory",
      'S100A8-Neu','Ccl4-Neu','ISG-Neu'
    ),
    targets.use = c(
      "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
      "Proliferative CD8+ T","CD8+ T memory",
      'S100A8-Neu','Ccl4-Neu','ISG-Neu'
    ),
    remove.isolate = TRUE,
    title.name = 'C5aRA_PD1 vs PD1',
    cells.level = c(
      "Cancer cell","Fibroblast",
      "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
      "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
      "Mac/Mono_Mki67",
      "S100A8-Neu","Ccl4-Neu","ISG-Neu",
      "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
      "Proliferative CD8+ T","CD8+ T memory",
      "CD4+ T memory","Treg","Naive T","Proliferative Treg",
      "cDC1_Xcr1","tDC","B cell","Mast_Mcpt2"
    ),
    weight.scale = T)
}

# Figs9 C -----------------------------------------------------------------
single_cluster_fun <- function(
    source_cell = source_cell,
    mfrow = c(2,2),
    target_cell = target_cell,
    cellchat_mock = NULL,
    cellchat_C5aRa = NULL,
    cellchat_PD1 = NULL,
    cellchat_C5aRa_PD1 = NULL){
  par(mfrow = mfrow, xpd=TRUE,oma = c(0, 0, 0, 0),mar = c(0, 0, 0, 2))
  if(!is.null(cellchat_mock)){
    cell_type_order <- c(source_cell,target_cell) %>% 
      .[. %in% rownames(cellchat_mock@net$count)]
    color.use <- my36colors[cell_type_order]
    if(source_cell %in% rownames(cellchat_mock@net$count)){
      netVisual_circle_self(
        net = cellchat_mock@net$count, 
        cell_type_order = cell_type_order,
        angle_coords = -45,
        color.use = color.use,
        vertex.weight = table(cellchat_mock@idents)[cell_type_order],
        vertex.label.cex = 1.8,
        vertex.size.max = 15,
        weight.scale = T, 
        label.edge= TRUE, 
        edge.label.cex = 2,
        sources.use = source_cell,
        targets.use = target_cell,
        remove.isolate = TRUE,
        arrow.width = 2,
        arrow.size = 1,
        title.name = 'Mock',
        margin = c(0,0,0.5,0)
      )
    } 
  }
  if(!is.null(cellchat_C5aRa)){
    cell_type_order <- c(source_cell,target_cell) %>% 
      .[. %in% rownames(cellchat_C5aRa@net$count)]
    color.use <- my36colors[cell_type_order]
    if (source_cell %in% rownames(cellchat_C5aRa@net$count)) {
      netVisual_circle_self(
        net = cellchat_C5aRa@net$count,
        cell_type_order = cell_type_order,
        angle_coords = -45,
        color.use = color.use,
        vertex.weight = table(cellchat_C5aRa@idents)[cell_type_order],
        vertex.label.cex = 1.8,
        vertex.size.max = 15,
        weight.scale = T,
        label.edge = TRUE,
        edge.label.cex = 2,
        sources.use = source_cell,
        targets.use = target_cell,
        remove.isolate = TRUE,
        arrow.width = 2,
        arrow.size = 1,
        title.name = 'C5aRa',
        margin = c(0, 0, 0.5, 0)
      )
    }
  }
  if(!is.null(cellchat_PD1)){
    cell_type_order <- c(source_cell,target_cell) %>% 
      .[. %in% rownames(cellchat_PD1@net$count)]
    color.use <- my36colors[cell_type_order]
    if (source_cell %in% rownames(cellchat_PD1@net$count)){
      netVisual_circle_self(
        net = cellchat_PD1@net$count, 
        cell_type_order = cell_type_order,
        angle_coords = -45,
        color.use = color.use,
        vertex.weight = table(cellchat_PD1@idents)[cell_type_order],
        vertex.label.cex = 1.8,
        vertex.size.max = 15,
        weight.scale = T, 
        label.edge= TRUE, 
        edge.label.cex = 2,
        sources.use = source_cell,
        targets.use = target_cell,
        remove.isolate = TRUE,
        arrow.width = 2,
        arrow.size = 1,
        title.name = 'PD1',
        margin = c(0,0,0.5,0)
      )
    }
  }
  if(!is.null(cellchat_C5aRa_PD1)){
    cell_type_order <- c(source_cell,target_cell) %>% 
      .[. %in% rownames(cellchat_C5aRa_PD1@net$count)]
    color.use <- my36colors[cell_type_order]
    if (source_cell %in% rownames(cellchat_C5aRa_PD1@net$count)){
      netVisual_circle_self(
        net = cellchat_C5aRa_PD1@net$count, 
        cell_type_order = cell_type_order,
        angle_coords = -45,
        color.use = color.use,
        vertex.weight = table(cellchat_C5aRa_PD1@idents)[cell_type_order],
        vertex.label.cex = 1.8,
        vertex.size.max = 15,
        weight.scale = T, 
        label.edge= TRUE, 
        edge.label.cex = 2,
        sources.use = source_cell,
        targets.use = target_cell,
        remove.isolate = TRUE,
        arrow.width = 1,
        arrow.size = 1,
        title.name = 'C5aRa_PD1',
        margin = c(0,0,0.5,0)
      )
    }
  }
}
target_cell <- c(
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "CD8+ T memory",
  "Naive T","Treg",
  "CD4+ T memory"
)
single_cluster_fun(
  source_cell = 'ISG-Neu',
  target_cell = target_cell,
  mfrow = c(1,1),
  cellchat_mock = NULL,
  cellchat_C5aRa = NULL,
  cellchat_PD1 = NULL,
  cellchat_C5aRa_PD1 = cellchat_C5aRa_PD1
)
# Figs9 D(left) -----------------------------------------------------------------
res_plot_self <- function(scRNA_11_sub,
                          diffgene_source = 'diffgroup',
                          celltype_diffgene = NULL,
                          scRNA_all = NULL,
                          diff_gene_from = 'double',
                          upstream_ligands_n = 20,
                          top_n_target_gene = 100,
                          target_gene = NULL,
                          target_ligands = NULL,
                          num_show_gene = 15,
                          min.pct = 0.1,
                          p_ligand_target_network_x_size = 20,
                          axis_text_y = 30,
                          width_fig = c(1,4),
                          axis_title_y = '',
                          axis_title_x = ''){
  organism = "mouse"
  if(organism == "human"){
    lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
    weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
  } else if(organism == "mouse"){
    lr_network = readRDS('E:/nicheneter/lr_network_mouse_21122021.rds')
    ligand_target_matrix = readRDS('E:/nicheneter/ligand_target_matrix_nsga2r_final_mouse.rds')
    weighted_networks = readRDS('E:/nicheneter/weighted_networks_nsga2r_final_mouse.rds')
    
  }
  lr_network = lr_network %>% distinct(from, to)
  head(lr_network)
  ligand_target_matrix[1:5,1:5] 
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
  head(weighted_networks$lr_sig) 
  
  head(weighted_networks$gr)
  
  
  seuratObj = alias_to_symbol_seurat(scRNA_11_sub, "mouse")
  Idents(seuratObj) <- seuratObj$anno_res
  expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  list_expressed_genes_sender = sender_celltypes %>% 
    unique() %>% 
    lapply(get_expressed_genes, seuratObj, 0.10) 
  expressed_genes_sender = list_expressed_genes_sender %>%
    unlist() %>% 
    unique()
  if(diffgene_source == 'diffgroup'){
    seurat_obj_receiver= subset(seuratObj, idents = receiver)
    seurat_obj_receiver = SetIdent(
      seurat_obj_receiver, 
      value = seurat_obj_receiver[["batch", drop = TRUE]])
    DE_table_receiver = FindMarkers(
      object = seurat_obj_receiver, 
      ident.1 = condition_oi, 
      ident.2 = condition_reference,
      min.cells.group = 1, 
      min.pct = min.pct) %>% 
      rownames_to_column("gene")
    
    if(diff_gene_from == 'double'){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }else if(diff_gene_from == condition_oi){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & avg_log2FC >= 0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }else if(diff_gene_from == condition_reference){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & avg_log2FC <= -0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }
  }else if(diffgene_source == 'diffcluster'){
    seurat_obj_receiver= subset(
      scRNA_all, subset = anno_res %in% celltype_diffgene
    ) %>% subset(subset = batch =='C5aRA_PD1')
    seurat_obj_receiver = SetIdent(
      seurat_obj_receiver, 
      value = seurat_obj_receiver[["anno_res", drop = TRUE]])
    
    DE_table_receiver = FindMarkers(
      object = seurat_obj_receiver, 
      ident.1 = receiver,
      min.cells.group = 1,
      min.pct = 0.10) %>% 
      rownames_to_column("gene")
    geneset_oi = DE_table_receiver %>% 
      filter(p_val <= 0.05 & avg_log2FC >= 0.25) %>%
      pull(gene)
    geneset_oi = geneset_oi %>% 
      .[. %in% rownames(ligand_target_matrix)]
  }
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% 
    filter(from %in% expressed_ligands & 
             to %in% expressed_receptors) %>% 
    pull(from) %>% unique()
  ligand_activities = predict_ligand_activities(
    geneset = geneset_oi, 
    background_expressed_genes = background_expressed_genes, 
    ligand_target_matrix = ligand_target_matrix, 
    potential_ligands = potential_ligands)
  
  ligand_activities = ligand_activities %>% 
    arrange(-aupr_corrected) %>% 
    mutate(rank = rank(plyr::desc(aupr_corrected)))
  best_upstream_ligands = ligand_activities %>%
    top_n(upstream_ligands_n, aupr_corrected) %>%
    arrange(-aupr_corrected) %>%
    pull(test_ligand) %>% unique()
  active_ligand_target_links_df = 
    best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi, 
           ligand_target_matrix = ligand_target_matrix, n = top_n_target_gene) %>% 
    bind_rows() %>% drop_na()
  
  active_ligand_target_links = 
    prepare_ligand_target_visualization(
      ligand_target_df = active_ligand_target_links_df, 
      ligand_target_matrix = ligand_target_matrix, 
      cutoff = 0.25)
  
  order_ligands = intersect(
    best_upstream_ligands, 
    colnames(active_ligand_target_links)) %>% 
    rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% 
    unique() %>% 
    intersect(rownames(active_ligand_target_links)) %>% 
    make.names()
  rownames(active_ligand_target_links) = 
    rownames(active_ligand_target_links) %>% 
    make.names() 
  colnames(active_ligand_target_links) = 
    colnames(active_ligand_target_links) %>% 
    make.names()
  
  vis_ligand_target = 
    active_ligand_target_links[order_targets,order_ligands] %>% t()
  if(ncol(vis_ligand_target)>num_show_gene){vis_ligand_target <- vis_ligand_target[,1:num_show_gene]}
  p_ligand_target_network = 
    make_heatmap_ggplot_self(
      matrix = vis_ligand_target,
      y_name ="Prioritized ligands",
      x_name = "Predicted target genes", 
      color = "purple",
      tile_color = 'white',
      legend_position = "top", 
      x_axis_position = "top",
      legend_title = "Regulatory potential")  + 
    theme(
      panel.border = element_rect(fill = NA),
      axis.text.y = element_text(face = "italic",size = 20),
      axis.title = element_text(face = "italic",size = 22),
      axis.title.x = element_blank(),
      legend.title = element_text(face = "italic",size = 22),
      legend.text = element_text(face = "italic",size = 12)
    )
  if(!is.null(target_gene)){
    p_ligand_target_network <- 
      p_ligand_target_network+
      theme(
        axis.text.x = element_text(
          face = "italic",
          size = p_ligand_target_network_x_size,
          hjust = 0.5,
          vjust = 0.5,
          color = ifelse(colnames(vis_ligand_target) %in% target_gene, "red", "black"))
      )
  }
  p_ligand_target_network
  lr_network_top = 
    lr_network %>% 
    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% 
    distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  
  rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
  colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()
  
  vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
  p_ligand_aupr = make_heatmap_ggplot_self(
    matrix = vis_ligand_aupr,
    y_name = paste("Prioritized ligands in\n",axis_title_y),
    x_name = "Ligand activity",
    color = "darkorange",
    color_low = '#FFCF95',
    tile_color = 'white',
    legend_position = "top",
    x_axis_position = "top",
    x_axis = FALSE,
    legend_title = "AUPR\n(target gene prediction ability)  "
  ) +
    theme(
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "italic",size = 32),
      legend.title = element_text(face = "italic",size = 22),
      axis.text.y = element_text(face = "italic",size = axis_text_y),
      legend.text = element_text(face = "italic",size = 18)
    )
  if(!is.null(target_ligands)){
    p_ligand_aupr <- p_ligand_aupr + 
      theme(
        axis.text.y = element_text(
          face = "italic",
          size = axis_text_y,
          color = ifelse(rownames(vis_ligand_aupr) %in% target_ligands, "red", "black"))
      )
  }
  order_ligands_adapted <- str_replace_all(order_ligands, "\\.", "-")
  
  figures_without_legend = cowplot::plot_grid(
    p_ligand_aupr +
      theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        plot.margin = margin(t = 5,r = 1,b = 5,l = 0,unit = "pt")
      ),
    p_ligand_target_network +
      theme(legend.position = "none",
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(t = 5,r = 5,b = 5,l = 5,unit = "pt")
      ),
    align = "h",
    axis = 'l',
    nrow = 1,
    rel_widths = width_fig
  )
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, 1.5))
  
  combined_plot = cowplot::plot_grid(
    figures_without_legend, 
    legends, 
    rel_heights = c(10,2), 
    nrow = 2, 
    align = "hv")
  if(diff_gene_from == 'double'){
    title = paste(condition_oi,' VS ',condition_reference,' differential gene of ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25)
  }else if(diff_gene_from == condition_oi){
    title = paste('Predicted target genes in ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 45)
  }else if(diff_gene_from == condition_reference){
    title = paste('Predicted target genes in ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25)
  }
  return(final_plot)
}
sender_celltypes = c(
  "interferon-stimulated NAN"
)
receiver = c(
  "Gzmb_T cell",
  "Dysfunction_T cell",
  "Cytolytic_T cell"
)

condition_oi = "MOCK"
condition_reference = "C5aRa"
scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% c(
  sender_celltypes,
  receiver
)) %>% 
  subset(subset = batch %in% c(condition_oi,condition_reference))
scRNA_11_sub@meta.data %>% group_by(batch,anno_res) %>% summarise(Counts = n())
target_gene <- c("Bhlhe40","Ccl3","Ccl5","Hif1a","Ccl4","H2.T23","Il10ra",
                 "Fos","Ier3","Stat1",'Jun')
target_ligands <- c('Icam1',"Cxcl10","Cd80")
p3 <- res_plot_self(
  scRNA_11_sub = scRNA_11_sub,
  diff_gene_from = condition_reference,
  upstream_ligands_n = 10,
  top_n_target_gene = 100,
  target_gene = target_gene,
  num_show_gene = 40,#9,
  p_ligand_target_network_x_size = 30,
  target_ligands = target_ligands,
  axis_text_y = 20,
  width_fig = c(1,4.5),
  axis_title_x = 'Gzmb+/Dysfunction/Cytolytic T cellS',
  axis_title_y = "Neu-Isg cell")
p3
# Figs9 D(right) -----------------------------------------------------------------
res_plot_self_2 <- function(scRNA_12_sub,
                            diffgene_source = 'diffgroup',
                            celltype_diffgene = NULL,
                            scRNA_all = NULL,
                            diff_gene_from = 'double',
                            upstream_ligands_n = 20,
                            top_n_target_gene = 100,
                            target_gene = NULL,
                            target_ligands = NULL,
                            num_show_gene = 15,
                            min.pct = 0.1,
                            p_ligand_target_network_x_size = 20,
                            axis_text_y = 30,
                            width_fig = c(1,4),
                            axis_title_y = '',
                            axis_title_x = ''){
  organism = "mouse"
  if(organism == "human"){
    lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
    weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
  } else if(organism == "mouse"){
    lr_network = readRDS('E:/nicheneter/lr_network_mouse_21122021.rds')
    ligand_target_matrix = readRDS('E:/nicheneter/ligand_target_matrix_nsga2r_final_mouse.rds')
    weighted_networks = readRDS('E:/nicheneter/weighted_networks_nsga2r_final_mouse.rds')
    
  }
  lr_network = lr_network %>% distinct(from, to)
  head(lr_network)
  ligand_target_matrix[1:5,1:5] 
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
  head(weighted_networks$lr_sig)  
  
  head(weighted_networks$gr) 
  
  
  seuratObj = alias_to_symbol_seurat(scRNA_12_sub, "mouse")
  Idents(seuratObj) <- seuratObj$anno_res_3
  expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  ## sender
  list_expressed_genes_sender = sender_celltypes %>% 
    unique() %>% 
    lapply(get_expressed_genes, seuratObj, 0.10) 
  expressed_genes_sender = list_expressed_genes_sender %>%
    unlist() %>% 
    unique()
  
  if(diffgene_source == 'diffgroup'){
    seurat_obj_receiver= subset(seuratObj, idents = receiver)
    seurat_obj_receiver = SetIdent(
      seurat_obj_receiver, 
      value = seurat_obj_receiver[["batch", drop = TRUE]])
    DE_table_receiver = FindMarkers(
      object = seurat_obj_receiver, 
      ident.1 = condition_oi, 
      ident.2 = condition_reference,
      min.cells.group = 1, 
      min.pct = min.pct) %>% 
      rownames_to_column("gene")
    
    if(diff_gene_from == 'double'){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }else if(diff_gene_from == condition_oi){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & avg_log2FC >= 0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }else if(diff_gene_from == condition_reference){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & avg_log2FC <= -0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }
  }else if(diffgene_source == 'diffcluster'){
    seurat_obj_receiver= subset(
      scRNA_all, subset = anno_res_3 %in% celltype_diffgene
    ) %>% subset(subset = batch =='C5aRA_PD1')
    seurat_obj_receiver = SetIdent(
      seurat_obj_receiver, 
      value = seurat_obj_receiver[["anno_res_3", drop = TRUE]])
    
    DE_table_receiver = FindMarkers(
      object = seurat_obj_receiver, 
      ident.1 = receiver,
      min.cells.group = 1,
      min.pct = min.pct) %>% 
      rownames_to_column("gene")
    geneset_oi = DE_table_receiver %>% 
      filter(p_val <= 0.05 & avg_log2FC >= 0.25) %>%
      pull(gene)
    geneset_oi = geneset_oi %>% 
      .[. %in% rownames(ligand_target_matrix)]
  }
  
  
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% 
    filter(from %in% expressed_ligands & 
             to %in% expressed_receptors) %>% 
    pull(from) %>% unique()
  ligand_activities = predict_ligand_activities(
    geneset = geneset_oi, 
    background_expressed_genes = background_expressed_genes, 
    ligand_target_matrix = ligand_target_matrix, 
    potential_ligands = potential_ligands)
  
  ligand_activities = ligand_activities %>% 
    arrange(-aupr_corrected) %>% 
    mutate(rank = rank(plyr::desc(aupr_corrected)))
  best_upstream_ligands = ligand_activities %>%
    top_n(upstream_ligands_n, aupr_corrected) %>%
    arrange(-aupr_corrected) %>%
    pull(test_ligand) %>% unique()
  active_ligand_target_links_df = 
    best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi, 
           ligand_target_matrix = ligand_target_matrix, n = top_n_target_gene) %>% 
    bind_rows() %>% drop_na()
  
  active_ligand_target_links = 
    prepare_ligand_target_visualization(
      ligand_target_df = active_ligand_target_links_df, 
      ligand_target_matrix = ligand_target_matrix, 
      cutoff = 0.25)
  
  order_ligands = intersect(
    best_upstream_ligands, 
    colnames(active_ligand_target_links)) %>% 
    rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% 
    unique() %>% 
    intersect(rownames(active_ligand_target_links)) %>% 
    make.names()
  rownames(active_ligand_target_links) = 
    rownames(active_ligand_target_links) %>% 
    make.names() 
  colnames(active_ligand_target_links) = 
    colnames(active_ligand_target_links) %>% 
    make.names()  
  
  vis_ligand_target = 
    active_ligand_target_links[order_targets,order_ligands] %>% t()
  if(ncol(vis_ligand_target)>num_show_gene){vis_ligand_target <- vis_ligand_target[,1:num_show_gene]}
  p_ligand_target_network = 
    make_heatmap_ggplot_self(
      matrix = vis_ligand_target,
      y_name ="Prioritized ligands",
      x_name = "Predicted target genes", 
      color = "purple",
      tile_color = 'white',
      legend_position = "top", 
      x_axis_position = "top",
      legend_title = "Regulatory potential")  + 
    theme(
      panel.border = element_rect(fill = NA),
      axis.text.y = element_text(face = "italic",size = 20),
      axis.title = element_text(face = "italic",size = 22),
      axis.title.x = element_blank(),
      legend.title = element_text(face = "italic",size = 22),
      legend.text = element_text(face = "italic",size = 12)
    )
  if(!is.null(target_gene)){
    p_ligand_target_network <- 
      p_ligand_target_network+
      theme(
        axis.text.x = element_text(
          face = "italic",
          size = p_ligand_target_network_x_size,
          hjust = 0.5,
          vjust = 0.5,
          color = ifelse(colnames(vis_ligand_target) %in% target_gene, "red", "black"))
      )
  }
  p_ligand_target_network
  lr_network_top = 
    lr_network %>% 
    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% 
    distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  ligand_aupr_matrix = ligand_activities %>% dplyr::select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  
  rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
  colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()
  
  vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
  p_ligand_aupr = make_heatmap_ggplot_self(
    matrix = vis_ligand_aupr,
    y_name = paste("Prioritized ligands in\n",axis_title_y),
    x_name = "Ligand activity",
    color = "darkorange",
    color_low = '#FFCF95',
    tile_color = 'white',
    legend_position = "top",
    x_axis_position = "top",
    x_axis = FALSE,
    legend_title = "AUPR\n(target gene prediction ability)  "
  ) +
    theme(
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "italic",size = 32),
      legend.title = element_text(face = "italic",size = 22),
      axis.text.y = element_text(face = "italic",size = axis_text_y),
      # legend_title = element_text(face = "italic",size = 22),
      legend.text = element_text(face = "italic",size = 18)
    )
  if(!is.null(target_ligands)){
    p_ligand_aupr <- p_ligand_aupr + 
      theme(
        axis.text.y = element_text(
          face = "italic",
          size = axis_text_y,
          color = ifelse(rownames(vis_ligand_aupr) %in% target_ligands, "red", "black"))
      )
  }
  order_ligands_adapted <- str_replace_all(order_ligands, "\\.", "-")
  
  figures_without_legend = cowplot::plot_grid(
    p_ligand_aupr +
      theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        plot.margin = margin(t = 5,r = 1,b = 5,l = 0,unit = "pt")
      ),
    p_ligand_target_network +
      theme(legend.position = "none",
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            # axis.text.x = element_text(size = p_ligand_target_network_x_size),
            plot.margin = margin(t = 5,r = 5,b = 5,l = 5,unit = "pt")
      ),
    align = "h",
    axis = 'l',
    nrow = 1,
    rel_widths = width_fig
    )
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    # ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, 1.5))
  
  combined_plot = cowplot::plot_grid(
    figures_without_legend, 
    legends, 
    rel_heights = c(10,2), 
    nrow = 2, 
    align = "hv")
  # 添加标题
  if(diff_gene_from == 'double'){
    title = paste(condition_oi,' VS ',condition_reference,' differential gene of ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25)
  }else if(diff_gene_from == condition_oi){
    title = paste('Predicted target genes in ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 45)
  }else if(diff_gene_from == condition_reference){
    title = paste('Predicted target genes in ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25)
  }
  return(final_plot)
}
sender_celltypes = c(
  "S100A8/ISG-Neu",
  "Ccl4-Neu"
)
receiver = c(
  "Cytolytic/Gzmb_T cell",
  "Dysfunction_T cell"
)
condition_oi = "PD1"
condition_reference = "C5aRA_PD1"
scRNA_12_sub <- subset(scRNA_12,subset = anno_res_3 %in% c(
  sender_celltypes,
  receiver
)) %>% 
  subset(subset = batch %in% c(condition_oi,condition_reference))
scRNA_12_sub@meta.data %>% group_by(batch,anno_res_3) %>% summarise(Counts = n())
target_gene <- c("Icam1","Ifng","Ccl5","Ifnar2","Il7r","Irf1",
                 "Bhlhe40","Icos","Hif1a")
p3 <- res_plot_self_2(
  scRNA_12_sub = scRNA_12_sub,
  diff_gene_from = condition_reference,
  upstream_ligands_n = 5,
  top_n_target_gene = 100,
  target_gene = target_gene,
  target_ligands = target_ligands,
  num_show_gene = 36,
  axis_title_y = 'Neutrophil cells',
  axis_title_x = 'CD8+ T cells',#Macrophage/Monocyte
  p_ligand_target_network_x_size = 30)
p3
# Figs9 E -----------------------------------------------------------------
single_cluster_fun <- function(
    source_cell = source_cell,
    target_cell = target_cell,
    cellchat_mock = cellchat_mock,
    cellchat_C5aRa = cellchat_C5aRa,
    cellchat_PD1 = cellchat_PD1,
    cellchat_C5aRa_PD1 = cellchat_C5aRa_PD1){
  par(mfrow = c(2,2), xpd=TRUE,oma = c(0, 0, 0, 0),mar = c(0, 0, 0, 2))
  cell_type_order <- c(source_cell,target_cell) %>% 
    .[. %in% rownames(cellchat_mock@net$count)]
  color.use <- my36colors[cell_type_order]
  if(source_cell %in% rownames(cellchat_mock@net$count)){
    netVisual_circle_self(
      net = cellchat_mock@net$count, 
      cell_type_order = cell_type_order,
      angle_coords = -45,
      color.use = color.use,
      vertex.weight = table(cellchat_mock@idents)[cell_type_order],
      vertex.label.cex = 1.8,
      vertex.size.max = 15,
      weight.scale = T, 
      label.edge= TRUE, 
      edge.label.cex = 2,
      sources.use = source_cell,
      targets.use = target_cell,
      remove.isolate = TRUE,
      # top = 0.25,
      arrow.width = 2,
      arrow.size = 1,
      title.name = 'Mock',
      margin = c(0,0,0.5,0)
    )
  } 
  
  cell_type_order <- c(source_cell,target_cell) %>% 
    .[. %in% rownames(cellchat_C5aRa@net$count)]
  color.use <- my36colors[cell_type_order]
  if (source_cell %in% rownames(cellchat_C5aRa@net$count)) {
    netVisual_circle_self(
      net = cellchat_C5aRa@net$count,
      cell_type_order = cell_type_order,
      angle_coords = -45,
      color.use = color.use,
      vertex.weight = table(cellchat_C5aRa@idents)[cell_type_order],
      vertex.label.cex = 1.8,
      vertex.size.max = 15,
      weight.scale = T,
      label.edge = TRUE,
      edge.label.cex = 2,
      sources.use = source_cell,
      targets.use = target_cell,
      remove.isolate = TRUE,
      # top = 0.25,
      arrow.width = 2,
      arrow.size = 1,
      title.name = 'C5aRa',
      margin = c(0, 0, 0.5, 0)
    )
  }
  cell_type_order <- c(source_cell,target_cell) %>% 
    .[. %in% rownames(cellchat_PD1@net$count)]
  color.use <- my36colors[cell_type_order]
  if (source_cell %in% rownames(cellchat_PD1@net$count)){
    netVisual_circle_self(
      net = cellchat_PD1@net$count, 
      cell_type_order = cell_type_order,
      angle_coords = -45,
      color.use = color.use,
      vertex.weight = table(cellchat_PD1@idents)[cell_type_order],
      vertex.label.cex = 1.8,
      vertex.size.max = 15,
      weight.scale = T, 
      label.edge= TRUE, 
      edge.label.cex = 2,
      sources.use = source_cell,
      targets.use = target_cell,
      remove.isolate = TRUE,
      # top = 0.25,
      arrow.width = 2,
      arrow.size = 1,
      title.name = 'PD1',
      margin = c(0,0,0.5,0)
    )
  }
  cell_type_order <- c(source_cell,target_cell) %>% 
    .[. %in% rownames(cellchat_C5aRa_PD1@net$count)]
  color.use <- my36colors[cell_type_order]
  if (source_cell %in% rownames(cellchat_C5aRa_PD1@net$count)){
    netVisual_circle_self(
      net = cellchat_C5aRa_PD1@net$count, 
      cell_type_order = cell_type_order,
      angle_coords = -45,
      color.use = color.use,
      vertex.weight = table(cellchat_C5aRa_PD1@idents)[cell_type_order],
      vertex.label.cex = 1.8,
      vertex.size.max = 15,
      weight.scale = T, 
      label.edge= TRUE, 
      edge.label.cex = 2,
      sources.use = source_cell,
      targets.use = target_cell,
      remove.isolate = TRUE,
      # top = 0.25,
      arrow.width = 1,
      arrow.size = 1,
      title.name = 'C5aRa_PD1',
      margin = c(0,0,0.5,0)
    )
  }
}
source_cell = 'interferon-stimulated NAN'
cell_type_order <- c(source_cell,target_cell) %>% 
  .[. %in% rownames(cellchat_C5aRa_PD1@net$count)]
color.use <- my36colors[cell_type_order]
netVisual_circle_self(
  net = cellchat_C5aRa_PD1@net$count, 
  cell_type_order = cell_type_order,
  angle_coords = -45,
  color.use = color.use,
  vertex.weight = table(cellchat_C5aRa_PD1@idents)[cell_type_order],
  vertex.label.cex = 1.8,
  vertex.size.max = 15,
  weight.scale = T, 
  label.edge= TRUE, 
  edge.label.cex = 2,
  sources.use = source_cell,
  targets.use = target_cell,
  remove.isolate = TRUE,
  # top = 0.25,
  arrow.width = 1,
  arrow.size = 1,
  title.name = 'C5aRa_PD1',
  margin = c(0,0,0.5,0)
)
# Figs9 F(left) -----------------------------------------------------------------
res_plot_self <- function(scRNA_11_sub,
                          diffgene_source = 'diffgroup',
                          celltype_diffgene = NULL,
                          scRNA_all = NULL,
                          diff_gene_from = 'double',
                          upstream_ligands_n = 20,
                          top_n_target_gene = 100,
                          target_gene = NULL,
                          target_ligands = NULL,
                          num_show_gene = 15,
                          min.pct = 0.1,
                          p_ligand_target_network_x_size = 20,
                          axis_text_y = 30,
                          width_fig = c(1,4),
                          axis_title_y = '',
                          axis_title_x = ''){
  organism = "mouse"
  if(organism == "human"){
    lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
    weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
  } else if(organism == "mouse"){
    lr_network = readRDS('E:/nicheneter/lr_network_mouse_21122021.rds')
    ligand_target_matrix = readRDS('E:/nicheneter/ligand_target_matrix_nsga2r_final_mouse.rds')
    weighted_networks = readRDS('E:/nicheneter/weighted_networks_nsga2r_final_mouse.rds')
    
  }
  lr_network = lr_network %>% distinct(from, to)
  head(lr_network)
  ligand_target_matrix[1:5,1:5]  
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
  head(weighted_networks$lr_sig) 
  
  head(weighted_networks$gr) 
  
  
  seuratObj = alias_to_symbol_seurat(scRNA_11_sub, "mouse")
  Idents(seuratObj) <- seuratObj$anno_res
  expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  list_expressed_genes_sender = sender_celltypes %>% 
    unique() %>% 
    lapply(get_expressed_genes, seuratObj, 0.10)  
  expressed_genes_sender = list_expressed_genes_sender %>%
    unlist() %>% 
    unique()
  if(diffgene_source == 'diffgroup'){
    seurat_obj_receiver= subset(seuratObj, idents = receiver)
    seurat_obj_receiver = SetIdent(
      seurat_obj_receiver, 
      value = seurat_obj_receiver[["batch", drop = TRUE]])
    DE_table_receiver = FindMarkers(
      object = seurat_obj_receiver, 
      ident.1 = condition_oi, 
      ident.2 = condition_reference,
      min.cells.group = 1, 
      min.pct = min.pct) %>% 
      rownames_to_column("gene")
    
    if(diff_gene_from == 'double'){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }else if(diff_gene_from == condition_oi){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & avg_log2FC >= 0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }else if(diff_gene_from == condition_reference){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & avg_log2FC <= -0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }
  }else if(diffgene_source == 'diffcluster'){
    seurat_obj_receiver= subset(
      scRNA_all, subset = anno_res %in% celltype_diffgene
    ) %>% subset(subset = batch =='C5aRA_PD1')
    seurat_obj_receiver = SetIdent(
      seurat_obj_receiver, 
      value = seurat_obj_receiver[["anno_res", drop = TRUE]])
    
    DE_table_receiver = FindMarkers(
      object = seurat_obj_receiver, 
      ident.1 = receiver,
      min.cells.group = 1,
      min.pct = 0.10) %>% 
      rownames_to_column("gene")
    geneset_oi = DE_table_receiver %>% 
      filter(p_val <= 0.05 & avg_log2FC >= 0.25) %>%
      pull(gene)
    geneset_oi = geneset_oi %>% 
      .[. %in% rownames(ligand_target_matrix)]
  }
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% 
    filter(from %in% expressed_ligands & 
             to %in% expressed_receptors) %>% 
    pull(from) %>% unique()
  ligand_activities = predict_ligand_activities(
    geneset = geneset_oi, 
    background_expressed_genes = background_expressed_genes, 
    ligand_target_matrix = ligand_target_matrix, 
    potential_ligands = potential_ligands)
  
  ligand_activities = ligand_activities %>% 
    arrange(-aupr_corrected) %>% 
    mutate(rank = rank(plyr::desc(aupr_corrected)))
  best_upstream_ligands = ligand_activities %>%
    top_n(upstream_ligands_n, aupr_corrected) %>%
    arrange(-aupr_corrected) %>%
    pull(test_ligand) %>% unique()
  active_ligand_target_links_df = 
    best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi, 
           ligand_target_matrix = ligand_target_matrix, n = top_n_target_gene) %>% 
    bind_rows() %>% drop_na()
  
  active_ligand_target_links = 
    prepare_ligand_target_visualization(
      ligand_target_df = active_ligand_target_links_df, 
      ligand_target_matrix = ligand_target_matrix, 
      cutoff = 0.25)
  
  order_ligands = intersect(
    best_upstream_ligands, 
    colnames(active_ligand_target_links)) %>% 
    rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% 
    unique() %>% 
    intersect(rownames(active_ligand_target_links)) %>% 
    make.names()
  rownames(active_ligand_target_links) = 
    rownames(active_ligand_target_links) %>% 
    make.names()  
  colnames(active_ligand_target_links) = 
    colnames(active_ligand_target_links) %>% 
    make.names()  
  
  vis_ligand_target = 
    active_ligand_target_links[order_targets,order_ligands] %>% t()
  if(ncol(vis_ligand_target)>num_show_gene){vis_ligand_target <- vis_ligand_target[,1:num_show_gene]}
  p_ligand_target_network = 
    make_heatmap_ggplot_self(
      matrix = vis_ligand_target,
      y_name ="Prioritized ligands",
      x_name = "Predicted target genes", 
      color = "purple",
      tile_color = 'white',
      legend_position = "top", 
      x_axis_position = "top",
      legend_title = "Regulatory potential")  + 
    theme(
      panel.border = element_rect(fill = NA),
      axis.text.y = element_text(face = "italic",size = 20),
      axis.title = element_text(face = "italic",size = 22),
      axis.title.x = element_blank(),
      legend.title = element_text(face = "italic",size = 22),
      legend.text = element_text(face = "italic",size = 12)
    )
  if(!is.null(target_gene)){
    p_ligand_target_network <- 
      p_ligand_target_network+
      theme(
        axis.text.x = element_text(
          face = "italic",
          # size = 12,
          size = p_ligand_target_network_x_size,
          hjust = 0.5,
          vjust = 0.5,
          color = ifelse(colnames(vis_ligand_target) %in% target_gene, "red", "black"))
      )
  }
  p_ligand_target_network
  lr_network_top = 
    lr_network %>% 
    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% 
    distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  ## 7) Summary visualizations of the NicheNet analysis 
  # ligand activity heatmap
  ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  
  rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
  colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()
  
  vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
  p_ligand_aupr = make_heatmap_ggplot_self(
    matrix = vis_ligand_aupr,
    y_name = paste("Prioritized ligands in\n",axis_title_y),
    x_name = "Ligand activity",
    color = "darkorange",
    color_low = '#FFCF95',
    tile_color = 'white',
    legend_position = "top",
    x_axis_position = "top",
    x_axis = FALSE,
    legend_title = "AUPR\n(target gene prediction ability)  "
  ) +
    theme(
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "italic",size = 32),
      legend.title = element_text(face = "italic",size = 22),
      axis.text.y = element_text(face = "italic",size = axis_text_y),
      # legend_title = element_text(face = "italic",size = 22),
      legend.text = element_text(face = "italic",size = 18)
    )
  if(!is.null(target_ligands)){
    p_ligand_aupr <- p_ligand_aupr + 
      theme(
        axis.text.y = element_text(
          face = "italic",
          size = axis_text_y,
          color = ifelse(rownames(vis_ligand_aupr) %in% target_ligands, "red", "black"))
      )
  }
  # p_ligand_aupr
  
  # ligand expression Seurat dotplot
  order_ligands_adapted <- str_replace_all(order_ligands, "\\.", "-")
  
  figures_without_legend = cowplot::plot_grid(
    p_ligand_aupr +
      # labs(title = 'Ligand activity') +
      theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        plot.margin = margin(t = 5,r = 1,b = 5,l = 0,unit = "pt")
      ),
    # p_ligand_lfc +
    #   theme(
    #     legend.position = "none",
    #     axis.ticks = element_blank(),
    #     axis.text.y = element_blank(),
    #     plot.margin = margin(t = 5,r = 1,b = 5,l = 1,unit = "pt")
    #   ),
    p_ligand_target_network +
      theme(legend.position = "none",
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            # axis.text.x = element_text(size = p_ligand_target_network_x_size),
            plot.margin = margin(t = 5,r = 5,b = 5,l = 5,unit = "pt")
      ),
    align = "h",
    axis = 'l',
    nrow = 1,
    rel_widths = width_fig
  )
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, 1.5))
  
  combined_plot = cowplot::plot_grid(
    figures_without_legend, 
    legends, 
    rel_heights = c(10,2), 
    nrow = 2, 
    align = "hv")
  if(diff_gene_from == 'double'){
    title = paste(condition_oi,' VS ',condition_reference,' differential gene of ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25)
  }else if(diff_gene_from == condition_oi){
    title = paste('Predicted target genes in ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 45)
  }else if(diff_gene_from == condition_reference){
    title = paste('Predicted target genes in ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25)
  }
  return(final_plot)
}
sender_celltypes = c(
  "Gzmb_T cell",
  "Dysfunction_T cell",
  "Cytolytic_T cell"
)
receiver = c(
  "NPBNs","Exhausted TAN",
  "interferon-stimulated NAN"
)
condition_oi = "PD1"
condition_reference = "C5aRA_PD1"
scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% c(
  sender_celltypes,
  receiver
)) %>% 
  subset(subset = batch %in% c(condition_oi,condition_reference))
scRNA_11_sub@meta.data %>% group_by(batch,anno_res) %>% summarise(Counts = n())
target_gene <- c("B2m","Cd47","Cxcl10","Ddx58","H2.D1","H2.K1",
                 "H2.Q4","H2.Q6","H2.T24","Icam1","Ifi35",
                 "Ifit2","Ifit3","Ifit3b","Ifitm3","Ifitm6",
                 "Irf1","Irf7")
target_ligands <- c("Ifng","Tnf","Icam1","Ccl4","Ccl5")
p3 <- res_plot_self(
  scRNA_11_sub = scRNA_11_sub,
  diff_gene_from = condition_reference,
  upstream_ligands_n = 15,
  top_n_target_gene = 100,
  num_show_gene = 20,
  target_gene = target_gene,
  target_ligands = target_ligands,
  p_ligand_target_network_x_size = 30,
  axis_text_y = 30,width_fig = c(1,4.5),
  axis_title_x = 'Neutrophil cells',
  axis_title_y = 'Gzmb+/Dysfunction/Cytolytic T cellS')
p3
# Figs9 F(right) -----------------------------------------------------------------
res_plot_self_2 <- function(scRNA_12_sub,
                            diffgene_source = 'diffgroup',
                            celltype_diffgene = NULL,
                            scRNA_all = NULL,
                            diff_gene_from = 'double',
                            upstream_ligands_n = 20,
                            top_n_target_gene = 100,
                            target_gene = NULL,
                            target_ligands = NULL,
                            num_show_gene = 15,
                            min.pct = 0.1,
                            p_ligand_target_network_x_size = 20,
                            axis_text_y = 30,
                            width_fig = c(1,4),
                            axis_title_y = '',
                            axis_title_x = ''){
  organism = "mouse"
  if(organism == "human"){
    lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
    weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
  } else if(organism == "mouse"){
    lr_network = readRDS('E:/nicheneter/lr_network_mouse_21122021.rds')
    ligand_target_matrix = readRDS('E:/nicheneter/ligand_target_matrix_nsga2r_final_mouse.rds')
    weighted_networks = readRDS('E:/nicheneter/weighted_networks_nsga2r_final_mouse.rds')
    
  }
  lr_network = lr_network %>% distinct(from, to)
  head(lr_network)
  ligand_target_matrix[1:5,1:5]  
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
  head(weighted_networks$lr_sig)  
  
  head(weighted_networks$gr)
  
  
  seuratObj = alias_to_symbol_seurat(scRNA_12_sub, "mouse")
  Idents(seuratObj) <- seuratObj$anno_res_3
  expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  ## sender
  list_expressed_genes_sender = sender_celltypes %>% 
    unique() %>% 
    lapply(get_expressed_genes, seuratObj, 0.10)  
  expressed_genes_sender = list_expressed_genes_sender %>%
    unlist() %>% 
    unique()
  
  if(diffgene_source == 'diffgroup'){
    seurat_obj_receiver= subset(seuratObj, idents = receiver)
    seurat_obj_receiver = SetIdent(
      seurat_obj_receiver, 
      value = seurat_obj_receiver[["batch", drop = TRUE]])
    DE_table_receiver = FindMarkers(
      object = seurat_obj_receiver, 
      ident.1 = condition_oi, 
      ident.2 = condition_reference,
      min.cells.group = 1, 
      min.pct = min.pct) %>% 
      rownames_to_column("gene")
    
    if(diff_gene_from == 'double'){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }else if(diff_gene_from == condition_oi){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & avg_log2FC >= 0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }else if(diff_gene_from == condition_reference){
      geneset_oi = DE_table_receiver %>% 
        filter(p_val <= 0.05 & avg_log2FC <= -0.25) %>%
        pull(gene)
      geneset_oi = geneset_oi %>% 
        .[. %in% rownames(ligand_target_matrix)]
    }
  }else if(diffgene_source == 'diffcluster'){
    seurat_obj_receiver= subset(
      scRNA_all, subset = anno_res_3 %in% celltype_diffgene
    ) %>% subset(subset = batch =='C5aRA_PD1')
    seurat_obj_receiver = SetIdent(
      seurat_obj_receiver, 
      value = seurat_obj_receiver[["anno_res_3", drop = TRUE]])
    
    DE_table_receiver = FindMarkers(
      object = seurat_obj_receiver, 
      ident.1 = receiver,
      min.cells.group = 1,
      min.pct = min.pct) %>% 
      rownames_to_column("gene")
    geneset_oi = DE_table_receiver %>% 
      filter(p_val <= 0.05 & avg_log2FC >= 0.25) %>%
      pull(gene)
    geneset_oi = geneset_oi %>% 
      .[. %in% rownames(ligand_target_matrix)]
  }
  
  
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  
  potential_ligands = lr_network %>% 
    filter(from %in% expressed_ligands & 
             to %in% expressed_receptors) %>% 
    pull(from) %>% unique()
  ligand_activities = predict_ligand_activities(
    geneset = geneset_oi, 
    background_expressed_genes = background_expressed_genes, 
    ligand_target_matrix = ligand_target_matrix, 
    potential_ligands = potential_ligands)
  
  ligand_activities = ligand_activities %>% 
    arrange(-aupr_corrected) %>% 
    mutate(rank = rank(plyr::desc(aupr_corrected)))
  best_upstream_ligands = ligand_activities %>%
    top_n(upstream_ligands_n, aupr_corrected) %>%
    arrange(-aupr_corrected) %>%
    pull(test_ligand) %>% unique()
  active_ligand_target_links_df = 
    best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi, 
           ligand_target_matrix = ligand_target_matrix, n = top_n_target_gene) %>% 
    bind_rows() %>% drop_na()
  
  active_ligand_target_links = 
    prepare_ligand_target_visualization(
      ligand_target_df = active_ligand_target_links_df, 
      ligand_target_matrix = ligand_target_matrix, 
      cutoff = 0.25)
  
  order_ligands = intersect(
    best_upstream_ligands, 
    colnames(active_ligand_target_links)) %>% 
    rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% 
    unique() %>% 
    intersect(rownames(active_ligand_target_links)) %>% 
    make.names()
  rownames(active_ligand_target_links) = 
    rownames(active_ligand_target_links) %>% 
    make.names()  
  colnames(active_ligand_target_links) = 
    colnames(active_ligand_target_links) %>% 
    make.names()  
  
  vis_ligand_target = 
    active_ligand_target_links[order_targets,order_ligands] %>% t()
  if(ncol(vis_ligand_target)>num_show_gene){vis_ligand_target <- vis_ligand_target[,1:num_show_gene]}
  p_ligand_target_network = 
    make_heatmap_ggplot_self(
      matrix = vis_ligand_target,
      y_name ="Prioritized ligands",
      x_name = "Predicted target genes", 
      color = "purple",
      tile_color = 'white',
      legend_position = "top", 
      x_axis_position = "top",
      legend_title = "Regulatory potential")  + 
    theme(
      panel.border = element_rect(fill = NA),
      axis.text.y = element_text(face = "italic",size = 20),
      axis.title = element_text(face = "italic",size = 22),
      axis.title.x = element_blank(),
      legend.title = element_text(face = "italic",size = 22),
      legend.text = element_text(face = "italic",size = 12)
    )
  if(!is.null(target_gene)){
    p_ligand_target_network <- 
      p_ligand_target_network+
      theme(
        axis.text.x = element_text(
          face = "italic",
          size = p_ligand_target_network_x_size,
          hjust = 0.5,
          vjust = 0.5,
          color = ifelse(colnames(vis_ligand_target) %in% target_gene, "red", "black"))
      )
  }
  p_ligand_target_network
  lr_network_top = 
    lr_network %>% 
    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% 
    distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  ligand_aupr_matrix = ligand_activities %>% dplyr::select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  
  rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
  colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()
  
  vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
  p_ligand_aupr = make_heatmap_ggplot_self(
    matrix = vis_ligand_aupr,
    y_name = paste("Prioritized ligands in\n",axis_title_y),
    x_name = "Ligand activity",
    color = "darkorange",
    color_low = '#FFCF95',
    tile_color = 'white',
    legend_position = "top",
    x_axis_position = "top",
    x_axis = FALSE,
    legend_title = "AUPR\n(target gene prediction ability)  "
  ) +
    theme(
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "italic",size = 32),
      legend.title = element_text(face = "italic",size = 22),
      axis.text.y = element_text(face = "italic",size = axis_text_y),
      legend.text = element_text(face = "italic",size = 18)
    )
  if(!is.null(target_ligands)){
    p_ligand_aupr <- p_ligand_aupr + 
      theme(
        axis.text.y = element_text(
          face = "italic",
          size = axis_text_y,
          color = ifelse(rownames(vis_ligand_aupr) %in% target_ligands, "red", "black"))
      )
  }
  order_ligands_adapted <- str_replace_all(order_ligands, "\\.", "-")
  
  figures_without_legend = cowplot::plot_grid(
    p_ligand_aupr +
      theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        plot.margin = margin(t = 5,r = 1,b = 5,l = 0,unit = "pt")
      ),
    p_ligand_target_network +
      theme(legend.position = "none",
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(t = 5,r = 5,b = 5,l = 5,unit = "pt")
      ),
    align = "h",
    axis = 'l',
    nrow = 1,
    rel_widths = width_fig
  )
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, 1.5))
  
  combined_plot = cowplot::plot_grid(
    figures_without_legend, 
    legends, 
    rel_heights = c(10,2), 
    nrow = 2, 
    align = "hv")
  if(diff_gene_from == 'double'){
    title = paste(condition_oi,' VS ',condition_reference,' differential gene of ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25)
  }else if(diff_gene_from == condition_oi){
    title = paste('Predicted target genes in ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 45)
  }else if(diff_gene_from == condition_reference){
    title = paste('Predicted target genes in ',axis_title_x,sep = '')
    final_plot = cowplot::ggdraw() +
      cowplot::draw_plot(combined_plot,height = 0.95) +
      cowplot::draw_label(title, x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25)
  }
  return(final_plot)
}
sender_celltypes = c(
  "Cytolytic/Gzmb_T cell",
  "Dysfunction_T cell"
)
receiver = c(
  "S100A8/ISG-Neu",
  "Ccl4-Neu"
)
condition_oi = "PD1"
condition_reference = "C5aRA_PD1"
scRNA_12_sub <- subset(scRNA_12,subset = anno_res_3 %in% c(
  sender_celltypes,
  receiver
)) %>% 
  subset(subset = batch %in% c(condition_oi,condition_reference))
scRNA_12_sub@meta.data %>% group_by(batch,anno_res_3) %>% summarise(Counts = n())
target_gene <- c("B2m","Cd47","Cxcl10","Ddx58","H2.D1","H2.K1",
                 "H2.Q4","H2.Q6","H2.T24","Icam1","Ifi35",
                 "Ifit2","Ifit3","Ifit3b","Ifitm3","Ifitm6",
                 "Irf1","Irf7")
target_ligands <- c("Ifng","Tnf","Icam1","Ccl4","Ccl5")
p3 <- res_plot_self_2(
  scRNA_12_sub = scRNA_12_sub,
  diff_gene_from = condition_reference,
  upstream_ligands_n = 5,
  top_n_target_gene = 250,
  target_gene = target_gene,
  target_ligands = target_ligands,
  num_show_gene = 36,
  axis_title_y = 'CD8+ T cells',
  axis_title_x = 'Neutrophil cells',
  p_ligand_target_network_x_size = 30)
p3