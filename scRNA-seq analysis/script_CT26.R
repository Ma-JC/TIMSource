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
# Function --------------------------------------------------------------------
stacked_bar_self <- function(meta_data_all, meta_data, cell_type,col_name = 'anno_res') {
  meta_data_all <- meta_data_all %>% 
    mutate(
      batch = case_when(
        batch == 'Mock' ~ 'Mock',
        batch == 'C5aRA' ~ 'C5aRA',
        batch == 'PD1' ~ 'PD1',
        batch == 'C5aRA_PD1' ~ 'C5aRA_PD1'
      )
    )
  meta_data <- meta_data %>% 
    mutate(
      batch = case_when(
        batch == 'Mock' ~ 'Mock',
        batch == 'C5aRA' ~ 'C5aRA',
        batch == 'PD1' ~ 'PD1',
        batch == 'C5aRA_PD1' ~ 'C5aRA_PD1'
      )
    )
  data_batch <- table(meta_data_all$batch) %>%
    as.data.frame() %>%
    rename_all( ~ c('batch', 'num'))
  data_plot <- meta_data %>%
    group_by(!!sym(col_name), batch) %>%
    summarise(Counts = n(),.groups = 'drop') %>%
    mutate(
      batch_counts = case_when(
        batch == "C5aRA" ~ data_batch$num[data_batch$batch == 'C5aRA'],
        batch == "C5aRA_PD1" ~ data_batch$num[data_batch$batch == 'C5aRA_PD1'],
        batch == "Mock" ~ data_batch$num[data_batch$batch == 'Mock'],
        batch == "PD1" ~ data_batch$num[data_batch$batch == 'PD1']
      )
    ) %>%
    group_by(!!sym(col_name)) %>% 
    mutate(percentage_inner = Counts / batch_counts * 100) %>%
    mutate(percentage = percentage_inner / sum(percentage_inner) * 100) %>%
    mutate(batch = factor(batch, levels = rev(c(
      "Mock", "C5aRA", "PD1", "C5aRA_PD1"
    )))) %>%
    mutate(anno_res = factor(!!sym(col_name), levels = cell_type)) %>%
    mutate(facet_group = case_when(
      batch %in% c("Mock", "C5aRA") ~ 'C5aRA VS Mock',
      batch %in% c("PD1", "C5aRA_PD1") ~ 'C5aRA_PD1 VS PD1'
    ))
  data_add <- data_plot %>%
    group_by(anno_res, facet_group) %>%
    summarize(n_rows = n(),.groups = 'drop') %>%
    left_join(data_plot, .) %>%
    mutate(add_rows = ifelse(
      n_rows == 1,
      yes = setdiff(
        c("Mock", "C5aRA", "PD1", "C5aRA_PD1"),
        unique(batch) %>% as.character()
      ),
      no = ''
    )) %>%
    ungroup() %>%
    filter(add_rows != '') %>% 
    mutate(
      batch = add_rows,
      Counts = 0,
      batch_counts = 0,
      percentage_inner = 0,
      percentage = 0
    ) %>%
    dplyr::select(-c(add_rows, n_rows))
  if(nrow(data_add)>0){
    data_plot <-  data_add %>%
      bind_rows(data_plot) %>%
      arrange(anno_res, batch) %>%
      mutate(batch = factor(batch, levels = c("Mock", "C5aRA", "PD1", "C5aRA_PD1")))
  }else{
    data_plot <-  data_plot %>%
      arrange(anno_res, batch) %>%
      mutate(batch = factor(batch, levels = c("Mock", "C5aRA", "PD1", "C5aRA_PD1")))
  }
  p <-
    ggplot(data_plot, mapping = aes(x = anno_res, y = percentage, fill = batch)) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.7),
      color = "black",
      width = 0.7,
      linewidth = 0.25
    ) +
    scale_fill_manual(values = c('#0099B4', '#42B540', '#ED0000', '#00468B') %>% rev()) +
    facet_grid(facet_group ~ .,
               scales = 'free_y') +
    theme_bw() +
    theme(
      legend.position = 'top',
      legend.justification = 'center',
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = -45,
        hjust = 0,
        size = 16,
        color = "black",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 16,
        color = "black",
        face = "bold",
        vjust = 0.5,
        hjust = 0.5
      ),
      axis.title.y = element_text(
        size = 25,
        color = "black",
        face = "bold"
      ),
      axis.text.y.right = element_blank(),
      plot.title = element_text(hjust = 0.5, vjust = 6),
      plot.margin = margin(30, 2, 2, 2),
      axis.line = element_line(colour = 'grey30', size = 0.2),
      panel.spacing = unit(2, "mm"),
      panel.border = element_rect(
        fill = NA,
        linetype = 'solid',
        linewidth = 1
      ),
      panel.grid.major = element_blank(),
      panel.grid = element_blank(),
      strip.text.y = element_text(
        size = 16,
        face = "bold",
        color = "#FFFFFF",
        vjust = 0.5,
        margin = margin(b = 3, t = 3,r = 3,l = 5)
      ),
      strip.background = element_rect(
        colour = 'black',
        fill = '#2878B5',
        size = 1
      ),
      strip.placement = "inside"
    )
  return(p)
}
stacked_bar_self_2 <- function(meta_data_all, meta_data, cell_type,col_name = 'anno_low') {
  data_batch <- table(meta_data_all$batch) %>%
    as.data.frame() %>%
    rename_all( ~ c('batch', 'num'))
  data_plot <- meta_data %>%
    group_by(!!sym(col_name), batch) %>%
    summarise(Counts = n(),.groups = 'drop') %>%
    mutate(
      batch_counts = case_when(
        batch == "C5aRA" ~ data_batch$num[data_batch$batch == 'C5aRA'],
        batch == "C5aRA_PD1" ~ data_batch$num[data_batch$batch == 'C5aRA_PD1'],
        batch == "Mock" ~ data_batch$num[data_batch$batch == 'Mock'],
        batch == "PD1" ~ data_batch$num[data_batch$batch == 'PD1']
      )
    ) %>%
    group_by(!!sym(col_name)) %>% 
    mutate(percentage_inner = Counts / batch_counts * 100) %>%
    mutate(percentage = percentage_inner / sum(percentage_inner) * 100) %>%
    mutate(batch = factor(batch, levels = rev(c(
      "Mock", "C5aRA", "PD1", "C5aRA_PD1"
    )))) %>%
    mutate(anno_res = factor(!!sym(col_name), levels = cell_type)) %>%
    mutate(facet_group = case_when(
      batch %in% c("Mock", "C5aRA") ~ 'C5aRA VS Mock',
      batch %in% c("PD1", "C5aRA_PD1") ~ 'C5aRA_PD1 VS PD1'
    ))
  data_add <- data_plot %>%
    group_by(anno_res, facet_group) %>%
    summarize(n_rows = n(),.groups = 'drop') %>%
    left_join(data_plot, .) %>%
    mutate(add_rows = ifelse(
      n_rows == 1,
      yes = setdiff(
        c("Mock", "C5aRA", "PD1", "C5aRA_PD1"),
        unique(batch) %>% as.character()
      ),
      no = ''
    )) %>%
    ungroup() %>%
    filter(add_rows != '') %>% 
    mutate(
      batch = add_rows,
      Counts = 0,
      batch_counts = 0,
      percentage_inner = 0,
      percentage = 0
    ) %>%
    dplyr::select(-c(add_rows, n_rows))
  if(nrow(data_add)>0){
    data_plot <-  data_add %>%
      bind_rows(data_plot) %>%
      arrange(anno_res, batch) %>%
      mutate(batch = factor(batch, levels = c("Mock", "C5aRA", "PD1", "C5aRA_PD1")))
  }else{
    data_plot <-  data_plot %>%
      arrange(anno_res, batch) %>%
      mutate(batch = factor(batch, levels = c("Mock", "C5aRA", "PD1", "C5aRA_PD1")))
  }
  p <-
    ggplot(data_plot, mapping = aes(x = anno_res, y = percentage, fill = batch)) +
    geom_bar(
      stat = "identity",
      position = 'stack',
      color = "black",
      width = 0.6,
      linewidth = 0.25
    ) +
    scale_fill_manual(values = c('#0099B4', '#42B540', '#ED0000', '#00468B') %>% rev()) +
    guides(fill = guide_legend(nrow = 1))+
    theme_bw() +
    theme(
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      legend.position = 'top',
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = -45,
        hjust = 0,
        size = 16,
        color = "black",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 16,
        color = "black",
        face = "bold",
        vjust = 0.5,
        hjust = 0.5
      ),
      axis.title.y = element_text(
        size = 25,
        color = "black",
        face = "bold"
      ),
      axis.text.y.right = element_blank(),
      plot.title = element_text(hjust = 0.5, vjust = 6),
      plot.margin = margin(30, 2, 2, 2),
      axis.line = element_line(colour = 'grey30', size = 0.2,linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      panel.spacing = unit(2, "mm"),
      panel.border = element_rect(
        fill = NA,
        linetype = 'solid',
        linewidth = 0.5
      ),
      panel.grid = element_line(),
      strip.text.y = element_text(
        size = 16,
        face = "bold",
        color = "#FFFFFF",
        vjust = 0.5,
        margin = margin(b = 3, t = 3,r = 3,l = 5)
      ),
      strip.background = element_rect(
        colour = 'black',
        fill = '#2878B5',
        size = 0.5
      ),
      strip.placement = "inside"
    )
  return(p)
}
stacked_bar_self_batch <- function(meta_data_all, meta_data, cell_type,col_name = 'anno_low',mycolors = my36colors) {
  data_batch <- table(meta_data_all$batch) %>%
    as.data.frame() %>%
    rename_all( ~ c('batch', 'num'))
  data_plot <- meta_data %>%
    group_by(!!sym(col_name), batch) %>%
    summarise(Counts = n(),.groups = 'drop') %>%
    mutate(
      batch_counts = case_when(
        batch == "C5aRA" ~ data_batch$num[data_batch$batch == 'C5aRA'],
        batch == "C5aRA_PD1" ~ data_batch$num[data_batch$batch == 'C5aRA_PD1'],
        batch == "Mock" ~ data_batch$num[data_batch$batch == 'Mock'],
        batch == "PD1" ~ data_batch$num[data_batch$batch == 'PD1']
      )
    ) %>%
    group_by(batch) %>% 
    mutate(percentage_inner = Counts / batch_counts * 100) %>%
    mutate(percentage = percentage_inner / sum(percentage_inner) * 100) %>%
    mutate(batch = factor(batch, levels = c(
      "Mock", "C5aRA", "PD1", "C5aRA_PD1"
    ))) %>%
    mutate(anno_res = factor(!!sym(col_name), levels = cell_type)) %>%
    mutate(facet_group = case_when(
      batch %in% c("Mock", "C5aRA") ~ 'C5aRA VS Mock',
      batch %in% c("PD1", "C5aRA_PD1") ~ 'C5aRA_PD1 VS PD1'
    )) 
  p <-
    ggplot(data_plot, mapping = aes(x = batch, y = percentage, fill = anno_res,stratum=anno_res, alluvium=anno_res)) +
    geom_bar(
      stat = "identity",
      position = 'stack',
      color = "black",
      width = 0.6,
      linewidth = 0.25
    ) +
    geom_flow(width=0.5,alpha=0.4, knot.pos=0.5) +
    scale_fill_manual(values = mycolors) +
    guides(fill = guide_legend(ncol = 1))+
    theme_bw() +
    theme(
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      legend.position = 'right',
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = -45,
        hjust = 0,
        size = 16,
        color = "black",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 16,
        color = "black",
        face = "bold",
        vjust = 0.5,
        hjust = 0.5
      ),
      axis.title.y = element_text(
        size = 25,
        color = "black",
        face = "bold"
      ),
      axis.text.y.right = element_blank(),
      plot.title = element_text(hjust = 0.5, vjust = 6),
      plot.margin = margin(30, 2, 2, 2),
      axis.line = element_line(colour = 'grey30', size = 0.2,linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      panel.spacing = unit(2, "mm"),
      panel.border = element_rect(
        fill = NA,
        linetype = 'solid',
        linewidth = 0.5
      ),
      panel.grid = element_blank(),
      panel.grid.major = element_blank()
    )
  return(p)
}
plot_gene <-
  function(scRNA_11_sub,gene_single,
           target = c("Mac_C1qa", "Mac_Vegfa", "MonoMac_Ccr2", "Mac_Mki67/C1qa","Mono_Isg15",'Mast_Mcpt2'),
           label_size = 4) {
    gene_exp <- FetchData(scRNA_11_sub,vars = gene_single)
    dat<- data.frame(scRNA_11_sub@meta.data, 
                     scRNA_11_sub@reductions$umap@cell.embeddings,
                     'gene_exp' = gene_exp[,1],
                     'anno_low' = scRNA_11_sub$anno_low) %>% 
      mutate(
        seurat_annotation = ifelse(
          test = anno_res%in%target,
          yes = anno_res %>% str_replace('_','-'),
          no = anno_low)
      )%>% 
      mutate(
        seurat_annotation = ifelse(
          test = anno_res%in%c('tDC','cDC1_Xcr1'),
          yes = 'DC',
          no = seurat_annotation)
      ) %>% 
      mutate(
        seurat_annotation = ifelse(
          test = anno_res%in%c('B cell'),
          yes = 'B cell',
          no = seurat_annotation)
      ) 
    class_avg <- dat %>%
      group_by(seurat_annotation) %>%
      summarise(
        umap_1 = median(umap_1),
        umap_2 = median(umap_2)
      )
    p <- ggplot(dat, aes(umap_1, umap_2))  +
      geom_point(aes(colour  = gene_exp)) +
      scale_color_gradient(low = 'grey80',high = "blue") +
      ggrepel::geom_label_repel(aes(label = seurat_annotation),
                                data = class_avg,
                                label.size = 0,
                                force = 30,
                                force_pull = 30,
                                alpha = 0.7,
                                fill = 'grey90',
                                size = label_size,
                                segment.color = NA)+
      labs(title = paste('The Expression of',gene_single,sep = ''),x = 'Umap 1',y= 'Umap 2')+
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5,size = 25),
        axis.text = element_blank(),
        axis.title = element_text(size = 20),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank()
      )
    return(p)
  }
ssgsea_violon_self <- function(data,ssgsea_i){
  title <- ssgsea_i %>% str_split(':') %>% unlist() %>% .[2]
  scRNA_use_ggsea_sub <- data[data$ssgsea == ssgsea_i,]
  scRNA_use_ggsea_sub$facet  <- ifelse(
    scRNA_use_ggsea_sub$group %in% c('Mock', 'C5aRA'),
    yes = 'MOCK vs C5aRa',
    no = 'PD1 vs C5aRA_PD1'
  )
  p <- ggplot(scRNA_use_ggsea_sub, aes(x = ssgsea, y = Expression,fill = group)) + 
    geom_violin(trim=FALSE,color="white") + 
    geom_boxplot(width=0.2,position=position_dodge(0.9),show.legend = FALSE,size = 0.5)+  
    scale_fill_lancet() +
    labs(title = title,y = 'Fraction') +
    stat_compare_means(aes(group =  group),
                       label = "p.format",
                       method = "wilcox.test",
                       hide.ns = FALSE,
                       size = 7,
                       symnum.args = list(
                         cutpoints = c(0, 0.01, 0.05, 1),
                         symbols = c("***", "*", "not")
                       )
    ) +
    facet_grid(~facet, scales="free_x",space = "free") +
    theme_bw() + 
    theme(plot.title = element_text(size = 12,color="black",hjust = 0.5), 
          axis.text.x = element_text(angle = 45, hjust = 1,size = 20 ),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.text = element_text(size= 12),
          legend.title= element_text(size= 12)) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=16, color="black",face="bold"),
      axis.title.y = element_text(size=22,colour = 'black',face="bold",vjust = 0,hjust = 0.5),
      axis.text.y.right = element_blank(),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 24),
      plot.margin = margin(2,2,2,2),
      axis.line = element_line(colour = 'grey30',size = 0.2), 
      panel.spacing=unit(2, "mm"), 
      panel.border = element_rect(fill = NA,linetype = 'solid',linewidth = 1),
      panel.grid.major.y = element_line(),
      strip.text.x = element_text(size=22, face="bold",color = "#FFFFFF",
                                  vjust = 0.5,margin = margin(b = 3,t=3)),
      strip.background = element_rect(colour="grey30", fill="grey60",size = 1)
    )
  p
  return(p)
}
ssgsea_violon_self_v2 <- function(data_plot,ssgsea_i){
  data_plot <- data_plot %>% 
    filter(ssgsea == ssgsea_i)
  head(data_plot)
  p <- ggplot(data_plot, aes(x = group, y = Expression,fill = group)) +
    # geom_violin() + 
    geom_violin(trim=FALSE,color="white",scale = 'width') +
    geom_boxplot(width=0.2,position=position_dodge(0.9),show.legend = FALSE,size = 0.5)+  
    scale_x_discrete(labels = function(x) stringr::str_wrap(x,width = 30,whitespace_only = FALSE)) +
    scale_fill_manual(values = c('lightslategrey','#0099B4', '#42B540', '#ED0000', '#00468B') %>% rev()) +
    stat_compare_means(aes(group =  as.factor(group),x = as.factor(group)),
                       # label = "p.format",
                       label = paste0("p = ", after_stat("p.format")),
                       method = "wilcox.test",
                       hide.ns = FALSE,
                       # symnum.args = list(
                       #   cutpoints = c(0, 0.01, 0.05, 1),
                       #   symbols = c("***", "*", "NS")
                       # ),
                       comparisons = 
                         list(
                           c('Mock','C5aRA'),
                           c('PD1','C5aRA_PD1')
                         ),
                       step.increase = 0,
                       size = 10
    ) +
    # facet_wrap(~GO_path,nrow = 3) +
    theme_bw()+ 
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 22,color="black",face="bold",angle =45,hjust = 1,vjust = 1),
      axis.text.y = element_text(size=16, color="black",face="bold"),
      axis.title.y = element_text(size=22,colour = 'black',face="bold",vjust = 0,hjust = 0.5),
      axis.text.y.right = element_blank(),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 24),
      plot.margin = margin(2,2,2,2),
      axis.line = element_line(colour = 'grey30',size = 0.2), 
    )
  return(p)
}
adjust_to_5_or_0 <- function(x) {
  last_digit <- (x * 10) %% 10
  if (last_digit < 5) {
    x <- x %/% 1
  } else {
    x <- x %/% 1 + 0.5
  }
  return(x)
}
gsea_plot_unselect_self <- function(
    gsea_results_df_sub,
    title_plot = "MisgDB_H terms of Mac(C5aRA)",
    title_vjust = -10,
    title_hjust = -0.1,
    y_size = 16,
    hjust_x_title = 0.3,
    color_name = 'hallmark'
){
  ylim_end_1 <- nrow(gsea_results_df_sub)
  xlim_left_1 <- 0#gsea_results_df_sub$NES %>% min()
  xlim_right_1 <- gsea_results_df_sub$NES %>% max()
  color_list = list(
    'hallmark'= c('#E69D9D','#c82423'),
    'kegg' = c('#C0E4E1','#8ECFC9'),
    'others' = c("grey90","#3131F2")
  )
  p_max <- gsea_results_df_sub$pvalue %>% min() %>% -log10(.) %>% floor()
  p_min <- gsea_results_df_sub$pvalue %>% max() %>% -log10(.) %>% ceiling()
  p_set <- signif((p_max-p_min)/3,digits = 3)
  p1 <- ggplot(data = gsea_results_df_sub,
               mapping = aes(x = NES, y = Description,fill = -log10(pvalue))) +
    geom_bar(stat = "identity",color = 'black',size = 0.3,width = 0.8) +
    scale_fill_gradient(
      low = color_list[[color_name]][1],
      high = color_list[[color_name]][2],
      breaks = c(p_min,adjust_to_5_or_0(p_min+p_set),adjust_to_5_or_0(p_min+p_set*2),p_max),
      labels = c(p_min,adjust_to_5_or_0(p_min+p_set),adjust_to_5_or_0(p_min+p_set*2),p_max)) +
    guides(fill = guide_colorbar(
      title.position = 'left',
      title.theme = element_text(angle = 90,hjust = 0.5,vjust = 0.5,face = 'bold',size = 16),
      barheight = 8,label.vjust = 1
    )) +
    scale_y_discrete(
      #expand = c(0, 0),
      limits = c(levels(gsea_results_df_sub$Description),
                 '', ''),
      breaks = levels(gsea_results_df_sub$Description),
      position = 'right'
    ) +
    scale_x_continuous(breaks = c(-1.5,-0.5,0.5,1.5)) +
    geom_segment(aes(
      x = xlim_left_1,
      xend = xlim_left_1,
      y = 0,
      yend = ylim_end_1 + 0.8,
    ), size = 0.3) +
    geom_segment(aes(
      x = min((gsea_results_df_sub$NES %>% min())-0.1,0),
      xend = xlim_right_1 + 0.1,
      y = 0,
      yend = 0
    ), size = 0.5) +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(
        hjust = title_hjust,
        vjust = title_vjust,
        size = 20,
        face = 'bold'
      ),
      plot.margin = margin(t = title_vjust +2, r = 2,b =  2,l =  2),
      axis.text = element_text(
        size = 16,
        face = 'bold',
        color = 'black'
      ),
      axis.text.y = element_text(
        vjust = 0.5,hjust = 0.5,face = 'bold',size = y_size),
      axis.title.x = element_text(
        size = 16,
        face = 'bold',
        hjust = hjust_x_title),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(face = 'bold',size = 16),
      legend.text = element_text(size = 16)
    ) +
    labs(title = title_plot,
         x = 'NES')
  # p1
  return(p1)
}


my36colors <-c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
  "#e377c2","#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78",
  "#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","#c7c7c7",
  "#dbdb8d","#9edae5","#7698b3","#d6616b","#a55194","#ce6dbd",
  "#756bb1","#8c6d31","#b5cf6b","#7b4173","#cedb9c","#6b6ecf",
  "#9c9ede","#bd9e39","#d9d9d9","#ad494a","#8ca252","#e7ba52"
) %>% sample(36,replace = FALSE)
my8colors <- c("#7698b3","#9c9ede","#2ca02c","#d62728",
               "#7f7f7f","#9edae5","#e7ba52","#e377c2")
scRNA_11 <- readRDS('D:/jiaoxi/scRNA_11_res.rds')
meta_data <- scRNA_11@meta.data
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
scRNA_11@meta.data <- meta_data
scRNA_11 <- subset(scRNA_11,cells = rownames(meta_data)[!(meta_data$anno_res %>% is.na())])
scRNA_11 <- subset(scRNA_11,cells = rownames(meta_data)[!(meta_data$anno_res %in% c(
  "low_qc cell","uncertain_1","uncertain_2"
))])
scRNA_11$anno_res %>% is.na() %>% table()
scRNA_11$anno_res %>% table()
scRNA_11$anno_res[scRNA_11$anno_res == 'NPBNs'] = 
  'S100A8-Neu'
scRNA_11$anno_res[scRNA_11$anno_res == 'Exhausted TAN'] = 
  'Ccl4-Neu'
scRNA_11$anno_res[scRNA_11$anno_res == 'interferon-stimulated NAN'] = 
  'ISG-Neu'
scRNA_T <- readRDS('D:/jiaoxi/scRNA_T_cd8_再注释.rds')
scRNA_T$anno_lv3 %>% unique()
meta_data <- scRNA_T@meta.data
meta_data$anno_lv3 %>% table()
meta_data$anno_lv3 <- 
  scRNA_11$anno_res[match(rownames(meta_data),colnames(scRNA_11))]
meta_data <- 
  meta_data %>% 
  mutate(
    batch = case_when(
      batch == 'MOCK' ~ 'Mock',
      batch == 'C5aRa' ~ 'C5aRA',
      TRUE ~ batch
    )
  ) %>% 
  mutate(
    batch = factor(batch,levels = c('Mock','C5aRA','PD1','C5aRA_PD1'))
  )
meta_data$batch %>% unique()
scRNA_T@meta.data <- meta_data
scRNA_T_sub <- subset(
  scRNA_T,
  subset = anno_lv3 != 'low_qc cell'
)
scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% c(
  "Cancer cell","Fibroblast","NK",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "Dysfunction_T cell",
  "Cytolytic_T cell",
  "Naive T","cDC1_Xcr1",
  "Proliferative CD8+ T","CD8+ T memory",
  "Treg",
  # "Neutrophil",
  # "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "S100A8-Neu","Ccl4-Neu","ISG-Neu",
  "Mast_Mcpt2","B cell",
  "Mono_Isg15",
  "Proliferative Treg","CD4+ T memory",
  "tDC","Gzmb_T cell"
))

scRNA_11_sub <-
  FindVariableFeatures(object = scRNA_11_sub, nfeatures = 3000)
scRNA_11_sub <- RunPCA(scRNA_11_sub,
                       features = VariableFeatures(object = scRNA_11_sub))
scRNA_11_sub <- RunUMAP(scRNA_11_sub, dims = 1:20, label = T)
DimPlot(
  object = scRNA_11_sub,
  group.by = c('anno_res'),
  label = TRUE
) + NoLegend()
# Fig3 I left----------------------------------------------------------------
scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% c(
  "Cancer cell","Fibroblast","NK",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "Dysfunction_T cell",
  "Cytolytic_T cell",
  "Naive T","cDC1_Xcr1",
  "Proliferative CD8+ T","CD8+ T memory",
  "Treg",
  # "Neutrophil",
  # "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "S100A8-Neu","Ccl4-Neu","ISG-Neu",
  "Mast_Mcpt2","B cell",
  "Mono_Isg15",
  "Proliferative Treg","CD4+ T memory",
  "tDC","Gzmb_T cell"
))

scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% c(
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "Mast_Mcpt2","cDC1_Xcr1","tDC","B cell","NK"
))
meta_data_sub <- scRNA_11_sub@meta.data
# meta_data_sub[,c('umap_1','umap_2')] <- 
# scRNA_11_sub@reductions$umap@cell.embeddings
meta_data_sub <- meta_data_sub %>%
  mutate(
    anno_low = case_when(
      anno_res %in% c(
        "Dysfunction_T cell",
        "Cytolytic_T cell",
        "Gzmb_T cell",
        "Proliferative CD8+ T",
        "CD8+ T memory",
        "Naive T",
        "Treg",
        "Proliferative Treg",
        "CD4+ T memory"
      ) ~ 'T Cells',
      anno_res %in% c("B cell") ~ 'B Cells',
      anno_res %in% c("NK") ~ 'Natural Killer (NK) Cells',
      anno_res %in% c(
        "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
        "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
        "Mac/Mono_Mki67"
      ) ~ 'Mac/Mono Cells',
      anno_res %in% c("Mast_Mcpt2") ~ 'Mast Cells',
      # anno_res %in% c("cDC1_Xcr1", "tDC") ~ 'Dendritic Cells (DCs)',
      anno_res %in% c("cDC1_Xcr1") ~ 'cDC Cells (cDCs)',
      anno_res %in% c("tDC") ~ 'tDC Cells (tDCs)',
      anno_res %in% c("NPBNs","Exhausted TAN","interferon-stimulated NAN") ~ 'Neutrophil Cells'
    )
  ) %>% 
  mutate(anno_low = factor(anno_low,levels = c(
    "T Cells","Natural Killer (NK) Cells","B Cells",
    "Mac/Mono Cells","Neutrophil Cells",
    "cDC Cells (cDCs)","tDC Cells (tDCs)",
    "Mast Cells"
  )))
# meta_data_sub$anno_res[is.na(meta_data_sub$anno_low)]
meta_data_sub$anno_low %>% unique()
scRNA_11_sub@meta.data <- meta_data_sub

### anno_low的cell type ------------------------------------------------------


scRNA_11_sub <-
  FindVariableFeatures(object = scRNA_11_sub, nfeatures = 4000)
scRNA_11_sub <- RunPCA(
  scRNA_11_sub,
  features = VariableFeatures(object = scRNA_11_sub,nfeatures = 3000))
scRNA_11_sub <- RunUMAP(scRNA_11_sub, dims = 1:20, label = T)#30
p <- DimPlot(
  object = scRNA_11_sub,
  group.by = c('anno_low'),
  label = TRUE,
  pt.size = 0.8,
  label.size = 6,
  label.box = TRUE,
  # cols = my36colors,
  repel = 30,
  alpha = 0.8
)
p
cell_leves <- c(
  "T Cells","Natural Killer (NK) Cells","B Cells",
  "Mac/Mono Cells","Neutrophil Cells",
  "cDC Cells (cDCs)","tDC Cells (tDCs)",
  "Mast Cells"
)
label_legend <- paste(
  1:length(cell_leves),
  cell_leves,sep = ': '
)
data_plot <- p$data %>% 
  mutate(anno_low = factor(anno_low,levels = cell_leves)) %>% 
  mutate(label_inner = as.numeric(anno_low) %>% as.factor()) %>% 
  mutate(label_legend = paste(label_inner,anno_low,sep = ': '))
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
ggplot() +
  geom_point(
    data = data_plot,
    aes(x = umap_1,y = umap_2,color = label_legend),
    size = 0.05) +
  geom_label_repel(data = cell_label_loc,
                   aes(
                     x = umap_1,
                     y = umap_2,
                     label = label_inner),
                   size = 20,
                   alpha = 0.8,
                   fill = NA,
                   label.padding = unit(0, "lines"),
                   label.r = unit(0,'lines'),
                   label.size = NA,
                   segment.color = NA,
                   show.legend = FALSE)+
  geom_segment(aes(
    x = data_plot$umap_1 %>% min(), 
    y = data_plot$umap_2 %>% min() , 
    xend = (data_plot$umap_1 %>% min())+
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2, 
    yend = (data_plot$umap_2 %>% min())
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.3) +
  geom_segment(aes(
    x = data_plot$umap_1 %>% min(), 
    y = data_plot$umap_2 %>% min(), 
    xend = (data_plot$umap_1 %>% min()), 
    yend = (data_plot$umap_2 %>% min()) + 
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2 
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.3) +
  scale_color_manual(
    values = my8colors,
    # breaks = 1:length(cell_leves),
    label = label_legend
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))+
  labs(x = 'Umap 1',y= 'Umap 2')+
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin =  margin(1, 1, 1,1, "lines"),
    legend.title = element_blank(),
    # legend.position = 'none',
    panel.border = element_blank(),
    legend.text = element_text(size = 35),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_text(size = 30,hjust = 0.1,vjust = 7),
    axis.title.y = element_text(size = 30,vjust = -6,hjust = 0.08)
  )
# Fig3 I right----------------------------------------------------------------
scRNA_11$anno_res %>% unique()
scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% c(
  "Cancer cell","Fibroblast","Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  # "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "S100A8-Neu","Ccl4-Neu","ISG-Neu",
  "Mast_Mcpt2","cDC1_Xcr1","tDC","B cell","NK"
))
meta_data_sub <- scRNA_11_sub@meta.data
# meta_data_sub[,c('umap_1','umap_2')] <- 
# scRNA_11_sub@reductions$umap@cell.embeddings
meta_data_sub <- meta_data_sub %>%
  mutate(
    anno_low = case_when(
      anno_res %in% c(
        "Dysfunction_T cell",
        "Cytolytic_T cell",
        "Gzmb_T cell",
        "Proliferative CD8+ T",
        "CD8+ T memory",
        "Naive T",
        "Treg",
        "Proliferative Treg",
        "CD4+ T memory"
      ) ~ 'T Cells',
      anno_res %in% c("B cell") ~ 'B Cells',
      anno_res %in% c("NK") ~ 'Natural Killer (NK) Cells',
      anno_res %in% c(
        "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
        "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
        "Mac/Mono_Mki67"
      ) ~ 'Mac/Mono Cells',
      anno_res %in% c("Mast_Mcpt2") ~ 'Mast Cells',
      anno_res %in% c("cDC1_Xcr1", "tDC") ~ 'Dendritic Cells (DCs)',
      # anno_res %in% c("cDC1_Xcr1") ~ 'cDC Cells (cDCs)',
      # anno_res %in% c("tDC") ~ 'tDC Cells (tDCs)',
      anno_res %in% c("S100A8-Neu","Ccl4-Neu","ISG-Neu") ~ 'Neutrophil Cells',
      TRUE ~ anno_res
    )
  ) %>% 
  mutate(anno_low = factor(anno_low,levels = c(
    "Cancer cell","Fibroblast","T Cells","Natural Killer (NK) Cells","B Cells",
    "Mac/Mono Cells","Neutrophil Cells",
    "Dendritic Cells (DCs)",
    "Mast Cells"
  ))) %>% 
  mutate(
    batch = factor(batch,levels = c('Mock','C5aRA','PD1','C5aRA_PD1'))
  )
# meta_data_sub$anno_res[is.na(meta_data_sub$anno_low)]
meta_data_sub$anno_low %>% unique()
# scRNA_11_sub@meta.data <- meta_data_sub
meta_data_sub$anno_res <- meta_data_sub$anno_low
tmp <- table(meta_data_sub$anno_res,meta_data_sub$batch) %>%
  as.data.frame() %>% 
  pivot_wider(id_cols = 'Var1',names_from = 'Var2',values_from = 'Freq') %>% 
  as.data.frame() %>% 
  column_to_rownames('Var1')
celltype_select <- c(
  "T Cells","Natural Killer (NK) Cells","B Cells",
  "Mac/Mono Cells","Neutrophil Cells",
  "Dendritic Cells (DCs)",
  "Mast Cells"
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
    # group = factor(group,levels = c('Depletion','Enrichment')),
  ) %>% 
  group_by(celltype) %>% 
  mutate(roe_scale = scale(`Ro/e`,center = FALSE)[,1]) 
# ungroup() %>% 
# mutate(roe_scale = roe_scale - min(roe_scale))
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

# Fig4 T cell--------------------------------------------------------------------
## Fig4 A ----------------------------------------------------------------
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
### Umap --------------------------------------------------------------------
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
# my10colors <- c('#7b4173','#80B1D3','#98df8a','#d6616b',
#                 "#7698b3","#FDB462","#2ca02c","#d62728",
#                 "#f7b6d2","#8DD3C7")
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

## Fig4 B ------------------------------------------------------------------

scRNA_11$anno_res %>% unique()
scRNA_11_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  "Cancer cell","Fibroblast",
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory","NK",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  # "NPBNs","Exhausted TAN","interferon-stimulated NAN",
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
  # "NPBNs","Exhausted TAN","interferon-stimulated NAN",
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
  # "NPBNs","Exhausted TAN","interferon-stimulated NAN",
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

# pheatmap(mat = data_plot,cluster_cols = FALSE,cluster_rows = FALSE,scale = 'row')

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
    # group = factor(group,levels = c('Depletion','Enrichment')),
  ) %>% 
  group_by(celltype) %>% 
  mutate(roe_scale = scale(`Ro/e`,center = FALSE)[,1]) 
# ungroup() %>% 
# mutate(roe_scale = roe_scale - min(roe_scale))
data_ggplot$roe_scale <- data_ggplot$`Ro/e`
limit_min <- data_ggplot$roe_scale %>% min()
limit_max <- data_ggplot$roe_scale %>% max()

ggplot(data = data_ggplot,
       mapping = aes(x = batch,y = celltype)) +
  geom_tile(aes(fill = roe_scale)) +
  geom_text(aes(label = round(roe_scale,digits = 2))) +
  guides(
    # label = guide_legend(title.theme = element_blank(),override.aes = list(size = 6)),
    fill = guide_colorbar(title = 'Ro/e',title.vjust = 1)
  ) +
  scale_fill_gradient(low = 'grey90',high =  '#3131F2',limit = c(limit_min,limit_max)) +
  # scale_y_discrete(position = 'right') +
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
    #分面标签样式
    strip.background = element_blank()
  )

## Fig4 C ------------------------------------------------------------------
scRNA_T$anno_lv3 %>% unique()
cell_type <- scRNA_T$anno_lv3 %>% unique() %>% .[c(1,2,4,5,8)]
cell_type
scRNA_T_sub <- subset(scRNA_T,subset = anno_lv3%in%cell_type)
Idents(scRNA_T_sub) <- scRNA_T_sub$batch
scRNA_T_sub$batch %>% table()
#### C5aRA_PD1 VS PD1 ------------------------------------------------------------------
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
##### Hallmark基因集 -------------------------------------------------------------
gmtfile ='./mh.all.v2023.1.Mm.symbols.gmt'
geneset <- read.gmt(gmtfile)
# fgseaRes<- fgsea(geneset, stats = ranks, nperm = 1000)
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
  gsea_results_df_sub$Description[-c(9,10,11,12)]
gsea_results_df_sub_des_list <- 
  gsea_results_df_sub_des_list[-c(5,6,7,9)]
gsea_results_df_sub_des_list <- 
  gsea_results_df_sub_des_list[1:5]
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

ggsave(
  'D:/jiaoxi/plot_res/hallmark_cd8t_C5aRAPD1_PD1.tiff',
  width = 9,
  height = 5)
ggsave(
  'D:/jiaoxi/plot_res/结果图/CD8T/hallmark_cd8t_C5aRAPD1_PD1.pdf',
  width = 9,
  height = 5)
#### Mock vs C5aRA -----------------------------------------------------------
scRNA_T_tmp <-
  subset(scRNA_T_sub, 
         subset = batch %in% c('C5aRa', 'MOCK')
  )
marker_T_sub <- FindMarkers(
  object = scRNA_T_tmp,
  ident.1 = 'C5aRa',
  ident.2 = 'MOCK',
  logfc.threshold = 0,
  min.pct = 0.25,
  only.pos = FALSE)
data <- 
  marker_T_sub %>% 
  mutate(change = as.factor(ifelse(p_val < 0.05 & abs(avg_log2FC) > 1,
                                   ifelse(avg_log2FC > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene') %>% 
  dplyr::rename(logFC = names(.)[3],adj.P.Val = names(.)[6])
head(data)
data_diff <- data[data$change != 'No change',] %>% 
  arrange(desc(change))
# write.csv(data_diff,'D:/jiaoxi/diffgene/select_CD8_T_c5aRa.csv')
# write.csv(data,'D:/jiaoxi/diffgene/unselect_CD8_T_c5aRa.csv')

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
##### Hallmark基因集 -------------------------------------------------------------
gmtfile ='./mh.all.v2023.1.Mm.symbols.gmt'
geneset <- read.gmt(gmtfile)
# fgseaRes<- fgsea(geneset, stats = ranks, nperm = 1000)
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
gsea_results_df_sub_des_list <- 
  gsea_results_df_sub_des_list[c(1,2,3,4,5,8,15,16,17,18,19)]
gsea_results_df_sub_des_list <- 
  gsea_results_df_sub_des_list[1:6]
gsea_results_df_sub <- gsea_results_df %>% 
  filter(group == 'Credible') %>% 
  arrange(desc(NES)) %>% 
  filter(Description %in% gsea_results_df_sub_des_list) %>% 
  mutate(Description = factor(gsea_results_df_sub_des_list))
# gsea_results_df_sub_2 <- gsea_results_df_sub
# gsea_results_df_sub <- gsea_results_df_sub_2
gsea_results_df_sub$Description
# gsea_results_df_sub <- 
#   gsea_results_df_sub[c(1,2,3,4,5,8,15,16,17,18,19),]

gsea_plot_self(
  gsea_results_df_sub = gsea_results_df_sub,
  title_plot = "MisgDB_H terms of CD8+ T(C5aRA)",
  title_vjust = -10,
  title_hjust = -0.1,
  color_name = 'hallmark'
)

## Fig4 D TCR clone proportion------------------------------------------------------------------
meta_data_T <- read.csv('D:/jiaoxi/meta_data_T.csv',row.names = 1)
meta_data_T$anno_lv3[meta_data_T$anno_lv3=='proliferative T'] = 'proliferative CD8+ T'
cell_type <- meta_data_T$anno_lv3 %>% unique()
cell_type
cell_type_select <- cell_type[c(1,2,3,4,5,8)]
level_group <- c('Hyperexpanded','Large','Medium','Small','Single')
data_plot <- meta_data_T %>% 
  filter(anno_lv3 %in% cell_type_select) %>% 
  group_by(batch,group) %>% 
  summarise(Counts = n()) %>%
  filter(!is.na(group)) %>% 
  group_by(batch) %>% 
  mutate(prop = Counts/sum(Counts)) %>% 
  mutate(batch_two = case_when(
    batch %in% c("MOCK","C5aRa") ~ 'C5aRA vs Mock',
    batch %in% c("PD1","C5aRA_PD1") ~ 'C5aRA_PD1 vs PD1'
  )) %>% 
  mutate(batch = factor(batch,levels = c("MOCK","C5aRa","PD1","C5aRA_PD1")))
data_plot_sub <- data_plot
ggplot(data = data_plot_sub,aes(x = batch,y = prop,fill = group)) +
  geom_bar(stat = 'identity',width = 0.5) +
  geom_alluvium(aes(x = batch,stratum = group,alluvium = group),
                alpha = 0.6,
                width = 0.3,knot.pos = 0.3,
                curve_type = "xspline") +
  scale_fill_jama() +
  facet_wrap( ~ batch_two,scales = 'free_x',ncol = 1,strip.position = 'top') +
  theme_bw() +
  # labs(title = 'CD8+ T/Naive T',x = 'Sample',y = '') +#'Percentage'
  guides(fill = guide_legend(
    ncol = 2
  )) +
  theme(
    legend.position = 'bottom',
    legend.justification = c(0,1),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    axis.text.x = element_text(size=22, angle=0, hjust=0.5, color="black",face="bold"),
    axis.text.y = element_blank(),
    axis.title = element_text(size=25,colour = 'black',vjust = -0.8,hjust = 0.5),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.2), 
    strip.text = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
  )


## Fig4 G T cell monocle3------------------------------------------------------------------
cds <- readRDS('monocle2_T.rds')
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
  color_cells_by = "anno_lv3",
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
## Fig4 E top10 TCR------------------------------------------------------------------
meta_cds <- cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] %>% 
  as.data.frame() %>% 
  rename_all(~c('Umap_1','Umap_2')) %>% 
  mutate(psuedotime = cds@principal_graph_aux@listData$UMAP$pseudotime) %>% 
  mutate(t_clonotype_id = cds@colData$t_clonotype_id) %>% 
  mutate(t_cdr3s_aa = cds@colData$t_cdr3s_aa) %>% 
  mutate(batch = cds@colData$batch) %>% 
  mutate(anno_lv3 = cds@colData$anno_lv3) %>% 
  mutate(t_clonotype_label = cds@colData$t_clonotype_label)
colnames(meta_cds)
top10_cdr <- meta_cds %>% 
  group_by(t_cdr3s_aa) %>% 
  summarise(
    Counts = n(),
    name = paste(unique(t_clonotype_label),collapse = ';'),
    batch_n = unique(t_clonotype_label) %>% length()
  ) %>% 
  arrange(desc(Counts)) %>% 
  filter(complete.cases(t_cdr3s_aa)) %>% 
  dplyr::slice(1:10) %>% 
  pull(t_cdr3s_aa)

meta_cds$tmp[meta_cds$t_cdr3s_aa %in% top10_cdr] <- 
  meta_cds$t_cdr3s_aa[meta_cds$t_cdr3s_aa %in% top10_cdr] 
meta_cds$tmp <- factor(meta_cds$tmp,levels = top10_cdr)


ggplot(data = meta_cds,aes(x=Umap_1,y = Umap_2,shape = anno_lv3)) +
  geom_point(fill ='grey30',show.legend = FALSE)+
  geom_point(aes(fill = tmp,color = tmp),na.rm = TRUE,size= 3)+
  scale_color_manual(values = c('#2878B5','#9AC9DB','#F8AC8C','#C82423','#FF8884')) +  
  scale_shape_manual(values = c(21:25),name = 'Cell type') +
  scale_color_discrete(
    name = 'CT26 | Top 10 clonal expanded TCRS',
    na.translate = FALSE
  )+
  guides(
    color = guide_legend(
      override.aes = list(size = 5),
      order = 1,
      ncol = 2,title.hjust = 0,title.vjust = 0,
      title.position = 'top'),
    fill = 'none',
    shape = guide_legend(
      override.aes = list(fill = 'grey30',size = 5),
      order = 0,
      nrow = 1,title.hjust = 0,title.vjust = 0,
      title.position = 'top')
  ) +
  labs(x = 'Umap 1',y= 'Umap 2')+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5,size = 25),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(face = 'bold',size = 22),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 20,face = 'bold'),
    legend.position = 'top',
    legend.box = "vertical",legend.box.just = 'left',
    legend.justification = 'left'
  )

## Fig4 F TCR colone score------------------------------------------------------------------
file_path <- 'D:/jiaoxi/TCR/immunarch_PD1'
immdata_10x <- repLoad(file_path)
tmp <- immdata_10x$data$C5aRa
meta_data_T <- read.csv('D:/jiaoxi/meta_data_T.csv',row.names = 1)
meta_data_T$anno_lv3[meta_data_T$anno_lv3=='proliferative T'] = 'proliferative CD8+ T'
colnames(meta_data_T)
data_tcr <- meta_data_T %>%
  filter(!is.na(t_cdr3s_aa)) %>%
  group_by(batch, anno_lv3, t_cdr3s_aa) %>%
  summarise(Clones = n()) %>%
  group_by(batch, anno_lv3) %>%
  mutate(p = Clones / sum(Clones),
         N = sum(Clones)) %>%
  summarise(clonality = 1 + sum(p * log(p)) / N[1]) %>%
  ungroup() %>%
  mutate(
    batch = factor(batch, levels = c('MOCK', 'C5aRa', 'PD1', 'C5aRA_PD1')),
    group = paste(batch,anno_lv3,sep = '_')
  ) %>%
  filter(group %in% tmp$group)

### C5aRa  vs MOCK ====
data_tcr$batch %>% unique()
data_tcr$anno_lv3 %>% unique()
data_tcr_sub <- data_tcr %>% 
  filter(batch %in% c('MOCK','C5aRa')) %>% 
  group_by(anno_lv3) %>% 
  filter(n() == 2) %>% 
  filter(!(anno_lv3 %in% c('Dysfunction_T cell')))
p1 <- ggplot(data_tcr_sub, aes(x = batch, y = clonality)) +
  geom_boxplot(
    aes(fill = batch), 
    show.legend = F,
    width = 0.6) +  
  scale_fill_manual(values = c('#00468B', '#ED0000')) + 
  geom_point(size = 3,color='#374E54') + 
  geom_point(size = 3,shape=21) +  
  geom_line(
    aes(group = anno_lv3), 
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
    method = 't.test',
    paired = TRUE,
    label.x.npc = 0.5,label = 'p.format'
  ) + x 
labs(title = 'C5aRa  vs Mock')
p1
### C5aRa_PD1  vs PD1 ====
data_tcr$batch %>% unique()
data_tcr$anno_lv3 %>% unique()
data_tcr_sub <- data_tcr %>% 
  filter(batch %in% c('PD1','C5aRA_PD1')) %>% 
  group_by(anno_lv3) %>% 
  filter(n() == 2) %>% 
  filter(!(anno_lv3 %in% c('Dysfunction_T cell')))
p2 <- ggplot(data_tcr_sub, aes(x = batch, y = clonality)) +
  geom_boxplot(
    aes(fill = batch), 
    show.legend = F,
    width = 0.6) + 
  scale_fill_manual(values =  c('#42B540','#0099B4')) +  
  geom_point(size = 3,color='#374E54') + 
  geom_point(size = 3,shape=21) +  
  geom_line(
    aes(group = anno_lv3), 
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
    method = "t.test",
    paired = TRUE,
    label.x.npc = 0.5
  ) + 
  labs(title = 'C5aRa_PD1  vs PD1')
p2
p <- p1+p2+theme(axis.title.y = element_blank())
p

## fig 4G T cell monocle2------------------------------------------------------------------
scRNA_T <- readRDS('D:/jiaoxi/scRNA_T_cd8_再注释.rds')
scRNA_sub <- subset(scRNA_T,
                    subset = anno_lv2 %in% 
                      c('CD8+ T',
                        'Gzmb/Havcr2 T',
                        'CD8+ T memory',
                        'proliferative T'
                      )#'naive T',
)
data <- GetAssayData(scRNA_sub, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNA_sub@meta.data
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
  umap.n_neighbors = 10)

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
  color_cells_by = "anno_lv3",
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

# Fig5 --------------------------------------------------------------------
## Fig5 A TAM umap ---------------------------------------------------------
scRNA_11_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "Mast_Mcpt2","cDC1_Xcr1","tDC"
))
cell_label_correct <- c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  'Neutrophil Cells','Neutrophil Cells','Neutrophil Cells',
  "Mast_Mcpt2","cDC1_Xcr1","tDC"
)
names(cell_label_correct) <- c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "Mast_Mcpt2","cDC1_Xcr1","tDC"
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

## Fig5 B ssGSEA -----------------------------------------------------------
data_sig <- read.csv('D:/jiaoxi/ssgsea/PMID_37832554.csv',header =  FALSE)
geneSet <- lapply(1:nrow(data_sig), function(x){
  gene_list <- data_sig[x,2:ncol(data_sig)] %>% as.character() %>% .[nchar(.)>0]
  return(gene_list)
})
names(geneSet) <- data_sig[,1]
human <- readRDS('D:/jiaoxi/ssgsea/gene_tran/human.rds')
mouse <- readRDS('D:/jiaoxi/ssgsea/gene_tran/mouse.rds')
for(signature_i in names(geneSet)){
  gene_list <- geneSet[[signature_i]] %>% toupper()
  geneMm <-
    getLDS(
      attributes = "hgnc_symbol",
      filters = "hgnc_symbol",
      values = gene_list,
      mart = human,
      attributesL = "mgi_symbol",
      martL = mouse,
      uniqueRows = TRUE
    )
  geneSet[[signature_i]] <- geneMm$MGI.symbol
}
geneSet_sub1 <- geneSet
names(geneSet_sub1)
geneSet_sub2 <- readRDS('D:/jiaoxi/ssgsea/gene_tran/geneSet.rds')
names(geneSet_sub2)
name_select <- c(
  "CancerCell:M1 signature",                     
  "CancerCell:Angiogenesis",
  "CancerDiscovery:Complement"
)
geneSet_sub3 <- list(
  'M1_PMID33545035' = c("IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6","CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7"),
  'M2_PMID33545035' = c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4"),
  'Angiogenesis_PMID33545035' = c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA"),
  'Phagocytosis_PMID33545035' = c("MRC1","CD163","MERTK","C1QB")
)
human <- readRDS('D:/jiaoxi/ssgsea/gene_tran/human.rds')
mouse <- readRDS('D:/jiaoxi/ssgsea/gene_tran/mouse.rds')
for(signature_i in names(geneSet_sub3)){
  gene_list <- geneSet_sub3[[signature_i]] %>% toupper()
  geneMm <-
    getLDS(
      attributes = "hgnc_symbol",
      filters = "hgnc_symbol",
      values = gene_list,
      mart = human,
      attributesL = "mgi_symbol",
      martL = mouse,
      uniqueRows = TRUE
    )
  geneSet_sub3[[signature_i]] <- geneMm$MGI.symbol
}
gene_list <- msigdbr(species = "Mus musculus",category = 'H')
gene_list$gs_cat %>% unique()
term_name <- c(
  "HALLMARK_GLYCOLYSIS","HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)
gene_list_sub <- gene_list %>% 
  filter( gs_name %in% term_name)%>%
  dplyr::select(gs_name,gene_symbol) %>% 
  rename_all(~c('term','gene'))
gene_list_sub$term %>% unique() %>% str_replace_all('HALLMARK_','') %>% str_to_lower()
geneSet_sub4 <- split(gene_list_sub$gene,f = gene_list_sub$term)
names(geneSet_sub4) <- names(geneSet_sub4)%>% str_replace_all('HALLMARK_','') %>% str_to_lower()
geneSet <- c(geneSet_sub1,geneSet_sub2[name_select],geneSet_sub3,geneSet_sub4)
names(geneSet)
geneSet_use <- c(
  # "M2 Signature", 
  "ANTIGEN PROCESSING AND PRESENTATION",
  "RESPONSE TO INTERFERON GAMMA", 
  "M1_PMID33545035",
  "M2_PMID33545035",
  "Angiogenesis_PMID33545035",
  "Phagocytosis_PMID33545035",
  "glycolysis",
  "oxidative_phosphorylation"
)
scRNA_11$anno_res %>% unique()
scRNA_11$anno_low %>% unique()
scRNA_use <- subset(
  scRNA_11,
  subset = anno_res %in% c(
    "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
    "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
    "Mac/Mono_Mki67"
  ) 
)
scRNA_use$anno_res %>% unique()
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
  mutate(group = scRNA_use$batch,anno_res = scRNA_use$anno_res) %>% 
  rownames_to_column("sample")
names(a)
scRNA_use_ggsea <- gather(a,key = ssgsea, value = Expression, -c(group,sample,anno_res)) %>% 
  mutate(
    group = factor(group,levels = c('Mock','C5aRA','PD1','C5aRA_PD1'))
  ) %>% 
  mutate(
    ssgsea = factor(ssgsea)
  )
scRNA_use_ggsea$ssgsea %>% unique()
scRNA_use_ggsea %>% colnames()
library(pheatmap)
data_plot <- scRNA_use_ggsea %>% 
  filter(ssgsea %in% geneSet_use) %>% 
  group_by(batch,ssgsea) %>% 
  summarise(Score = mean(Expression)) %>% 
  pivot_wider(id_cols = 'anno_res',names_from = 'ssgsea',values_from = 'Score') %>% 
  column_to_rownames('anno_res') %>% t()
data_plot <- data_plot[,c(
  "Mono_Isg","Mono_Ly6i","Mono_Arhgap26",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67"
)]
colnames(data_plot)
data_plot_tmp <- data_plot[c(
  "RESPONSE TO INTERFERON GAMMA", 
  "M1_PMID33545035",
  "ANTIGEN PROCESSING AND PRESENTATION",
  "Angiogenesis_PMID33545035",
  "glycolysis",
  "M2_PMID33545035",
  "Phagocytosis_PMID33545035",
  "oxidative_phosphorylation"
) %>% rev(),]
rownames(data_plot_tmp) <- c(
  "RESPONSE TO INTERFERON GAMMA", 
  "M1 Signature",
  "ANTIGEN PROCESSING AND PRESENTATION",
  "Angiogenesis",
  "glycolysis",
  "M2 Signature",
  "Phagocytosis",
  "oxidative phosphorylation"
) %>% rev() %>% str_to_title()
pdf('D:/jiaoxi_res/20231206/ssgsea_TAM.pdf',width = 6,height = 5)
ComplexHeatmap::pheatmap(
  mat = data_plot_tmp,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  scale = "row",
  cluster_rows = F,
  treeheight_row = 0,
  cellwidth = 25,
  cellheight = 25,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  fontsize_col = 13,
  fontsize_row = 12,
  legend = TRUE,
  angle_col = '315',
  heatmap_legend_param = list(title = 'ssGSEA Score',
                              title_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5),
                              title_position = "lefttop-rot",
                              legend_direction = 'vertical',
                              legend_position = "topright",
                              legend_height = unit(30,units = "mm"),
                              labels_gp = grid::gpar(fontsize = 12, fontface = "bold",hjust = 0.5,vjust = 0.5))
)
dev.off()

## Fig5 C cell proportion --------------------------------------------------

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
meta_data_sub <- meta_data_sub %>%
  mutate(
    anno_low = case_when(
      anno_res %in% c("S100A8-Neu","Ccl4-Neu","ISG-Neu") ~ 'Neutrophil Cells',
      TRUE ~ anno_res
    )
  ) %>% 
  mutate(anno_low = factor(anno_low,levels = c(
    "Cancer cell","Fibroblast",
    "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
    "Proliferative CD8+ T","CD8+ T memory",
    "Naive T","Treg","Proliferative Treg","CD4+ T memory","NK",
    "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
    "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
    "Mac/Mono_Mki67",
    'Neutrophil Cells',
    "Mast_Mcpt2","B cell","cDC1_Xcr1","tDC"
  ))) %>% 
  mutate(
    batch = factor(batch,levels = c('Mock','C5aRA','PD1','C5aRA_PD1'))
  )
meta_data_sub$anno_low %>% unique()
meta_data_sub$anno_res  <- meta_data_sub$anno_low
celltype_select <- c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  'Neutrophil Cells',
  "Mast_Mcpt2","cDC1_Xcr1","tDC"
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

# Fig5 D M1/M2 signature --------------------------------------------------

scRNA_11_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67"
))
cell_label_correct <- c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67"
)
names(cell_label_correct) <- c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67"
)
scRNA_11_tmp$anno_tmp <- cell_label_correct[scRNA_11_tmp$anno_res] %>% as.character()
scRNA_11_tmp$anno_tmp %>% unique()
scRNA_11_tmp$batch <- factor(scRNA_11_tmp$batch,levels = c(
  "Mock", "C5aRA", "PD1", "C5aRA_PD1"
))
gene_select <-list(
  'M2 gene' = c("Tnfsf8","Fn1","Ctsa","Tgfb1","Spp1"),
  'M1 gene' = c("Tnf","Nos2","Il15","Ccl4","Fcgr4","Irf1","Cd40","Cd86","Cxcl10","Cxcl9"),
  'antigen presentation'= c('H2-Ab1', 'H2-Aa', 'H2-Eb1')
)
DotPlot(object = scRNA_11_tmp,
        features = gene_select %>% unlist() %>% unique(),
        group.by = 'batch') +
  coord_flip() +
  scale_color_gradient(
    name = 'Mean Exp',
    low = 'grey90',
    high = '#3131F2',
    breaks = c(-1,0,1),
    labels = c(-1,0,1)) +
  scale_size_continuous(name = 'Fraction') +
  guides(
    color = guide_colourbar(nrow = 1,title.position = 'top'),
    size = guide_legend(ncol = 1,title.position = 'top')
  ) +
  labs(x = '', y= 'Features') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0),
    legend.direction = 'vertical',
    axis.title.x = element_blank()
  )

# Fig5 D GSEA -------------------------------------------------------------

scRNA_11_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67"
))
cell_label_correct <- c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67"
)
names(cell_label_correct) <- c(
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67"
)
scRNA_11_tmp$anno_tmp <- cell_label_correct[scRNA_11_tmp$anno_res] %>% as.character()
scRNA_11_tmp$anno_tmp %>% unique()
scRNA_myeloid_tmp <-
  subset(scRNA_11_tmp, 
         subset = batch %in% c('Mock', 'C5aRA')
  )

marker_myeloid_sub <- FindMarkers(
  object = scRNA_myeloid_tmp,
  ident.1 = 'C5aRA',
  ident.2 = 'Mock',
  group.by = 'batch',
  logfc.threshold = 0,
  min.pct = 0.25,
  only.pos = FALSE)
cluster0.genes<- marker_myeloid_sub %>%
  arrange(desc(avg_log2FC)) %>%
  rownames_to_column() %>% 
  dplyr::rename(feature = names(.)[1]) %>% 
  dplyr::select(feature,avg_log2FC) %>% 
  dplyr::rename(logFC = avg_log2FC)
ranks<- cluster0.genes$logFC
names(ranks) <- cluster0.genes$feature
head(ranks)
geneList <- ranks
gene_list <- msigdbr(species = "Mus musculus")
gene_list$gs_subcat %>% unique()
gene_list_sub <- gene_list %>% 
  filter((gs_subcat == 'CP:KEGG') | (gs_cat == 'H'))%>% 
  dplyr::select(gs_name,gene_symbol) %>% 
  rename_all(~c('term','gene'))
egmt <- GSEA(geneList = geneList, 
             TERM2GENE=gene_list_sub, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
gsea_results_df <- egmt@result %>% 
  arrange(desc(NES)) %>%
  mutate(group = ifelse(test = pvalue<0.05,yes = 'Credible','NS')) %>% 
  mutate(Description = ifelse(
    test = grepl('KEGG_',Description),
    yes = str_replace(Description,'KEGG_','') %>% 
      str_split("_") %>% 
      map_chr(~ str_to_title(.x) %>% paste(collapse = " ")),
    no = str_replace(Description,'HALLMARK_','') %>% 
      str_split("_") %>% 
      map_chr(~ str_to_title(.x) %>% paste(collapse = " "))
  )) %>% 
  mutate(Description = factor(Description,levels = Description %>% unique() %>% rev()))
gsea_results_df$group %>% table()
gsea_results_df_sub <- gsea_results_df %>%
  arrange(desc(NES)) %>% 
  filter(group == 'Credible')

gsea_results_df_sub$Description
gsea_results_df_sub_des_list <- c(
  "Tnfa Signaling Via Nfkb",
  "Interferon Gamma Response",
  "Inflammatory Response",
  "Kras Signaling Up",
  "Interferon Alpha Response",
  "Oxidative Phosphorylation",
  "Myc Targets V1",
  "E2f Targets",
  "Cytosolic Dna Sensing Pathway",
  "Toll Like Receptor Signaling Pathway",
  "Antigen Processing And Presentation",
  "Natural Killer Cell Mediated Cytotoxicity"
)
gsea_results_df_sub <- gsea_results_df %>% 
  arrange(desc(NES)) %>% 
  filter(Description %in% gsea_results_df_sub_des_list) %>% 
  filter(!(Description == 'Oxidative Phosphorylation' & grepl('KEGG_',ID))) %>% 
  mutate(Description = factor(Description))

gsea_plot_unselect_self(
  gsea_results_df_sub = gsea_results_df_sub,
  title_plot = "MisgDB_H&KEGG terms of Mac(C5aRA)",
  title_vjust = -8,
  title_hjust = -0.1,
  y_size = 20,
  hjust_x_title = 0.5,
  color_name = 'hallmark'
) +
  theme(legend.position = c(2.8,0.9))

# Fig5 F cellchat------------------------------------------------------------------
library(CellChat)
library(patchwork)
library(pbmc3k.SeuratData)
library(Seurat)
options(stringsAsFactors = FALSE)

scRNA_11 <- readRDS('D:/jiaoxi/scRNA_11_res.rds')
cell_type_T <- c("Dysfunction_T cell","Cytolytic_T cell",
                 "CD8+ T memory","Gzmb_T cell",
                 "Proliferative CD8+ T")
cell_type_TAM <- c("Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
                   "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
                   "Mac/Mono_Mki67","cDC1_Xcr1","tDC")
cell_type_others <- c(
  "Cancer cell","Fibroblast",
  "Naive T",
  "Treg","NK",
  "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "Mast_Mcpt2","B cell",
  "Proliferative Treg","CD4+ T memory"
)
celltype <- c(cell_type_TAM,cell_type_T,cell_type_others)
celltype
scRNA_use <- subset(scRNA_11,subset = anno_res %in% celltype)
scRNA_use$res_use <- scRNA_use$anno_res
scRNA_use$res_use %>% unique()
scRNA_use$batch %>% unique()
scRNA_use$res_use <- factor(scRNA_use$res_use,levels = celltype)

group.cellType <- c(
  'Mono','Mono','Mono',
  'Mac/Mono','Mac/Mono','Mac/Mono','Mac/Mono',
  'DC','DC',"Dysfunction_T cell",'Cytolytic_T cell',
  'CD8+ T memory',"Gzmb_T cell",'Proliferative CD8+ T',
  "Cancer cell","Fibroblast",
  "Naive T",
  "Treg","NK",
  "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "Mast_Mcpt2","B cell",
  "Proliferative Treg","CD4+ T memory"
) 
names(group.cellType) <- unique(scRNA_use$res_use) %>% levels()
scRNA_use$res_low <- group.cellType[scRNA_use$res_use] %>% as.character()

scRNA_mock <- subset(scRNA_use,subset = batch ==  'MOCK')
Idents(scRNA_mock) <- scRNA_mock$res_low
scRNA_C5aRa <- subset(scRNA_use,subset = batch ==  'C5aRa')
Idents(scRNA_C5aRa) <- scRNA_C5aRa$res_low
scRNA_C5aRA_PD1 <- subset(scRNA_use,subset = batch ==  'C5aRA_PD1')
Idents(scRNA_C5aRA_PD1) <- scRNA_C5aRA_PD1$res_low
scRNA_PD1 <- subset(scRNA_use,subset = batch ==  'PD1')
Idents(scRNA_PD1) <- scRNA_PD1$res_low
cell_type_myeloid <- c('Mac/Mono','Mono')


netVisual_diffInteraction_self <- function(object, comparison = c(1, 2), measure = c("count", 
                                                                                     "weight", "count.merged", "weight.merged"), color.use = NULL, 
                                           color.edge = c("#b2182b", "#2166ac"), title.name = NULL, 
                                           sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
                                           top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
                                           vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black", 
                                           edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
                                           label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
                                           edge.curved = 0.2, shape = "circle", layout = in_circle(), 
                                           margin = 0.2, arrow.width = 1, arrow.size = 0.2) 
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
  }
  else if (measure %in% c("weight", "weight.merged")) {
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
    cells.level <- rownames(net.diff)
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
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
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
  }
  else {
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
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 2)
  }
  gg <- recordPlot()
  return(gg)
}

### sRNA_mock  ------------------------------------------------------------
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
cellchat <- computeCommunProbPathway(cellchat,thresh = 1)
cellchat <- aggregateNet(cellchat,thresh = 1) 
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(
  cellchat@net$count, 
  vertex.weight = groupSize,
  weight.scale = T, 
  label.edge= F, 
  title.name = "Number of interactions")

netVisual_circle(
  cellchat@net$weight, 
  vertex.weight = groupSize, 
  weight.scale = T, 
  label.edge= F, 
  title.name = "Interaction weights/strength")
netVisual_circle(
  cellchat@net$count, 
  sources.use = cell_type_myeloid,
  targets.use = cell_type_T,
  vertex.weight = groupSize,
  weight.scale = T, 
  label.edge= F, 
  title.name = "Number of interactions")
netVisual_circle(
  cellchat@net$weight, 
  sources.use = cell_type_myeloid,
  targets.use = cell_type_T,
  vertex.weight = groupSize, 
  weight.scale = T, 
  label.edge= F, 
  title.name = "Interaction weights/strength")
cellchat_mock <- cellchat
df.net_mock <- subsetCommunication(cellchat_mock,thresh = 0.05)
# write.csv(df.net_mock,'D:/jiaoxi/CT26_net_mock.csv')
# df.net_mock_tmp <- subsetCommunication(cellchat_mock,thresh = 1)

### sRNA_C5aRa ------------------------------------------------------------
data.input <- GetAssayData(scRNA_C5aRa, assay = "RNA", slot = "data")  
scRNA_C5aRa$res_use %>% table()
labels <- Idents(scRNA_C5aRa)
identity <- data.frame(group = labels, row.names = names(labels)) 
cellchat <- createCellChat(object = data.input, meta = scRNA_C5aRa@meta.data,group.by = "res_low")
cellchat <- setIdent(cellchat, ident.use = "res_low")  
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))  
CellChatDB <- CellChatDB.mouse 
colnames(CellChatDB$interaction) # 查看一下数据库信息
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
cellchat <- aggregateNet(cellchat,thresh = 1) 
groupSize <- as.numeric(table(cellchat@idents))
cellchat_C5aRa <- cellchat
df.net_C5aRa <- subsetCommunication(cellchat_C5aRa,thresh = 1)

### scRNA_C5aRA_PD1 ------------------------------------------------------------
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
colnames(CellChatDB$interaction) 息
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
cellchat <- aggregateNet(cellchat,thresh = 1) 
cellchat_C5aRa_PD1 <- cellchat
df.net_C5aRa_PD1 <- subsetCommunication(cellchat_C5aRa_PD1,thresh = 1)
## scRNA_PD1  ------------------------------------------------------------
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
cellchat <- computeCommunProbPathway(cellchat,thresh = 1)
cellchat <- aggregateNet(cellchat,thresh = 1)
groupSize <- as.numeric(table(cellchat@idents))
cellchat_PD1 <- cellchat
df.net_PD1 <- subsetCommunication(cellchat_PD1,thresh = 1)

## plot --------------------------------------------------------------------


pairLR.use.C5aRa_PD1 <-
  data.frame(
    interaction_name = c(
      'H2-Q6_CD8A',
      'H2-Q1_CD8A','CD86_CD28',
      'H2-Q7_CD8A','GM7030_CD8A') %>% rev()
  )
cell_type_myeloid_sub1 <- c(
  "Mac","Mono",#"MonoMac",
  'DC',"Neutrophil"#,"Mast"
)
library(ggsci)
p <- netVisual_bubble_self(
  cellchat,
  pairLR.use = pairLR.use.C5aRa_PD1,
  sources.use = cell_type_myeloid_sub1,#[c(1)],
  targets.use = cell_type_T[c(4,1)],#cell_type_T[c(1,2,4)],
  comparison = c(1, 2,3,4),
  color.text = pal_npg()(4),
  angle.x = 90,
  color.heatmap = 'viridis',
  remove.isolate = FALSE,
  font.size = 20,
  font.size.title = 24,
  line.on = TRUE,
  line.size = 1,
  grid.on = TRUE,
  sort.by.source = TRUE,
  sort.by.target = TRUE,
  title.name = paste(
    "Up-regulated signaling",
    sep = '')
)
p

data_plot <- p$data %>% 
  mutate(
    target = target %>% 
      str_remove(' cell') %>% 
      str_replace('_',' ')
  ) %>% 
  mutate(
    interaction_name_2 = interaction_name_2 %>%
      str_replace_all(' - ','->') %>% 
      str_replace_all(' ','') %>% 
      factor(levels = c(
        "Gm7030->Cd8a","H2-q7->Cd8a","Cd86->Cd28",
        "H2-q1->Cd8a","H2-q6->Cd8a"
      ))
  ) %>% 
  mutate(
    dataset = case_when(
      dataset == 'mock' ~ 'Mock',
      dataset == 'C5aRa' ~ 'C5aRA',
      dataset == 'PD1' ~ 'PD1',
      dataset == 'C5aRa_PD1' ~ 'C5aRA_PD1'
    )
  )
colnames(data_plot)
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
n.colors = 10;color.heatmap='viridis'
color.use <- tryCatch({
  RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
}, error = function(e) {
  (scales::viridis_pal(option = color.heatmap, direction = 1))(n.colors)
})
# scales::viridis_pal(option = color.heatmap, direction = -1)(n.colors)
data_plot$dataset %>% unique()
data_plot$source %>% unique()
data_plot$target %>% unique()
# data_plot$dataset %in% c('mock','C5aRa')

#### Mock VS C5aRA -----------------------------------------------------------

data_plot_sub1 <- data_plot %>% 
  filter(
    dataset %in% c('Mock','C5aRA'),
    source %in% c('Mac','Mono','DC')
  ) %>% 
  mutate(
    dataset = dataset %>% factor(levels = c('Mock','C5aRA'))
  )
p1 <- ggplot(data = data_plot_sub1,
             aes(
               x = dataset,
               y = interaction_name_2,
               color = prob,
               size = pval
             )) +
  geom_point() +
  facet_grid(~ source+target) +
  scale_radius(
    range = c(4, 6),
    breaks = sort(unique(data_plot$pval)),
    labels = names(values)[values %in% sort(unique(data_plot$pval))],
    name = "p-value"
  ) +
  scale_colour_gradientn(
    colors = colorRampPalette(color.use)(99),
    na.value = "white",
    limits = c(quantile(data_plot$prob,
                        0, na.rm = T), quantile(data_plot$prob, 1, na.rm = T)),
    breaks = c(quantile(data_plot$prob, 0, na.rm = T), quantile(data_plot$prob,
                                                                1, na.rm = T)),
    labels = c("min", "max")
  ) +
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
  theme_classic() +
  theme(
    legend.text = element_text(size= 12),
    legend.title= element_text(size= 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = -45,hjust = 0,size=16, color="black",face="bold"),
    axis.text.y = element_text(
      size = 16,
      color = "black",
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.title.y = element_blank(),
    axis.text.y.right = element_blank(),
    legend.justification = c(0, 1),
    plot.title = element_text(hjust = 0.5,vjust = 6),
    plot.margin = margin(30,2,2,2),
    axis.line = element_line(colour = 'grey30',size = 0.2),
    panel.spacing=unit(2, "mm"),
    panel.border = element_rect(fill = NA,linetype = 'solid',linewidth = 1),
    panel.grid.major.y = element_line(),
    panel.grid.major = element_line(size = 3,linetype = 1),
    strip.text.x = element_text(size=12, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),
    strip.background = element_rect(
      colour='black',
      fill=c('#2878B5','#2878B5'),
      size = 1),
    strip.placement = "inside"
  )

strip_back_color <- c(
  '#2878B5','#9AC9DB','#F8AC8C',
  '#C82423','#FF8884','#A9B8C6',
  '#96C37D','#F3D266','#C497B2',
  '#14517c',"#BD956A", '#585658'
)

#### PD1 VS C5aRA_PD1 --------------------------------------------------------
data_plot$dataset %>% unique()
data_plot_sub1 <- data_plot %>% 
  filter(
    dataset %in% c('PD1','C5aRA_PD1'),
    source %in% c('Mac','Mono','DC')
  ) %>% 
  mutate(
    dataset = dataset %>% factor(levels = c('PD1','C5aRA_PD1'))
  )
p2 <- ggplot(data = data_plot_sub1,
             aes(
               x = dataset,
               y = interaction_name_2,
               color = prob,
               size = pval
             )) +
  geom_point() +
  facet_grid(~ source+target) +
  scale_radius(
    range = c(4, 6),
    breaks = sort(unique(data_plot$pval)),
    labels = names(values)[values %in% sort(unique(data_plot$pval))],
    name = "p-value"
  ) +
  scale_colour_gradientn(
    colors = colorRampPalette(color.use)(99),
    na.value = "white",
    limits = c(quantile(data_plot$prob,
                        0, na.rm = T), quantile(data_plot$prob, 1, na.rm = T)),
    breaks = c(quantile(data_plot$prob, 0, na.rm = T), quantile(data_plot$prob,
                                                                1, na.rm = T)),
    labels = c("min", "max")
  ) +
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
  theme_classic() +
  theme(
    legend.text = element_text(size= 12),
    legend.title= element_text(size= 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = -45,hjust = 0,size=16, color="black",face="bold"),
    axis.text.y = element_text(
      size = 16,
      color = "black",
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.title.y = element_blank(),
    axis.text.y.right = element_blank(),
    legend.justification = c(0, 1),
    plot.title = element_text(hjust = 0.5,vjust = 6),
    plot.margin = margin(30,2,2,2),
    axis.line = element_line(colour = 'grey30',size = 0.2), 
    panel.spacing=unit(2, "mm"), 
    panel.border = element_rect(fill = NA,linetype = 'solid',linewidth = 1),
    panel.grid.major.y = element_line(),
    panel.grid.major = element_line(size = 3,linetype = 1),
    strip.text.x = element_text(size=12, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),
    strip.background = element_rect(
      colour='black',
      fill=c('#2878B5','#2878B5'),
      size = 1),
    strip.placement = "inside"
  )
p <- p1/p2
p

# Fig5 G cellchat between Mono and T --------------------------------------
scRNA_12 <- readRDS('D:/jiaoxi_12/scRNA_12_anno.rds')
scRNA_use <- scRNA_12
rm(scRNA_12)
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
  # "Neutrophil",
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
#### scRNA_PD1  ------------------------------------------------------------

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
if(FALSE){
  ##### TAM ===> T ====
  source_cell <- c(
    "Mac_C1qa","Mac_Vegfa","MonoMac_Ccr2","Mac_Mki67/C1qa"#,"Mono_Isg15"
  )
  target_cell <- c(
    "Dysfunction_T cell","Cytolytic_T cell"#,"Gzmb_T cell"
    # "Proliferative CD8+ T","CD8+ T memory"
  )
  cell_type_order <- c(source_cell,target_cell) %>% 
    .[. %in% rownames(cellchat_PD1@net$count)]
  color.use <- my36colors[cell_type_order]
  netVisual_circle_self(
    net = cellchat_PD1@net$count, 
    cell_type_order = cell_type_order,
    angle_coords = -90,
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
    margin = c(0,0,0,0)
  )
  ##### T ===> TAM====
  par(mfrow = c(1,2), xpd=TRUE,oma = c(0, 0, 0, 0),mar = c(0, 0, 0, 2))
  source_cell <- c(
    # "Dysfunction_T cell","Proliferative CD8+ T"
    "Cytolytic_T cell","CD8+ T memory","Gzmb_T cell"
  )
  target_cell <- c("Mac_C1qa",
                   "Mac_Vegfa",
                   "MonoMac_Ccr2",
                   "Mac_Mki67/C1qa",
                   "Mono_Isg15"
  )
  cell_type_order <- c(source_cell,target_cell) %>% 
    .[. %in% rownames(cellchat_PD1@net$count)]
  color.use <- my36colors[cell_type_order]
  netVisual_circle_self(
    net = cellchat_PD1@net$count, 
    cell_type_order = cell_type_order,
    angle_coords = -210,
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
    title.name = NULL,
    margin = c(0,0,0,0)
  )
}
#### scRNA_C5aRA_PD1  ------------------------------------------------------------

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
colnames(CellChatDB$interaction) # 查看一下数据库信息
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
if(FALSE){
  ##### TAM ===> T ====
  source_cell <- c(
    "Mac_C1qa","Mac_Vegfa","MonoMac_Ccr2","Mac_Mki67/C1qa"#,"Mono_Isg15"
  )
  target_cell <- c(
    "Dysfunction_T cell","Cytolytic_T cell"#,"Gzmb_T cell"
    # "Proliferative CD8+ T","CD8+ T memory"
  )
  cell_type_order <- c(source_cell,target_cell) %>% 
    .[. %in% rownames(cellchat_C5aRa_PD1@net$count)]
  color.use <- my36colors[cell_type_order]
  netVisual_circle_self(
    net = cellchat_C5aRa_PD1@net$count, 
    cell_type_order = cell_type_order,
    angle_coords = -90,
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
    arrow.width = 2,
    arrow.size = 1,
    title.name = 'C5aRa_PD1',
    margin = c(0,0,0,0)
  )
  #### T ===> TAM====
  source_cell <- c(
    # "Dysfunction_T cell","Proliferative CD8+ T"
    "Cytolytic_T cell","CD8+ T memory","Gzmb_T cell"
  )
  target_cell <- c("Mac_C1qa",
                   "Mac_Vegfa",
                   "MonoMac_Ccr2",
                   "Mac_Mki67/C1qa",
                   "Mono_Isg15"
  )
  cell_type_order <- c(source_cell,target_cell) %>% 
    .[. %in% rownames(cellchat_C5aRa_PD1@net$count)]
  color.use <- my36colors[cell_type_order]
  netVisual_circle_self(
    net = cellchat_C5aRa_PD1@net$count, 
    cell_type_order = cell_type_order,
    angle_coords = -210,#-210,
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
    arrow.width = 2,
    arrow.size = 1,
    title.name = NULL,
    margin = c(0,0,0,0)
  )
}

# Fig5 H Nichenet between Mono_Isg and Gzmb+ CD8+ T------------------------------------------------------------------
get_lfc_batch <-
  function (seurat_obj,
            condition_colname,
            condition_oi,
            condition_reference,
            celltype_col = "celltype",
            expression_pct = 0.1){
    requireNamespace("Seurat")
    requireNamespace("dplyr")
    seuratObj_sender = seurat_obj
    seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[[condition_colname,
                                                                            drop = TRUE]])
    DE_table_sender = FindMarkers(
      object = seuratObj_sender,
      ident.1 = condition_oi,
      ident.2 = condition_reference,
      min.pct = expression_pct,
      logfc.threshold = 0.05
    ) %>%
      rownames_to_column("gene")
    SeuratV4 = c("avg_log2FC") %in% colnames(DE_table_sender)
    if (SeuratV4 == TRUE) {
      DE_table_sender = DE_table_sender %>% as_tibble() %>%
        select(-p_val) %>% select(gene, avg_log2FC)
    }
    else {
      DE_table_sender = DE_table_sender %>% as_tibble() %>%
        select(-p_val) %>% select(gene, avg_logFC)
    }
    colnames(DE_table_sender) = c("gene", condition_oi)
    return(DE_table_sender)
  }

make_heatmap_ggplot_self <-
  function (matrix,
            y_name,
            x_name,
            y_axis = TRUE,
            x_axis = TRUE,
            x_axis_position = "top",
            legend_position = "top",
            color = "blue",
            color_low = 'whitesmoke',
            tile_color = "black",
            legend_title = "score",
            ...){
    if (!is.matrix(matrix))
      stop("matrix should be a matrix")
    if (!is.character(y_name) | length(y_name) != 1)
      stop("y_name should be a character vector of length 1")
    if (!is.character(x_name) | length(x_name) != 1)
      stop("x_name should be a character vector of length 1")
    if (!is.logical(y_axis) | length(y_axis) != 1)
      stop("y_axis should be a TRUE or FALSE")
    if (!is.logical(x_axis) | length(x_axis) != 1)
      stop("x_axis should be a TRUE or FALSE")
    if ((x_axis_position %in% c("top", "bottom")) == FALSE)
      stop("x_axis_position should be top or bottom")
    if ((legend_position %in% c("top", "bottom", "left", "right",
                                "none")) == FALSE)
      stop("legend_position should be top, bottom, left, right or none")
    if (!is.character(color) | length(color) != 1)
      stop("color should be character vector of length 1")
    requireNamespace("dplyr")
    requireNamespace("ggplot2")
    matrix_df_vis = matrix %>% data.frame() %>% rownames_to_column("y") %>%
      as_tibble() %>% gather(x, "score",-y) %>% mutate(
        y = factor(y,
                   levels = rownames(matrix), ordered = TRUE),
        x = factor(x,
                   levels = colnames(matrix), ordered = TRUE)
      )
    plot_object = matrix_df_vis %>% 
      ggplot(aes(x, y, fill = score)) +
      geom_tile(color = tile_color, size = 0.2) + 
      scale_fill_gradient(
        low = color_low,
        high = color,
        breaks = c(min(matrix_df_vis$score), max(matrix_df_vis$score)),
        labels = c(round(min(matrix_df_vis$score), 3),
                   round(max(matrix_df_vis$score), 3))
      ) + 
      theme_minimal()
    if (x_axis == FALSE) {
      if (y_axis == TRUE) {
        plot_object = plot_object + theme(
          panel.grid.minor = element_line(color = "transparent"),
          panel.grid.major = element_line(color = "transparent"),
          legend.position = legend_position,
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(...),
          axis.text.y = element_text(...)
        )
        plot_object = plot_object + ylab(paste0(y_name))
      }
      else if (y_axis == FALSE) {
        plot_object = plot_object + theme(
          panel.grid.minor = element_line(color = "transparent"),
          panel.grid.major = element_line(color = "transparent"),
          legend.position = legend_position,
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()
        )
        plot_object = plot_object
      }
    }
    else if (x_axis == TRUE) {
      if (y_axis == TRUE) {
        plot_object = plot_object + theme(
          panel.grid.minor = element_line(color = "transparent"),
          panel.grid.major = element_line(color = "transparent"),
          legend.position = legend_position,
          axis.ticks = element_line(size = 0),
          axis.text.x.top = element_text(angle = 90, hjust = 0,vjust = 0.5,
                                         ...),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1, ...),
          axis.title = element_text(...),
          axis.text.y = element_text(...)
        )
        plot_object = plot_object + scale_x_discrete(position = x_axis_position) +
          xlab(paste0(x_name)) + ylab(paste0(y_name))
      }
      else if (y_axis == FALSE) {
        plot_object = plot_object + theme(
          panel.grid.minor = element_line(color = "transparent"),
          panel.grid.major = element_line(color = "transparent"),
          legend.position = legend_position,
          axis.ticks = element_line(size = 0),
          axis.text.x.top = element_text(angle = 90, hjust = 0,vjust = 0.5,
                                         ...),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1, ...),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()
        )
        plot_object = plot_object + scale_x_discrete(position = x_axis_position) +
          xlab(paste0(x_name))
      }
    }
    plot_object = plot_object + labs(fill = legend_title)
  }
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
  ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
  head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
  
  head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
  
  
  seuratObj = alias_to_symbol_seurat(scRNA_11_sub, "mouse")
  Idents(seuratObj) <- seuratObj$anno_res
  expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = 0.10)
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  list_expressed_genes_sender = sender_celltypes %>% 
    unique() %>% 
    lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
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
    make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = 
    colnames(active_ligand_target_links) %>% 
    make.names() # make.names() for heatmap visualization of genes like H2-T23
  
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
    # labels = c("A", "B",'C')
    # rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target))
  )
  # figures_without_legend
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    # ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
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
  # 显示最终图表
  # final_plot
  return(final_plot)
}
scRNA_11 <- readRDS('D:/jiaoxi/scRNA_11_res.rds')
scRNA_11$batch %>% unique()
scRNA_11$anno_res %>% unique()
#### Mono_Isg ==> Gzmb_T cell -----------------------------------------------------------
sender_celltypes = c(
  "Mono_Isg"
)
receiver = c(
  "Gzmb_T cell"
)
condition_oi = "PD1"
condition_reference = "C5aRA_PD1"
scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% c(
  sender_celltypes,
  receiver
)) %>% 
  subset(subset = batch %in% c(condition_oi,condition_reference))
scRNA_11_sub@meta.data %>% group_by(batch,anno_res) %>% summarise(Counts = n())
target_gene <- c('Gzmb')
p3 <- res_plot_self(
  scRNA_11_sub = scRNA_11_sub,
  diff_gene_from = condition_reference,
  upstream_ligands_n = 6,
  top_n_target_gene = 250,
  target_gene = target_gene,
  target_ligands = c('Ccl3','Ccl4','Il15','Il18'),
  p_ligand_target_network_x_size = 30,
  axis_text_y = 30,
  num_show_gene = 50,
  min.pct = 0.01,
  axis_title_y = 'Mono-Isg cell',
  axis_title_x = 'Gzmb T cells')
p3
#### CD8 T ==> TAM -----------------------------------------------------------
sender_celltypes = c(
  # "Cytolytic_T cell"
  # "Proliferative CD8+ T",
  # "CD8+ T memory",
  "Gzmb_T cell",
  # "Dysfunction_T cell",
  "Cytolytic_T cell"
)
receiver = c(
  "Mono_Arhgap26",
  "Mono_Ly6i",
  "Mono_Isg",
  "Mac/Mono_Arg1",
  "Mac/Mono_Cx3cr1",
  "Mac/Mono_Ccl8",
  "Mac/Mono_Mki67"
)
condition_oi = "PD1"
condition_reference = "C5aRA_PD1"
scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% c(
  sender_celltypes,
  receiver
)) %>% 
  subset(subset = batch %in% c(condition_oi,condition_reference))
scRNA_11_sub@meta.data %>% group_by(batch,anno_res) %>% summarise(Counts = n())
target_gene <- c("B2m","Cd38","Cd40","Cd74","Cxcl1","Cxcl10","Cxcl9","H2.Ab1","H2.D1",
                 "H2.K1","H2.M3","H2.Q1","H2.Q10","H2.Q2","H2.Q4","H2.Q7","H2.T23",
                 "H2.T24","Ifit2","Irf1","Irf8","Tap1","Tap2")
target_ligands <- c("Tnf","Ifng")
p3 <- res_plot_self(
  scRNA_11_sub = scRNA_11_sub,
  diff_gene_from = condition_reference,
  upstream_ligands_n = 6,
  top_n_target_gene = 200,
  target_gene = target_gene,
  target_ligands = target_ligands,
  num_show_gene = 36,
  axis_title_y = 'Cytolytic & Gzmb T',
  axis_title_x = 'Macrophage/Monocyte',
  p_ligand_target_network_x_size = 30)
p3
# Fig6 --------------------------------------------------------------------

# Fig6 A umap and Ro/e score------------------------------------------------------------------
## umap ====
scRNA_11_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  "NPBNs","Exhausted TAN","interferon-stimulated NAN"
))
cell_label_correct <- c(
  "S100A8-Neu","Ccl4-Neu","ISG-Neu"
)
names(cell_label_correct) <- c(
  "NPBNs","Exhausted TAN","interferon-stimulated NAN"
)
scRNA_11_tmp$anno_tmp <- cell_label_correct[scRNA_11_tmp$anno_res] %>% as.character()
scRNA_11_tmp$anno_tmp %>% unique()
scRNA_11_tmp <-
  FindVariableFeatures(object = scRNA_11_tmp, nfeatures = 3000)
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
  mutate(anno_tmp = factor(anno_tmp,levels = cell_leves)) %>% 
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
my3colors <- c('#d6616b',"#7698b3","#2ca02c")
ggplot() +
  geom_point(
    data = data_plot,
    aes(x = umap_1,y = umap_2,color = label_inner),
    size = 3.8) +
  geom_label_repel(data = cell_label_loc,
                   aes(
                     x = umap_1,
                     y = umap_2,
                     label = anno_tmp),
                   label.size = 0,
                   size = 12,
                   alpha = 0.8,
                   fill = 'grey90',
                   segment.color = NA,
                   show.legend = FALSE)+
  #x
  geom_segment(aes(
    x = data_plot$umap_1 %>% min() * 1.2, 
    y = data_plot$umap_2 %>% min() * 1.1 , 
    xend = (data_plot$umap_1 %>% min())* 1.2+
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2, 
    yend = (data_plot$umap_2 %>% min())* 1.1
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.5) +
  #Y
  geom_segment(aes(
    x = data_plot$umap_1 %>% min()* 1.2, 
    y = data_plot$umap_2 %>% min() * 1.1, 
    xend = (data_plot$umap_1 %>% min())* 1.2, 
    yend = (data_plot$umap_2 %>% min()) + 
      (max(data_plot$umap_1) - min(data_plot$umap_1))*0.2 
  ),
  arrow = arrow(length = unit(0.3, "cm")),
  size = 0.5) +
  scale_color_manual(
    values = my3colors
    # breaks = 1:length(cell_leves),
    # label = label_legend
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

## Ro/e --------------------------------------------------------------------
scRNA_11$anno_res %>% unique()
scRNA_11_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  "Cancer cell","Fibroblast",
  "Dysfunction_T cell","Cytolytic_T cell","Gzmb_T cell",
  "Proliferative CD8+ T","CD8+ T memory",
  "Naive T","Treg","Proliferative Treg","CD4+ T memory","NK",
  "Mono_Arhgap26","Mono_Ly6i","Mono_Isg",
  "Mac/Mono_Arg1","Mac/Mono_Cx3cr1","Mac/Mono_Ccl8",
  "Mac/Mono_Mki67",
  # "NPBNs","Exhausted TAN","interferon-stimulated NAN",
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
  # "NPBNs","Exhausted TAN","interferon-stimulated NAN",
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
  # "NPBNs","Exhausted TAN","interferon-stimulated NAN",
  "S100A8-Neu","Ccl4-Neu","ISG-Neu",
  "Mast_Mcpt2","B cell","cDC1_Xcr1","tDC"
)
scRNA_11_tmp$anno_tmp <- cell_label_correct[scRNA_11_tmp$anno_res] %>% as.character()
scRNA_11_tmp$anno_tmp %>% unique()
meta_data_sub <- scRNA_11_tmp@meta.data
meta_data_sub$anno_res  <- meta_data_sub$anno_tmp
celltype_select <- c(
  "S100A8-Neu","Ccl4-Neu","ISG-Neu"
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
    # group = factor(group,levels = c('Depletion','Enrichment')),
  ) %>% 
  group_by(celltype) %>% 
  mutate(roe_scale = scale(`Ro/e`,center = FALSE)[,1]) 
# ungroup() %>% 
# mutate(roe_scale = roe_scale - min(roe_scale))
data_ggplot$roe_scale <- data_ggplot$`Ro/e`
limit_min <- data_ggplot$roe_scale %>% min()
limit_max <- data_ggplot$roe_scale %>% max()

ggplot(data = data_ggplot,
       mapping = aes(x = batch,y = celltype)) +
  geom_tile(aes(fill = roe_scale)) +
  geom_text(aes(label = round(roe_scale,digits = 2))) +
  guides(
    # label = guide_legend(title.theme = element_blank(),override.aes = list(size = 6)),
    fill = guide_colorbar(title = 'Ro/e',title.vjust = 1)
  ) +
  scale_fill_gradient(low = 'grey90',high =  '#3131F2',limit = c(limit_min,limit_max)) +
  # scale_y_discrete(position = 'right') +
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
    #分面标签样式
    strip.background = element_blank()
  )

# Fig6 B Dotplot-------------------------------------------------------------------

scRNA_11_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  "NPBNs","Exhausted TAN","interferon-stimulated NAN"
))
cell_label_correct <- c(
  "S100A8-Neu","Ccl4-Neu","ISG-Neu"
)
names(cell_label_correct) <- c(
  "NPBNs","Exhausted TAN","interferon-stimulated NAN"
)
scRNA_11_tmp$anno_tmp <- cell_label_correct[scRNA_11_tmp$anno_res] %>% as.character()
scRNA_11_tmp$anno_tmp %>% unique()
scRNA_11_tmp$anno_tmp <- factor(scRNA_11_tmp$anno_tmp,levels = cell_label_correct)
gene_select <- c(
  "S100a8","S100a9","S100a11","Spi1","Retnlg","Ccl3",
  "Ccl4","Ddit3","Spp1","Ctsb","Il1rn","Slc2a1","Cxcl2",
  "Ifit1","Ifit2","Ifit3","Cxcl10","Irf1","Ly6e","Gbp5",
  "H2-D1","H2-K1","Isg15","Isg20"
)
DotPlot(object = scRNA_11_tmp,
        features = gene_select,
        group.by = 'anno_tmp') +
  coord_flip() +
  scale_color_gradient(
    name = 'Mean Exp',
    low = 'grey90',
    high = '#3131F2',
    breaks = c(-1,0,1),
    labels = c(-1,0,1)) +
  scale_size_continuous(name = 'Fraction') +
  guides(
    color = guide_colourbar(nrow = 1,title.position = 'left',direction = 'vertical'),
    size = guide_legend(ncol = 1,title.position = 'left',direction = 'vertical')
  ) +
  labs(x = '', y= 'Features') +
  theme(
    axis.text.x = element_text(angle = -45,hjust = 0,size = 24),
    axis.text.y = element_text(size = 18),
    legend.direction = 'vertical',
    legend.title = element_text(angle = 90,hjust = 0.5,vjust = 0.5),
    axis.title.x = element_blank()
  )

# Fig6 C GSEA -------------------------------------------------------------

options(stringsAsFactors = F)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)  #install.packages("msigdbr")
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
library(pheatmap)
gene.expr <-  as.matrix(scRNA_neu[["RNA"]]@data)
dim(gene.expr)
GO_df_all <- msigdbr(species = "Mus musculus",
                     category = "C5")  
GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
GO_df$tmp <- GO_df$gs_name %>% 
  str_remove('GOBP_') %>% 
  str_remove('GOCC_') %>% 
  str_remove('GOMF_') %>% 
  str_replace_all('_',' ') %>% 
  str_to_title()
GO_df_use <- GO_df %>% 
  filter(
    tmp %in% signature_list
  )
tmp <- GO_df %>% 
  filter(
    tmp %in% signature_list
  ) %>% 
  dplyr::select(gs_exact_source,gs_name) %>% 
  filter(!duplicated(gs_exact_source)) %>% 
  mutate(
    gs_name = gs_name %>% 
      str_remove('GOBP_') %>% 
      str_remove('GOCC_') %>% 
      str_remove('GOMF_') %>%
      str_replace_all('_',' ') %>% 
      str_to_title()
  )
go_use <- split(GO_df_use$gene_symbol,GO_df_use$gs_name) %>% 
  lapply(.,function(x){paste(x,collapse = ',')}) %>% 
  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column('ge')
go_list <- split(GO_df$gene_symbol, GO_df$gs_name)
geneset <- go_list
gsva_mat <- gsva(expr=gene.expr, 
                 gset.idx.list=geneset, 
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())#调用所有核

meta <- scRNA_neu@meta.data %>% arrange(anno_res)
data <- gsva_mat[,rownames(meta)]
group <- factor(meta[,"anno_res"],ordered = F)
levels(group)
data1 <-NULL
for(i in levels(group)){
  # i= levels(group)[1]
  ind <-which(group==i)
  dat <- apply(data[,ind], 1, mean)
  data1 <-cbind(data1,dat)
}
colnames(data1) <-levels(group)
result<- t(scale(t(data1)))
tmp <- rownames(result) %>% 
  str_remove('GOBP_') %>% 
  str_remove('GOCC_') %>% 
  str_remove('GOMF_') %>%
  str_replace_all('_',' ')
pathway_select <- c(
  "Chemoattractant activity","Response to hyperoxia",
  "Chemokine activity","Pyruvate metabolic process",
  "Response to type I interferon",
  "Defense response to virus","Innate immune response",
  "Myotube cell development","Muscle fiber development",
  "Kidney vasculature development","Ribosome",
  "Ncrna metabolic process","Regulation of cell development",
  "Myeloid cell development"
) %>% toupper()
table(pathway_select %in% tmp)
tmp <- pathway_select[is.na(match(pathway_select,tmp))] %>% 
  str_replace_all(' ','_') %>% 
  str_to_upper()
lapply(tmp, function(x){return(
  rownames(result)[grepl(x,rownames(result))]
)}) %>% unlist() %>% unique()
rownames(result)[grepl('DEFENSE_RESPONSE_TO_VIRUS',rownames(result))]
rownames(result)[grepl('VASCULATURE',rownames(result))]
result_sub <- result
rownames(result_sub) <- rownames(result_sub) %>% 
  str_remove('GOBP_') %>% 
  str_remove('GOCC_') %>% 
  str_remove('GOMF_') %>% 
  str_replace_all('_',' ') %>% 
  str_to_title()
signature_list <- c(
  "Response To Type I Interferon",  
  "Innate Immune Response", 
  "Regulation Of Cell Development", 
  "Myotube Cell Development",   
  "Myeloid Cell Development",   
  "Chemokine Activity", 
  "Chemoattractant Activity",   
  "Pyruvate Metabolic Process",
  "Response To Hyperoxia",  
  "Regulation Of Cellular Response To Hypoxia",
  "Ncrna Metabolic Process",
  "Regulation Of Glycolytic Process", 
  "Ribosome"
)
signature_list <- rownames(result_sub)[grepl('Angiogenesis',rownames(result_sub),ignore.case = TRUE)]
cell_type_order <- c("S100A8-Neu","Ccl4-Neu","ISG-Neu")
result_sub <- result_sub[signature_list,cell_type_order]
pdf('./GSVA_Neu.pdf', width = 6, height = 7)
ComplexHeatmap::pheatmap(
  mat = result_sub,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  scale = "row",
  cluster_rows = F,cellwidth = 25,cellheight = 25,
  cluster_cols = F,
  show_rownames = T,
  show_colnames = T,
  fontsize_col = 18,
  fontsize_row = 18,
  angle_col = '315',
  heatmap_legend_param = list(title = 'GSVA Score',
                              title_gp = grid::gpar(fontsize = 18, fontface = "bold"),
                              title_position = "lefttop-rot",
                              legend_direction = 'vertical',
                              legend_position = "topright",
                              labels_gp = grid::gpar(fontsize = 16, fontface = "bold"))
)
dev.off()

# Fig6D AUCell ------------------------------------------------------------
library(AUCell) 
library(clusterProfiler)
library(pheatmap)
library(reshape2)
library(tidyr)
library(ggpubr)
# 定义函数 --------------------------------------------------------------------

violin_aucell_cluster <- function(data_plot){
  p <- ggplot(data_plot, aes(x = anno_res, y = AUCell_score,fill = anno_res)) +
    # geom_violin() + 
    geom_violin(trim=FALSE,color="white",scale = 'width') +
    geom_boxplot(width=0.2,position=position_dodge2(0.9),show.legend = FALSE,size = 1)+  
    scale_x_discrete(labels = function(x) stringr::str_wrap(x,width = 30,whitespace_only = FALSE)) +
    scale_fill_manual(values = c('lightslategrey','#0099B4', '#42B540', '#ED0000', '#00468B') %>% rev()) +
    stat_compare_means(aes(group =  as.factor(anno_res),x = as.factor(anno_res)),
                       # label = "p.format",
                       label = paste0("p = ", after_stat("p.format")),
                       method = "wilcox.test",
                       hide.ns = FALSE,
                       # symnum.args = list(
                       #   cutpoints = c(0, 0.01, 0.05, 1),
                       #   symbols = c("***", "*", "NS")
                       # ),
                       comparisons =
                         list(
                           c('0','1'),
                           c('0','2'),
                           c('1','2')
                         ),
                       step.increase = 0.1
    ) +
    facet_wrap(~GO_path,nrow = 1) +
    theme_bw()+
    theme(
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = -45,
        hjust = 0,
        size = 16,
        color = "black",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 16,
        color = "black",
        face = "bold",
        vjust = 0.5,
        hjust = 0.5
      ),
      axis.title.y = element_text(
        size = 25,
        color = "black",
        face = "bold"
      ),
      axis.text.y.right = element_blank(),
      legend.justification = c(0, 1),
      plot.title = element_text(hjust = 0.5, vjust = 6),
      plot.margin = margin(30, 2, 2, 2),
      axis.line = element_line(colour = 'grey30', size = 0.2),
      panel.spacing = unit(2, "mm"),
      panel.border = element_rect(
        fill = NA,
        linetype = 'solid',
        linewidth = 1
      ),
      panel.grid.major.y = element_line(),
      panel.grid.major = element_line(size = 3, linetype = 1),
      strip.text.x = element_text(
        size = 16,
        face = "bold",
        color = "#FFFFFF",
        vjust = 0.5,
        margin = margin(b = 3, t = 3,r = 3,l = 5)
      ),
      #分面标签样式
      strip.background = element_rect(
        colour = 'black',
        fill = '#2878B5',
        size = 1
      ),
      strip.placement = "inside"
    )
  return(p)
}
violin_aucell_batch <- function(data_plot){
  p <- ggplot(data_plot, aes(x = anno_res, y = AUCell_score,fill = anno_res)) +
    # geom_violin() + 
    geom_violin(trim=FALSE,color="white",scale = 'width') +
    geom_boxplot(width=0.2,position=position_dodge2(0.9),show.legend = FALSE,size = 1)+  
    scale_x_discrete(labels = function(x) stringr::str_wrap(x,width = 30,whitespace_only = FALSE)) +
    scale_fill_manual(values = c('lightslategrey','#0099B4', '#42B540', '#ED0000', '#00468B') %>% rev()) +
    stat_compare_means(aes(group =  as.factor(anno_res),x = as.factor(anno_res)),
                       # label = "p.format",
                       label = paste0("p = ", after_stat("p.format")),
                       method = "wilcox.test",
                       hide.ns = FALSE,
                       # symnum.args = list(
                       #   cutpoints = c(0, 0.01, 0.05, 1),
                       #   symbols = c("***", "*", "NS")
                       # ),
                       comparisons = 
                         list(
                           c('MOCK','C5aRa'),
                           c('PD1','C5aRA_PD1')
                         ),
                       step.increase = 0.1
    ) +
    facet_wrap(~GO_path,nrow = 3) +
    theme_bw()+
    theme(
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = -45,
        hjust = 0,
        size = 16,
        color = "black",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 16,
        color = "black",
        face = "bold",
        vjust = 0.5,
        hjust = 0.5
      ),
      axis.title.y = element_text(
        size = 25,
        color = "black",
        face = "bold"
      ),
      axis.text.y.right = element_blank(),
      legend.justification = c(0, 1),
      plot.title = element_text(hjust = 0.5, vjust = 6),
      plot.margin = margin(30, 2, 2, 2),
      axis.line = element_line(colour = 'grey30', size = 0.2),
      panel.spacing = unit(2, "mm"),
      panel.border = element_rect(
        fill = NA,
        linetype = 'solid',
        linewidth = 1
      ),
      panel.grid.major.y = element_line(),
      panel.grid.major = element_line(size = 3, linetype = 1),
      strip.text.x = element_text(
        size = 9,
        face = "bold",
        color = "#FFFFFF",
        vjust = 0.5,
        margin = margin(b = 3, t = 3,r = 3,l = 5)
      ),
      #分面标签样式
      strip.background = element_rect(
        colour = 'black',
        fill = '#2878B5',
        size = 1
      ),
      strip.placement = "inside"
    )
  return(p)
}
violin_aucell_batch_single <- function(data_plot,GO_path){
  data_plot <- data_plot %>% 
    filter(GO_path ==GO_path) %>% 
    mutate(facet_group = ifelse(
      anno_res %in% c('MOCK', 'C5aRa'),
      yes = 'MOCK vs C5aRa',
      no = 'PD1 vs C5aRA_PD1'
    ))
  p <- ggplot(data_plot, aes(x = GO_path, y = AUCell_score,fill = anno_res)) +
    # geom_violin() + 
    geom_violin(trim=FALSE,color="white",scale = 'width') +
    geom_boxplot(width=0.2,position=position_dodge(0.9),show.legend = FALSE,size = 1)+  
    scale_x_discrete(labels = function(x) stringr::str_wrap(x,width = 30,whitespace_only = FALSE)) +
    scale_fill_manual(values = c('lightslategrey','#0099B4', '#42B540', '#ED0000', '#00468B') %>% rev()) +
    stat_compare_means(aes(group = anno_res),
                       label = "p.format",
                       size = 5,
                       # label = paste0("p = ", after_stat("p.format")),
                       method = "wilcox.test",
                       hide.ns = FALSE,
                       # symnum.args = list(
                       #   cutpoints = c(0, 0.01, 0.05, 1),
                       #   symbols = c("***", "*", "NS")
                       # ),
                       # comparisons =
                       #   list(
                       #     c('MOCK','C5aRa'),
                       #     c('PD1','C5aRA_PD1')
                       #   ),
                       step.increase = 0.1
    ) +
    facet_wrap(~facet_group,scales="free_x") +
    theme_bw() + 
    theme(plot.title = element_text(size = 12,color="black",hjust = 0.5), 
          axis.text.x = element_text(angle = 45, hjust = 1,size = 20 ),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.text = element_text(size= 12),
          legend.title= element_text(size= 12)) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=16, color="black",face="bold"),
      axis.title.y = element_text(size=22,colour = 'black',face="bold",vjust = 0,hjust = 0.5),
      axis.text.y.right = element_blank(),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 24),
      plot.margin = margin(2,2,2,2),
      axis.line = element_line(colour = 'grey30',size = 0.2),
      panel.spacing=unit(2, "mm"), 
      panel.border = element_rect(fill = NA,linetype = 'solid',linewidth = 1),
      panel.grid.major.y = element_line(),
      strip.text.x = element_text(size=22, face="bold",color = "#FFFFFF",
                                  vjust = 0.5,margin = margin(b = 3,t=3)),
      strip.background = element_rect(colour="grey30", fill="grey60",size = 1)
    )
  return(p)
}
violin_aucell_batch_single_v2 <- function(data_plot,GO_path){
  data_plot <- data_plot %>% 
    filter(GO_path ==GO_path)
  p <- ggplot(data_plot, aes(x = anno_res, y = AUCell_score,fill = anno_res)) +
    # geom_violin() + 
    geom_violin(trim=FALSE,color="white",scale = 'width') +
    geom_boxplot(width=0.2,position=position_dodge2(0.9),show.legend = FALSE,size = 0.5)+  
    scale_x_discrete(labels = function(x) stringr::str_wrap(x,width = 30,whitespace_only = FALSE)) +
    scale_fill_manual(values = c('lightslategrey','#0099B4', '#42B540', '#ED0000', '#00468B') %>% rev()) +
    stat_compare_means(aes(group =  as.factor(anno_res),x = as.factor(anno_res)),
                       # label = "p.format",
                       label = paste0("p = ", after_stat("p.format")),
                       method = "wilcox.test",
                       hide.ns = FALSE,
                       # symnum.args = list(
                       #   cutpoints = c(0, 0.01, 0.05, 1),
                       #   symbols = c("***", "*", "NS")
                       # ),
                       comparisons = 
                         list(
                           c('MOCK','C5aRa'),
                           c('PD1','C5aRA_PD1')
                         ),
                       step.increase = 0,
                       size = 10
    ) +
    # facet_wrap(~GO_path,nrow = 3) +
    theme_bw()+ 
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 22,color="black",face="bold",angle =45,hjust = 1,vjust = 1),
      axis.text.y = element_text(size=16, color="black",face="bold"),
      axis.title.y = element_text(size=22,colour = 'black',face="bold",vjust = 0,hjust = 0.5),
      axis.text.y.right = element_blank(),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5,vjust = 0.5,size = 24),
      plot.margin = margin(2,2,2,2),
      axis.line = element_line(colour = 'grey30',size = 0.2), 
    )
  p
  return(p)
}
res_list <- list()
data_1 <- read.csv('D:/jiaoxi_res/20231201/go_path/azurophil_granule.txt',row.names = NULL,sep = '\t')
res_list[['azurophil_granule']] <- data_1$MGI.Gene.Marker.ID
data_3 <- read.csv('D:/jiaoxi_res/20231201/go_path/neutrophil_chemotaxis.txt',row.names = NULL,sep = '\t')
res_list[['neutrophil_chemotaxis']] <- data_3$MGI.Gene.Marker.ID
data_4 <- read.csv('D:/jiaoxi_res/20231201/go_path/Neutrophil_killing.txt',row.names = NULL,sep = '\t')
res_list[['Neutrophil_killing']] <- data_4$MGI.Gene.Marker.ID %>% unique()
data_5 <- read.csv('D:/jiaoxi_res/20231201/go_path/positive_regulation_of_apoptotic_process.txt',row.names = NULL,sep = '\t')
res_list[['positive_regulation_of_apoptotic_process']] <- data_5$MGI.Gene.Marker.ID
data_6 <- read.csv('D:/jiaoxi_res/20231201/go_path/tertiary_granule.txt',row.names = NULL,sep = '\t')
res_list[['tertiary_granule']] <- data_6$MGI.Gene.Marker.ID
geneSet <- readRDS('D:/jiaoxi/ssgsea/gene_tran/geneSet.rds')
res_list[['Interferon signaling']] <- geneSet[['Neutrophil:Interferon signaling']]
res_list[['Neutrophil degranulation']] <- geneSet[['Neutrophil:Neutrophil']]
human <- readRDS('D:/jiaoxi/ssgsea/gene_tran/human.rds')
mouse <- readRDS('D:/jiaoxi/ssgsea/gene_tran/mouse.rds')
for(signature_i in names(res_list)){
  gene_list <- res_list[[signature_i]] %>% toupper()
  geneMm <-
    getLDS(
      attributes = "hgnc_symbol",
      filters = "hgnc_symbol",
      values = gene_list,
      mart = human,
      attributesL = "mgi_symbol",
      martL = mouse,
      uniqueRows = TRUE
    )
  res_list[[signature_i]] <- geneMm$MGI.symbol
}
res_list <- lapply(res_list, function(x){
  x <- x[x %in% rownames(scRNA_neu_tmp)]
  tmp <- FetchData(scRNA_neu_tmp,vars = x)
  x[colSums(tmp)>0]
})
scRNA_neu_tmp <- subset(scRNA_11,subset = anno_res %in% c(
  'NPBNs','Exhausted TAN','interferon-stimulated NAN'
))

cells_rankings <- AUCell_buildRankings(scRNA_neu_tmp@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(
  res_list,
  cells_rankings, 
  aucMaxRank=nrow(cells_rankings)*0.1
)
data_plot <- cells_AUC@assays@data@listData[["AUC"]] %>%
  t() %>% as.data.frame() %>%
  mutate(anno_res = scRNA_neu_tmp$batch) %>%
  melt(id.vars = 'anno_res',
       value.name = 'AUCell_score',
       variable.name = 'GO_path') %>% 
  mutate(GO_path =  str_replace_all(GO_path,'_',' ')) %>% 
  mutate(anno_res =factor(anno_res,levels = c('MOCK','C5aRa','PD1','C5aRA_PD1'))) %>% 
  filter(anno_res %in% c('MOCK','C5aRa','C5aRA_PD1')) %>% 
  filter(GO_path %in% c('Neutrophil killing','Interferon signaling','Neutrophil degranulation'))
ggplot(data_plot, aes(x = anno_res, y = AUCell_score,fill = anno_res)) +
  geom_violin(trim=FALSE,color="white",scale = 'width') +
  geom_boxplot(width=0.2,position=position_dodge(0.9),show.legend = FALSE,size = 0.5)+ 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x,width = 30,whitespace_only = FALSE)) +
  scale_fill_manual(values = c('lightslategrey','#0099B4', '#ED0000', '#00468B') %>% rev()) +
  stat_compare_means(
    aes(x = anno_res),
    method = "t.test",
    step.increase = 0.15,
    label = "p.adj.signif",  
    comparisons = list(c('MOCK', 'C5aRa'), c('MOCK', 'C5aRA_PD1'))
  ) +
  facet_wrap(~GO_path,scale = 'free_y',ncol = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NULL,colour = 'black'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=16, color="black",face="bold"),
    axis.title.y = element_text(size=22,colour = 'black',face="bold",vjust = 0,hjust = 0.5),
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

# Fig6 E monocle3 of Nue --------------------------------------------------

library(Seurat)
library(tibble)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(monocle3)
scRNA_11 <- readRDS('D:/jiaoxi/scRNA_11_res.rds')
scRNA_neu <- subset(scRNA_11, subset = anno_res %in% c(
  'NPBNs','Exhausted TAN','interferon-stimulated NAN'
)
)
cell_label_correct <- c(
  "S100A8-Neu","Ccl4-Neu","ISG-Neu"
)
names(cell_label_correct) <- c(
  "NPBNs","Exhausted TAN","interferon-stimulated NAN"
)
scRNA_neu$anno_res <- cell_label_correct[scRNA_neu$anno_res] %>% as.character()
scRNA_neu$anno_res %>% unique()
scRNA_neu <-
  FindVariableFeatures(object = scRNA_neu, nfeatures = 3000)
scRNA_neu <- RunPCA(scRNA_neu,
                    features = VariableFeatures(object = scRNA_neu))
scRNA_neu <- FindNeighbors(scRNA_neu, dims = 1:16)
scRNA_neu <- FindClusters(scRNA_neu, resolution = 1)
scRNA_neu <- FindClusters(scRNA_neu, resolution = 3)
scRNA_neu <- RunUMAP(scRNA_neu, dims = 1:16, label = T)
DimPlot(
  object = scRNA_neu,
  group.by = c('anno_res'),
  label = TRUE
) + NoLegend()
## monocle3 ----------------------------------------------------------------
data <- GetAssayData(scRNA_neu, assay = 'RNA', slot = 'counts')
cell_metadata <- scRNA_neu@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds,num_dim = 16)
cds <- reduce_dimension(
  cds,
  preprocess_method = 'PCA',
  max_components = 2,
  umap.n_neighbors = 10)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA_neu, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
colnames(int.embed) <- NULL
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(
  cds = cds,
  color_cells_by = 'anno_res',
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
my36colors <-c('#14517c', "#96C37D", 
               "#D8383A","#BD956A", '#585658')
cds <- learn_graph(
  cds = cds,
  use_partition = TRUE,
  learn_graph_control = list(
    minimal_branch_len = 10
  )
)
cds <- order_cells(cds = cds)
p1 <- plot_cells(
  cds,
  color_cells_by = "anno_res",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 6,
  cell_size = 1.4,
  trajectory_graph_color = 'black',
  label_cell_groups = FALSE,
  trajectory_graph_segment_size = 1
) +
  scale_color_manual(values = my36colors) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 25),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.position = c(0.83,0.1),
    legend.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 20)
  )
p1
p2 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 1.4,
  trajectory_graph_segment_size = 1,
  trajectory_graph_color = 'black',
) +
  guides(
    color = guide_colorbar(
      title.position = 'top',
      label.position = 'bottom')) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 25),
    axis.text = element_blank(),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 20,hjust = 0.5,vjust = 0.5),
    legend.text = element_text(size = 20),
    legend.direction = 'horizontal',
    legend.position = c(0.85,0.1),
    legend.background = element_blank()
  )
p2
p1+p2
# Fig6 F anti-gene --------------------------------------------------------
data_geneset <- read.csv('D:/jiaoxi/ant-tumor_pro-tumor.tsv',sep = '\t')
data_geneset$member <- data_geneset$member %>% str_to_title()
data_geneset <- data_geneset[data_geneset$member %in% rownames(cds),]
gene_list <- split(data_geneset$member %>% str_to_title(),data_geneset$group)


gene.expr <-  as.matrix(scRNA_neu[["RNA"]]@data)
gsva_mat <- gsva(expr=gene.expr, 
                 gset.idx.list=gene_list, 
                 kcdf="Gaussian" ,
                 verbose=T, 
                 parallel.sz = parallel::detectCores())
meta <- scRNA_neu@meta.data
meta$anno_res %>% unique()
data_tmp <- gsva_mat[, rownames(meta)] %>% t()


pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(scRNA_neu@meta.data)]
data <- data_tmp %>% as.data.frame() %>%
  mutate(anno_res = meta$anno_res[match(rownames(data_tmp), rownames(meta))],
         batch = meta$batch[match(rownames(data_tmp), rownames(meta))],
         pseudotime = pseudotime[match(rownames(data_tmp),names(pseudotime))] %>% as.numeric()
  ) %>%
  filter(anno_res %in% c("S100A8-Neu",'ISG-Neu')) %>% 
  mutate(anno_res = factor(anno_res,levels = c("S100A8-Neu",'ISG-Neu'))) %>% 
  arrange(pseudotime) %>% 
  mutate(x_rank = 1:n())
names(data)
tmp <- lm(anti_tumor ~ poly(x_rank, 1), data = data) %>% summary(.)
p_value <- 1-pf(tmp$fstatistic[1],tmp$fstatistic[2],tmp$fstatistic[3])
p3 <- ggplot(data = data, aes(x = x_rank, y = (anti_tumor))) +
  geom_point(aes(color = anno_res)) +
  geom_smooth(
    method = "lm",
    formula = y ~ poly(x, 1),
    se = TRUE,
    color = 'grey30',
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.1,
    vjust = 1.1,
    size = 6,
    label = sprintf(
      "r = %.2f, p = %s",
      cor(data$x_rank, data$anti_tumor),
      formatC(
        p_value,
        format = "e",
        digits = 2
      )
    )
  ) +
  scale_color_manual(values = c("#D8383A","#96C37D")) +
  guides(color = guide_legend(title = 'Cell label', override.aes = list(size = 5))) +
  labs(x = 'pseudotime', 
       y = 'anti tumor Signature Score \n(Tyler Keeley et al.)') +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(),
    legend.position = 'top',
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 20)
  )
p3

# Fig6 G cellchat of Neu --------------------------------------------------
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
  # "Neutrophil",
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
  # "Neutrophil",
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


scRNA_C5aRA_PD1 <- subset(scRNA_use,subset = batch ==  'C5aRA_PD1')
Idents(scRNA_C5aRA_PD1) <- scRNA_C5aRA_PD1$res_low
scRNA_PD1 <- subset(scRNA_use,subset = batch ==  'PD1')
Idents(scRNA_PD1) <- scRNA_PD1$res_low


## scRNA_PD1  ------------------------------------------------------------

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
## scRNA_C5aRA_PD1  ------------------------------------------------------------
input <-
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

pdf(
  file = './cellchat_diff_Neu_anno_res.pdf',
  width = 9,
  height = 7,
  bg = 'white'
)
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
dev.off()


# Fig6 H nichenet ---------------------------------------------------------

sender_celltypes = c(
  "interferon-stimulated NAN"
)
receiver = c(
  "Gzmb_T cell"
)
condition_oi = "PD1"
condition_reference = "C5aRA_PD1"
scRNA_11_sub <- subset(scRNA_11,subset = anno_res %in% c(
  sender_celltypes,
  receiver
)) %>% 
  subset(subset = batch %in% c(condition_oi,condition_reference))
scRNA_11_sub@meta.data %>% group_by(batch,anno_res) %>% summarise(Counts = n())
target_gene <- c("Gzmb","H2-Q6","Irf7","Ccl4","Ccr5","H2-T22")
target_ligands <- c('Il15','Il18','Inf')
p3 <- res_plot_self(
  scRNA_11_sub = scRNA_11_sub,
  diffgene_source = 'diffgroup',
  diff_gene_from = condition_reference,
  upstream_ligands_n = 14,
  top_n_target_gene = 250,
  target_gene = target_gene,
  target_ligands = target_ligands,
  p_ligand_target_network_x_size = 30,
  axis_text_y = 20,
  num_show_gene = 50,
  axis_title_y = 'Neu-Isg cell',
  axis_title_x = 'Gzmb_T cells')
p3
# endline -----------------------------------------------------------------


