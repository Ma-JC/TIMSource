library(maftools)
library(dplyr)
library(stringr)
setwd('/mnt/data3/jiaoxi/新辅助治疗队列/vcf/vcf_anno_company/SNP/')
dir()

createOncoMatrix <- function (m, g = NULL, chatty = TRUE, add_missing = FALSE, cbio = FALSE) 
{
  if (is.null(g)) {
    stop("Please provde atleast two genes!")
  }
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, 
                     mafObj = FALSE)
  if (nrow(subMaf) == 0) {
    if (add_missing) {
      numericMatrix = matrix(data = 0, nrow = length(g), 
                             ncol = length(levels(getSampleSummary(x = m)[, 
                                                                          Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[, 
                                                               Tumor_Sample_Barcode])
      oncoMatrix = matrix(data = "", nrow = length(g), 
                          ncol = length(levels(getSampleSummary(x = m)[, 
                                                                       Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[, 
                                                            Tumor_Sample_Barcode])
      vc = c("")
      names(vc) = 0
      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, 
                  vc = vc))
    }
    else {
      return(NULL)
    }
  }
  if (add_missing) {
    subMaf[, `:=`(Hugo_Symbol, factor(x = Hugo_Symbol, levels = g))]
  }
  cnv_events = c(c("Amp", "Del"), as.character(subMaf[Variant_Type == 
                                                        "CNV"][, .N, Variant_Classification][, Variant_Classification]))
  cnv_events = unique(cnv_events)
  if (cbio) {
    vc = c("Nonstop_Mutation", "Frame_Shift_Del", "Missense_Mutation", 
           "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", 
           "In_Frame_Del", "In_Frame_Ins")
    vc.cbio = c("Truncating", "Truncating", "Missense", 
                "Truncating", "Truncating", "Truncating", "In-frame", 
                "In-frame")
    names(vc.cbio) = vc
    subMaf[, `:=`(Variant_Classification_temp, vc.cbio[as.character(subMaf$Variant_Classification)])]
    subMaf$Variant_Classification_temp = ifelse(test = is.na(subMaf$Variant_Classification_temp), 
                                                yes = as.character(subMaf$Variant_Classification), 
                                                no = subMaf$Variant_Classification_temp)
    subMaf[, `:=`(Variant_Classification, as.factor(as.character(Variant_Classification_temp)))]
    subMaf[, `:=`(Variant_Classification_temp, NULL)]
  }
  oncomat = data.table::dcast(data = subMaf[, .(Hugo_Symbol, 
                                                Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ 
                                Tumor_Sample_Barcode, fun.aggregate = function(x, cnv = cnv_events) {
                                  x = as.character(x)
                                  xad = x[x %in% cnv]
                                  xvc = x[!x %in% cnv]
                                  if (length(xvc) > 0) {
                                    xvc = ifelse(test = length(xvc) > 1, yes = "Multi_Hit", 
                                                 no = xvc)
                                  }
                                  x = ifelse(test = length(xad) > 0, yes = paste(xad, 
                                                                                 xvc, sep = ";"), no = xvc)
                                  x = gsub(pattern = ";$", replacement = "", x = x)
                                  x = gsub(pattern = "^;", replacement = "", x = x)
                                  return(x)
                                }, value.var = "Variant_Classification", fill = "", drop = FALSE)
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[, -1, drop = FALSE])
  variant.classes = as.character(unique(subMaf[, Variant_Classification]))
  variant.classes = c("", variant.classes, "Multi_Hit")
  names(variant.classes) = 0:(length(variant.classes) - 1)
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)]) + 
                                     1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  oncomat.copy <- oncomat
  for (i in 1:length(variant.classes2)) {
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  if (nrow(oncomat) == 1) {
    mdf = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, 
                vc = variant.classes))
  }
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  if (ncol(mdf) == 1) {
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE), ])
    colnames(mdf) = sampleId
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf), 
    ])
    colnames(oncomat.copy) = sampleId
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, 
                vc = variant.classes))
  }
  else {
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), 
    ]
    nMut = mdf[, ncol(mdf)]
    mdf = mdf[, -ncol(mdf)]
    mdf.temp.copy = mdf
    mdf[mdf != 0] = 1
    tmdf = t(mdf)
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), 
                                  decreasing = TRUE)), ])
    mdf.temp.copy = mdf.temp.copy[rownames(mdf), ]
    mdf.temp.copy = mdf.temp.copy[, colnames(mdf)]
    mdf = mdf.temp.copy
    oncomat.copy <- oncomat.copy[, colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf), ]
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, 
                vc = variant.classes, cnvc = cnv_events))
  }
}
get_vcColors <- function (alpha = 1, websafe = FALSE, named = TRUE) 
{
  if (websafe) {
    col = c("#F44336", "#E91E63", "#9C27B0", "#673AB7", 
            "#3F51B5", "#2196F3", "#03A9F4", "#00BCD4", "#009688", 
            "#4CAF50", "#8BC34A", "#CDDC39", "#FFEB3B", "#FFC107", 
            "#FF9800", "#FF5722", "#795548", "#9E9E9E", "#607D8B")
  }
  else {
    col = c(RColorBrewer::brewer.pal(11, name = "Paired"), 
            RColorBrewer::brewer.pal(11, name = "Spectral")[1:3], 
            "black", "violet", "royalblue", "#7b7060", "#535c68")
    col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  }
  if (named) {
    names(col) = names = c("Nonstop_Mutation", "Frame_Shift_Del", 
                           "IGR", "Missense_Mutation", "Silent", "Nonsense_Mutation", 
                           "RNA", "Splice_Site", "Intron", "Frame_Shift_Ins", 
                           "In_Frame_Del", "ITD", "In_Frame_Ins", "Translation_Start_Site", 
                           "Multi_Hit", "Amp", "Del", "Complex_Event", "pathway")
  }
  col
}
update_vc_codes <- function (om_op) 
{
  uniq_vc = as.character(unique(unlist(as.numeric(unlist(apply(om_op$numericMatrix, 
                                                               2, unique))))))
  missing_vc = uniq_vc[!uniq_vc %in% names(om_op$vc)]
  temp_names = names(om_op$vc)
  om_op$vc = c(om_op$vc, rep("Complex_Event", length(missing_vc)))
  names(om_op$vc) = c(temp_names, rep(missing_vc, length(missing_vc)))
  om_op$vc
}
update_colors <- function (x, y) 
{
  x = as.character(x)
  avail_colors = as.character(y[!names(y) %in% x])
  missing_entries = as.character(x)[!as.character(x) %in% 
                                      names(y)]
  missing_entries = missing_entries[!missing_entries %in% 
                                      ""]
  if (length(missing_entries) > 0) {
    if (length(missing_entries) > length(avail_colors)) {
      avail_colors = sample(x = colors(distinct = TRUE), 
                            size = length(missing_entries), replace = FALSE)
      names(avail_colors) = missing_entries
      y = c(y, avail_colors)
    }
    else {
      avail_colors = avail_colors[1:length(missing_entries)]
      names(avail_colors) = missing_entries
      y = c(y, avail_colors)
    }
  }
  y
}
plot_layout <- function (clinicalFeatures = NULL, drawRowBar = TRUE, drawColBar = TRUE, 
                         draw_titv = FALSE, exprsTbl = NULL, legend_height = 4, anno_height = 1) 
{
  if (is.null(clinicalFeatures)) {
    if (draw_titv) {
      if (!drawRowBar & !drawColBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, rep(0, 3)), 
                          nrow = 3, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          4, legend_height), widths = c(4, 0.5))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 5), 
                          nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          4, legend_height), widths = c(1, 4))
        }
      }
      else if (!drawRowBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, rep(0, 
                                                   4)), nrow = 4, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 
                                                          12, 4, legend_height), widths = c(4, 0.5))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 
                                                          12, 4, legend_height), widths = c(1, 4))
        }
      }
      else if (!drawColBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 5), 
                          nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          4, 4), widths = c(4, 1))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 7, 7), nrow = 3, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          4, legend_height), widths = c(1, 4, 1))
        }
      }
      else {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(4, 
                                                         1), heights = c(4, 12, 4, legend_height))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 8, 9, 10, 10, 10), nrow = 4, ncol = 3, 
                          byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(1, 
                                                         4, 1), heights = c(4, 12, 4, legend_height))
        }
      }
    }
    else {
      if (!drawRowBar & !drawColBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, rep(0, 2)), 
                          nrow = 2, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          legend_height), widths = c(4, 0.5))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 3), nrow = 2, 
                          ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          legend_height), widths = c(1, 4))
        }
      }
      else if (!drawRowBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, rep(0, 3)), 
                          nrow = 3, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 
                                                          12, legend_height), widths = c(4, 0.5))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 5), 
                          nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 
                                                          12, legend_height), widths = c(1, 4))
        }
      }
      else if (!drawColBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 3), nrow = 2, 
                          ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          legend_height), widths = c(4, 1))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 4, 4), 
                          nrow = 2, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          legend_height), widths = c(1, 4, 1))
        }
      }
      else {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 5), 
                          nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(4, 
                                                         1), heights = c(4, 12, legend_height))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 7, 7), nrow = 3, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(1, 
                                                         4, 1), heights = c(4, 12, legend_height))
        }
      }
    }
  }
  else {
    if (draw_titv) {
      if (!drawRowBar & !drawColBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, rep(0, 
                                                   4)), nrow = 4, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          anno_height, 4, legend_height), widths = c(4, 
                                                                                                     0.5))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          anno_height, 4, legend_height), widths = c(1, 
                                                                                                     4))
        }
      }
      else if (!drawRowBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, rep(0, 
                                                      5)), nrow = 5, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 
                                                          12, anno_height, 4, legend_height), widths = c(4, 
                                                                                                         0.5))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 8, 9, 9), nrow = 5, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 
                                                          12, anno_height, 4, legend_height), widths = c(1, 
                                                                                                         4))
        }
      }
      else if (!drawColBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          anno_height, 4, legend_height), widths = c(4, 
                                                                                                     1))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 8, 9, 10, 10, 10), nrow = 4, ncol = 3, 
                          byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          anno_height, 4, legend_height), widths = c(1, 
                                                                                                     4, 1))
        }
      }
      else {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 8, 9, 9), nrow = 5, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(4, 
                                                         1), heights = c(4, 12, anno_height, 4, legend_height))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 8, 9, 10, 11, 12, 13, 13, 13), nrow = 5, 
                          ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(1, 
                                                         4, 1), heights = c(4, 12, anno_height, 4, 
                                                                            legend_height))
        }
      }
    }
    else {
      if (!drawRowBar & !drawColBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, rep(0, 3)), 
                          nrow = 3, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          anno_height, legend_height), widths = c(4, 
                                                                                                  0.5))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 5), 
                          nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          anno_height, legend_height), widths = c(1, 
                                                                                                  4))
        }
      }
      else if (!drawRowBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, rep(0, 
                                                   4)), nrow = 4, ncol = 2, byrow = FALSE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 
                                                          12, anno_height, legend_height), widths = c(4, 
                                                                                                      0.5))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(4, 
                                                          12, anno_height, legend_height), widths = c(1, 
                                                                                                      4))
        }
      }
      else if (!drawColBar) {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 5), 
                          nrow = 3, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          anno_height, legend_height), widths = c(4, 
                                                                                                  1))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 7, 7), nrow = 3, ncol = 3, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, heights = c(12, 
                                                          anno_height, legend_height), widths = c(1, 
                                                                                                  4, 1))
        }
      }
      else {
        if (is.null(exprsTbl)) {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(4, 
                                                         1), heights = c(4, 12, anno_height, legend_height))
        }
        else {
          mat_lo = matrix(data = c(1, 2, 3, 4, 5, 6, 
                                   7, 8, 9, 10, 10, 10), nrow = 4, ncol = 3, 
                          byrow = TRUE)
          lo = graphics::layout(mat = mat_lo, widths = c(1, 
                                                         4, 1), heights = c(4, 12, anno_height, legend_height))
        }
      }
    }
  }
  lo
}

# Clinical information arrangement ------------------------------------------------------------------

data_clinical <- read.csv(
  '/mnt/data3/jiaoxi/新辅助治疗队列/vcf/meta_data/clinical_data.csv',
  sep = ',',row.names = NULL
) %>% 
  rename_all(~c('ID','PPN')) 
data_clinical_pcr <- read.csv(
  '/mnt/data3/jiaoxi/新辅助治疗队列/vcf/meta_data/clinical_data_pCR.csv',
  sep = ',',row.names = NULL) %>% 
  select("ID","sample_id") %>% 
  left_join(data_clinical,by = 'ID') %>% 
  filter(!is.na(PPN)) %>% 
  rename(Tumor_Sample_Barcode = sample_id) %>% 
  select(Tumor_Sample_Barcode,ID,PPN)
colnames(data_clinical_pcr)
# Complement related gene screening ----------------------------------------------------------------

all_pathway <- readRDS('/mnt/data3/jiaoxi/新辅助治疗队列/vcf/all_pathway.rds')
genes_complement <- all_pathway[['GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION']] %>% unique()
# var_maf= annovarToMaf(annovar = "all_sample.txt", 
#                       Center = 'NA', 
#                       refBuild = 'hg19', 
#                       tsbCol = 'Tumor_Sample_Barcode', 
#                       table = 'ensGene',MAFobj =FALSE,
#                       sep = "\t")
var_maf <-   annovarToMaf(
  annovar = "all_sample.txt",
  # annovar = 'all_sample_unfilter.txt',
  Center = 'CSI-NUS',
  refBuild = 'hg19',
  tsbCol = 'Tumor_Sample_Barcode',
  table = 'ensGene'
) 
# colnames(tmp)
colnames(var_maf)
var_maf$AAChange.refGene %>% unique()
var_maf <- var_maf %>% 
  # filter(
  #   AAChange.refGene %in% c('exonic','exonic;splicing')
  # ) %>% 
  select(
    Chromosome,Start_Position,End_Position,
    Reference_Allele,Tumor_Seq_Allele2,Gene.refGene,
    GeneDetail.refGene,ExonicFunc.refGene,Func.refGene,
    Tumor_Sample_Barcode,AAChange.refGene
  ) %>% 
  rename_all(~c(
    "Chr","Start","End", 
    "Ref","Alt","Gene.ensGene",
    "GeneDetail.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene",
    "Tumor_Sample_Barcode","Func.ensGene" 
  )) %>% #筛选使用到的样本
  filter(
    Tumor_Sample_Barcode %in% Sample_Select
  )
# 将转换后的 MAF 数据保存为临时文件
temp_maf_file <- tempfile(fileext = ".txt")
write.table(var_maf, file = temp_maf_file, sep = "\t", quote = FALSE, row.names = FALSE)
var_maf <-   annovarToMaf(
  annovar = temp_maf_file,
  Center = 'CSI-NUS',
  refBuild = 'hg19',
  tsbCol = 'Tumor_Sample_Barcode',
  table = 'ensGene'
) %>% 
  mutate(
    Hugo_Symbol = Gene.refGene
  )
temp_maf_file <- tempfile(fileext = ".maf")
write.table(var_maf, file = temp_maf_file, sep = "\t", quote = FALSE, row.names = FALSE)
# Read temporary files and convert them to MAF objects for maftools
maf <- read.maf(maf = temp_maf_file)
# Obtain mutation statistics of genes ====
gene_summary <- getGeneSummary(maf)
gene_summary %>% colnames()
gene_Mut <- gene_summary %>% 
  filter(total>0.01) %>% 
  pull(Hugo_Symbol)
genes_complement_Mut <- gene_Mut[gene_Mut %in% genes_complement]
genes_complement_Mut
# plotmafSummary(maf = maf, addStat = 'median')
plotmafSummary(
  maf = maf,
  rmOutlier = TRUE,
  addStat = 'median',
  dashboard = TRUE,showBarcodes = TRUE,
  titvRaw = FALSE,top = 15
)#Draw a summary file of the maf
oncoplot(
  maf = maf,
  genes = c("C8B","CPN1","CD5L","CPB2","CFP","VTN"),
  draw_titv = FALSE,
  drawRowBar = FALSE,
  fontSize = 2,
  legendFontSize = 3,
)
png("/mnt/data3/jiaoxi/新辅助治疗队列/result_figs_0814/wes_maftools.png", ,width = 260,height = 600)
oncoplot_self(
  maf = maf,
  genes = c("C8B","CPN1","CD5L","CPB2","CFP","VTN"),
  draw_titv = FALSE,
  drawRowBar = FALSE,
  fontSize = 1.3,
  bgCol  = "grey90",
  # sepwd_genes = 3,
  # sepwd_samples = 3,
  legendFontSize = 2.5
)
dev.off()
pdf("/mnt/data3/jiaoxi/新辅助治疗队列/result_figs_0814/wes_maftools.pdf", ,width = 2.6,height = 6)
oncoplot_self(
  maf = maf,
  genes = c("C8B","CPN1","CD5L","CPB2","CFP","VTN"),
  draw_titv = FALSE,
  drawRowBar = FALSE,
  fontSize = 1.3,
  bgCol  = "grey90",
  # sepwd_genes = 3,
  # sepwd_samples = 3,
  legendFontSize = 2.5
)
dev.off()
# Screening mutant samples ------------------------------------------------------------------


sample_maf <- maf@maf.silent %>% 
  as.data.frame() %>% 
  select(
    c(
      # 'Hugo_Symbol',
      'Tumor_Sample_Barcode',
      # 'Variant_Classification',
      # 'Variant_Type'
    )
  ) %>% 
  distinct_all() %>% 
  pull(Tumor_Sample_Barcode) %>% 
  str_remove('/sample/') %>% unique()
sample_maf
meta_sample <- data_clinical_pcr$PPN %>% unique()
maf_use <- p %>% str_remove('/sample/') %>% unique()
`%noin%` <- Negate(`%in%`)
# maf_use[maf_use %noin% tmp_1]
maf_use[maf_use %noin% meta_sample]
meta_sample[meta_sample %noin% sample_maf]
sample_maf[sample_maf %noin% meta_sample]
sample_maf[sample_maf %noin% maf_use]
data_meta <- data_clinical_pcr %>% 
  # filter(
  #   PPN %in% sample_maf
  # ) %>% 
  mutate(
    group_maf = ifelse(
      test = PPN %in% (
        p %>% str_remove('/sample/') %>% unique()
      ),
      yes = 'mutation',
      no = 'normal'
    ),
    group_maf = ifelse(
      test = PPN %in% sample_maf,
      yes = group_maf,
      no = 'without wes'
    )
  )
write.csv(data_meta,'/mnt/data3/jiaoxi/新辅助治疗队列/vcf/data_meta.csv')
# endline -----------------------------------------------------------------


