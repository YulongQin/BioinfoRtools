### ---------------
###
### Create: qyl
### Date: 2023-10-11 17:34:00
### Email: 524583919@qq.com
### Function：
### - GSVA and ssGSEA analysis
### -
###
### Update Log: 2023-10-11
###
### ---------------

#### ----- Function ------ ####
#' GSVA and ssGSEA analysis
#'
#' @title GSVA and ssGSEA analysis
#' @description GSVA and ssGSEA analysis
#' @param exp_df TPM matrix
#' @param gsva_mat_go_scale the result of GSVA analysis, gsva_mat_go_scale
#' @param gsva_mat_kegg_scale the result of GSVA analysis, gsva_mat_kegg_scale
#' @param ssgsea_mat_go_scale the result of ssGSEA analysis, ssgsea_mat_go_scale
#' @param ssgsea_mat_kegg_scale the result of ssGSEA analysis, ssgsea_mat_kegg_scale
#' @param group the group of samples
#' @param case the name of case
#' @param ctrl the name of ctrl
#' @param path_gmt_GO the path of GO gmt file
#' @param path_gmt_KEGG the path of KEGG gmt file
#' @param grp_nm the name of group
#' @param dir_nm the name of directory
#' @param OUT_df_index whether to output the index of the result
#' @param updown_index whether to output the index of up and down
#' @param updown_thres the threshold of up and down
#' @param p_value the threshold of p value
#' @param padj_value the threshold of padj value
#'
#' @return NULL
#' @author Yulong Qin
#' @seealso \code{\link{GO_KEGG}},\code{\link{GSEA_GO_KEGG}}, \code{\link{GSVA_ssGSEA}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/GSVA_ssGSEA.R")
#' GSVA_ssGSEA(
#'   exp_df = exp_df,
#'   group = group,
#'   case = "Case", ctrl = "Ctrl",
#'   grp_nm = "GSVA_ssGSEA", dir_nm = "GSVA_ssGSEA",
#'   OUT_df_index = FALSE, updown_index = FALSE, updown_thres = 1,
#'   p_value = 0.05, padj_value = 0.05
#' )
#' }
#' @import GSVA
#' @import GSEABase
#' @import clusterProfiler
#' @importFrom ggplotify as.ggplot
#' @import ggplot2
#' @import limma
#' @importFrom pheatmap pheatmap
#' @import patchwork
#' @importFrom magrittr add
#' @importFrom dplyr select mutate filter arrange group_by summarise slice
#' @importFrom tibble column_to_rownames
#' @export GSVA_ssGSEA
GSVA_ssGSEA <- function(exp_df = NULL, # 输入的是原始的TPM矩阵，
                        gsva_mat_go_scale = NULL,
                        gsva_mat_kegg_scale = NULL,
                        ssgsea_mat_go_scale = NULL,
                        ssgsea_mat_kegg_scale = NULL,
                        group = NULL,
                        case = "Case", ctrl = "Ctrl",
                        path_gmt_GO = "./inputdata/c5.go.v2023.1.Hs.symbols_BP.gmt",
                        path_gmt_KEGG = "./inputdata/c5.go.v2023.1.Hs.symbols_BP.gmt",
                        grp_nm = "GSVA_ssGSEA", dir_nm = "GSVA_ssGSEA",
                        OUT_df_index = FALSE, updown_index = FALSE, updown_thres = 1,
                        p_value = 0.05, padj_value = 0.05) {
  ### 1.library ####
  cat(paste0("1.library start at ", Sys.time(), "\n"))
  # suppressMessages({
  #   library(limma)
  #   library(magrittr)
  #   library(ggplotify)
  #   library(tidyverse)
  #   library(clusterProfiler)
  #   library(GSVA)
  #   library(GSEABase) # getGmt,与read.gmt不同
  #   library(pheatmap)
  #   library(patchwork)
  #   # library(BiocParallel)
  # })

  output_dir <- paste0("./outputdata/", dir_nm, "/", grp_nm)
  photo_dir <- paste0("./photo/", dir_nm, "/", grp_nm)
  dir.create(output_dir, recursive = T)
  dir.create(photo_dir, recursive = T)
  patchwork_index <- F

  if (!is.null(gsva_mat_go_scale) | !is.null(gsva_mat_kegg_scale) |
    !is.null(ssgsea_mat_go_scale) | !is.null(ssgsea_mat_kegg_scale)) {
    exp_df <- NULL
  }
  if (!is.null(exp_df) & is.null(gsva_mat_go_scale) & is.null(gsva_mat_kegg_scale) &
    is.null(ssgsea_mat_go_scale) & is.null(ssgsea_mat_kegg_scale)) {
    ### 2.exp and gmt ####
    cat(paste0("2.exp and gmt start at ", Sys.time(), "\n"))
    exp_df <- exp_df %>%
      na.omit() %>%
      filter(rowSums(.) > 0) %>%
      add(1) %>%
      log2() %>%
      as.matrix()
    go_list <- getGmt(path_gmt_GO)
    kegg_list <- getGmt(path_gmt_KEGG)

    ### 3.GSVA ####
    ## 3.1 GO_BP
    cat(paste0("3.1 GSVA GO_BP start at ", Sys.time(), "\n"))
    gsva_mat_go <- gsva(
      expr = exp_df, # 行为基因，列为样本
      gset.idx.list = go_list,
      method = "gsva",
      kcdf = "Gaussian",
      verbose = T,
      parallel.sz = 16 # 调用所有核
    )
    # saveRDS(gsva_mat_go,file = "./outputdata/gsva_mat_go.rds")
    write.table(gsva_mat_go, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/GSVA_score_GO_BP_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )

    gsva_mat_go_scale <- gsva_mat_go %>%
      t() %>%
      scale(center = T, scale = T) %>%
      t() %>%
      as.data.frame()
    write.table(gsva_mat_go_scale, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/GSVA_score_scale_GO_BP_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )

    ## 3.2 KEGG
    cat(paste0("3.2 GSVA KEGG start at ", Sys.time(), "\n"))
    gsva_mat_kegg <- gsva(
      expr = exp_df, # 行为基因，列为样本
      gset.idx.list = kegg_list,
      method = "gsva",
      kcdf = "Gaussian",
      verbose = T,
      parallel.sz = 16 # 调用所有核
    )
    # saveRDS(gsva_mat_kegg,file = "./outputdata/gsva_mat_kegg.rds")
    write.table(gsva_mat_kegg, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/GSVA_score_KEGG_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )

    gsva_mat_kegg_scale <- gsva_mat_kegg %>%
      t() %>%
      scale(center = T, scale = T) %>%
      t() %>%
      as.data.frame()
    write.table(gsva_mat_kegg_scale, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/GSVA_score_scale_KEGG_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )

    ### 4.ssGSEA ####
    ## 4.1 GO_BP
    cat(paste0("4.1 ssGSEA GO_BP start at ", Sys.time(), "\n"))
    ssgsea_mat_go <- gsva(
      expr = exp_df, # 行为基因，列为样本
      gset.idx.list = go_list,
      method = "ssgsea", #
      kcdf = "Gaussian",
      verbose = T,
      parallel.sz = 16 # 调用所有核
    )
    # saveRDS(gsva_mat_go,file = "./outputdata/gsva_mat_go.rds")
    write.table(ssgsea_mat_go, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/ssGSEA_score_GO_BP_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )

    ssgsea_mat_go_scale <- ssgsea_mat_go %>%
      t() %>%
      scale(center = T, scale = T) %>%
      t() %>%
      as.data.frame()
    write.table(ssgsea_mat_go_scale, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/ssGSEA_score_scale_GO_BP_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )

    ## 4.2 KEGG
    cat(paste0("4.1 ssGSEA KEGG start at ", Sys.time(), "\n"))
    ssgsea_mat_kegg <- gsva(
      expr = exp_df, # 行为基因，列为样本
      gset.idx.list = kegg_list,
      method = "ssgsea",
      kcdf = "Gaussian",
      verbose = T,
      parallel.sz = 16 # 调用所有核
    )
    # saveRDS(gsva_mat_kegg,file = "./outputdata/gsva_mat_kegg.rds")
    write.table(ssgsea_mat_kegg, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/ssGSEA_score_KEGG_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )

    ssgsea_mat_kegg_scale <- ssgsea_mat_kegg %>%
      t() %>%
      scale(center = T, scale = T) %>%
      t() %>%
      as.data.frame()
    write.table(ssgsea_mat_kegg_scale, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/ssGSAE_score_scale_KEGG_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )
  }

  ### 5.limma_DEGs ####
  ## 5.1 GSVA_GO_BP
  if (!is.null(exp_df) | !is.null(gsva_mat_go_scale) | !is.null(group)) {
    cat(paste0("5.1 limma GSVA_GO_BP start at ", Sys.time(), "\n"))
    group <- group %>%
      factor(levels = c(ctrl, case))
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    rownames(design) <- colnames(gsva_mat_go)
    vs <- paste0(case, "-", ctrl)
    contrast_matrix <- makeContrasts(
      contrasts = vs,
      levels = design
    )
    fit1 <- lmFit(gsva_mat_go, design)
    fit2 <- contrasts.fit( # 两种方法的区别：更改了设计矩阵fit，这影响下面的coef作用于fit1还是fit2
      fit = fit1,
      contrasts = contrast_matrix,
      coefficients = NULL
    )
    fit2 <- eBayes(fit2)
    DEGs_df <- topTable(
      fit2,
      coef = 1, # 只看一种则指定contrasts中哪种组合数字,若只有一种NULL等价于1
      number = Inf
    )
    DEGs_df <- DEGs_df %>%
      as.data.frame() %>%
      mutate(gene = rownames(.), .before = 1) %>%
      mutate(group = ifelse(P.Value < p_value & adj.P.Val < padj_value & abs(logFC) >= updown_thres,
        ifelse(logFC >= updown_thres, "Up", "Down"), "No-Change"
      )) %>%
      mutate(FC_4 = ifelse(P.Value < p_value & adj.P.Val < padj_value &
        abs(logFC) >= 2, gene, ""))
    if (OUT_df_index) {
      OUT_df <- DEGs_df[c(1, 2, 6, 8, 9)]
      colnames(OUT_df) <- c("gene", "log2FC", "p_adj", "group", "FC_4")
    } else {
      OUT_df <- DEGs_df
    }
    if (updown_index) {
      OUT2_df <- OUT_df %>% # 仅有updown
        filter(group == "Up" | group == "Down") %>%
        arrange(p_adj)
    } else {
      OUT2_df <- OUT_df
    }
    DEGs_gsva_go <- OUT2_df
    write.table(DEGs_gsva_go, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/DEGs_GSVA_GO_BP_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )
  }

  ## 5.2 GSVA_KEGG
  if (!is.null(exp_df) | !is.null(gsva_mat_kegg_scale) | !is.null(group)) {
    cat(paste0("5.2 limma GSVA_KEGG start at ", Sys.time(), "\n"))
    group <- group %>%
      factor(levels = c(ctrl, case))
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    rownames(design) <- colnames(gsva_mat_kegg)
    vs <- paste0(case, "-", ctrl)
    contrast_matrix <- makeContrasts(
      contrasts = vs,
      levels = design
    )
    fit1 <- lmFit(gsva_mat_kegg, design)
    fit2 <- contrasts.fit( # 两种方法的区别：更改了设计矩阵fit，这影响下面的coef作用于fit1还是fit2
      fit = fit1,
      contrasts = contrast_matrix,
      coefficients = NULL
    )
    fit2 <- eBayes(fit2)
    DEGs_df <- topTable(
      fit2,
      coef = 1, # 只看一种则指定contrasts中哪种组合数字,若只有一种NULL等价于1
      number = Inf
    )
    DEGs_df <- DEGs_df %>%
      as.data.frame() %>%
      mutate(gene = rownames(.), .before = 1) %>%
      mutate(group = ifelse(P.Value < p_value & adj.P.Val < padj_value & abs(logFC) >= updown_thres,
        ifelse(logFC >= updown_thres, "Up", "Down"), "No-Change"
      )) %>%
      mutate(FC_4 = ifelse(P.Value < p_value & adj.P.Val < padj_value &
        abs(logFC) >= 2, gene, ""))
    if (OUT_df_index) {
      OUT_df <- DEGs_df[c(1, 2, 6, 8, 9)]
      colnames(OUT_df) <- c("gene", "log2FC", "p_adj", "group", "FC_4")
    } else {
      OUT_df <- DEGs_df
    }
    if (updown_index) {
      OUT2_df <- OUT_df %>% # 仅有updown
        filter(group == "Up" | group == "Down") %>%
        arrange(p_adj)
    } else {
      OUT2_df <- OUT_df
    }
    DEGs_gsva_kegg <- OUT2_df
    write.table(DEGs_gsva_kegg, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/DEGs_GSVA_KEGG_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )
  }

  ## 5.3 ssGSEA_GO_BP
  if (!is.null(exp_df) | !is.null(ssgsea_mat_go_scale) | !is.null(group)) {
    cat(paste0("5.3 limma ssGSEA_GO_BP start at ", Sys.time(), "\n"))
    group <- group %>%
      factor(levels = c(ctrl, case))
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    rownames(design) <- colnames(ssgsea_mat_go)
    vs <- paste0(case, "-", ctrl)
    contrast_matrix <- makeContrasts(
      contrasts = vs,
      levels = design
    )
    fit1 <- lmFit(ssgsea_mat_go, design)
    fit2 <- contrasts.fit( # 两种方法的区别：更改了设计矩阵fit，这影响下面的coef作用于fit1还是fit2
      fit = fit1,
      contrasts = contrast_matrix,
      coefficients = NULL
    )
    fit2 <- eBayes(fit2)
    DEGs_df <- topTable(
      fit2,
      coef = 1, # 只看一种则指定contrasts中哪种组合数字,若只有一种NULL等价于1
      number = Inf
    )
    DEGs_df <- DEGs_df %>%
      as.data.frame() %>%
      mutate(gene = rownames(.), .before = 1) %>%
      mutate(group = ifelse(P.Value < p_value & adj.P.Val < padj_value & abs(logFC) >= updown_thres,
        ifelse(logFC >= updown_thres, "Up", "Down"), "No-Change"
      )) %>%
      mutate(FC_4 = ifelse(P.Value < p_value & adj.P.Val < padj_value &
        abs(logFC) >= 2, gene, ""))
    if (OUT_df_index) {
      OUT_df <- DEGs_df[c(1, 2, 6, 8, 9)]
      colnames(OUT_df) <- c("gene", "log2FC", "p_adj", "group", "FC_4")
    } else {
      OUT_df <- DEGs_df
    }
    if (updown_index) {
      OUT2_df <- OUT_df %>% # 仅有updown
        filter(group == "Up" | group == "Down") %>%
        arrange(p_adj)
    } else {
      OUT2_df <- OUT_df
    }
    DEGs_ssgsea_go <- OUT2_df
    write.table(DEGs_ssgsea_go, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/DEGs_ssGSEA_GO_BP_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )
  }

  ## 5.4 ssGSAE_KEGG
  if (!is.null(exp_df) | !is.null(ssgsea_mat_kegg_scale) | !is.null(group)) {
    cat(paste0("5.4 limma ssGSEA_KEGG start at ", Sys.time(), "\n"))
    group <- group %>%
      factor(levels = c(ctrl, case))
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    rownames(design) <- colnames(ssgsea_mat_kegg)
    vs <- paste0(case, "-", ctrl)
    contrast_matrix <- makeContrasts(
      contrasts = vs,
      levels = design
    )
    fit1 <- lmFit(ssgsea_mat_kegg, design)
    fit2 <- contrasts.fit( # 两种方法的区别：更改了设计矩阵fit，这影响下面的coef作用于fit1还是fit2
      fit = fit1,
      contrasts = contrast_matrix,
      coefficients = NULL
    )
    fit2 <- eBayes(fit2)
    DEGs_df <- topTable(
      fit2,
      coef = 1, # 只看一种则指定contrasts中哪种组合数字,若只有一种NULL等价于1
      number = Inf
    )
    DEGs_df <- DEGs_df %>%
      as.data.frame() %>%
      mutate(gene = rownames(.), .before = 1) %>%
      mutate(group = ifelse(P.Value < p_value & adj.P.Val < padj_value & abs(logFC) >= updown_thres,
        ifelse(logFC >= updown_thres, "Up", "Down"), "No-Change"
      )) %>%
      mutate(FC_4 = ifelse(P.Value < p_value & adj.P.Val < padj_value &
        abs(logFC) >= 2, gene, ""))
    if (OUT_df_index) {
      OUT_df <- DEGs_df[c(1, 2, 5, 6, 8, 9)]
      colnames(OUT_df) <- c("gene", "log2FC", "pvalue", "p_adj", "group", "FC_4")
    } else {
      OUT_df <- DEGs_df
    }
    if (updown_index) {
      OUT2_df <- OUT_df %>% # 仅有updown
        filter(group == "Up" | group == "Down") %>%
        arrange(p_adj)
    } else {
      OUT2_df <- OUT_df
    }
    DEGs_ssgsea_kegg <- OUT2_df
    write.table(DEGs_ssgsea_kegg, # 行为通路，列为样本，值为此样本的此通路值
      file = paste0(output_dir, "/DEGs_ssGSEA_KEGG_", grp_nm, ".txt"),
      sep = "\t", quote = F, row.names = T, col.names = NA
    )
  }

  ### 6.pathway heatmap ####
  ## 6.1 GSVA_GO_BP
  if (!is.null(exp_df) | !is.null(gsva_mat_go_scale) | !is.null(group)) {
    cat(paste0("6.1 GSVA_GO_BP start at ", Sys.time(), "\n"))
    DEGs_gsva_go$Des_short <- sapply(DEGs_gsva_go$gene, shorten_names)
    DEGs_gsva_go <- DEGs_gsva_go %>%
      arrange(adj.P.Val) %>%
      mutate(Des_short_unique = ifelse(!duplicated(Des_short), Des_short, paste0(Des_short, 2))) %>%
      mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique))) %>%
      mutate(case_mean = as.numeric(apply(gsva_mat_go_scale[!!(group) == !!(case)], 1, mean))) %>%
      mutate(ctrl_mean = as.numeric(apply(gsva_mat_go_scale[!!(group) == !!(ctrl)], 1, mean)))

    hm_in <- DEGs_gsva_go[1:20, ] %>%
      data.frame(row.names = NULL) %>%
      column_to_rownames("Des_short_unique") %>%
      .[12:11]
    hm_in <- hm_in[apply(hm_in, 1, var) != 0, ] %>%
      na.omit()
    hm_ann_col <- data.frame(
      Group = c(case, ctrl), # grps[,1]
      row.names = c("case_mean", "ctrl_mean") # 行名为样本,与分组对应
    )
    color <- c("#0000CD", "#DC143C")
    names(color) <- c(ctrl, case)
    hm_ann_colors <- list(
      Group = alpha(color, 0.8) # 含有名字的向量
    )
    name <- paste0(grp_nm, " GSVA_GO_BP")
    star_df <- cbind("", DEGs_gsva_go$adj.P.Val[1:20]) %>%
      as.data.frame() %>%
      mutate(V2 = as.numeric(V2)) %>%
      mutate(V2 = ifelse(V2 >= 0.05, "", ifelse(V2 < 0.001, "***", ifelse(V2 < 0.01, "**", "*"))))

    heatmap1 <- pheatmap(hm_in,
      ##
      cluster_cols = F,
      cluster_rows = F,
      scale = "none",
      ##
      angle_col = 45,
      ##
      cellwidth = 12,
      color = colorRampPalette(c("#00468BFF", "white", "#ED0000FF"))(1000),
      ##
      border = F,
      border_color = NA,
      main = name,
      annotation_col = hm_ann_col,
      annotation_colors = hm_ann_colors,
      ##
      show_colnames = T,
      show_rownames = T,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 8,
      display_numbers = star_df,
      silent = T
    )
    p1 <- as.ggplot(heatmap1) +
      theme(plot.margin = margin(l = 100, r = 10)) # 调整边距
    pdf(
      file = paste0(photo_dir, "/HTmap_GSVA_GO_BP_", grp_nm, ".pdf"),
      width = 6, height = 5
    )
    print(p1)
    dev.off()
  }

  ## 6.2 GSVA_KEGG
  if (!is.null(exp_df) | !is.null(gsva_mat_kegg_scale) | !is.null(group)) {
    cat(paste0("6.2 GSVA_KEGG start at ", Sys.time(), "\n"))
    DEGs_gsva_kegg$Des_short <- sapply(DEGs_gsva_kegg$gene, shorten_names)
    DEGs_gsva_kegg <- DEGs_gsva_kegg %>%
      arrange(adj.P.Val) %>%
      mutate(Des_short_unique = ifelse(!duplicated(Des_short), Des_short, paste0(Des_short, 2))) %>%
      mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique))) %>%
      mutate(case_mean = as.numeric(apply(gsva_mat_kegg_scale[!!(group) == !!(case)], 1, mean))) %>%
      mutate(ctrl_mean = as.numeric(apply(gsva_mat_kegg_scale[!!(group) == !!(ctrl)], 1, mean)))

    hm_in <- DEGs_gsva_kegg[1:20, ] %>%
      data.frame(row.names = NULL) %>%
      column_to_rownames("Des_short_unique") %>%
      .[12:11]
    hm_in <- hm_in[apply(hm_in, 1, var) != 0, ] %>%
      na.omit()
    hm_ann_col <- data.frame(
      Group = c(case, ctrl), # grps[,1]
      row.names = c("case_mean", "ctrl_mean") # 行名为样本,与分组对应
    )
    color <- c("#0000CD", "#DC143C")
    names(color) <- c(ctrl, case)
    hm_ann_colors <- list(
      Group = alpha(color, 0.8) # 含有名字的向量
    )
    name <- paste0(grp_nm, " GSVA_KEGG")
    star_df <- cbind("", DEGs_gsva_kegg$adj.P.Val[1:20]) %>%
      as.data.frame() %>%
      mutate(V2 = as.numeric(V2)) %>%
      mutate(V2 = ifelse(V2 >= 0.05, "", ifelse(V2 < 0.001, "***", ifelse(V2 < 0.01, "**", "*"))))

    heatmap2 <- pheatmap(hm_in,
      ##
      cluster_cols = F,
      cluster_rows = F,
      scale = "none",
      ##
      angle_col = 45,
      ##
      cellwidth = 12,
      color = colorRampPalette(c("#00468BFF", "white", "#ED0000FF"))(1000),
      ##
      border = F,
      border_color = NA,
      main = name,
      annotation_col = hm_ann_col,
      annotation_colors = hm_ann_colors,
      ##
      show_colnames = T,
      show_rownames = T,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 8,
      display_numbers = star_df,
      silent = T
    )
    p2 <- as.ggplot(heatmap2) +
      theme(plot.margin = margin(l = 100, r = 10)) # 调整边距
    pdf(
      file = paste0(photo_dir, "/HTmap_GSVA_KEGG_", grp_nm, ".pdf"),
      width = 6, height = 5
    )
    print(p2)
    dev.off()
  }

  ## 6.3 ssGSEA_GO_BP
  if (!is.null(exp_df) | !is.null(ssgsea_mat_go_scale) | !is.null(group)) {
    cat(paste0("6.3 ssGSEA_GO_BP start at ", Sys.time(), "\n"))
    DEGs_ssgsea_go$Des_short <- sapply(DEGs_ssgsea_go$gene, shorten_names)
    DEGs_ssgsea_go <- DEGs_ssgsea_go %>%
      arrange(adj.P.Val) %>%
      mutate(Des_short_unique = ifelse(!duplicated(Des_short), Des_short, paste0(Des_short, 2))) %>%
      mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique))) %>%
      mutate(case_mean = as.numeric(apply(ssgsea_mat_go_scale[!!(group) == !!(case)], 1, mean))) %>%
      mutate(ctrl_mean = as.numeric(apply(ssgsea_mat_go_scale[!!(group) == !!(ctrl)], 1, mean)))

    hm_in <- DEGs_ssgsea_go[1:20, ] %>%
      data.frame(row.names = NULL) %>%
      column_to_rownames("Des_short_unique") %>%
      .[12:11]
    hm_in <- hm_in[apply(hm_in, 1, var) != 0, ] %>%
      na.omit()
    hm_ann_col <- data.frame(
      Group = c(case, ctrl), # grps[,1]
      row.names = c("case_mean", "ctrl_mean") # 行名为样本,与分组对应
    )
    color <- c("#0000CD", "#DC143C")
    names(color) <- c(ctrl, case)
    hm_ann_colors <- list(
      Group = alpha(color, 0.8) # 含有名字的向量
    )
    name <- paste0(grp_nm, " ssGSEA_GO_BP")
    star_df <- cbind("", DEGs_ssgsea_go$adj.P.Val[1:20]) %>%
      as.data.frame() %>%
      mutate(V2 = as.numeric(V2)) %>%
      mutate(V2 = ifelse(V2 >= 0.05, "", ifelse(V2 < 0.001, "***", ifelse(V2 < 0.01, "**", "*"))))

    heatmap3 <- pheatmap(hm_in,
      ##
      cluster_cols = F,
      cluster_rows = F,
      scale = "none",
      ##
      angle_col = 45,
      ##
      cellwidth = 12,
      color = colorRampPalette(c("#00468BFF", "white", "#ED0000FF"))(1000),
      ##
      border = F,
      border_color = NA,
      main = name,
      annotation_col = hm_ann_col,
      annotation_colors = hm_ann_colors,
      ##
      show_colnames = T,
      show_rownames = T,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 8,
      display_numbers = star_df,
      silent = T
    )
    p3 <- as.ggplot(heatmap3) +
      theme(plot.margin = margin(l = 100, r = 10)) # 调整边距
    pdf(
      file = paste0(photo_dir, "/HTmap_ssGSEA_GO_BP_", grp_nm, ".pdf"),
      width = 6, height = 5
    )
    print(p3)
    dev.off()
  }

  ## 6.4 ssGSEA_KEGG
  if (!is.null(exp_df) | !is.null(ssgsea_mat_kegg_scale) | !is.null(group)) {
    cat(paste0("6.4 ssGSEA_KEGG start at ", Sys.time(), "\n"))
    DEGs_ssgsea_kegg$Des_short <- sapply(DEGs_ssgsea_kegg$gene, shorten_names)
    DEGs_ssgsea_kegg <- DEGs_ssgsea_kegg %>%
      arrange(adj.P.Val) %>%
      mutate(Des_short_unique = ifelse(!duplicated(Des_short), Des_short, paste0(Des_short, 2))) %>%
      mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique))) %>%
      mutate(case_mean = as.numeric(apply(ssgsea_mat_kegg_scale[!!(group) == !!(case)], 1, mean))) %>%
      mutate(ctrl_mean = as.numeric(apply(ssgsea_mat_kegg_scale[!!(group) == !!(ctrl)], 1, mean)))

    hm_in <- DEGs_ssgsea_kegg[1:20, ] %>%
      data.frame(row.names = NULL) %>%
      column_to_rownames("Des_short_unique") %>%
      .[12:11]
    hm_in <- hm_in[apply(hm_in, 1, var) != 0, ] %>%
      na.omit()
    hm_ann_col <- data.frame(
      Group = c(case, ctrl), # grps[,1]
      row.names = c("case_mean", "ctrl_mean") # 行名为样本,与分组对应
    )
    color <- c("#0000CD", "#DC143C")
    names(color) <- c(ctrl, case)
    hm_ann_colors <- list(
      Group = alpha(color, 0.8) # 含有名字的向量
    )
    name <- paste0(grp_nm, " ssGSEA_KEGG")
    star_df <- cbind("", DEGs_ssgsea_kegg$adj.P.Val[1:20]) %>%
      as.data.frame() %>%
      mutate(V2 = as.numeric(V2)) %>%
      mutate(V2 = ifelse(V2 >= 0.05, "", ifelse(V2 < 0.001, "***", ifelse(V2 < 0.01, "**", "*"))))

    heatmap4 <- pheatmap(hm_in,
      ##
      cluster_cols = F,
      cluster_rows = F,
      scale = "none",
      ##
      angle_col = 45,
      ##
      cellwidth = 12,
      color = colorRampPalette(c("#00468BFF", "white", "#ED0000FF"))(1000),
      ##
      border = F,
      border_color = NA,
      main = name,
      annotation_col = hm_ann_col,
      annotation_colors = hm_ann_colors,
      ##
      show_colnames = T,
      show_rownames = T,
      fontsize = 10,
      fontsize_row = 8,
      fontsize_col = 8,
      display_numbers = star_df,
      silent = T
    )
    p4 <- as.ggplot(heatmap4) +
      theme(plot.margin = margin(l = 100, r = 10)) # 调整边距
    pdf(
      file = paste0(photo_dir, "/HTmap_ssGSEA_KEGG_", grp_nm, ".pdf"),
      width = 6, height = 5
    )
    print(p4)
    dev.off()
  }

  ### 7.patchwork ####
  if (exists("p1") | exists("p2") | exists("p3") | exists("p4")) {
    cat(paste0("7.patchwork start at ", Sys.time(), "\n"))
    pdf(
      file = paste0(photo_dir, "/GSVA_ssGSEA_GO_KEGG.pdf"),
      width = 12, height = 10, family = "serif"
    )
    print(
      p1 + p2 + p3 + p4 +
        plot_layout(nrow = 2, guides = "auto")
    )
    dev.off()
  }
}

# shorten_names <- function(x, n_word=6, n_char=40){ # 单词数>n_word,单词长度>n_char
#   if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char))
#   {
#     if (nchar(x) > n_char) x <- substr(x, 1, n_char)
#     x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
#                      collapse=" "), "...", sep="")
#     return(x)
#   }else{
#     return(x)
#   }
# }

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/GSVA_ssGSEA.R")
  GSVA_ssGSEA(
    exp_df = NULL, # 输入的是原始的TPM矩阵，
    gsva_mat_go_scale = NULL,
    gsva_mat_kegg_scale = NULL,
    ssgsea_mat_go_scale = NULL,
    ssgsea_mat_kegg_scale = NULL,
    group, case = "Case", ctrl = "Ctrl",
    path_gmt_GO = "./inputdata/c5.go.v2023.1.Hs.symbols_BP.gmt",
    path_gmt_KEGG = "./inputdata/c5.go.v2023.1.Hs.symbols_BP.gmt",
    grp_nm = "GSVA_ssGSEA", dir_nm = "GSVA_ssGSEA",
    OUT_df_index = FALSE, updown_index = FALSE, updown_thres = 1,
    p_value = 0.05, padj_value = 0.05
  )
}
