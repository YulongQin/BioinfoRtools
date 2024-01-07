### ---------------
###
### Create: qyl
### Date: 2023-08-15 14:58:00
### Email: 524583919@qq.com
### Function：
### - Using DESeq2 for DEGs analysis
### -
###
### Update Log: 2023-08-15
###
### ---------------

#### ----- Function ------ ####
#' Using DESeq2 for DEGs analysis
#'
#' @title Using DESeq2 for DEGs analysis
#' @description
#' This function uses DESeq2 package to identify differentially expressed genes (DEGs)
#' between two groups of samples. It takes in gene count data and returns a list of
#' DEGs with adjusted p-values and log2 fold changes.
#' @param countdata expression matrix
#' @param group the group of samples
#' @param case the case of samples
#' @param ctrl the control of samples
#' @param covariate the covariate of samples
#' @param updown_thres the threshold of updown
#' @param rate_sap the rate of sap
#' @param OUT_df_index the index of output df
#' @param updown_index the index of updown
#' @param p_value the p value
#' @param padj_value the padj value
#' @param grp_nm the name of group
#' @param dir_nm the name of dir
#'
#' @return data.frame
#' @author Yulong Qin
#' @seealso \code{\link{DEGs_Wilcoxon}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/DEGs_DESeq2.R")
#' df <- DEGs_DESeq2(countdata, group,
#'   case = "case", ctrl = "ctrl", rate_sap = 0.1,
#'   OUT_df_index = T, updown_index = T
#' )
#' }
#' @importFrom utils read.table  write.table
#' @importFrom stats na.omit as.formula model.matrix lm p.adjust residuals var
#' @importFrom grDevices png pdf dev.off colorRampPalette
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results counts
#' @importFrom dplyr select mutate filter arrange group_by summarise slice
#' @export
#'
DEGs_DESeq2 <- function(countdata, group, case = "case", ctrl = "ctrl", covariate = NULL,
                        updown_thres = 1, rate_sap = 0.1,
                        OUT_df_index = FALSE, updown_index = FALSE,
                        p_value = 0.05, padj_value = 0.05,
                        grp_nm = "DEGs_DESeq2", dir_nm = "DEGs_DESeq2") {
  ### 1.library ####
  # suppressMessages({
  #   library(DESeq2)
  #   library(tidyverse)
  # })
  output_dir <- paste0("./outputdata/", dir_nm, "/", grp_nm)
  # photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm)
  dir.create(output_dir, recursive = T)
  # dir.create(photo_dir,recursive = T)

  ### 2.countdata ####
  countdata <- countdata %>%
    round(digits = 0) %>% # ！！！ 取整，Deseq2才要
    na.omit()
  countdata <- countdata[rowSums(countdata) > 0, ] # 删除全0行
  countdata <- countdata[rowSums(countdata > 0) >= rate_sap * ncol(countdata), ]

  ### 3.coldata/clinic ####
  coldata <- data.frame(
    row.names = colnames(countdata), # countdata的列名==coldata的行名
    condition = factor(group, levels = c(case, ctrl))
  )
  if (!is.null(covariate)) {
    coldata <- cbind(covariate, coldata) # 关注分组放在最后一列
  }
  coldata <- coldata[, apply(coldata, 2, function(x) length(unique(x))) > 1, drop = F] # 去除只有一种值的列
  if (ncol(coldata) == 0) {
    stop("coldata has no columns,please check your coldata!")
  }
  all(rownames(coldata) == colnames(countdata))

  ### 4.DESeq2 ####
  str1 <- paste0("~", paste(colnames(coldata), collapse = "+"))
  str1 <- as.formula(str1)
  # deign_matrix <- model.matrix(str1,data=coldata) # 不需要自己构建设计矩阵
  dds <- DESeqDataSetFromMatrix(
    countData = countdata,
    colData = coldata,
    design = str1
  )
  dds <- dds[rowSums(counts(dds)) > 1, ] # 二次过滤，其实作用和前面rate_sap*ncol(countdata)一样
  dds <- DESeq(dds)

  ### 5.results ####
  res <- results(dds, contrast = c("condition", case, ctrl))
  # res <- results(dds,contrast = list(c(paste("condition",case,"vs",ctrl,sep="_"))))
  DEGs_df <- res %>%
    as.data.frame() %>%
    mutate(gene = rownames(.), .before = 1) %>%
    mutate(group = ifelse(pvalue < p_value & padj < padj_value & abs(log2FoldChange) >= updown_thres,
      ifelse(log2FoldChange >= updown_thres, "Up", "Down"), "No-Change"
    )) %>%
    mutate(FC_4 = ifelse(pvalue < p_value & padj < padj_value &
      abs(log2FoldChange) >= 2, gene, ""))
  if (OUT_df_index) {
    OUT_df <- DEGs_df[c(1, 3, 6, 7, 8, 9)]
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

  write.table(OUT2_df,
    file = paste0(output_dir, "/DEGs_", grp_nm, ".txt"),
    row.names = T, col.names = NA, quote = F, sep = "\t"
  )
  invisible(return(OUT2_df))
}

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/DEGs_DESeq2.R")
  df <- DEGs_DESeq2(countdata, group,
    case = "case", ctrl = "ctrl", rate_sap = 0.1,
    OUT_df_index = FALSE, updown_index = FALSE, updown_thres = 1,
    p_value = 0.05, padj_value = 0.05,
    grp_nm = "DEGs_DESeq2", dir_nm = "DEGs_DESeq2"
  )
}
