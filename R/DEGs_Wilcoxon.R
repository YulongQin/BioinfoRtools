### ---------------
###
### Create: qyl
### Date: 2023-11-12 14:24:00
### Email: 524583919@qq.com
### Function：
### - Using Wilcoxon for DEGs analysis
### -
###
### Update Log: 2023-11-12
###
### ---------------

#### ----- Function ------ ####
#' Using Wilcoxon for DEGs analysis
#'
#' @title Using Wilcoxon for DEGs analysis
#' @description This function uses Wilcoxon to identify differentially expressed genes (DEGs)
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
#' @seealso \code{\link{DEGs_DESeq2}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/DEGs_Wilcoxon.R")
#' df <- DEGs_Wilcoxon(countdata, group,
#'   case = "case", ctrl = "ctrl", rate_sap = 0.1,
#'   OUT_df_index = T, updown_index = T
#' )
#' }
#' @importFrom utils write.table
#' @importFrom stats lm p.adjust residuals wilcox.test
#' @importFrom dplyr select mutate filter arrange group_by summarise slice
#' @export
#'
DEGs_Wilcoxon <- function(countdata, group, case = "case", ctrl = "ctrl", covariate = NULL,
                          updown_thres = 1, rate_sap = 0.3,
                          OUT_df_index = FALSE, updown_index = FALSE,
                          p_value = 0.05, padj_value = 0.05,
                          grp_nm = "DEGs_Wilcoxon", dir_nm = "DEGs_Wilcoxon") {
  ### 1.library ####
  # suppressMessages({
  #   library(tidyverse)
  # })
  output_dir <- paste0("./outputdata/", dir_nm, "/", grp_nm)
  # photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm)
  dir.create(output_dir, recursive = T)
  # dir.create(photo_dir,recursive = T)

  ### 2.countdata ####
  countdata <- countdata %>%
    # round(digits = 0) %>% # ！！！ 取整，Deseq2才要
    na.omit()
  countdata <- countdata[rowSums(countdata) > 0, ] # 删除全0行
  countdata <- countdata[rowSums(countdata > 0) >= rate_sap * ncol(countdata), ]

  ### 3.coldata/clinic ####
  coldata <- data.frame(
    row.names = colnames(countdata), # countdata的列名==coldata的行名
    condition = factor(group, levels = c(case, ctrl))
  )
  all(rownames(coldata) == colnames(countdata))

  ### 4.covariate ####
  if (!is.null(covariate)) {
    print(paste0("Covariate correction start at ", Sys.time()))
    countdata <- t(countdata)
    combined_data <- cbind(countdata, covariate) %>%
      as.data.frame()
    str1 <- paste0("countdata ~", paste(colnames(covariate), collapse = "+"))
    str1 <- as.formula(str1)
    model <- lm(str1, # 拟合线性模型，使用协变量来矫正RNA表达矩阵
      data = combined_data
    )
    # 获取线性模型的残差，即矫正后的矩阵。相当于距离直线的距离，直线就是矫正了批次的基准线
    countdata <- residuals(model) %>%
      t() %>%
      as.data.frame()
    write.table(coldata,
      file = paste0(output_dir, "/coldata_lm_covariate_", grp_nm, ".txt"),
      row.names = T, col.names = NA, quote = F, sep = "\t"
    )
    write.table(countdata,
      file = paste0(output_dir, "/countdata_lm_covariate_", grp_nm, ".txt"),
      row.names = T, col.names = NA, quote = F, sep = "\t"
    )
  }

  ### 5.Wilcoxon ####
  print(paste0("Wilcoxon test start at ", Sys.time()))
  res_df <- as.data.frame(matrix(numeric(0), ncol = 5))
  for (i_num in 1:nrow(countdata)) {
    i_gene <- rownames(countdata)[i_num]
    i_df <- data.frame(
      gene = as.numeric(countdata[i_num, ]),
      group = group
    )
    case_vec <- i_df$gene[i_df$group == case]
    case_mean <- mean(case_vec)
    ctrl_vec <- i_df$gene[i_df$group == ctrl]
    ctrl_mean <- mean(ctrl_vec)
    log2FC <- log2((case_mean + 10^-6) / (ctrl_mean + 10^-6)) # 防止出现全0值

    res <- wilcox.test(
      x = case_vec,
      y = ctrl_vec,
      alternative = "two.sided",
      mu = 0, paired = F, exact = F, correct = T
    )
    res_df[i_num, 1] <- case_mean
    res_df[i_num, 2] <- ctrl_mean
    res_df[i_num, 3] <- log2FC
    res_df[i_num, 4] <- res$statistic
    res_df[i_num, 5] <- res$p.value
    rownames(res_df)[i_num] <- i_gene
  }
  colnames(res_df) <- c("Case_mean", "Ctrl_mean", "log2FC", "statistic", "pvalue")
  res_df$padj <- p.adjust(res_df$pvalue, method = "BH")
  print(paste0("Wilcoxon test end at ", Sys.time()))

  ### 5.results ####
  DEGs_df <- res_df %>%
    as.data.frame() %>%
    mutate(gene = rownames(.), .before = 1) %>%
    mutate(group = ifelse(pvalue < p_value & padj < padj_value & abs(log2FC) >= updown_thres,
      ifelse(log2FC >= updown_thres, "Up", "Down"), "No-Change"
    )) %>%
    mutate(FC_4 = ifelse(pvalue < p_value & padj < padj_value &
      abs(log2FC) >= 2, gene, ""))
  if (OUT_df_index) {
    OUT_df <- DEGs_df[c(1, 4, 6, 7, 8, 9)]
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
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/DEGs_Wilcoxon.R")
  df <- DEGs_Wilcoxon(countdata, group,
    case = "case", ctrl = "ctrl", covariate = NULL,
    updown_thres = 1, rate_sap = 0.3,
    OUT_df_index = FALSE, updown_index = FALSE,
    p_value = 0.05, padj_value = 0.05,
    grp_nm = "DEGs_Wilcoxon", dir_nm = "DEGs_Wilcoxon"
  )
}
