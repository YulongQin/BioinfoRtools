### ---------------
###
### Create: qyl
### Date: 2023-09-04 16:40:00
### Email: 524583919@qq.com
### Function：
### - GSEA analysis
### -
###
### Update Log: 2023-09-04
###
### ---------------

#### ----- Function ------ ####
#' GSEA analysis
#'
#' @title GSEA analysis
#' @description GSEA analysis
#' @param DEGs the input DEGs must be a df with only two columns of c('SYMBOL','LogFC')
#' @param simplify_index whether to simplify the results
#' @param maxGSSize the maximum number of genes in the gene set
#' @param grp_nm the name of the group
#' @param dir_nm the name of the directory
#'
#' @return NULL
#' @author Yulong Qin
#' @seealso \code{\link{GO_KEGG}},\code{\link{GSEA_GO_KEGG}}, \code{\link{GSVA_ssGSEA}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/GSEA_GO_KEGG.R")
#' GSEA_GO_KEGG(
#'   DEGs = NULL, simplify_index = F,
#'   grp_nm = "GSEA_GO_KEGG", dir_nm = "GSEA_GO_KEGG"
#' )
#' }
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import ggplot2
#' @import DOSE
#' @import patchwork
#' @export GSEA_GO_KEGG
GSEA_GO_KEGG <- function(DEGs = NULL, # 输入的DEGs必须是经过处理的仅有c('SYMBOL','LogFC')两列的df
                         simplify_index = F, maxGSSize = 500,
                         grp_nm = "GSEA_GO_KEGG", dir_nm = "GSEA_GO_KEGG") {
  ### 1.library ####
  # suppressMessages({
  #   library(tidyverse)
  #   library(magrittr)
  #   library(clusterProfiler) # enrichGo/enrichKEGG/bitr
  #   library(org.Hs.eg.db) # 人类gene ID
  #   library(DOSE) # setReadable
  #   library(patchwork)
  # })
  output_dir <- paste0("./outputdata/", dir_nm, "/", grp_nm)
  photo_dir <- paste0("./photo/", dir_nm, "/", grp_nm)
  dir.create(output_dir, recursive = T)
  dir.create(photo_dir, recursive = T)
  # source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/shorten_names.R")
  patchwork_index <- F

  ### 2.LogFC and gmt
  if (ncol(DEGs) > 2) {
    stop("Error: The number of columns in the DEGs exceeds the maximum allowed(2).")
  }
  gsea_df <- DEGs
  colnames(gsea_df) <- c("SYMBOL", "LogFC")
  gsea_df <- gsea_df[order(gsea_df$LogFC, decreasing = T), ] # 按照logFC从大到小排序
  # go_gmtfile <- read.gmt("F:/Bioinformatic_database/MsigDB/c5.go.v2023.1.Hs.entrez.gmt")
  # kegg_gmtfile <- read.gmt("F:/Bioinformatic_database/MsigDB/c2.cp.kegg.v2023.1.Hs.entrez.gmt")

  ### 3.ID coversion
  gsea_SYMBOL <- toupper(gsea_df$SYMBOL) # LogFC排序后基因名
  gsea_ENTREZID <- bitr(gsea_SYMBOL,
    fromType = "SYMBOL",
    toType = c("ENSEMBL", "ENTREZID"),
    OrgDb = org.Hs.eg.db
  )
  gsea_df <- merge(gsea_df, gsea_ENTREZID, by = "SYMBOL")

  ### 4.named logFC
  gsea_LogFC <- gsea_df$LogFC
  names(gsea_LogFC) <- as.character(gsea_df$ENTREZID) # 有时也用gsea_df$SYMBOL
  gsea_LogFC <- sort(gsea_LogFC, decreasing = T)

  ### 5.gseGO
  go_ont <- c("BP", "CC", "MF")
  for (go_i in go_ont) {
    ##
    gsea_GO <- gseGO( # 离线GSEA_GO分析
      gsea_LogFC,
      ont = go_i, # "BP"/"MF"/"CC"/"ALL"
      OrgDb = org.Hs.eg.db, # 该物种对应的org包
      keyType = "ENTREZID",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = maxGSSize,
      eps = 1e-16,
      pvalueCutoff = ifelse(simplify_index, 0.2, 1), # 这个必须为1，影响后面的gsea_GO@result
      pAdjustMethod = "BH",
      verbose = TRUE,
      seed = T, # GSEA存在置换检验，所以每次的结果都不太一样，可以设置seed保持一样
      by = "fgsea" # fast gsea方法
    )

    ## simplify
    if (simplify_index) {
      gsea_GO <- simplify(gsea_GO,
        cutoff = 0.7, by = "p.adjust",
        select_fun = min, measure = "Wang"
      )
    }

    ##
    gsea_GO <- setReadable(gsea_GO, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    gsea_GO_res <- gsea_GO@result
    if (nrow(gsea_GO_res) == 0) {
      patchwork_index <- T
      next
    }
    gsea_GO_res$Des_short <- sapply(gsea_GO_res$Description, shorten_names)
    gsea_GO_res$group <- grp_nm
    write.table(gsea_GO_res,
      file = paste0(output_dir, "/GSEA_GO_", go_i, "_res.txt"),
      row.names = F, col.names = T, quote = F, sep = "\t"
    )

    # go_data_up
    go_data_up <- gsea_GO_res %>%
      filter(NES > 0) %>%
      mutate(new_NES = round(NES, digits = 2))
    if (nrow(go_data_up) >= 20) {
      go_data <- go_data_up[1:20, ]
    } else if (nrow(go_data_up) >= 0) {
      go_data <- go_data_up[1:nrow(go_data_up), ]
    } else {
      # patchwork_index=T
      # next
    }
    go_data <- go_data %>%
      arrange(new_NES, -p.adjust) %>%
      mutate(Des_short_unique = paste0(Des_short, 1:nrow(go_data))) %>%
      mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

    pdf(
      file = paste0(photo_dir, "/GSEA_GO_", go_i, "_up_bar.pdf"),
      width = 6.5, height = 4
    )
    p1 <- ggplot(go_data, aes(x = new_NES, y = Des_short_unique, fill = p.adjust)) +
      geom_col(position = "identity") +
      geom_text(aes(label = new_NES), hjust = -0.15) +
      theme_test() +
      labs(
        x = "NES(Abs)", y = "",
        title = paste0("GSEA_GO_up of ", go_i, " of ", grp_nm)
      ) +
      scale_fill_gradient(low = "#ED0000FF", high = "#00468BFF") +
      scale_x_continuous(expand = c(0, 0, 0.2, 0)) +
      scale_y_discrete(expand = c(0, 1), labels = go_data$Des_short) +
      theme(
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
        axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11, face = "bold", hjust = 0),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    print(p1)
    assign(paste0("GO_", go_i, "_up_bar"), p1)
    dev.off()

    # go_data_down
    go_data_down <- gsea_GO_res %>%
      filter(NES < 0) %>%
      mutate(
        new_NES = -round(NES, digits = 2),
        new_enrichmentScore = -enrichmentScore
      )
    if (nrow(go_data_down) >= 20) {
      go_data <- go_data_down[1:20, ]
    } else if (nrow(go_data_down) >= 0) {
      go_data <- go_data_down[1:nrow(go_data_down), ]
    } else {
      # patchwork_index=T
      # next
    }
    go_data <- go_data %>%
      arrange(new_NES, -p.adjust) %>%
      mutate(Des_short_unique = paste0(Des_short, 1:nrow(go_data))) %>%
      mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

    pdf(
      file = paste0(photo_dir, "/GSEA_GO_", go_i, "_down_bar.pdf"),
      width = 6.5, height = 4
    )
    p2 <- ggplot(go_data, aes(x = new_NES, y = Des_short_unique, fill = p.adjust)) +
      geom_col(position = "identity") +
      geom_text(aes(label = new_NES), hjust = -0.15) +
      theme_test() +
      labs(
        x = "NES(Abs)", y = "",
        title = paste0("GSEA_GO_down of ", go_i, " of ", grp_nm)
      ) +
      scale_fill_gradient(low = "#ED0000FF", high = "#00468BFF") +
      scale_x_continuous(expand = c(0, 0, 0.2, 0)) +
      scale_y_discrete(expand = c(0, 1), labels = go_data$Des_short) +
      theme(
        plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
        axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 11, face = "bold", hjust = 0),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    print(p2)
    assign(paste0("GO_", go_i, "_down_bar"), p2)
    dev.off()
  }

  ### 6.gseKEGG
  gsea_KEGG <- gseKEGG( # 在线/离线GSEA_KEGG，容易网络出错
    gsea_LogFC,
    organism = "hsa",
    keyType = "kegg",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = maxGSSize,
    eps = 1e-16,
    pvalueCutoff = ifelse(simplify_index, 0.2, 1),
    pAdjustMethod = "BH",
    verbose = TRUE,
    use_internal_data = T, # F：在线
    seed = T,
    by = "fgsea"
  )
  gsea_KEGG <- setReadable(gsea_KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  gsea_KEGG_res <- gsea_KEGG@result
  if (nrow(gsea_KEGG_res) == 0) {
    patchwork_index <- T
    return()
  }
  gsea_KEGG_res$Des_short <- sapply(gsea_KEGG_res$Description, shorten_names)
  gsea_KEGG_res$group <- grp_nm
  write.table(gsea_KEGG_res,
    file = paste0(output_dir, "/GSEA_KEGG_res.txt"),
    row.names = F, col.names = T, quote = F, sep = "\t"
  )

  # kegg_data_up
  kegg_data_up <- gsea_KEGG_res %>%
    filter(NES > 0) %>%
    mutate(new_NES = round(NES, digits = 2))
  if (nrow(kegg_data_up) >= 20) {
    kegg_data <- kegg_data_up[1:20, ]
  } else if (nrow(kegg_data_up) >= 0) {
    kegg_data <- kegg_data_up[1:nrow(kegg_data_up), ]
  } else {
    # patchwork_index=T # 当kegg_data为空时，mutate不会报错
    # return() # no loop, can not use next
  }
  kegg_data <- kegg_data %>%
    arrange(new_NES, -p.adjust) %>%
    mutate(Des_short_unique = paste0(Des_short, 1:nrow(kegg_data))) %>%
    mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

  pdf(
    file = paste0(photo_dir, "/GSEA_KEGG_up_bar.pdf"),
    width = 6.5, height = 4
  )
  p3 <- ggplot(kegg_data, aes(x = new_NES, y = Des_short_unique, fill = p.adjust)) +
    geom_col(position = "identity") +
    geom_text(aes(label = new_NES), hjust = -0.15) +
    theme_test() +
    labs(
      x = "NES(Abs)", y = "",
      title = paste0("GSEA_KEGG_up of ", grp_nm)
    ) +
    scale_fill_gradient(low = "#ED0000FF", high = "#00468BFF") +
    scale_x_continuous(expand = c(0, 0, 0.2, 0)) +
    scale_y_discrete(expand = c(0, 1), labels = kegg_data$Des_short) +
    theme(
      plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
      axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 11, face = "bold", hjust = 0),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  print(p3)
  assign(paste0("KEGG_up_bar"), p3)
  dev.off()

  # kegg_data_down
  kegg_data_down <- gsea_KEGG_res %>%
    filter(NES < 0) %>%
    mutate(
      new_NES = -round(NES, digits = 2),
      new_enrichmentScore = -enrichmentScore
    )
  if (nrow(kegg_data_down) >= 20) {
    kegg_data <- kegg_data_down[1:20, ]
  } else if (nrow(kegg_data_down) >= 0) {
    kegg_data <- kegg_data_down[1:nrow(kegg_data_down), ]
  } else {
    # patchwork_index=T
    # return() # no loop, can not use next
  }
  kegg_data <- kegg_data %>%
    arrange(new_NES, -p.adjust) %>%
    mutate(Des_short_unique = paste0(Des_short, 1:nrow(kegg_data))) %>%
    mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

  pdf(
    file = paste0(photo_dir, "/GSEA_KEGG_down_bar.pdf"),
    width = 6.5, height = 4
  )
  p4 <- ggplot(kegg_data, aes(x = new_NES, y = Des_short_unique, fill = p.adjust)) +
    geom_col(position = "identity") +
    geom_text(aes(label = new_NES), hjust = -0.15) +
    theme_test() +
    labs(
      x = "NES(Abs)", y = "",
      title = paste0("GSEA_KEGG_down of ", grp_nm)
    ) +
    scale_fill_gradient(low = "#ED0000FF", high = "#00468BFF") +
    scale_x_continuous(expand = c(0, 0, 0.2, 0)) +
    scale_y_discrete(expand = c(0, 1), labels = kegg_data$Des_short) +
    theme(
      plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
      axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 11, face = "bold", hjust = 0),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  print(p4)
  assign(paste0("KEGG_down_bar"), p4)
  dev.off()

  ### 7.combine pictures
  if (patchwork_index) {
    return()
  }

  pdf(
    file = paste0(photo_dir, "/GSEA_GO_KEGG_updown.pdf"),
    width = 26, height = 8, family = "serif"
  )
  print(
    GO_BP_up_bar + GO_MF_up_bar + GO_CC_up_bar + KEGG_up_bar +
      GO_BP_down_bar + GO_MF_down_bar + GO_CC_down_bar + KEGG_down_bar +
      plot_layout(nrow = 2, guides = "auto")
  )
  dev.off()
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
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/GSEA_GO_KEGG.R")
  GSEA_GO_KEGG(
    DEGs = NULL, simplify_index = F,
    grp_nm = "GSEA_GO_KEGG", dir_nm = "GSEA_GO_KEGG"
  )
}
