### ---------------
###
### Create: qyl
### Date: 2023-05-26 20:09:00
### Email: 524583919@qq.com
### Function：
### - GO and KEGG analysis
### -
###
### Update Log: 2023-05-26
###
### ---------------

#### ----- Function ------ ####
#' GO and KEGG analysis
#'
#' @title GO and KEGG analysis
#' @description GO and KEGG analysis
#' @param updown updown genes
#' @param up up genes
#' @param down down genes
#' @param simplify_index whether simplify
#' @param maxGSSize max gene set size
#' @param grp_nm the name of group
#' @param dir_nm the name of dir
#'
#' @return NULL
#' @author Yulong Qin
#' @seealso \code{\link{GO_KEGG}},\code{\link{GSEA_GO_KEGG}}, \code{\link{GSVA_ssGSEA}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/GO_KEGG.R")
#' GO_KEGG(
#'   updown = gene_bf_go_updown, up = gene_bf_go_up, down = gene_bf_go_down,
#'   simplify_index = F, grp_nm = "GO_KEGG", dir_nm = "GO_KEGG"
#' )
#' }
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import DOSE
#' @import ggplot2
#' @import patchwork
#' @importFrom rlang is_empty
#' @export GO_KEGG
GO_KEGG <- function(updown = NULL, up = NULL, down = NULL, simplify_index = F,
                    maxGSSize = 500, grp_nm = "GO_KEGG", dir_nm = "GO_KEGG") {
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

  ### 2.ID conversion ####
  ent_updown <- bitr(updown, # 基因名大写
    fromType = "SYMBOL",
    toType = c("ENSEMBL", "ENTREZID"), # 'ENTREZID'最适合富集
    OrgDb = org.Hs.eg.db
  )$ENTREZID
  ent_up <- bitr(up,
    fromType = "SYMBOL",
    toType = c("ENSEMBL", "ENTREZID"),
    OrgDb = org.Hs.eg.db
  )$ENTREZID
  ent_down <- bitr(down,
    fromType = "SYMBOL",
    toType = c("ENSEMBL", "ENTREZID"),
    OrgDb = org.Hs.eg.db
  )$ENTREZID

  ### 3.GO and KEGG ####
  vector_ent <- c()
  if (!is_empty(ent_updown)) {
    vector_ent <- c(vector_ent, as.character(quote(ent_updown))) # "ent_updown"
  }
  if (!is_empty(ent_up)) {
    vector_ent <- c(vector_ent, as.character(quote(ent_up)))
  }
  if (!is_empty(ent_down)) {
    vector_ent <- c(vector_ent, as.character(quote(ent_down)))
  }
  go_ont <- c("BP", "CC", "MF")
  for (go_i in vector_ent) { # "ent_updown"
    tryCatch(
      {
        ## 3.1 go
        for (go_j in go_ont) { # "BP","CC","MF"
          tryCatch(
            {
              go <- enrichGO(
                gene = get(go_i), # entrezgene_id
                OrgDb = "org.Hs.eg.db", # 物种注释包
                keyType = "ENTREZID", # ID类型
                ont = go_j, # ALL,MF:molecular function ,BP:biological process ,CC:cellular compotent
                pAdjustMethod = "BH", # 多重假设检验矫正方法
                pvalueCutoff = ifelse(simplify_index, 0.2, 1), # 设置过滤标准，影响go的矩阵数据，影响后续许多分析
                qvalueCutoff = 1, # 同上
                maxGSSize = maxGSSize,
                pool = T, # ont=all,同时输出3个通路
                readable = T
              ) # 将go结果中的ENTREZID转为gene symbol
              # assign(paste0('go_',go_i,'_',go_j),go)

              ## simplify
              if (simplify_index) {
                go <- simplify(go,
                  cutoff = 0.7, by = "p.adjust",
                  select_fun = min, measure = "Wang"
                )
              }

              ## go_res
              # go <- setReadable(go,OrgDb=org.Hs.eg.db, keyType = "ENTREZID") # enrichGO中已经有readable = T
              go_res <- go@result # go@result是原始的所有结果，go保存的是按照指标过滤后的结果
              if (nrow(go_res) == 0) {
                patchwork_index <- T
                next
              }
              go_res$Des_short <- sapply(go_res$Description, shorten_names)
              go_res$group <- grp_nm
              # assign(paste0('go_',go_i,'_',go_j,'_res'),go_res)
              write.table(go_res,
                file = paste0(output_dir, "/GO_", go_i, "_", go_j, "_res.txt"),
                row.names = F, col.names = T, quote = F, sep = "\t"
              )

              ## barplot
              if (nrow(go_res) >= 20) {
                go_data <- go_res[1:20, ] # 自动挑选感兴趣的通路
              } else if (nrow(go_res) >= 0) {
                go_data <- go_res[1:nrow(go_res), ]
              } else {
                # patchwork_index=T
                # next
              }
              if (F) {
                go_data <- read.table(
                  file = "./inputdata/go_data.txt", # 手动挑选，手动读取
                  sep = "\t",
                  header = T,
                  row.names = NULL
                )
              }
              go_data <- go_data %>%
                arrange(Count, -p.adjust) %>%
                mutate(Des_short_unique = paste0(Des_short, 1:nrow(go_data))) %>%
                mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

              pdf(
                file = paste0(photo_dir, "/GO_", go_i, "_", go_j, "_bar.pdf"),
                width = 6.5, height = 4.5
              )
              p1 <- ggplot(go_data, aes(x = Count, y = Des_short_unique, fill = p.adjust)) +
                geom_col(position = "identity") +
                geom_text(aes(label = Count), hjust = -0.3) +
                theme_test() +
                labs(
                  x = "Counts", y = "",
                  title = paste0("GO of ", go_i, "_", go_j, " of ", grp_nm)
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
              assign(paste0("GO_", go_i, "_", go_j, "_bar"), p1)
              dev.off()

              ## bubble diagram
              GeneRatio_all <- go_res$GeneRatio[1] %>%
                gsub(".*/", "", .) %>%
                as.numeric()
              go_data <- go_data %>%
                mutate(GeneRatio_all = GeneRatio_all, .before = 10)

              pdf(
                file = paste0(photo_dir, "/GO_", go_i, "_", go_j, "_bubble.pdf"),
                width = 6.5, height = 4.5
              )
              p2 <- ggplot(go_data, aes(x = Count / GeneRatio_all, y = Des_short_unique, size = Count, color = p.adjust)) +
                geom_point() +
                theme_bw() +
                scale_color_gradient(low = "#ED0000FF", high = "#00468BFF") +
                scale_x_continuous(expand = c(0.1, 0, 0.1, 0)) +
                scale_y_discrete(expand = c(0, 1), labels = go_data$Des_short) +
                labs(
                  x = "Gene Ratio", y = "", color = "P.adjust",
                  title = paste0("GO of ", go_i, "_", go_j, " of ", grp_nm)
                ) +
                guides(color = guide_colorbar(order = 0), size = guide_legend(order = 1)) +
                theme(
                  plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
                  axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
                  axis.text = element_text(size = 10, face = "bold", hjust = 0.5),
                  legend.title = element_text(size = 11, face = "bold", hjust = 0),
                  plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
                )
              print(p2)
              assign(paste0("GO_", go_i, "_", go_j, "_bubble"), p2)
              dev.off()

              rm(go, go_res, go_data, GeneRatio_all)
            },
            error = function(e) print(e)
          )
        }
      },
      error = function(e) print(e)
    )

    ## 3.2 KEGG
    tryCatch(
      {
        kegg <- enrichKEGG(get(go_i),
          organism = "hsa", # KEGG中物种缩写
          pvalueCutoff = ifelse(simplify_index, 0.2, 1), # 卡的p.adjust的值
          keyType = "kegg", # KEGG中编号格式
          pAdjustMethod = "BH", # 多重假设检验矫正
          qvalueCutoff = 1, # 卡的qvalue值
          maxGSSize = maxGSSize,
          use_internal_data = T
        ) # T:本地的KEGG.db数据库时,F:使用在线数据库
        # assign(paste0('kegg_',go_i),kegg)

        ## kegg_res
        kegg <- setReadable(kegg, # 自动ent→sym
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID"
        )
        kegg_res <- kegg@result
        kegg_res$Des_short <- sapply(kegg_res$Description, shorten_names)
        kegg_res$group <- grp_nm
        # assign(paste0('kegg_',go_i,'_res'),kegg_res)
        write.table(kegg_res,
          file = paste0(output_dir, "/KEGG_", go_i, "_res.txt"),
          row.names = F, col.names = T, quote = F, sep = "\t"
        )

        ## barplot
        if (nrow(kegg_res) >= 20) {
          kegg_data <- kegg_res[1:20, ] # 挑选感兴趣的通路
        } else if (nrow(kegg_res) >= 0) {
          kegg_data <- kegg_res[1:nrow(kegg_res), ]
        } else {
          # patchwork_index=T
          # next
        }
        if (F) {
          kegg_data <- read.table(
            file = "./inputdata/kegg_data.txt", # 手动挑选，手动读取
            sep = "\t",
            header = T,
            row.names = NULL
          )
        }
        kegg_data <- kegg_data %>%
          arrange(Count, -p.adjust) %>%
          mutate(Des_short_unique = paste0(Des_short, 1:nrow(kegg_data))) %>%
          mutate(Des_short_unique = factor(Des_short_unique, levels = unique(Des_short_unique)))

        pdf(
          file = paste0(photo_dir, "/KEGG_", go_i, "_bar.pdf"),
          width = 6.5, height = 4.5
        )
        p3 <- ggplot(kegg_data, aes(x = Count, y = Des_short_unique, fill = p.adjust)) +
          geom_col(position = "identity") +
          geom_text(aes(label = Count), hjust = -0.3) +
          theme_test() +
          labs(
            x = "Counts", y = "",
            title = paste0("KEGG of ", go_i, " of ", grp_nm)
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
        assign(paste0("KEGG_", go_i, "_bar"), p3)
        dev.off()

        ## bubble diagram
        GeneRatio_all <- kegg_res$GeneRatio[1] %>%
          gsub(".*/", "", .) %>%
          as.numeric()
        kegg_data <- kegg_data %>%
          mutate(GeneRatio_all = GeneRatio_all, .before = 10)

        pdf(
          file = paste0(photo_dir, "/KEGG_", go_i, "_bubble.pdf"),
          width = 6.5, height = 4.5
        )
        p4 <- ggplot(kegg_data, aes(x = Count / GeneRatio_all, y = Des_short_unique, size = Count, color = p.adjust)) +
          geom_point() +
          theme_bw() +
          scale_color_gradient(low = "#ED0000FF", high = "#00468BFF") +
          scale_x_continuous(expand = c(0.1, 0, 0.1, 0)) +
          scale_y_discrete(expand = c(0, 1), labels = kegg_data$Des_short) +
          labs(
            x = "Gene Ratio", y = "", color = "P.adjust",
            title = paste0("KEGG of ", go_i, " of ", grp_nm)
          ) +
          guides(color = guide_colorbar(order = 0), size = guide_legend(order = 1)) +
          theme(
            plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
            axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
            axis.text = element_text(size = 10, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 11, face = "bold", hjust = 0),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
        print(p4)
        assign(paste0("KEGG_", go_i, "_bubble"), p4)
        dev.off()

        rm(kegg, kegg_res, kegg_data)
      },
      error = function(e) print(e)
    )
  }

  ### 4. combine pictures
  if (patchwork_index) {
    return()
  }

  ## up
  tryCatch(
    {
      if (!is_empty(ent_up)) {
        pdf(
          file = paste0(photo_dir, "/GO_KEGG_up.pdf"),
          width = 26, height = 8, family = "serif"
        )
        print(
          GO_ent_up_BP_bar + GO_ent_up_MF_bar + GO_ent_up_CC_bar + KEGG_ent_up_bar +
            GO_ent_up_BP_bubble + GO_ent_up_MF_bubble + GO_ent_up_CC_bubble + KEGG_ent_up_bubble +
            plot_layout(nrow = 2, guides = "auto")
        )
        dev.off()
      }
    },
    error = function(e) print(e)
  )

  ## down
  tryCatch(
    {
      if (!is_empty(ent_down)) {
        pdf(
          file = paste0(photo_dir, "/GO_KEGG_down.pdf"),
          width = 26, height = 8, family = "serif"
        )
        print(
          GO_ent_down_BP_bar + GO_ent_down_MF_bar + GO_ent_down_CC_bar + KEGG_ent_down_bar +
            GO_ent_down_BP_bubble + GO_ent_down_MF_bubble + GO_ent_down_CC_bubble + KEGG_ent_down_bubble +
            plot_layout(nrow = 2, guides = "auto")
        )
        dev.off()
      }
    },
    error = function(e) print(e)
  )

  ## updpwn
  tryCatch(
    {
      if (!is_empty(ent_updown)) {
        pdf(
          file = paste0(photo_dir, "/GO_KEGG_updown.pdf"),
          width = 26, height = 8, family = "serif"
        )
        print(
          GO_ent_updown_BP_bar + GO_ent_updown_MF_bar + GO_ent_updown_CC_bar + KEGG_ent_updown_bar +
            GO_ent_updown_BP_bubble + GO_ent_updown_MF_bubble + GO_ent_updown_CC_bubble + KEGG_ent_updown_bubble +
            plot_layout(nrow = 2, guides = "auto")
        )
        dev.off()
      }
    },
    error = function(e) print(e)
  )
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
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/GO_KEGG.R")
  GO_KEGG(
    updown = gene_bf_go_updown, up = gene_bf_go_up, down = gene_bf_go_down,
    simplify_index = F, grp_nm = "GO_KEGG", dir_nm = "GO_KEGG"
  )
}
