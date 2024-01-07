### ---------------
###
### Create: qyl
### Date: 2022-06-30 17:46:00
### Email: 524583919@qq.com
### Function：
### - Draw heat maps according to groups
### -
###
### Update Log: 2022-06-30
###
### ---------------

#### ----- Function ------ ####
#' Draw heat maps according to groups
#'
#' @title Draw heat maps according to groups
#' @description Draw heat maps according to groups
#' @param countdata the count data matrix
#' @param group the group information
#' @param log2_index whether to log2 the data
#' @param scale_index whether to scale the data
#'
#' @return NULL
#' @author Yulong Qin
#' @seealso \code{\link{PCA}}, \code{\link{HTmap}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/HTmap.R")
#' df <- iris[1:4] %>%
#'   t() %>%
#'   as.data.frame()
#' HTmap(df, iris$Species)
#' }
#' @import ggplot2
#' @importFrom ggsci scale_color_nejm scale_fill_nejm scale_color_igv scale_fill_igv
#' @importFrom scales alpha
#' @importFrom ggplotify as.ggplot
#' @importFrom pheatmap pheatmap
#' @importFrom dplyr across everything mutate if_else
#' @export
#'
HTmap <- function(countdata, group, log2_index = T, scale_index = "row") {
  ### 1.library ####
  # suppressMessages({
  #   library(pheatmap)
  #   library(magrittr)
  #   library(scales)
  #   library(ggsci)
  #   library(ggplot2)
  #   library(ggplotify)
  # })

  ### 2.热图矩阵 ####
  hm_in <- countdata %>%
    mutate(across(everything(), ~ if_else(is.na(.), 0, .))) %>% # 不应直接删除有缺失行na.omit()
    # filter(rowSums(.)>0) %>% # 不需要过滤行
    {
      if (log2_index) add(., 1) else .
    } %>% # 理解+1与+10^6的区别
    {
      if (log2_index) log2(.) else .
    } %>%
    filter(apply(., 1, var) != 0) # filter中加apply可以对每行遍历

  # countdata <- countdata[rowSums(countdata)>0,] # 删除全0行
  # hm_in <- log2(countdata+1)
  # hm_in <- hm_in[apply(hm_in, 1, var) != 0,] %>%
  #   na.omit()

  ### 3.行列注释 ####
  group <- factor(group, levels = as.character(unique(group)))
  hm_ann_col <- data.frame(
    Group = group,
    row.names = colnames(hm_in) # 行名为样本,与分组对应
  )
  # color <- c('#0000CD','#DC143C') # pal_lancet()(2) %>% rev()
  color <- pal_lancet()(length(unique(group)))
  names(color) <- levels(group) # ctrl在前，case在后
  hm_ann_colors <- list(
    Group = alpha(color, 0.8) # 含有名字的向量
  )

  ### 4.热图 ####
  p1 <- pheatmap(hm_in,
    cluster_cols = F, # 行列聚类
    cluster_rows = T,
    scale = scale_index, # 行归一化z-score
    # cellwidth = 6, #设定cellwidth，不设定cellheigh(跟随图片)
    color = c(
      colorRampPalette(colors = c("blue", "white"))(100),
      colorRampPalette(colors = c("white", "red"))(100)
    ),
    border_color = NA, # "grey60"
    main = "",
    annotation_row = NULL,
    annotation_col = hm_ann_col,
    annotation_colors = hm_ann_colors,
    show_colnames = F, # 行列名
    show_rownames = T,
    fontsize = 8,
    fontsize_row = 6, # 行列名字体大小
    fontsize_col = 6,
    display_numbers = F, # 添加数字
    silent = TRUE
  )
  as.ggplot(p1) +
    theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4)) # 调整边距
}

#### ----- Examples ------ ####
if (F) {
  # 数据不用提前log2，可以绘制两组或多组，需要指定向量
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/HTmap.R")
  df <- iris[1:4] %>%
    t() %>%
    as.data.frame()
  HTmap(df, iris$Species)
}
