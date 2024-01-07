### ---------------
###
### Create: qyl
### Date: 2023-04-01 22:29:00
### Email: 524583919@qq.com
### Function：
### - PCA plots are drawn according to groups
### -
###
### Update Log: 2023-04-01
###
### ---------------

#### ----- Function ------ ####
#' PCA plots are drawn according to groups
#'
#' this function use `ggplot()` to draw PCA plots according to groups
#'
#' @name PCA
#' @title PCA plots are drawn according to groups
#' @description PCA plots are drawn according to groups
#' @param countdata the count data matrix
#' @param group the group information
#' @param grp_nm the name of the group
#' @param scale whether to scale the data
#'
#' @return NULL
#' @author Yulong Qin
#' @seealso \code{\link{PCA}}, \code{\link{HTmap}}
#' @seealso [ggplot()], [ggplot2::ggplot()] which this function wraps.
#'
#' @section Useful mutate functions:
#'
#' * [`+`], [`-`], [log()], etc., for their usual mathematical meanings
#'
#' ...
#'
#' @examples \dontrun{
#'
#' # this is a example
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/PCA.R")
#' df <- iris[1:4] %>%
#'   t() %>%
#'   as.data.frame()
#' PCA(df, iris$Species, scale = F)
#' }
#'
#' @import ggplot2
#' @importFrom stats prcomp
#' @importFrom ggsci scale_color_nejm scale_fill_nejm scale_color_igv scale_fill_igv pal_lancet
#' @export
#'
PCA <- function(countdata, group, grp_nm = "PCA", scale = T) {
  UseMethod("PCA",countdata)
}

#' @export
#' @rdname PCA
PCA.data.frame <- function(countdata, group, grp_nm = "PCA", scale = T) {
  ### 1.library ####
  # suppressMessages({
  #   library(tidyverse)
  #   library(magrittr)
  #   library(ggsci)
  # })

  ### 2.data ####
  countdata <- countdata %>% # countdata[updown,saps[,1]]
    # round(digits = 0) %>% # ！！！ 取整，Deseq2才要
    na.omit()
  countdata <- countdata[rowSums(countdata) > 0, ] # 删除全0行
  countdata <- countdata[rowSums(countdata > 0) > 0.3 * ncol(countdata), ]
  countdata <- countdata[apply(countdata, 1, var) != 0, ] # ！！！去除方差为0的基因，PCA特有

  ### 3.scale ####
  countdata <- countdata[, !is.na(group)]
  group <- na.omit(group)
  if (scale == T) {
    pca_info_sum <- countdata %>%
      t() %>%
      prcomp(
        center = TRUE, # Z-score标准化，center中心化，scale标准化
        scale = TRUE
      ) %>%
      summary()
  } else {
    pca_info_sum <- countdata %>%
      t() %>%
      prcomp(
        center = F,
        scale = F
      ) %>%
      summary()
  }
  pca_data <- data.frame(
    Sample = rownames(pca_info_sum$x),
    Group = as.factor(group),
    pca_info_sum$x # PC对各个样本的解释度
  )

  ### 4.plot ####
  p <- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
    stat_ellipse(
      geom = "polygon", type = "norm", level = 0.95,
      alpha = 0.2, color = NA
    ) +
    geom_point(size = 1.5, shape = 16) +
    theme_test() +
    labs(
      x = paste0("PC1(", round(pca_info_sum$importance[2, 1] * 100, 2), "%)"),
      y = paste0("PC2(", round(pca_info_sum$importance[2, 2] * 100, 2), "%)"),
      title = grp_nm
    ) +
    guides(fill = "none") + # 去除图例
    theme(
      axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  if (length(unique(group)) < 9) {
    p <- p + scale_color_nejm() +
      scale_fill_nejm()
  } else {
    p <- p + scale_color_igv() +
      scale_fill_igv()
  }
  if (scale == T) {
    p <- p + geom_hline(yintercept = 0, linetype = 2, color = "#9C9C9C", linewidth = 0.5) +
      geom_vline(xintercept = 0, linetype = 2, color = "#9C9C9C", linewidth = 0.5)
  } else {
    p
  }

  return(p)
}

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/PCA.R")
  df <- iris[1:4] %>%
    t() %>%
    as.data.frame()
  PCA(df, iris$Species, scale = F)
}
