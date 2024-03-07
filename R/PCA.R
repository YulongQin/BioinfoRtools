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
#' @param scale_index whether to scale the data
#' @param method the method of dimensionality reduction
#' @param grp_nm the name of the group
#'
#' @return plot
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
#' df <- iris[1:4] %>% t() %>% as.data.frame()
#' PCA(df,iris$Species,scale_index = T
#' PCA(df,iris$Species,method = "UMAP",scale_index = T)
#' PCA(df,iris$Species,method = "tSNE",scale_index = T)
#' PCA(df,iris$Species,method = "LDA",scale_index = T)
#' PCA(df,iris$Species,method = "lda",scale_index = T)
#' }
#'
#' @import ggplot2
#' @importFrom stats prcomp
#' @importFrom ggsci scale_color_nejm scale_fill_nejm scale_color_igv scale_fill_igv pal_lancet
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @importFrom MASS lda
#' @export
#'
PCA <- function(countdata,group,scale_index=T,
                method="PCA",grp_nm="PCA") {
  UseMethod("PCA",countdata)
}

#' @export
#' @rdname PCA
PCA.data.frame <- function(countdata,group,scale_index=T,
                           method="PCA",grp_nm="PCA"){

  ### 1.library ####
  # suppressMessages({
  #   library(tidyverse)
  #   library(magrittr)
  #   library(ggsci)
  #   library(Rtsne)
  #   library(umap)
  #   library(MASS)
  # })

  ### 2.data ####
  countdata <- countdata %>% # countdata[updown,saps[,1]]
    # round(digits = 0) %>% # ！！！ 取整，Deseq2才要
    na.omit()
  countdata <- countdata[rowSums(countdata)>0,] # 删除全0行
  countdata <- countdata[rowSums(countdata>0)>0.3*ncol(countdata),]
  countdata <- countdata[apply(countdata, 1, var)!=0 ,] # ！！！去除方差为0的基因，PCA特有

  ### 3.scale ####
  countdata <- countdata[,!is.na(group)]
  group <- na.omit(group)
  if(scale_index==T){
    count_scale <- countdata %>%
      t() %>%
      scale(center = TRUE, scale = TRUE) %>%
      as.data.frame()
  }else{
    count_scale <- countdata %>%
      t() %>% as.data.frame()
  }

  if(method=="PCA"){
    pca_info_sum <- count_scale %>%
      prcomp(center = F, # Z-score标准化，center中心化，scale标准化
             scale= F) %>%
      summary()
    pca_data <- data.frame(Sample = rownames(pca_info_sum$x),
                           Group = as.factor(group),
                           pca_info_sum$x) # PC对各个样本的解释度
  }else if(method=="UMAP"){
    umap_info <- umap::umap(
      d=count_scale,
      random_state=123,
      preserve.seed = TRUE # 使用固定seed
    )
    umap_data <- data.frame(umap_info$layout,group)
    colnames(umap_data) <- c("UMAP1","UMAP2","Group")
    head(umap_data)
  }else if(method=="tSNE"){
    rtsne_info <- Rtsne(
      count_scale,
      perplexity = 30,
      check_duplicates = FALSE
    )
    tsne_data <- data.frame(rtsne_info$Y,group)
    colnames(tsne_data) <- c("tSNE_1","tSNE_2","Group")
  }else if(method=="LDA"){
    lda_info <- lda(count_scale,group)
    lda_data <- data.frame(predict(lda_info)$x,group)
    colnames(lda_data) <- c("LD1","LD2","Group")
  }else{
    stop("method error")
  }

  ### 4.plot ####
  if(method=="PCA"){
    if(scale_index==T){
      p <- ggplot(data = pca_data,aes(x=PC1,y=PC2,color=Group,fill=Group))+
        stat_ellipse(geom ="polygon",type = "norm",level = 0.95,
                     alpha=0.2,color=NA)+
        geom_hline(yintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
        geom_vline(xintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
        geom_point(size=1.5,shape=16)+
        theme_test()+
        labs(x=paste0("PC1(",round(pca_info_sum$importance[2,1]*100,2),"%)"),
             y=paste0("PC2(",round(pca_info_sum$importance[2,2]*100,2),"%)"),
             title = grp_nm )+
        guides(fill='none')+ #去除图例
        theme(axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
              plot.title = element_text(size = 16,face = 'bold',hjust = 0.5))
    }else{
      p <- ggplot(data = pca_data,aes(x=PC1,y=PC2,color=Group,fill=Group))+
        stat_ellipse(geom ="polygon",type = "norm",level = 0.95,
                     alpha=0.2,color=NA)+
        # geom_hline(yintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
        # geom_vline(xintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
        geom_point(size=1.5,shape=16)+
        theme_test()+
        labs(x=paste0("PC1(",round(pca_info_sum$importance[2,1]*100,2),"%)"),
             y=paste0("PC2(",round(pca_info_sum$importance[2,2]*100,2),"%)"),
             title = grp_nm )+
        guides(fill='none')+ #去除图例
        theme(axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
              plot.title = element_text(size = 16,face = 'bold',hjust = 0.5))
    }
    if(length(unique(group))<9){
      p <- p + scale_color_nejm()+
        scale_fill_nejm()
    }else if(length(unique(group))<=50){
      p <- p + scale_color_igv()+
        scale_fill_igv()
    }else{
      stop("group too many")
    }
  }else if(method=="UMAP"){
    p <- ggplot(data = umap_data,aes(x=UMAP1,y=UMAP2,color=Group,fill=Group))+
      stat_ellipse(geom ="polygon",type = "norm",level = 0.95,
                   alpha=0.1,color=NA)+
      geom_hline(yintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
      geom_vline(xintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
      geom_point(size=1.5,shape=16)+
      theme_test()+
      labs(x="UMAP1",y="UMAP2",title = "UMAP" )+
      guides(fill='none')+ #去除图例
      theme(axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
            plot.title = element_text(size = 16,face = 'bold',hjust = 0.5))
    if(length(unique(group))<9){
      p <- p + scale_color_nejm()+
        scale_fill_nejm()
    }else if(length(unique(group))<=50){
      p <- p + scale_color_igv()+
        scale_fill_igv()
    }else{
      stop("group too many")
    }
  }else if(method=="tSNE"){
    p <- ggplot(data = tsne_data,aes(x=tSNE_1,y=tSNE_2,color=Group,fill=Group))+
      stat_ellipse(geom ="polygon",type = "norm",level = 0.95,
                   alpha=0.1,color=NA)+
      geom_hline(yintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
      geom_vline(xintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
      geom_point(size=1.5,shape=16)+
      theme_test()+
      labs(x="tSNE_1",y="tSNE_2",title = "tSNE" )+
      guides(fill='none')+ #去除图例
      theme(axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
            plot.title = element_text(size = 16,face = 'bold',hjust = 0.5))
    if(length(unique(group))<9){
      p <- p + scale_color_nejm()+
        scale_fill_nejm()
    }else if(length(unique(group))<=50){
      p <- p + scale_color_igv()+
        scale_fill_igv()
    }else{
      stop("group too many")
    }
  }else if(method=="LDA"){
    p <- ggplot(data = lda_data,aes(x=LD1,y=LD2,color=Group,fill=Group))+
      stat_ellipse(geom ="polygon",type = "norm",level = 0.95,
                   alpha=0.1,color=NA)+
      geom_hline(yintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
      geom_vline(xintercept= 0,linetype=2,color="#9C9C9C",linewidth=0.5)+
      geom_point(size=1.5,shape=16)+
      theme_test()+
      labs(x="LD1",y="LD2",title = "LDA" )+
      guides(fill='none')+ #去除图例
      theme(axis.title = element_text(size = 12,face = 'bold',hjust = 0.5),
            plot.title = element_text(size = 16,face = 'bold',hjust = 0.5))
    if(length(unique(group))<9){
      p <- p + scale_color_nejm()+
        scale_fill_nejm()
    }else if(length(unique(group))<=50){
      p <- p + scale_color_igv()+
        scale_fill_igv()
    }else{
      stop("group too many")
    }
  }else{}
  p

  return(p)
}

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/PCA.R")
  df <- iris[1:4] %>% t() %>% as.data.frame()
  PCA(df,iris$Species,scale_index = T)
  PCA(df,iris$Species,method = "UMAP",scale_index = T)
  PCA(df,iris$Species,method = "tSNE",scale_index = T)
  PCA(df,iris$Species,method = "LDA",scale_index = T)
  PCA(df,iris$Species,method = "lda",scale_index = T)
}
