### ---------------
###
### Create: qyl
### Date: 2023-12-13 15:06:00
### Email: 524583919@qq.com
### Function：
### - Mfuzz analysis
### -
###
### Update Log: 2023-12-13
###
### ---------------

#### ----- Function ------ ####
#' Mfuzz analysis
#'
#' @title Mfuzz analysis
#' @description Mfuzz analysis
#' @param exp the input expression matrix
#' @param grp the input group
#' @param cluster_num the number of cluster
#' @param min_acore the min acore of cluster
#' @param GO_KEGG whether to do GO and KEGG analysis
#' @param grp_nm the name of group
#' @param dir_nm the name of directory
#'
#' @return NULL
#' @author Yulong Qin
#' @seealso \code{\link{WGCNA_pipeline}}, \code{\link{miRNA_mRNA}}, \code{\link{Mfuzz_pipeline}}
#'
#' @examples \dontrun{
#'   source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/Mfuzz_pipeline.R")
#'   Mfuzz_pipeline(exp = exp,grp = grp,
#'               cluster_num = 20,min_acore = 0.3,
#'               GO_KEGG = F,
#'               grp_nm = "Mfuzz",dir_nm = "Mfuzz")
#'               }
#' @import Mfuzz
#' @export
#'
Mfuzz_pipeline <- function(exp = NULL,grp = NULL,
                           cluster_num = 20,min_acore = 0.3,
                           GO_KEGG = F,
                           grp_nm = "Mfuzz",dir_nm = "Mfuzz"){

  ### 1.library ####
  cat(paste0("1.library start at ",Sys.time(),"\n"))
  # suppressMessages({
  #   library(tidyverse)
  #   library(Mfuzz)
  # })
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm)
  photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm)
  unlink(output_dir,recursive = T)
  unlink(photo_dir,recursive = T)
  dir.create(output_dir,recursive = T)
  dir.create(photo_dir,recursive = T)

  ### 2.process data ####
  cat(paste0("2.process data start at ",Sys.time(),"\n"))
  exp_mean <- exp %>%
    t() %>%
    as.data.frame() %>%
    mutate(group=grp) %>%
    group_by(group) %>%
    summarise(across(everything(), mean)) %>%
    t() %>%
    as.data.frame()
  colnames(exp_mean) <- exp_mean[1,]
  exp_mean <- exp_mean[-1,]
  exp_mean <- exp_mean %>%
    mutate(across(everything(),as.numeric))

  ### 3.Mfuzz ####
  cat(paste0("3.Mfuzz start at ",Sys.time(),"\n"))
  protein <- as.matrix(exp_mean)
  mfuzz_class <- new('ExpressionSet',exprs = protein)
  {
    mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
    mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
    mfuzz_class <- filter.std(mfuzz_class, min.std = 0,visu=F)
  }
  mfuzz_class <- standardise(mfuzz_class)
  set.seed(123)
  # cluster_num <- 20
  m <- mestimate(mfuzz_class)
  mfuzz_cluster <- mfuzz(mfuzz_class,
                         centers = cluster_num,
                         m = m)
  ### 4.plot ####
  cat(paste0("4.plot start at ",Sys.time(),"\n"))
  ifelse(cluster_num<=10,mfow_index <- c(5, 2),
         ifelse(cluster_num<=20,mfow_index <- c(5, 4),
                ifelse(cluster_num<=30,mfow_index <- c(5, 6),mfow_index <- c(5, 8))))
  width_index <- ifelse(cluster_num<=10,6,
                        ifelse(cluster_num<=20,12,
                               ifelse(cluster_num<=30,18,24)))
  pdf(file = paste0(photo_dir,"/mfuzz_merge_all_",cluster_num,".pdf"),
      width = width_index,height = 15)
  mfuzz.plot2(mfuzz_class,
              cl = mfuzz_cluster,
              mfrow = mfow_index,
              time.labels = colnames(protein),
              Xwidth=12,Xheight=10,
              min.mem = 0,
              x11=F)
  dev.off()
  pdf(file = paste0(photo_dir,"/mfuzz_merge_fil_",cluster_num,".pdf"),
      width = width_index,height = 15)
  mfuzz.plot2(mfuzz_class,
              cl = mfuzz_cluster,
              mfrow = mfow_index,
              time.labels = colnames(protein),
              Xwidth=12,Xheight=10,
              min.mem = 0.3,
              x11=F)
  dev.off()
  for (i in 1:cluster_num) {
    pdf(file = paste0(photo_dir,"/mfuzz_fil_",i,".pdf"),
        width = 10,height = 10)
    mat <- matrix(1:2,ncol=2,nrow=1,byrow=TRUE)
    layout(mat,widths=c(5,1))

    mfuzz.plot2(mfuzz_class,
                cl=mfuzz_cluster,
                mfrow=NA,
                colo="fancy",
                single=i,
                min.mem = 0.3,
                x11=F)
    mfuzzColorBar(col = "fancy",
                  main = "Membership",
                  cex.main = 1) # 添加颜色图例
    dev.off()
  }

  ### 5.cluster of genes ####
  cat(paste0("5.cluster of genes start at ",Sys.time(),"\n"))
  {
    cluster_size <- mfuzz_cluster$size
    names(cluster_size) <- paste0("cluster_",1:cluster_num)
    cluster_size # 未过滤的
    write.table(cluster_size,
                file = paste0(output_dir,"/cluster_size_all.txt"),
                sep = "\t",col.names = F,row.names = T,quote = F)

    protein_mem <- mfuzz_cluster$membership %>%
      as.data.frame() %>%
      dplyr::rename_with(~paste0("cluster_",.x),everything())
  }
  {
    protein_cluster <- mfuzz_cluster$cluster
    protein_oriexp_cluster <- cbind(protein[names(protein_cluster), ],
                                    protein_cluster,
                                    protein_mem)
    write.table(protein_oriexp_cluster,
                file = paste0(output_dir,"/protein_oriexp_cluster.txt"),
                sep = '\t', col.names = NA, quote = FALSE)
  }
  {
    protein_cluster <- mfuzz_cluster$cluster
    protein_standard <- mfuzz_class@assayData$exprs
    protein_standard_cluster <- cbind(protein_standard[names(protein_cluster), ],
                                      protein_cluster,
                                      protein_mem) %>%
      as.data.frame()
    write.table(protein_standard_cluster,
                file = paste0(output_dir,"/protein_standard_cluster.txt"),
                sep = '\t', col.names = NA, quote = FALSE)
  }

  ### 6.min_acore filtering ####
  cat(paste0("6.min_acore filtering start at ",Sys.time(),"\n"))
  acore_fil <- acore(
    mfuzz_class,mfuzz_cluster,
    min.acore = min_acore # 根据membership值过滤基因
  )
  saveRDS(acore_fil,file = paste0(output_dir,"/acore_genes_fil.Rds"))
  names(acore_fil) <- paste0("cluster_",1:cluster_num)
  cluster_size_fil <- sapply(1:length(acore_fil),function(i){
    a <- nrow(acore_fil[[i]])
    names(a) <- names(acore_fil[i])
    return(a)
  })
  write.table(cluster_size_fil,
              file = paste0(output_dir,"/cluster_size_fil.txt"),
              sep = "\t",col.names = F,row.names = T,quote = F)

  ### 7.GO_KEGG ####
  cat(paste0("7.GO_KEGG start at ",Sys.time(),"\n"))
  for (i_cluster in 1:cluster_num) {
    cluster_id <- paste0("cluster_",i_cluster)
    cat(paste0(cluster_id," start at ",Sys.time(),"\n"))

    dir.create(paste0(output_dir,"/",cluster_id))
    # dir.create(paste0(photo_dir,"/",cluster_id))
    acore_cluster <- acore_fil[[cluster_id]]$NAME
    write.table(acore_fil[[cluster_id]],
                file = paste0(output_dir,"/",cluster_id,"/fil_gene_",cluster_id,".txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)
    df_cluster <- protein_standard_cluster[acore_cluster,]
    write.table(df_cluster,
                file = paste0(output_dir,"/",cluster_id,"/protein_standard_cluster_",cluster_id,".txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)
    if(GO_KEGG){
      source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/GO_KEGG.R")
      GO_KEGG(updown = acore_cluster ,up = NULL,down = NULL,maxGSSize = 100,
              grp_nm = cluster_id,dir_nm = paste0("../",output_dir))
    }
  }

  return(protein_oriexp_cluster)
}

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/Mfuzz_pipeline.R")
  Mfuzz_pipeline(exp = NULL,grp = NULL,
                 cluster_num = 10,min_acore = 0.3,
                 grp_nm = "Mfuzz",dir_nm = "Mfuzz")
}
