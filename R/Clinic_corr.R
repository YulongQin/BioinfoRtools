### ---------------
###
### Create: qyl
### Date: 2023-10-27 15:48:00
### Email: 524583919@qq.com
### Function：
### - Calculate correlations between continuous and discrete variables
### -
###
### Update Log: 2023-10-27
###
### ---------------

#### ----- Function ------ ####
#' Calculate correlations between continuous and discrete variables
#'
#' @title Calculate correlations between continuous and discrete variables
#' @description Calculate correlations between continuous and discrete variables
#'
#' @param continuous the continuous data.frame
#' @param continuous2 the continuous2 data.frame
#' @param binary the binary data.frame
#' @param grp_nm the name of the group
#' @param dir_nm the name of the directory
#' @param skip_cont whether to skip continuous to continuous
#' @param skip_bina whether to skip binary to binary
#' @param skip_cont_bina whether to skip continuous to binary
#'
#' @return NULL
#' @author Yulong Qin
#' @seealso \code{\link{GO_KEGG}}, \code{\link{GSEA_GO_KEGG}}, \code{\link{GSVA_ssGSEA}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/Clinic_corr.R")
#' Clinic_corr(
#'   continuous = NULL, binary = NULL,
#'   grp_nm = "Clinic_corr1", dir_nm = "Clinic_corr",
#'   skip_cont = F, skip_bina = F
#' )
#' }
#' @importFrom psych corr.test
#' @importFrom stats fisher.test median wilcox.test
#' @export
#'
Clinic_corr <- function(continuous=NULL,continuous2=NULL,binary=NULL, # 行为样本，列为临床变量
                        grp_nm = "Clinic_corr1",dir_nm = "Clinic_corr",
                        skip_cont=F,skip_bina=F,skip_cont_bina=F){

  ### 0.check ####
  if(!is.null(continuous) & !is.null(binary)){
    if(nrow(continuous)!=nrow(binary)){
      stop("The continuous and binary need to have the same number of rows")
    }
  }

  ### 1.library ####
  # suppressMessages({
  #   library(tidyverse)
  #   library(psych)
  # })
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm) # dir_nm下有子目录，grp_nm一般为gene或其他指定值
  # photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm)
  dir.create(output_dir,recursive = T)
  # dir.create(photo_dir,recursive = T)

  ### 2.continuous to continuous ####
  ## 2.1 continuous to continuous
  if(!is.null(continuous) & !skip_cont){
    print("### continuous to continuous ####")

    n=1
    p_all_spearman <- as.data.frame(matrix(numeric(0),ncol = 5))
    colnames(p_all_spearman) <- c("var1","var2","length","r","pvalue")
    for (i in 1:(ncol(continuous)-1)) {
      for (j in (i+1):ncol(continuous)) {
        vec_1 <- continuous[,i] %>% as.numeric()
        vec_2 <- continuous[,j] %>% as.numeric()
        tryCatch({
          df <- corr.test(vec_1,vec_2,
                          use = "complete",method = "spearman",
                          adjust = "BH")
          if(is.na(df$p.adj)){next}
          p_all_spearman[n,1] <- colnames(continuous)[i]
          p_all_spearman[n,2] <- colnames(continuous)[j]
          p_all_spearman[n,3] <- df$n
          p_all_spearman[n,4] <- df$r
          p_all_spearman[n,5] <- df$p.adj
          n=n+1
        },
        error=function(e)e)
      }
    }
    p_all_spearman$padj <- p.adjust(p_all_spearman$pvalue,method = "BH")
    p_fil_spearman <- p_all_spearman %>%
      filter(padj<0.05) %>%
      filter(r>0.4 | r < -0.4) %>%
      filter(length>100)
    write.table(p_all_spearman,
                file = paste0(output_dir,'/p_all_spearman.txt'),
                sep = "\t",row.names = F,col.names = T,quote = F)
    write.table(p_fil_spearman,
                file = paste0(output_dir,'/p_fil_spearman.txt'),
                sep = "\t",row.names = F,col.names = T,quote = F)
  }

  ## 2.2 continuous to continuous2
  if(!is.null(continuous) & !is.null(continuous2)){
    print("### continuous to continuous2 ####")

    df <- corr.test(continuous,continuous2,
                    use = "pairwise",method = "spearman",
                    adjust = "BH")
    res_r <- df$r
    res_p <- df$p
    res_padj <- df$p.adj
    res_n <- df$n
    write.table(res_r,
                file = paste0(output_dir,'/continuous_continuous2_r.txt'),
                sep = "\t",row.names = T,col.names = NA,quote = F)
    write.table(res_p,
                file = paste0(output_dir,'/continuous_continuous2_p.txt'),
                sep = "\t",row.names = T,col.names = NA,quote = F)
    write.table(res_padj,
                file = paste0(output_dir,'/continuous_continuous2_padj.txt'),
                sep = "\t",row.names = T,col.names = NA,quote = F)
    write.table(res_n,
                file = paste0(output_dir,'/continuous_continuous2_n.txt'),
                sep = "\t",row.names = T,col.names = NA,quote = F)
  }

  ### 3.continuous to binary ####
  if(!is.null(continuous) & !is.null(binary) & !skip_cont_bina){
    print("### continuous to binary ####")

    n=1
    p_all_wilcox <- as.data.frame(matrix(numeric(0),ncol = 8))
    colnames(p_all_wilcox) <- c("var1","var2","length","median_0","median_1",
                                "Change","W","pvalue")
    for (i in 1:ncol(continuous)) {
      for (j in 1:ncol(binary)) {
        vec_1 <- continuous[,i] %>% as.numeric()
        vec_2 <- binary[,j] %>% as.numeric()
        tryCatch({
          df <- wilcox.test(vec_1~vec_2)
          mean_vec <- tapply(vec_1,vec_2,median,na.rm=T)
          if(is.na(df$p.value)){next}
          p_all_wilcox[n,1] <- colnames(continuous)[i]
          p_all_wilcox[n,2] <- colnames(binary)[j]
          p_all_wilcox[n,3] <- nrow(na.omit(cbind(vec_1,vec_2)))
          p_all_wilcox[n,4:5] <- mean_vec
          p_all_wilcox[n,6] <- ifelse(mean_vec[2]>mean_vec[1],"Up",
                                      ifelse(mean_vec[2]<mean_vec[1],"Down","Equal"))
          p_all_wilcox[n,7] <- df[["statistic"]]
          p_all_wilcox[n,8] <- df$p.value
          n=n+1
        },
        error=function(e)e)
      }
    }
    p_all_wilcox$padj <- p.adjust(p_all_wilcox$pvalue,method = "BH")
    p_fil_wilcox <- p_all_wilcox %>%
      filter(padj<0.05) %>%
      filter(length>100)
    write.table(p_all_wilcox,
                file = paste0(output_dir,'/p_all_wilcox.txt'),
                sep = "\t",row.names = F,col.names = T,quote = F)
    write.table(p_fil_wilcox,
                file = paste0(output_dir,'/p_fil_wilcox.txt'),
                sep = "\t",row.names = F,col.names = T,quote = F)
  }

  ### 4.binary to binary ####
  if(!is.null(binary) & !skip_bina){
    print("### binary to binary ####")

    n=1
    p_all_fisher <- as.data.frame(matrix(numeric(0),ncol = 10))
    colnames(p_all_fisher) <- c("var1","var2","length","table_00","table_01","table_10","table_11",
                                "Change","odds_ratio","pvalue")
    for (i in 1:(ncol(binary)-1)) {
      for (j in (i+1):ncol(binary)) {
        vec_1 <- binary[,i] %>% as.numeric()
        vec_2 <- binary[,j] %>% as.numeric()
        tryCatch({
          df <- fisher.test(vec_1,vec_2)
          tab_vec <- table(vec_1,vec_2)
          if(is.na(df$p.value)){next}
          p_all_fisher[n,1] <- colnames(binary)[i]
          p_all_fisher[n,2] <- colnames(binary)[j]
          p_all_fisher[n,3] <- nrow(na.omit(cbind(vec_1,vec_2)))
          p_all_fisher[n,4] <- tab_vec[1,1]
          p_all_fisher[n,5] <- tab_vec[1,2]
          p_all_fisher[n,6] <- tab_vec[2,1]
          p_all_fisher[n,7] <- tab_vec[2,2]
          p_all_fisher[n,8] <- ifelse(df$estimate>1,"Positive",
                                      ifelse(df$estimate<1,"Negative","Equal"))
          p_all_fisher[n,9] <- df$estimate
          p_all_fisher[n,10] <- df$p.value
          n=n+1
        },
        error=function(e)e)
      }
    }
    p_all_fisher$padj <- p.adjust(p_all_fisher$pvalue,method = "BH")
    p_fil_fisher <- p_all_fisher %>%
      filter(padj<0.05) %>%
      filter(length>100)
    write.table(p_all_fisher,
                file = paste0(output_dir,'/p_all_fisher.txt'),
                sep = "\t",row.names = F,col.names = T,quote = F)
    write.table(p_fil_fisher,
                file = paste0(output_dir,'/p_fil_fisher.txt'),
                sep = "\t",row.names = F,col.names = T,quote = F)
  }
}

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/Clinic_corr.R")
  Clinic_corr(continuous=NULL,continuous2=NULL,binary=NULL,
              grp_nm = "Clinic_corr1",dir_nm = "Clinic_corr",
              skip_cont=F,skip_bina=F)
}
