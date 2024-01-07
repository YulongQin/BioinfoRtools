### ---------------
###
### Create: qyl
### Date: 2023-08-15 14:58:00
### Email: 524583919@qq.com
### Function：
### - miRNA and mRNA predict
### -
###
### Update Log: 2023-08-15
###
### ---------------

#### ----- Function ------ ####
#' miRNA and mRNA predict
#'
#' @title miRNA and mRNA predict
#' @description miRNA and mRNA predict
#' @param mRNA_genes the input mRNA genes
#' @param miRNA_genes the input miRNA genes
#' @param table_nm the table name
#' @param type_nm the type name
#' @param search_index whether to search pubmed
#' @param local_merge_index whether to merge local data
#' @param search_nm the search name
#' @param grp_nm the group name
#' @param dir_nm the dir name
#'
#' @return the predict result
#' @author Yulong Qin
#' @seealso \code{\link{WGCNA_pipeline}}, \code{\link{miRNA_mRNA}}, \code{\link{Mfuzz_pipeline}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/miRNA_mRNA.R")
#' miRNA_mRNA(mRNA_genes = NULL, # 输入up/down,软件不能区分mRNA-miRNA间的对应关系，只能通过输入的updown控制
#'           miRNA_genes = NULL, # 输入up/down
#'           table_nm = "validated", type_nm="miup_mdown",
#'           search_index = F,local_merge_index=F,
#'           search_nm = "[TIAB] AND (sperm[TIAB] OR semen[TIAB])",
#'           grp_nm = "group",dir_nm = "miRNA")
#'           }
#' @importFrom multiMiR get_multimir columns
#' @importFrom easyPubMed get_pubmed_ids fetch_pubmed_data
#' @export
#'
miRNA_mRNA <- function(mRNA_genes = NULL, # 输入up/down,软件不能区分mRNA-miRNA间的对应关系，只能通过输入的updown控制
                       miRNA_genes = NULL, # 输入up/down
                       table_nm = "validated", type_nm="miup_mdown",
                       search_index = F,local_merge_index=F,
                       search_nm = "[TIAB] AND (sperm[TIAB] OR semen[TIAB])",
                       grp_nm = "group",dir_nm = "miRNA"){
  ### 1.library ####
  # suppressMessages({
  #   library(magrittr)
  #   library(tidyverse)
  #   library(multiMiR) # get_multimir
  #   library(easyPubMed)
  # })
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm)
  # photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm)
  dir.create(output_dir,recursive = T)
  # dir.create(photo_dir,recursive = T)

  ### 2.mRNA predict ####
  if(!is.null(mRNA_genes) & is.null(miRNA_genes)){
    example1 <- get_multimir(
      org = "hsa",
      mirna = NULL,
      target = mRNA_genes,
      table = table_nm,
      summary = TRUE,
      predicted.cutoff.type = "p",
      predicted.cutoff = 10,
      use.tibble = TRUE
    )
    example1_result <- example1@data
    example1_sum <- example1@summary
    example1_res_fil <- AnnotationDbi::select(
      example1,
      keytype = "type",
      keys = "validated",
      columns = columns(example1)
    )
    {
      dup_index <- duplicated(example1_res_fil[, c("mature_mirna_id", "target_entrez")])
      unique_pairs <- example1_res_fil[!dup_index, ]
      unique_miRNA <- unique(example1_res_fil$mature_mirna_id)
      unique_target <- unique(example1_res_fil$target_symbol)
    }
    write.table(unique_pairs,
                file = paste0(output_dir,grp_nm,"_",type_nm,".txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)
    unique_pairs_fil <- unique_pairs
    unique_pairs_fil <- unique_pairs_fil[which(unique_pairs_fil$support_type=="Functional MTI"),]
    write.table(unique_pairs_fil,
                file = paste0(output_dir,grp_nm,"_",type_nm,"_valid.txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)
    tmp_gene <- unique_pairs_fil$mature_mirna_id # mRNA相关的miRNA
  }

  ### 3.miRNA predict ####
  if(!is.null(miRNA_genes) & is.null(mRNA_genes)){
    example1 <- get_multimir(
      org = "hsa",
      mirna = miRNA_genes,
      target = NULL,
      table = table_nm,
      summary = TRUE,
      predicted.cutoff.type = "p",
      predicted.cutoff = 10,
      use.tibble = TRUE
    )
    example1_result <- example1@data
    example1_sum <- example1@summary
    example1_res_fil <- AnnotationDbi::select(
      example1,
      keytype = "type",
      keys = "validated",
      columns = columns(example1)
    )
    {
      dup_index <- duplicated(example1_res_fil[, c("mature_mirna_id", "target_entrez")])
      unique_pairs <- example1_res_fil[!dup_index, ]
      unique_miRNA <- unique(example1_res_fil$mature_mirna_id)
      unique_target <- unique(example1_res_fil$target_symbol)
    }
    write.table(unique_pairs,
                file = paste0(output_dir,grp_nm,"_",type_nm,".txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)
    unique_pairs_fil <- unique_pairs
    unique_pairs_fil <- unique_pairs_fil[which(unique_pairs_fil$support_type=="Functional MTI"),]
    write.table(unique_pairs_fil,
                file = paste0(output_dir,grp_nm,"_",type_nm,"_valid.txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)
    tmp_gene <- unique_pairs_fil$target_symbol # miRNA相关的mRNA
  }

  ### 4.mRNA and miRNA predict ####
  if(!is.null(miRNA_genes) & !is.null(mRNA_genes)){
    ## 4.1 multimir ####
    example1 <- get_multimir(
      org = "hsa",
      mirna = miRNA_genes,
      target = mRNA_genes,
      table = table_nm,
      summary = TRUE,
      predicted.cutoff.type = "p",
      predicted.cutoff = 10,
      use.tibble = TRUE
    )
    example1_result <- example1@data
    example1_sum <- example1@summary
    example1_res_fil <- AnnotationDbi::select(
      example1,
      keytype = "type",
      keys = "validated",
      columns = columns(example1)
    )
    {
      dup_index <- duplicated(example1_res_fil[, c("mature_mirna_id", "target_entrez")])
      unique_pairs <- example1_res_fil[!dup_index, ]
      unique_miRNA <- unique(example1_res_fil$mature_mirna_id)
      unique_target <- unique(example1_res_fil$target_symbol)
    }
    write.table(unique_pairs,
                file = paste0(output_dir,grp_nm,"_",type_nm,".txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)
    unique_pairs_fil <- unique_pairs
    unique_pairs_fil <- unique_pairs_fil[which(unique_pairs_fil$support_type=="Functional MTI"),]
    write.table(unique_pairs_fil,
                file = paste0(output_dir,grp_nm,"_",type_nm,"_valid.txt"),
                sep = "\t",col.names = NA,row.names = T,quote = F)
    tmp_gene <- unique_pairs_fil$target_symbol # miRNA相关的mRNA

    ## 4.2 local merge ####
    if(local_merge_index){
      ## ref data
      mirtar <- data.table::fread(
        file = "F:/Bioinformatic_database/Homo_sapiens/miRTarBase_2021/miRTarBase_MTI.txt",
        header = T,sep = "\t") %>%
        as.data.frame()

      NcPath_mi <- data.table::fread(
        file = "F:/Bioinformatic_database/Homo_sapiens/NcPath_20230907/MTIs.txt",
        header = T,sep = "\t") %>%
        as.data.frame()

      ## local mirtarbase file
      miRNA_index <- (mirtar$miRNA %in% miRNA_genes & mirtar$`Target Gene` %in% mRNA_genes)
      miRNA_df <- mirtar[miRNA_index,] %>%
        filter(`Support Type`=="Functional MTI")
      write.table(miRNA_df,
                  file = paste0(output_dir,grp_nm,"_",type_nm,"_mirtar.txt"),
                  sep = "\t",col.names = NA,row.names = T,quote = F)

      ## local NcPath file
      miRNA_index <- (NcPath_mi$miRNA %in% miRNA_genes & NcPath_mi$`Target Gene`%in% mRNA_genes)
      miRNA_df <- NcPath_mi[miRNA_index,] %>%
        filter(`Support Type`=="Functional MTI")
      write.table(miRNA_df,
                  file = paste0(output_dir,grp_nm,"_",type_nm,"_NcPath.txt"),
                  sep = "\t",col.names = NA,row.names = T,quote = F)
    }
  }

  ### 5.pubmedabstract ####
  if(!is.null(miRNA_genes)){
    if(search_index){
      dir.create(paste0(output_dir,grp_nm,"_",type_nm,"_pubmed/"))
      if(!is_empty(tmp_gene)){
        for(i in 1:length(tmp_gene)){
          tryCatch({
            i_gene <- tmp_gene[i]
            my_query <- paste0(i_gene,search_nm)
            my_entrez_id <- get_pubmed_ids(
              my_query,
              api_key = "54d30eb11bfb539cfbab47558615ccf9e809"
            )
            my_abstracts_txt <- fetch_pubmed_data(my_entrez_id,
                                                  retstart = 0,
                                                  retmax = 20,
                                                  format = "abstract", # abstract格式
                                                  encoding = "UTF8")
            cat(my_abstracts_txt,
                file = paste0(output_dir,grp_nm,"_",type_nm,"_pubmed/",i_gene,".txt"),
                sep = "\n")
          },error=function(e){print(e)})
        }
      }
    }
  }
  return(tmp_gene)
}

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/miRNA_mRNA.R")
  miRNA_mRNA(mRNA_genes = NULL, miRNA_genes = NULL,
             table_nm = "validated", type_nm="miup_mdown",
             search_index = F,local_merge_index=F,
             search_nm = "[TIAB] AND (sperm[TIAB] OR semen[TIAB])",
             grp_nm = "group",dir_nm = "miRNA")
}
