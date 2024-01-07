### ---------------
###
### Create: qyl
### Date: 2023-08-15 14:58:00
### Email: 524583919@qq.com
### Function：
### - serach pubmed by keywords
### -
###
### Update Log: 2023-08-15
###
### ---------------

#### ----- Function ------ ####
#' serach pubmed by keywords
#'
#' @title serach pubmed by keywords
#' @description serach pubmed by keywords
#' @param search_words keywords
#' @param fixed_words fixed words
#' @param sep the sep of keywords and fixed words
#' @param api_key the api_key of pubmed
#' @param grp_nm the group name
#' @param dir_nm the dir name
#'
#' @return my_abstracts_txt
#' @author Yulong Qin
#' @seealso \code{\link{shorten_names}}, \code{\link{split_paste}},
#' \code{\link{pubmed_search}}, \code{\link{manipulate_set}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/pubmed_search.R")
#' pubmed_search(
#'   search_words = NULL, fixed_words = NULL, sep = "AND",
#'   api_key = "54d30eb11bfb539cfbab47558615ccf9e809",
#'   grp_nm = "pubmed_search", dir_nm = "pubmed_search"
#' )
#' }
#'
#' @importFrom easyPubMed get_pubmed_ids fetch_pubmed_data
#' @export
#'
pubmed_search <- function(search_words = NULL, fixed_words = NULL, sep = "AND",
                          api_key = "54d30eb11bfb539cfbab47558615ccf9e809",
                          grp_nm = "pubmed_search", dir_nm = "pubmed_search") {
  ### 1.library ####
  # suppressMessages({
  #   library(easyPubMed)
  # })
  output_dir <- paste0("./outputdata/", dir_nm, "/", grp_nm)
  # photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm)
  dir.create(output_dir, recursive = T)
  # dir.create(photo_dir,recursive = T)

  ### 2.pubmed_search ####
  cat(paste0("pubmed_search start at ", Sys.time(), "\n"))
  for (i_search in 1:length(search_words)) {
    cat(paste0("\n", "i_search ", i_search, " start at ", Sys.time(), "\n"))
    search_nm <- search_words[i_search]

    my_query <- paste(search_nm, fixed_words, sep = paste0(" ", sep, " "))
    cat("Search strategy: ", my_query, "\n", sep = "")

    # 核心代码，经常卡顿
    cat(paste0("get_pubmed_ids start at ", Sys.time(), "\n"))
    my_entrez_id <- get_pubmed_ids(
      my_query,
      api_key = api_key
    )

    cat(paste0("fetch_pubmed_data start at ", Sys.time(), "\n"))
    my_abstracts_txt <- fetch_pubmed_data(
      my_entrez_id,
      retstart = 0,
      retmax = 20,
      format = "abstract", # abstract格式
      encoding = "UTF8"
    )
    cat(my_abstracts_txt,
      file = paste0(output_dir, "/", search_nm, ".txt"),
      sep = "\n"
    )
  }
  cat(paste0("\n", "pubmed_search end at ", Sys.time(), "\n"))

  return(invisible(my_abstracts_txt))
}

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/pubmed_search.R")
  pubmed_search(
    search_words = NULL, fixed_words = NULL, sep = "AND",
    api_key = "54d30eb11bfb539cfbab47558615ccf9e809",
    grp_nm = "pubmed_search", dir_nm = "pubmed_search"
  )
}
