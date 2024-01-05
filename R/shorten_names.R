### ---------------
###
### Create: qyl
### Date: 2023-04-01 22:29:00
### Email: 524583919@qq.com
### Function：
### - Shorten the description of excessively long GO enrichment
### -
###
### Update Log: 2023-04-01
###
### ---------------

#### ----- Function ------ ####
#' Shorten the description of excessively long GO enrichment
#'
#' @param x character vector
#' @param n_word number of words
#' @param n_char number of characters
#'
#' @return character vector
#' @author Yulong Qin
#' @seealso \code{\link{shorten_names}},\code{\link{split-paste}},
#' \code{\link{pubmed_search}},\code{\link{manipulate-set}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/shorten_names.R")
#' go_res$Des_short <- sapply(go_res$Description, shorten_names)
#' }
#'
#' @export
#'
shorten_names <- function(x, n_word = 6, n_char = 40) { # 单词数>n_word,单词长度>n_char
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char)) {
    if (nchar(x) > n_char) x <- substr(x, 1, n_char)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x, " ")[[1]]), n_word)],
      collapse = " "
    ), "...", sep = "")
    return(x)
  } else {
    return(x)
  }
}

#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/shorten_names.R")
  go_res$Des_short <- sapply(go_res$Description, shorten_names)
}
