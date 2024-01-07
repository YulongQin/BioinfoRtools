### ---------------
###
### Create: qyl
### Date: 2023-04-01 22:29:00
### Email: 524583919@qq.com
### Function：
### - Convert the string separator
### -
###
### Update Log: 2023-04-01
###
### ---------------

#### ----- Function ------ ####
#' convert the string separator
#'
#' @name split_paste
#' @title convert the string separator
#' @description convert the string separator
#' @param string the string
#' @param split the split
#' @param collapse the collapse
#'
#' @return a string
#' @author Yulong Qin
#' @seealso \code{\link{shorten_names}}, \code{\link{split_paste}},
#' \code{\link{pubmed_search}}, \code{\link{manipulate_set}}
#'
#' @examples \dontrun{
#' a <- "a_b_c_d_e"
#' split_paste(a, split = "_", collapse = "-")
#' }
#'
#' @importFrom stringr str_split
#' @export
#'
split_paste <- function(string, split = " ", collapse = "\t") {
  str_list <- stringr::str_split(string, pattern = split)
  str <- paste(str_list[[1]], collapse = collapse)
  return(str)
}

#### ----- Examples ------ ####
if (F) {
  a <- "a_b_c_d_e"
  split_paste(a, split = "_", collapse = "-")
}
