### ---------------
###
### Create: qyl
### Date: 2023-12-13 14:28:00
### Email: 524583919@qq.com
### Function：
### - Computes the intersection, union, and difference sets between sets
### -
###
### Update Log: 2023-12-13
###
### ---------------

#### ----- Function ------ ####
#' computes the intersection, union, and difference sets between sets
#'
#' @name manipulate_set
#' @title computes the intersection, union, and difference sets between sets
#' @description computes the intersection, union, and difference sets between sets
#' @param ... the input vectors
#'
#' @family manipulate_set
#'
#' @return a list of vectors
#' @author Yulong Qin
#' @seealso \code{\link{shorten_names}}, \code{\link{split_paste}},
#' \code{\link{pubmed_search}}, \code{\link{manipulate_set}}
#'
#' @examples \dontrun{
#' source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/manipulate_set.R")
#' get_intersection(a, b, c, d)
#' }
#'
NULL

#' @export
#' @rdname manipulate_set
get_intersection <- function(...) {
  vectors <- list(...)
  intersection <- vectors[[1]]

  for (i in 2:length(vectors)) {
    intersection <- intersect(intersection, vectors[[i]])
  }

  return(intersection)
}

#' @export
#' @rdname manipulate_set
get_union <- function(...) {
  vectors <- list(...)
  union_set <- unique(unlist(vectors))

  return(union_set)
}

#' @export
#' @rdname manipulate_set
get_intersection_list <- function(...) {
  vectors <- list(...)
  n <- length(vectors)
  combinations <- list()

  # Calculate intersection for each pair of vectors
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      index <- length(combinations) + 1
      intersection <- intersect(vectors[[i]], vectors[[j]])
      combinations[[index]] <- intersection
      names(combinations)[index] <- paste0("index_", i, "_", j)
    }
  }

  return(combinations)
}

#' @export
#' @rdname manipulate_set
get_union_list <- function(...) {
  vectors <- list(...)
  n <- length(vectors)
  combinations <- list()

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      index <- length(combinations) + 1
      intersection <- union(vectors[[i]], vectors[[j]])
      combinations[[index]] <- intersection
      names(combinations)[index] <- paste0("index_", i, "_", j)
    }
  }

  return(combinations)
}


#### ----- Examples ------ ####
if (F) {
  source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/manipulate_set.R")
  get_intersection(a, b, c, d)
}
