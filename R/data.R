#' Meta-data of samples
#'
#' Meta-data of a sample of 7 samples and 2 clinical information samples.
#'
#' @format A data frame with 7 rows and 2 variables:
#' \describe{
#'   \item{condition}{condition}
#'   \item{type}{type}
#' }
"coldata"

#' A transcriptome expression
#'
#' A transcriptome expression profile of 14599 genes in 7 samples.
#'
#' @format A data frame with 14599 rows and 7 variables:
#' \describe{
#'   \item{treated1}{treated1}
#'   \item{treated2}{treated2}
#'   \item{treated3}{treated3}
#'   \item{untreated1}{untreated1}
#'   \item{untreated2}{untreated2}
#'   \item{untreated3}{untreated3}
#'   \item{untreated4}{untreated4}
#' }
"countdata"

#' A gene set
#'
#' A gene set of 80 DEGs.
#'
"genes"
