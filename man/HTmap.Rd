% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HTmap.R
\name{HTmap}
\alias{HTmap}
\title{Draw heat maps according to groups}
\usage{
HTmap(
  countdata,
  group,
  log2_index = T,
  scale_index = "row",
  cluster_cols = F,
  cluster_rows = T,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_colnames = T,
  show_rownames = F
)
}
\arguments{
\item{countdata}{the count data matrix}

\item{group}{the group information}

\item{log2_index}{whether to log2 the data}

\item{scale_index}{whether to scale the data}

\item{cluster_cols}{whether to cluster the columns}

\item{cluster_rows}{whether to cluster the rows}

\item{clustering_distance_rows}{the distance method for clustering rows}

\item{clustering_distance_cols}{the distance method for clustering columns}

\item{clustering_method}{the clustering method}

\item{show_colnames}{whether to show the column names}

\item{show_rownames}{whether to show the row names}
}
\value{
plot
}
\description{
Draw heat maps according to groups
}
\details{
Draw heat maps according to groups
}
\examples{
\dontrun{
source(file = "F:/Bioinformatic_repository/02_R/code_R/A_Script_Function/HTmap.R")
df <- iris[1:4] \%>\%
  t() \%>\%
  as.data.frame()
HTmap(df, iris$Species)
}
}
\seealso{
\code{\link{PCA}}, \code{\link{HTmap}}
}
\author{
Yulong Qin
}
