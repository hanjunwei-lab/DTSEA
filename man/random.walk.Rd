% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random.walk.R
\name{random.walk}
\alias{random.walk}
\title{Function to implement Random Walk with Restart (RwR) algorithm on the
input graph}
\usage{
random.walk(
  network,
  p0,
  edge_weight = FALSE,
  gamma = 0.7,
  threshold = 1e-10,
  pt.post.processing = "log",
  pt.align = "median",
  verbose = FALSE
)
}
\arguments{
\item{network}{The input graph object. It should be either an igraph object
or an edge list matrix / data frame.}

\item{p0}{The starting vector on time t0.}

\item{edge_weight}{Logical to indicate whether the input graph contains
weight information.}

\item{gamma}{The restart probability used for RwR. The \code{gamma} takes the
value from \strong{0} to \strong{1}, controlling the probability that a node would go
back to its starting node.}

\item{threshold}{The threshold used for RwR. The \code{threshold} indicates the
stabilization status, which is a stopping criterion of RwR.}

\item{pt.post.processing}{The way to scale the \code{pt} vector. It can be
'none', 'zscore', and 'log'.}

\item{pt.align}{The way to normalize the output \code{pt} vector. It can be
\strong{mean} to manually cut the up- and down-regulated genes, \strong{median} to
avoid the influence of the distribution shape, or \strong{none} for no
normalization.}

\item{verbose}{Show the progress of the calculation.}
}
\value{
\code{pt} vector
}
\description{
Function \code{random.walk} is supposed to implement the original
Random Walk with Restart (RwR) on the input graph. If the seeds (i.e.,  a set
of starting nodes) are given, it intends to calculate the affinity score of
all nodes in the graph to the seeds.
}
\examples{
library(DTSEA)

# Load the data
data("example_disease_list", package = "DTSEA")
data("example_drug_target_list", package = "DTSEA")
data("example_ppi", package = "DTSEA")

# Perform random walk
p0 <- calculate_p0(nodes = example_ppi, disease = example_disease_list)
pt <- random.walk(network = example_ppi, p0 = p0)

# Perform GSEA analysis
# ....

# If you have obtained the supplemental data, then you can do random walk
# with restart in the real data set

# supp_data <- get_data(c("graph", "disease_related", "example_ppi"))
# p0 <- calculate_p0(nodes = supp_data[["graph"]],
#                    disease = supp_data[["disease_related"]])
# pt <- random.walk(network = supp_data[["example_ppi"]],
#                   p0 = p0)

}
