#' @title Function to implement Random Walk with Restart (RwR) algorithm on the
#' input graph
#' @description Function `random.walk` is supposed to implement the original
#' Random Walk with Restart (RwR) on the input graph. If the seeds (i.e.,  a set
#' of starting nodes) are given, it intends to calculate the affinity score of
#' all nodes in the graph to the seeds.
#' @param network The input graph object. It should be either an igraph object
#' or an edge list matrix / data frame.
#' @param p0 The starting vector on time t0.
#' @param edge_weight Logical to indicate whether the input graph contains
#' weight information.
#' @param gamma The restart probability used for RwR. The `gamma` takes the
#' value from **0** to **1**, controlling the probability that a node would go
#' back to its starting node.
#' @param threshold The threshold used for RwR. The `threshold` indicates the
#' stabilization status, which is a stopping criterion of RwR.
#' @param pt.post.processing The way to scale the `pt` vector. It can be
#' 'none', 'zscore', and 'log'.
#' @param pt.align The way to normalize the output `pt` vector. It can be
#' **mean** to manually cut the up- and down-regulated genes, **median** to
#' avoid the influence of the distribution shape, or **none** for no
#' normalization.
#' @param verbose Show the progress of the calculation.
#'
#' @importFrom dplyr %>% between arrange pull left_join mutate
#' @importFrom tibble tibble
#' @importFrom igraph graph_from_data_frame graph_from_edgelist is.igraph simplify induced_subgraph V components as_adjacency_matrix
#' @importFrom stats median
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @return `pt` vector
#' @export
#'
#' @examples
#' library(DTSEA)
#'
#' # Load the data
#' data("example_disease_list", package = "DTSEA")
#' data("example_drug_target_list", package = "DTSEA")
#' data("example_ppi", package = "DTSEA")
#'
#' # Perform random walk
#' p0 <- calculate_p0(nodes = example_ppi, disease = example_disease_list)
#' pt <- random.walk(network = example_ppi, p0 = p0)
#'
#' # Perform GSEA analysis
#' # ....
#'
#' # If you have obtained the supplemental data, then you can do random walk
#' # with restart in the real data set
#'
#' # supp_data <- get_data(c("graph", "disease_related", "example_ppi"))
#' # p0 <- calculate_p0(nodes = supp_data[["graph"]],
#' #                    disease = supp_data[["disease_related"]])
#' # pt <- random.walk(network = supp_data[["example_ppi"]],
#' #                   p0 = p0)
#'
random.walk <- function(network, p0,
                        edge_weight = FALSE, gamma = 0.7, threshold = 1e-10,
                        pt.post.processing = "log", pt.align = "median",
                        verbose = FALSE) {

  # add some nonsense statement to pass the test
  gene <- pt <- rpt <- NULL

  # First, we need to simplify the whole graph
  if (is.data.frame(network)) {
    graph <- graph_from_data_frame(network, directed = FALSE)
  } else if (is.matrix(network)) {
    graph <- graph_from_edgelist(network, directed = FALSE)
  } else if (!is.igraph(network)) {
    stop("Not a valid igraph object. ")
  } else {
    graph <- network
  }

  nodes <- V(graph)$name

  # The weight is the disease gene in the human metabolic network

  # Calculate similarity score using random walk with restart (RWR) algorithm

  # RwR begins

  W <- graph %>%
    as_adjacency_matrix(attr = switch(
      edge_weight,
      "weight",
      NULL
    ), type = "both") %>%
    as.matrix()

  # p0 must be named vector
  if (names(p0) %>% is.null() %>% sum()) {
    stop("The given vector p0 must be a named vector. ")
  } else if (names(p0) %>% is.na() %>% sum()) {
    stop("The names of the given vector p0 contain more than one NAs. ")
  } else if (sum(names(p0) == "")) {
    stop("The names of the given vector p0 contain more than one blank. ")
  } else if (table(unname(p0)) %>%
             names() %>%
             length() %>%
             { . < 2 }) {
    # p0 must have more than 1 levels
    stop("The given vector p0 is homogeneity. ")
  } else if (max(p0) > 1 || min(p0) < 0 || sum(p0) <= 0) {
    # the initial probabilities must be in the range of 0 to 1
    stop("The given vector p0 is corrupted. ")
  } else if (length(p0) != nrow(W)) {
    stop("The length of p0 is mismatched with the matrix W. ")
  } else if (!between(gamma, 0, 1)) {
    # the gamma must be in the range of 0 to 1
    stop("The restart probability gamma must be in the range of 0 to 1. ")
  } else if (threshold > 1e-5) {
    warning(paste(
      "The stop threshold is too large. Consider a smaller threshold",
      threshold / 1e5, ". "
    ))
  } else if (length(nodes) < 5) {
    warning(paste("The network size is too small. "))
  } else if (!(pt.post.processing %in% c("none", "log", "median"))) {
    warning("The `pt.post.processing` is set to be `none`. ")
    pt.post.processing <- "none"
  } else if (!(pt.align %in% c("none", "median", "mean"))) {
    warning("The `pt.align` is set to be `none`. ")
    pt.align <- "none"
  }

  invisible(gc())

  # Get the order of the matrix W
  p0 <- tibble(gene = names(p0), p0 = unname(p0)) %>%
    left_join(tibble(gene = rownames(W), order = seq_along(rownames(W))),
      by = "gene"
    ) %>%
    arrange(order) %>%
    pull(p0, gene)

  p0 <- t(p0) / sum(p0)
  PT <- p0

  delta <- 1

  Ng <- dim(W)[2]
  for (i in seq_len(Ng)) {
    sumr <- sum(W[i, ])
    if (sumr == 0) {
      W[i, ] <- numeric(length = length(W[i, ]))
    } else if (sumr > 0) {
      W[i, ] <- W[i, ] / sum(W[i, ])
    }
  }
  W <- t(W)

  if (verbose)
    pb <- txtProgressBar(min = 0, max = -log10(threshold), style = 3)
  while (delta > threshold) {
    PT1 <- (1 - gamma) * W
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma * p0)
    PT4 <- t(PT2) + PT3
    delta <- sum(abs(PT4 - PT))
    PT <- PT4
    if (verbose)
      setTxtProgressBar(pb, -log10(delta))
  }
  rwr.pt <- t(PT) %>%
    drop() %>%
    sort(decreasing = TRUE) %>%
    tibble(gene = names(.), pt = .) %>%
    mutate(
      pt = switch(pt.post.processing,
        "none" = pt,
        "log" = log(pt),
        "zscore" = scale(pt)
      ),
      rpt = switch(pt.align,
        "none" = pt,
        "mean" = scale(pt, center = TRUE, scale = FALSE),
        "median" = (pt - median(pt))
      )
    ) %>%
    pull(rpt, gene)

  return(rwr.pt)
}
