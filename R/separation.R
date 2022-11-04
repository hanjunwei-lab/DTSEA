#' @title A measure of network separation
#' @description Calculates the separation of two sets of nodes on a network. The
#' metric is calculated as in Menche et al. (2015).
#' @param graph The input graph object. It should be either an igraph object or
#' an edge list matrix/data frame.
#' @param set_a The first gene set
#' @param set_b The second gene set
#'
#' @importFrom igraph graph_from_data_frame graph_from_edgelist is.igraph
#' simplify delete.vertices degree V
#' @importFrom dplyr %>%
#' @return The separation and distance measurement of the specified two modules.
#' @export
#'
#' @examples
#' library(DTSEA)
#'
#' # Load the data
#' data("random_graph", package = "DTSEA")
#'
#' # Compute the separation metric
#' separation <- separation(
#'   graph = random_graph,
#'   set_a = c("4", "6", "8", "13"),
#'   set_b = c("8", "9", "10", "15", "18")
#' )
#' cat(separation, "\n")
#'
separation <- function(graph, set_a, set_b) {
  # first, check the input `graph` object
  if (is.data.frame(graph)) {
    graph <- graph_from_data_frame(graph, directed = FALSE)
  } else if (is.matrix(graph)) {
    graph <- graph_from_edgelist(graph, directed = FALSE)
  } else if (!(is.igraph(graph) || length(V(graph)$name))) {
    stop("The input graph is not a valid named igraph object. ")
  }

  # slim the graph object
  graph <- simplify(graph, remove.loops = TRUE, remove.multiple = TRUE) %>%
    delete.vertices(which(degree(.) == 0))

  # check the existence of set a and set b in the graph

  between <- calculate_between(graph, set_a, set_b)
  within.a <- calculate_within(graph, set_a)
  within.b <- calculate_within(graph, set_b)
  return(between - (within.a + within.b) / 2)
}

# calculate between variance
#' @title Calculate between variance in network
#' @description No description
#' @param graph The input graph object. It should be either an igraph object or
#' an edge list matrix/data frame.
#' @param set_a The first gene set
#' @param set_b The second gene set
#'
#' @importFrom igraph shortest.paths
#' @importFrom dplyr %>%
#' @keywords internal
#' @return a positive number
calculate_between <- function(graph, set_a, set_b) {
  groups <- list(set_a = set_a, set_b = set_b)

  lapply(seq_along(groups) %>% rev(), function(group.i) {
    sapply(groups[[3 - group.i]], function(i) {
      shortest.paths(graph, v = groups[[group.i]], to = i) %>% min()
    })
  }) %>%
    unlist() %>%
    mean()
}

# calculate within variance
#' @title Calculate within variance
#' @description No description
#' @param graph The input graph object. It should be either an igraph object or
#' an edge list matrix / data frame.
#' @param given_set The first gene set
#'
#' @importFrom igraph shortest.paths
#' @importFrom dplyr %>%
#' @keywords internal
#' @return a positive number
calculate_within <- function(graph, given_set) {
  within <- sapply(given_set, function(i) {
    node_set <- given_set[!(given_set %in% i)]
    shortest.paths(graph, v = node_set, to = i) %>% min()
  }) %>%
    mean()
  return(within)
}
