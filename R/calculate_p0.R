#' @title Function to calculate the p0 vector used in Random Walk with Restart
#' (RwR)
#' @description  The function provides a reliable approach to generating a p0
#' vector.
#'
#' @param nodes The `nodes` variable can either accept the igraph object or the
#' nodes vector.
#' @param disease The `disease` variable must specify the disease-affected
#' nodes in a short vector.
#'
#' @importFrom tidyr replace_na
#' @importFrom dplyr %>% mutate distinct filter pull
#' @importFrom igraph is.igraph V
#'
#' @return The resulting p0 vector.
#' @export
#'
#' @examples
#' library(DTSEA)
#' library(dplyr)
#'
#' # Load the data
#' data("example_disease_list", package = "DTSEA")
#' data("example_drug_target_list", package = "DTSEA")
#' data("example_ppi", package = "DTSEA")
#'
#' # Compute the p0 vector
#' p0 <- calculate_p0(nodes = example_ppi, disease = example_disease_list)
#'
#' # You can decrease the order of the p0 to get the most affected nodes.
#' p0 <- sort(p0, decreasing = TRUE) %>%
#'   names() %>%
#'   head(10)
#'
#' # If you have obtained the supplemental data, then you can compute the p0
#' # in the real data set
#'
#' # supp_data <- get_data(c("graph", "disease_related"))
#' # p0 <- calculate_p0(nodes = supp_data[["graph"]],
#' #                    disease = supp_data[["disease_related"]])
#'
calculate_p0 <- function(nodes, disease) {
  # add some nonsense statement to pass the test
  symbol <- NULL

  # check the node
  if (is.igraph(nodes)) {
    nodes <- V(nodes)$name
  }
  if (!is.vector(nodes)) {
    stop("Not a valid graph. ")
  } else if (length(nodes) < 4) {
    stop("Not an appropriate graph. ")
  } else if (!(length(nodes) > length(disease))) {
    stop("The disease vector is too long for DTSEA to proceed. ")
  } else if (!sum(disease %in% nodes)) {
    stop("Maybe there is a wrong `disease` vector provided. ")
  }

  weight <- data.frame(
    symbol = nodes,
    weight = as.numeric(nodes %in% disease)
  ) %>%
    mutate(weight = replace_na(weight, 0)) %>%
    distinct() %>%
    filter(symbol %in% nodes) %>%
    pull(weight, symbol)

  return(weight)
}
