#' @title Cronbach's alpha
#' @description Computes Cronbach's alpha
#' @param data A data frame or matrix contains n subjects * m raters.
#'
#' @return The Cronbach's alpha (unstandardized)
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom tidyr drop_na
#' @importFrom magrittr raise_to_power multiply_by
#' @importFrom stats na.omit var sd
#'
#' @examples
#' library(DTSEA)
#' library(tibble)
#'
#' # Load the data
#' data <- tribble(~x, ~y, ~z, 1, 1, 2, 5, 6, 5, 7, 8, 4, 2, 3, 2, 8, 6, 5)
#'
#' # Run Cronbach's alpha
#' cat(cronbach.alpha(data))
#'
cronbach.alpha <- function(data) {
  # data should be in either data frame or matrix
  if (is.data.frame(data)) {
    data <- drop_na(data) %>%
      as.matrix()
  } else if (is.matrix(data)) {
    data <- na.omit(data)
  } else {
    stop("The input data should be in either a data frame or a matrix.")
  }

  # check the column of data
  if (ncol(data) < 2) {
    stop("The input data should have more than two columns. ")
  }

  variance.total <- apply(data, 1, sum) %>%
    var()
  variance.byrow <- apply(data, 2, sd) %>%
    raise_to_power(2) %>%
    sum()
  alpha <- multiply_by(
    e1 = 1 + 1 / (ncol(data) - 1),
    e2 = (1 - variance.byrow / variance.total)
  )

  return(alpha)
}
