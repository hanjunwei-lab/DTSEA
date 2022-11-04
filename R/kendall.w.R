#' @title Kendall's coefficient of concordance W
#' @description Computes the Kendall's coefficient of concordance.
#' @param raw A data frame or matrix contains n subjects * m raters.
#' @param correct Logical. Indicates whether the W should be corrected for ties
#' within raters.
#'
#' @return The resulting list consists of `title`, `kendall.w`, `chisq`, `df`,
#' `pval`, `report`.
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom tidyr drop_na
#' @importFrom magrittr divide_by multiply_by subtract
#' @importFrom stats pchisq na.omit var
#'
#' @examples
#' library(DTSEA)
#' library(tibble)
#'
#' # Load the data
#' data <- tribble(~x, ~y, ~z, 1,1,2,  5,6,5,  7,8,4, 2,3,2, 8,6,5)
#'
#' # Run Kendall's W
#' print(kendall.w(data)$report)
kendall.w <- function(raw, correct = TRUE) {
  # if the `raw` is dataframe: convert to matrix
  if (is.data.frame(raw)) {
    ratings <- drop_na(raw) %>%
      as.matrix()
  } else if (is.matrix(raw)) {
    ratings <- na.omit(raw)
  } else {
    stop("The input data should be in either a data frame or a matrix.")
  }

  if (ncol(ratings) < 2)
    stop("The input data should have more than two columns. ")

  ratings.rank <- apply(ratings, 2, rank)
  if (correct) {
    # correct for ties
    ties.adjust <- apply(ratings.rank, 2, function(percol) {
      ties <- table(percol) %>%
        { .[. > 1] } %>%
        as.numeric()
      return(sum(ties^3 - ties))
    }) %>% sum()
  } else {
    ties.adjust <- 0
  }
  coeff <- divide_by(
    e1 = 12 * multiply_by(var(apply(ratings.rank, 1, sum)), nrow(ratings) - 1),
    e2 = subtract(ncol(ratings)^2 * (nrow(ratings)^3 - nrow(ratings)),
                  ncol(ratings) * ties.adjust)
  )

  chisq <- multiply_by(
    e1 = coeff,
    e2 = nrow(ratings) * ncol(ratings) - nrow(ratings)
  )
  df <- nrow(ratings) - 1
  pval <- pchisq(chisq, df, lower.tail = FALSE)

  return(list(
    title = paste("Kendall's coefficient of concordance",
                  ifelse(correct, "Wt", "W")),
    kendall.w = coeff,
    chisq = chisq,
    df = df,
    pval = pval,
    report = paste0("Kendall's coefficient W = ", round(coeff, 3), ", p = ", signif(pval, 3))
  ))
}
