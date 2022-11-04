#' Get extra data
#'
#' @param name Data name
#'
#' @return A list with the wanted data
#' @export
#' @keywords internal
#' @importFrom stats setNames
#' @importFrom utils data installed.packages
#' @importFrom dplyr %>%
#' @examples
#' # Do some stuff
#'
#' data <- get_data("ncbi_list")
#'
get_data <- function(name) {
  if (!"DTSEAdata" %in% as.data.frame(installed.packages())$Package) {
    name <- name[name %in% data(package = "DTSEAdata")$results[, "Item"]] %>%
      setNames(object = ., nm = .)

    return_list <- lapply(name, function(x) {
      e <- new.env()
      x <- data(list = x, envir = e)
      return(e[[x]])
    })

    if (length(return_list) == 1) {
      return_list <- return_list[[1]]
    }
  } else {
    return_list <- list()
  }

  return(return_list)
}
