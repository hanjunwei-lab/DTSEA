## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(dplyr)
library(magrittr)
library(DTSEA)

## -----------------------------------------------------------------------------
# Load the data
data("example_disease_list", package = "DTSEA")
data("example_drug_target_list", package = "DTSEA")
data("example_ppi", package = "DTSEA")

# Perform a simple DTSEA analysis using default optional parameters then sort
# the result dataframe by normalized enrichment scores (NES)
result <- DTSEA(network = example_ppi,
                disease = example_disease_list,
                drugs = example_drug_target_list, verbose = FALSE
) %>%
  arrange(desc(NES))
head(result)

## -----------------------------------------------------------------------------
select(result, -leadingEdge) %>%
  arrange(desc(NES)) %>%
  filter(NES > 0 & pval < .05) %>% head()

## -----------------------------------------------------------------------------
fgsea::plotEnrichment(
  pathway = example_drug_target_list %>%
    filter(drug_id == slice(result, 1)$drug_id) %>%
    pull(gene_target),
  stats = random.walk(network = example_ppi,
                      p0 = calculate_p0(nodes = example_ppi,
                                        disease = example_disease_list)
                      )
)

## -----------------------------------------------------------------------------
# Calculate p0
p0 <- calculate_p0(nodes = example_ppi, disease = example_disease_list)

# Then perform random walk
random.walk(network = example_ppi, p0 = p0, verbose = FALSE) %>%
  head()

## ----echo = FALSE-------------------------------------------------------------
# Imagine there are three prediction results with ten samples
x <- runif(10, min = 0, max = 5)
y <- runif(10, min = 0, max = 5)
z <- sqrt(x + y) + runif(10, min = -1, max = 1)
data <- data.frame(x, y, z)

## -----------------------------------------------------------------------------
# Just report the results
kendall.w(data)$report

# Or just report the alpha
cronbach.alpha(data)

## ---- eval=FALSE--------------------------------------------------------------
#  # Load the data
#  data("example_disease_list", package = "DTSEA")
#  data("example_drug_target_list", package = "DTSEA")
#  data("example_ppi", package = "DTSEA")
#  
#  # set up environment
#  
#  single.core <- function() {
#   suppressWarnings(capture.output(DTSEA(network = example_ppi,
#                                         disease = example_disease_list,
#                                         drugs = example_drug_target_list,
#                                         nproc = 0)))
#    NULL
#  }
#  
#  dual.core <- function() {
#   suppressWarnings(capture.output(DTSEA(network = example_ppi,
#                                         disease = example_disease_list,
#                                         drugs = example_drug_target_list,
#                                         nproc = 10)))
#    NULL
#  }
#  
#  system.time(single.core()) - system.time(dual.core())
#  

## ---- eval=FALSE--------------------------------------------------------------
#  if (!"devtools" %in% as.data.frame(installed.packages())$Package)
#    install.packages("devtools")
#  devtools::install_github("hanjunwei-lab/DTSEAdata")

