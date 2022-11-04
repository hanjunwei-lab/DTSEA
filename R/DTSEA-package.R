#' DTSEA
#' @title The Drug target set enrichment analysis (DTSEA)
#' @name DTSEA-package
#' @aliases DTSEA-package
#' @description The DTSEA implements a novel application to GSEA and extends the
#' adoption of GSEA.
#'
#' The Drug Target Set Enrichment Analysis (DTSEA) is a novel tool used to
#' identify the most effective drug set against a particular disease based on
#' the Gene Set Enrichment Analysis (GSEA).
#'
#' The central hypothesis of DTSEA is that the targets of potential candidates
#' for a specific disease (e.g., COVID-19) ought to be close to each other, or
#' at least not so far away from the disease. The DTSEA algorithm determines
#' whether a drug is potent for the chosen disease by the proximity between drug
#' targets and the disease-related genes. Under the central hypothesis of DTSEA,
#' the DTSEA consists of two main parts:
#'
#' 1. Evaluate the influence of the specific disease in the PPI network by the
#' random walk with restart algorithm. \cr
#' To evaluate the influence, we compute the disease-node distance by using the
#' random walk with restart (RwR) algorithm, then rank the nodes reversely.
#' 2. Evaluate the drug-disease associations based on GSEA. \cr
#' The GSEA approach is adopted in this part to identify whether candidate drug
#' targets are disease-related (top) or disease-unrelated (bottom) on the human
#' PPI list. The specific disease gene list is normalized by the median and is
#' set zero as the arbitrary cutoff point to classify the relations manually.
#'
#' In this package, we provide the example data, which is a small set of data to
#' demonstrate the usage and the main idea behind DTSEA.
#' We provide some extra data files, the real data we used in the DTSEA paper.
#' The supplementary package is now on the
#' [GitHub](https://github.com/hanjunwei-lab/DTSEAdata). Anyone can obtain this
#' package by the example code.
#'
#' @examples
#' # if (!"devtools" %in% as.data.frame(installed.packages())$Package)
#' #   install.packages("devtools")
#' # devtools::install_github("hanjunwei-lab/DTSEAdata")
#'
#' @importFrom utils data
#' @importFrom dplyr %>%
#
NULL

utils::globalVariables(c("."))

.onAttach <- function(libname, pkgname) {
#.onLoad <- function(libname, pkgname) {
  packageStartupMessage(
    "This package SHOULD NOT BE USED UNDER INTEL MATH KERNEL LIBRARY ON ANY OCCASION.
There is an avoidable but critical bug with Intel Math Kernel Library (MKL) on various operating systems.

======================================
For better performance, we recommend not using RStudio on Windows because RStudio cannot take advantage of the multi-core capabilities available on modern computers. ")
}
