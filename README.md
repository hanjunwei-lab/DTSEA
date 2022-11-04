
# Drug set enrichment analysis (DTSEA)

<!-- badges: start -->
<!-- badges: end -->
**DTSEA** contains the basic R functions and sample data for running the DTSEA algorithm. After installing and loading the package, users will be able to explore the framework of DTSEA. 

## More about DTSEA

The Drug Target Set Enrichment Analysis (DTSEA) is a novel tool used to identify the most effective drug set against a particular disease based on the Gene Set Enrichment Analysis (GSEA). It assumes the most effective drugs are those with a closer affinity to the specified disease. 

## Getting started

### Step 1. Prerequisites for installation

You should ensure that you have the necessary system dependencies configured. 

For Debian-based Linux system (Debian buster or bullseye), build-essential should be installed at first: 

```
sudo apt install -y build-essential
```

For Windows (8.1 / 10 / 11): [Rtools](https://cran.r-project.org/bin/windows/Rtools/) should be installed to the system path. 

The latest base R is recommended. The compatibility of the earlier version (v4.0.x) is under evaluation. 

### Step 2. Install the package

The dependency `fgsea` and `BiocParallel` are unavailable on the CRAN but available on [BioConductor](https://www.bioconductor.org/). So we need to install the BiocManager manually. 

``` r
if (!"BiocManager" %in% as.data.frame(installed.packages())$Package)
  install.packages("BiocManager")
BiocManager::install(c("fgsea", "BiocParallel"))
```
Then you can install the development version of DTSEA from [GitHub](https://github.com/) with:

``` r
if (!"devtools" %in% as.data.frame(installed.packages())$Package)
  install.packages("devtools")
devtools::install_github("hanjunwei-lab/DTSEA")
```

### Examples

Below is a basic example that shows how to solve a common problem:

``` r
library(DTSEA)

# Load the data
data(example_disease_list)
data(example_drug_target_list)
data(example_ppi)

# Perform a simple DTSEA analysis

result <- DTSEA(
    network = example_ppi,
    disease = example_disease_list,
    drugs = example_drug_target_list
)
```

You can enable the multicore feature to utilize the multicore advantages. Here is the benchmark. 

``` r
# set up environment

single.core <- function() {
  suppressWarnings(capture.output(DTSEA(network = example_ppi, disease = example_disease_list, drugs = example_drug_target_list, nproc = 0)))
  NULL
}

dual.core <- function() {
  suppressWarnings(capture.output(DTSEA(network = example_ppi, disease = example_disease_list, drugs = example_drug_target_list, nproc = 10)))
  NULL
}

single.core()
dual.core()
```

### Supplementary data files

In this package, we provide the example data, which is a small set of data to demonstrate the usage and the main idea behind DTSEA. 
We provide some extra data files, the real data we used in the DTSEA paper. 
The supplementary package is now on the [GitHub](https://github.com/hanjunwei-lab/DTSEAdata). Anyone can obtain this package by:

```{r echo=FALSE}
if (!"devtools" %in% as.data.frame(installed.packages())$Package)
  install.packages("devtools")
devtools::install_github("hanjunwei-lab/DTSEAdata")
```

## Known bugs

The Intel Math Kernel Library (MKL) performs poorly with this package when dealing with linear algebra operations. If you use MKL-based BLAS or MKL-based R distribution (like [Microsoft R Open](https://mran.microsoft.com/)), you will get unexpected or zero results in certain circumstances. Please install the Automatically Tuned Linear Algebra package (libatlas) or the multi-threaded OpenBlas library in order to get higher performance and reliable results with:

```
sudo apt install libatlas3-base -y
```
or
```
sudo apt install libopenblas-base -y
```
