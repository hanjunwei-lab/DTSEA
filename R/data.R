#' @title An example vector of disease nodes
#' @description The list was integrated the significantly differentially
#' expressed genes (DEGs) of GEO dataset GSE183071 and the work from Feng,
#' Song, Guo, and et al.
#'
#' @references
#' Gómez-Carballa A, Rivero-Calle I, Pardo-Seco J, Gómez-Rial J, Rivero-Velasco
#' C, Rodríguez-Núñez N, Barbeito-Castiñeiras G, Pérez-Freixo H, Cebey-López M,
#' Barral-Arca R, Rodriguez-Tenreiro C, Dacosta-Urbieta A, Bello X, Pischedda
#' S, Currás-Tuala MJ, Viz-Lasheras S, Martinón-Torres F, Salas A; GEN-COVID
#' study group. A multi-tissue study of immune gene expression profiling
#' highlights the key role of the nasal epithelium in COVID-19 severity.
#' Environ Res. 2022 Jul;210:112890. doi: 10.1016/j.envres.2022.112890. Epub
#' 2022 Feb 22. PMID: 35202626; PMCID: PMC8861187.
#'
#' Feng S, Song F, Guo W, Tan J, Zhang X, Qiao F, Guo J, Zhang L, Jia X.
#' Potential Genes Associated with COVID-19 and Comorbidity. Int J Med Sci.
#' 2022 Jan 24;19(2):402-415. doi: 10.7150/ijms.67815. PMID: 35165525; PMCID:
#' PMC8795808.
#'
#' @examples
#' library(DTSEA)
#'
#' data("example_disease_list", package = "DTSEA")
"example_disease_list"


#' @title An example data frame of drug target lists
#' @description Drug-target interactions were downloaded and integrated from
#' DrugBank and ChEMBL.
#' @format A data frame with 970 rows and 3 variables:
#' - `drug_id`: the DrugBank ID
#' - `drug_name`: the name of each drug
#' - `gene_target`: the targets of drugs
#' @references
#' Wishart DS, Feunang YD, Guo AC, Lo EJ, Marcu A, Grant JR, Sajed T, Johnson D,
#' Li C, Sayeeda Z, Assempour N, Iynkkaran I, Liu Y, Maciejewski A, Gale N,
#' Wilson A, Chin L, Cummings R, Le D, Pon A, Knox C, Wilson M. DrugBank 5.0: a
#' major update to the DrugBank database for 2018. Nucleic Acids Res. 2018 Jan
#' 4;46(D1):D1074-D1082. doi: 10.1093/nar/gkx1037. PMID: 29126136; PMCID:
#' PMC5753335.
#'
#' Gaulton A, Hersey A, Nowotka M, Bento AP, Chambers J, Mendez D, Mutowo P,
#' Atkinson F, Bellis LJ, Cibrián-Uhalte E, Davies M, Dedman N, Karlsson A,
#' Magariños MP, Overington JP, Papadatos G, Smit I, Leach AR. The ChEMBL
#' database in 2017. Nucleic Acids Res. 2017 Jan 4;45(D1):D945-D954. doi:
#' 10.1093/nar/gkw1074. Epub 2016 Nov 28. PMID: 27899562; PMCID: PMC5210557.
#'
#' @examples
#' library(DTSEA)
#' data("example_drug_target_list", package = "DTSEA")
"example_drug_target_list"


#' @title An example human gene functional interaction network object
#' @description We extracted the gene functional interaction network from
#' multiple sources with experimental evidence and then integrated them.
#' @format An igraph object
#' @references
#' Kanehisa M, Furumichi M, Sato Y, Ishiguro-Watanabe M, Tanabe M. KEGG:
#' integrating viruses and cellular organisms. Nucleic Acids Res. 2021 Jan 8;49
#' (D1):D545-D551. doi: 10.1093/nar/gkaa970. PMID: 33125081; PMCID: PMC7779016.
#'
#' @examples
#' library(DTSEA)
#' data("example_ppi", package = "DTSEA")
"example_ppi"



#' @title A random graph for the computation of the separation measure
#' @description The random graph was retrieved from Menche et al (2015).
#' @format An igraph object
#' @references
#' Menche J, Sharma A, Kitsak M, Ghiassian SD, Vidal M, Loscalzo J, Barabási AL.
#' Disease networks. Uncovering disease-disease relationships through the
#' incomplete interactome. Science. 2015 Feb 20;347(6224):1257601. doi:
#' 10.1126/science.1257601. PMID: 25700523; PMCID: PMC4435741.
#'
#' @examples
#' library(DTSEA)
#' data("random_graph", package = "DTSEA")
"random_graph"
