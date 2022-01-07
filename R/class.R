#' @title Definition of 'scCATCH' class
#'
#' @description An S4 class containing the data, meta, and results of inferred cell types.
#' @slot data A list containing normalized data. See \code{\link{demo_data}}.
#' @slot meta A data frame containing the meta data.
#' @slot para A list containing the parameters.
#' @slot markergene A data frame containing the identified markers for each cluster.
#' @slot celltype A data frame containing the cell types for each cluster.
#' @slot marker A data frame containing the known markers. See \code{\link{demo_marker}}.
#' @import methods
#' @name scCATCH
#' @rdname scCATCH
#' @aliases scCATCH-class
#' @exportClass scCATCH

setClass("scCATCH", representation(data = "list", meta = "data.frame", para = "list", markergene = "data.frame", celltype = "data.frame",
    marker = "data.frame"), prototype(data = list(), meta = data.frame(), para = list(), markergene = data.frame(), celltype = data.frame(),
    marker = data.frame()))

#' Find potential marker genes for each cluster
#'
#' @description Identify potential marker genes for each cluster.
#' @param object scCATCH object generated from \code{\link{createscCATCH}}.
#' @param species The specie of cells. The species must be defined. 'Human' or 'Mouse'. When if_use_custom_marker is set TRUE, no need to define the species.
#' @param cluster Select which clusters for potential marker genes identification. e.g. '1', '2', etc. The default is 'All' to find potential makrer genes for each cluster.
#' @param if_use_custom_marker Whether to use custom markers data.frame.
#' @param marker A data.frame containing marker genes. See \code{\link{demo_marker}}. Default is to use the system \code{\link{cellmatch}} data.frame.
#' @param cancer If the sample is from cancer tissue, then the cancer type may be defined. When if_use_custom_marker is set TRUE, no need to define the species.
#' @param tissue Tissue origin of cells must be defined. Select one or more related tissue types. When if_use_custom_marker is set TRUE, no need to define the species.
#' @param use_method '1' is to compare with other every cluster. '2' is to compare with other clusters together.
#' @param comp_cluster Number of clusters to compare. Default is to compare all other cluster for each cluster. Set it between 1 and length of unique clusters. More marker genes will be obtained for smaller comp_cluster.
#' @param cell_min_pct Include the gene detected in at least this many cells in each cluster.
#' @param logfc Include the gene with at least this fold change of average gene expression compared to every other clusters.
#' @param pvalue Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters.
#' @param verbose Show progress messages.
#' @return scCATCH object
#' @details Details of available tissues see \url{https://github.com/ZJUFanLab/scCATCH/wiki}
#' @import methods
#' @export

setGeneric("findmarkergene", def = function(object, species = NULL, cluster = "All", if_use_custom_marker = FALSE, marker = NULL,
    cancer = "Normal", tissue = NULL, use_method = "1", comp_cluster = NULL, cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05, verbose = TRUE) {
    standardGeneric("findmarkergene")
})

#' Evidence-based score and annotation for each cluster
#'
#' @description Evidence-based score and annotation for each cluster.
#' @param object scCATCH object generated from \code{\link{findmarkergene}}.
#' @param verbose Show progress messages.
#' @return scCATCH object containing the results of predicted cell types for each cluster.
#' @import methods
#' @export

setGeneric("findcelltype", def = function(object, verbose = TRUE) {
    standardGeneric("findcelltype")
})
