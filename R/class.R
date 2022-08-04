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
