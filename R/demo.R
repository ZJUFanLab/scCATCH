#' @title Demo data of single-cell RNA-seq data
#'
#' @description Demo data of single-cell RNA-seq data
#' @details \code{data} used in \code{\link{createscCATCH}} must be a \code{matrix} object, each column representing a cell, each row representing a gene.
#' @return A matrix.
#' @export demo_data
#' @examples data_demo <- demo_data()

demo_data <- function() {
    cellname <- paste0("cell", 1:6)
    genename <- c("A1BG", "A2M", "A2MP", "NAT1", "NAT2", "NAT20")
    countdata <- sample(x = c(rep(0, 50), 1:50), size = 36, replace = T)
    countdata <- matrix(countdata, nrow = 6, ncol = 6)
    rownames(countdata) <- genename
    colnames(countdata) <- cellname
    return(countdata)
}

#' @title Demo data of geneinfo
#'
#' @description Demo data of geneinfo
#' @details \code{geneinfo} used in \code{\link{rev_gene}} must be a \code{data.frame} object with three columns, namely \code{'symbol'}, \code{'synonyms'}, \code{'species'}.
#' @export demo_geneinfo
#' @examples geneinfo_demo <- demo_geneinfo()

demo_geneinfo <- function() {
    gene1 <- c("A1BG", "A1BG", "A2MP1", "Aco1")
    gene2 <- c("A1B", "ABG", "A2MP", "Aco")
    species <- c("Human", "Human", "Human", "Mouse")
    geneinfo_demo <- data.frame(symbol = gene1, synonyms = gene2, species = species, stringsAsFactors = F)
    return(geneinfo_demo)
}

#' @title Demo data of markers
#'
#' @description Demo data of markers
#' @details \code{markers} used in \code{\link{findmarkergene}} must be a \code{data.frame} object with eleven columns.
#' @export demo_marker
#' @examples markers_demo <- demo_marker()

demo_marker <- function() {
    species <- c("Human", "Human", "Human", "Human")
    tissues <- c("Liver", "Liver", "Liver", "Liver")
    cancers <- c("Normal", "Normal", "Hepatocellular Cancer", "Hepatocellular Cancer")
    conditions <- c("Normal cell", "Normal cell", "Cancer cell", "Cancer cell")
    subtype1 <- c("NA", "NA", "NA", "Regulatory")
    subtype2 <- c("CD4+", "CD8+", "NA", "NA")
    subtype3 <- c("NA", "NA", "Exhausted", "NA")
    celltype <- c("T Cell", "T Cell", "T Cell", "T Cell")
    genes <- c("CD4", "CD8A", "ABCG1", "ACP5")
    resources <- c("Experiment", "Experiment", "Single-cell sequencing", "Single-cell sequencing")
    pmid <- c("27781378", "27781378", "28622514", "28622514")
    markers <- data.frame(species = species, tissue = tissues, cancer = cancers, condition = conditions, subtype1 = subtype1,
        subtype2 = subtype2, subtype3 = subtype3, celltype = celltype, gene = genes, resource = resources, pmid = pmid, stringsAsFactors = F)
    return(markers)
}
