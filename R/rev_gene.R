#' @title Pre-processing step: revising gene symbols
#'
#' @description Revise genes according to NCBI Gene symbols updated in Jan. 1, 2022 for count matrix, user-custom cell marker data.frame.
#' @param data A matrix or dgCMatrix containing count or normalized data, each column representing a spot or a cell, each row representing a gene; Or a data.frame containing cell markers, use \code{\link{demo_marker}}.
#' @param data_type A character to define the type of \code{data}, select \code{'data'} for the data matrix, \code{'marker'} for the data.frame containing cell markers.
#' @param species Species of the data.\code{'Human'} or \code{'Mouse'}.
#' @param geneinfo A data.frame of the system data containing gene symbols of \code{'Human'} and \code{'Mouse'} updated on Jan. 1, 2022.
#' @return A new matrix or data.frame.
#' @importFrom crayon red cyan green
#' @export rev_gene

rev_gene <- function(data = NULL, data_type = NULL, species = NULL, geneinfo = NULL) {
    if (is.null(data)) {
        stop("Please provide the data for revsing gene symbols!")
    }
    if (is.null(data_type) | !is.character(data_type)) {
        stop("Please provide a correct data_type, i.e., 'data', or 'marker'!")
    }
    if (is.null(geneinfo)) {
        stop("Please provide geneinfo for revsing gene symbols, or use the system data like 'geneinfo = geneinfo'")
    }
    if (length(data_type) > 1 | !data_type %in% c("data", "marker")) {
        stop("Please provide a correct data_type, i.e., 'data', or 'marker'!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse")) {
        stop("Please provide a correct species, i.e., 'Human', or 'Mouse'!")
    }
    # define the species
    if (species == "Human") {
        geneinfo <- geneinfo[geneinfo$species == "Human", ]
    }
    if (species == "Mouse") {
        geneinfo <- geneinfo[geneinfo$species == "Mouse", ]
    }
    # revise matrix
    if (data_type == "data") {
        if (!is.matrix(data)) {
            if (!is(data, "dgCMatrix")) {
                stop("data must be a matrix or dgCMatrix when data_type is 'data'!")
            }
        } else {
            data <- as(data, Class = "dgCMatrix")
        }
        cat(crayon::cyan("Revising gene symbols for data matrix", "\n"))
        Sys.sleep(1)
        # revise gene symbols
        genename <- rownames(data)
        genename1 <- genename[genename %in% geneinfo$symbol]
        if (length(genename1) == 0) {
            stop("Please ensure the rownames of data are gene symbols! See demo_count()!")
        }
        genename2 <- genename[!genename %in% geneinfo$symbol]
        if (length(genename2) > 0) {
            genename3 <- genename2[genename2 %in% geneinfo$synonyms]
            if (length(genename3) > 0) {
                genename4 <- rep("NA", length(genename3))
                for (i in 1:length(genename3)) {
                  d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                  if (length(d1) == 1) {
                    genename4[i] <- d1
                  }
                }
                genename3 <- c(genename1, genename3)
                genename4 <- c(genename1, genename4)
                genedata <- data.frame(raw_name = genename3, new_name = genename4, stringsAsFactors = F)
                genedata <- genedata[!genedata$new_name == "NA", ]
                genedata1 <- as.data.frame(table(genedata$new_name), stringsAsFactors = F)
                genedata1 <- genedata1[genedata1$Freq == 1, ]
                genedata <- genedata[genedata$new_name %in% genedata1$Var1, ]
                data <- data[genedata$raw_name, ]
                rownames(data) <- genedata$new_name
            }
        } else {
            data <- data[rownames(data) %in% geneinfo$symbol, ]
        }
        cat(crayon::green("***Done***", "\n"))
        Sys.sleep(2)
    }
    # revise lrpairs
    if (data_type == "marker") {
        if (!is.data.frame(data)) {
            stop("data must be a data.frame when data_type is 'marker'!")
        }
        cat(crayon::cyan("Revising gene symbols for marker data.frame", "\n"))
        Sys.sleep(1)
        markers <- demo_marker()
        if (all(colnames(markers) %in% colnames(data))) {
            # genes
            genename <- unique(data$gene)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the gene of data are gene symbols! See demo_marker()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                  for (i in 1:length(genename3)) {
                    d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                    if (length(d1) == 1) {
                      data[data$ligand == genename3[i], ]$ligand <- d1
                    }
                  }
                }
            }
            data <- data[data$gene %in% geneinfo$symbol, ]
        } else {
            stop("Please provide a correct marker data.frame! See demo_marker()!")
        }
        cat(crayon::green("***Done***", "\n"))
        Sys.sleep(2)
    }
    return(data)
}
