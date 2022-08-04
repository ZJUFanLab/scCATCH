#' @title scCATCH object
#'
#' @description create scCATCH object using single-cell count data and cluster information.
#' @param data A matrix or dgCMatrix containing normalized single-cell RNA-seq data, each column representing a cell, each row representing a gene. See \code{\link{demo_data}}.
#' @param cluster A character containing the cluster information for each cell. The length of it must be equal to the ncol of the data.
#' @return scCATCH object
#' @importFrom methods as
#' @export createscCATCH

createscCATCH <- function(data, cluster) {
    if (!is.matrix(data)) {
        if (!is(data, "dgCMatrix")) {
            stop("count is not a matrix or dgCMatrix! See demo_data()")
        }
    } else {
        data <- as(data, Class = "dgCMatrix")
    }
    if (!is.character(cluster)) {
        stop("cluster is not a character!")
    }
    if (length(cluster) != ncol(data)) {
        stop("Length of cluster is not equal to the ncol of the data!")
    }
    meta <- data.frame(cell = colnames(data), cluster = cluster, stringsAsFactors = FALSE)
    object <- new("scCATCH", data = list(ndata = data), meta = meta)
    return(object)
}

#' Find potential marker genes for each cluster
#'
#' @description Identify potential marker genes for each cluster.
#' @param object scCATCH object generated from \code{\link{createscCATCH}}.
#' @param species The specie of cells. The species must be defined. 'Human' or 'Mouse'. When if_use_custom_marker is set TRUE, no need to define the species.
#' @param cluster Select which clusters for potential marker genes identification. e.g. '1', '2', etc. The default is 'All' to find potential makrer genes for each cluster.
#' @param if_use_custom_marker Whether to use custom markers data.frame.
#' @param marker A data.frame containing marker genes. See \code{\link{demo_marker}}. Default is to use the system \code{\link{cellmatch}} data.frame.
#' @param cancer If the sample is from cancer tissue, then the cancer type may be defined. When if_use_custom_marker is set TRUE, no need to define the cancer.
#' @param tissue Tissue origin of cells must be defined. Select one or more related tissue types. When if_use_custom_marker is set TRUE, no need to define the tissue.
#' @param use_method '1' is to compare with other every cluster. '2' is to compare with other clusters together.
#' @param comp_cluster Number of clusters to compare. Default is to compare all other cluster for each cluster. Set it between 1 and length of unique clusters. More marker genes will be obtained for smaller comp_cluster.
#' @param cell_min_pct Include the gene detected in at least this many cells in each cluster.
#' @param logfc Include the gene with at least this fold change of average gene expression compared to every other clusters.
#' @param pvalue Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters.
#' @param verbose Show progress messages.
#' @return scCATCH object
#' @details Details of available tissues see \url{https://github.com/ZJUFanLab/scCATCH/wiki}
#' @import methods progress
#' @importFrom stats wilcox.test
#' @export

findmarkergene <- function(object, species = NULL, cluster = "All", if_use_custom_marker = FALSE,
    marker = NULL, cancer = "Normal", tissue = NULL, use_method = "1", comp_cluster = NULL, cell_min_pct = 0.25, logfc = 0.25,
    pvalue = 0.05, verbose = TRUE) {
    if (!is(object, "scCATCH")) {
        stop("object must be scCATCH obect generated from createscCATCH()!")
    }
    if (is.null(marker)) {
        stop("Please provide the marker, either the system cellmatch or custom marker data.frame!")
    }
    ndata <- object@data$ndata
    ndata <- .filter_ndata(ndata)
    meta <- object@meta
    if (!if_use_custom_marker) {
        marker <- .filter_marker(marker, species, cancer, tissue)
    }
    marker <- marker[marker$gene %in% rownames(ndata),]
    if (length(unique(marker$gene)) < 2) {
        stop("No matched potential marker genes in the matrix!")
    }
    ndata <- ndata[rownames(ndata) %in% unique(marker$gene), ]
    # generating pair-wise clusters
    clu_pair <- .get_clu_pair(meta, cluster)
    clu_num <- clu_pair[[2]]
    clu_pair <- clu_pair[[1]]
    cluster <- unique(meta$cluster)
    if (use_method == "1") {
        clu_marker <- .get_marker_scCATCH1(ndata, meta, cluster, clu_num, clu_pair, cell_min_pct, logfc, pvalue, comp_cluster, verbose)
    }
    if (use_method == "2") {
        clu_marker <- .get_marker_scCATCH2(ndata, meta, cluster, clu_num, clu_pair, cell_min_pct, logfc, pvalue, verbose)
    }
    if (is.null(clu_marker)) {
        warning("No potential marker genes found in the matrix! You may adjust the cell clusters by applying other clustering algorithm and try again.")
    }
    if (nrow(clu_marker) == 0) {
        warning("No potential marker genes found in the matrix! You may adjust the cell clusters by applying other clustering algorithm and try again.")
    }
    if (nrow(clu_marker) > 0) {
        object@markergene <- clu_marker
    }
    para <- list(species = species, cluster = cluster, if_use_custom_marker = if_use_custom_marker, cancer = cancer, tissue = tissue,
        use_method = use_method, cell_min_pct = cell_min_pct, logfc = logfc, pvalue = pvalue)
    object@para <- para
    object@marker <- marker
    return(object)
}

#' Evidence-based score and annotation for each cluster
#'
#' @description Evidence-based score and annotation for each cluster.
#' @param object scCATCH object generated from \code{\link{findmarkergene}}.
#' @param verbose Show progress messages.
#' @return scCATCH object containing the results of predicted cell types for each cluster.
#' @import methods progress
#' @importFrom reshape2 melt
#' @export

findcelltype <- function(object, verbose = TRUE) {
    if (!is(object, "scCATCH")) {
        stop("object must be scCATCH obect generated from createscCATCH()!")
    }
    markergene <- object@markergene
    marker <- object@marker
    if (nrow(markergene) == 0) {
        stop("No marker genes found!")
    }
    # Evidence-based scoring and annotation
    Sys.sleep(2)
    clu_num <- unique(markergene$cluster)
    pb <- progress::progress_bar$new(format = "[:bar] Finished::percent Time :elapsedfull", total = length(clu_num), clear = FALSE, width = 60,
        complete = "+", incomplete = "-")
    clu_ann_res <- NULL
    for (i in 1:length(clu_num)) {
        cellsubtype1 <- "NA"
        cellsubtype2 <- "NA"
        cellsubtype3 <- "NA"
        celltype <- "NA"
        celltype_score <- "NA"
        clu_marker <- "NA"
        PMID <- "NA"
        clu_markergene <- unique(markergene[markergene$cluster == clu_num[i], ]$gene)
        if (length(clu_markergene) > 0) {
            res_clu <- .get_clu_ann(clu_markergene, marker)
            cellsubtype1 <- res_clu[[1]]
            cellsubtype2 <- res_clu[[2]]
            cellsubtype3 <- res_clu[[3]]
            celltype <- res_clu[[4]]
            celltype_score <- res_clu[[5]]
            clu_marker <- res_clu[[6]]
            PMID <- res_clu[[7]]
        } else {
            clu_markergene <- "NA"
        }
        # processing the format of cell type
        d1 <- clu_markergene[1]
        if (length(clu_markergene) > 1) {
            for (j in 2:length(clu_markergene)) {
                d1 <- paste(d1, clu_markergene[j], sep = ", ")
            }
        }
        clu_markergene <- d1
        d1 <- celltype[1]
        if (length(celltype) > 1) {
            for (j in 2:length(celltype)) {
                d1 <- paste(d1, celltype[j], sep = ", ")
            }
        }
        celltype <- d1
        d1 <- celltype_score[1]
        if (length(celltype_score) > 1) {
            for (j in 2:length(celltype_score)) {
                d1 <- paste(d1, celltype_score[j], sep = ", ")
            }
        }
        celltype_score <- d1
        d1 <- cellsubtype1[1]
        if (length(cellsubtype1) > 1) {
            for (j in 2:length(cellsubtype1)) {
                d1 <- paste(d1, cellsubtype1[j], sep = "; ")
            }
        }
        cellsubtype1 <- d1
        d1 <- cellsubtype2[1]
        if (length(cellsubtype2) > 1) {
            for (j in 2:length(cellsubtype2)) {
                d1 <- paste(d1, cellsubtype2[j], sep = "; ")
            }
        }
        cellsubtype2 <- d1
        d1 <- cellsubtype3[1]
        if (length(cellsubtype3) > 1) {
            for (j in 2:length(cellsubtype3)) {
                d1 <- paste(d1, cellsubtype3[j], sep = "; ")
            }
        }
        cellsubtype3 <- d1
        d1 <- clu_marker[1]
        if (length(clu_marker) > 1) {
            for (j in 2:length(clu_marker)) {
                d1 <- paste(d1, clu_marker[j], sep = ", ")
            }
        }
        clu_marker <- d1
        d1 <- PMID[1]
        if (length(PMID) > 1) {
            for (j in 2:length(PMID)) {
                d1 <- paste(d1, PMID[j], sep = ", ")
            }
        }
        PMID <- d1
        clu_ann <- data.frame(cluster = clu_num[i], cluster_marker = clu_markergene, cellsubtype3 = cellsubtype3, cellsubtype2 = cellsubtype2,
            cellsubtype1 = cellsubtype1, cell_type = celltype, celltype_score = celltype_score, celltype_related_marker = clu_marker,
            PMID = PMID, stringsAsFactors = FALSE)
        clu_ann_res <- rbind(clu_ann_res, clu_ann)
        if (verbose) {
            pb$tick()
        }
    }
    for (i in 1:nrow(clu_ann_res)) {
        d1 <- as.character(clu_ann_res[i, c("cellsubtype3", "cellsubtype2", "cellsubtype1", "cell_type")])
        d1 <- d1[!is.na(d1)]
        d1 <- d1[!d1 == "NA"]
        d2 <- d1[1]
        if (length(d1) > 1) {
            for (j in 2:length(d1)) {
                d2 <- paste(d2, d1[j], sep = " ")
            }
        }
        clu_ann_res[i, "cell_type"] <- d2
    }
    clu_ann_res <- clu_ann_res[, -c(3, 4, 5)]
    object@celltype <- clu_ann_res
    return(object)
}
