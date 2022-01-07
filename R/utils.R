.filter_ndata <- function(ndata) {
    gene_sum <- as.numeric(apply(ndata, 1, sum))
    ndata <- ndata[which(gene_sum > 0), ]
    return(ndata)
}

.percent_cell <- function(x) {
    return(length(x[x > 0])/length(x))
}

.filter_marker <- function(marker, species, cancer, tissue) {
    if (is.null(species)) {
        stop("Please provide the correct species, 'Human' or 'Mouse'!")
    }
    if (length(species) != 1) {
        stop("Please provide the correct species, 'Human' or 'Mouse'!")
    }
    if (!species %in% c("Human", "Mouse")) {
        stop("Please provide the correct species, 'Human' or 'Mouse'!")
    }
    marker <- marker[marker$species == species, ]
    if (cancer == "Normal") {
        marker <- marker[marker$cancer == "Normal", ]
        # check tissue
        if (is.null(tissue)) {
            stop("Please define the origin tissue of cells! Select one or more related tissue types.")
        }
        tissue_match <- NULL
        # check the tissue
        tissue_match <- tissue %in% marker$tissue
        tissue_match <- which(tissue_match == FALSE)
        if (length(tissue_match) > 0) {
            stop(paste(tissue[tissue_match], ", not matched with the tissue types in CellMatch database!", sep = ""))
        }
        marker <- marker[marker$tissue %in% tissue, ]
        cellmarkers <- marker$gene
        tissue1 <- tissue[1]
        if (length(tissue) > 1) {
            for (i in 2:length(tissue)) {
                tissue1 <- paste(tissue1, tissue[i], sep = ", ")
            }
        }
        if (length(cellmarkers) < 100) {
            warning(paste("There are only ", length(cellmarkers), " potential marker genes in CellMatch database for ",
                species, " on ", tissue1, "!", sep = ""))
        }
        if (length(cellmarkers) >= 100) {
            message(paste("There are ", length(cellmarkers), " potential marker genes in CellMatch database for ", species,
                " on ", tissue1, ".", sep = ""))
        }
    } else {
        # check cancer
        cancer_match <- NULL
        # check the tissue
        cancer_match <- cancer %in% marker$cancer
        cancer_match <- which(cancer_match == FALSE)
        if (length(cancer_match) > 0) {
            stop(paste(cancer[cancer_match], ", not matched with the cancer types in CellMatch database!", sep = ""))
        }
        marker <- marker[marker$cancer %in% cancer, ]
        # check tissue
        if (is.null(tissue)) {
            stop("Please define the origin tissue of cells! Select one or more related tissue types.")
        }
        tissue_match <- NULL
        # check the tissue
        tissue_match <- tissue %in% marker$tissue
        tissue_match <- which(tissue_match == FALSE)
        if (length(tissue_match) > 0) {
            stop(paste(tissue[tissue_match], ", not matched with the tissue types in CellMatch database!", sep = ""))
        }
        marker <- marker[marker$tissue %in% tissue, ]
        cellmarkers <- marker$gene
        cancer1 <- cancer[1]
        if (length(cancer) > 1) {
            for (i in 2:length(cancer)) {
                cancer1 <- paste(cancer1, cancer[i], sep = ", ")
            }
        }
        tissue1 <- tissue[1]
        if (length(tissue) > 1) {
            for (i in 2:length(tissue)) {
                tissue1 <- paste(tissue1, tissue[i], sep = ", ")
            }
        }
        if (length(cellmarkers) < 100) {
            warning(paste("There are only ", length(cellmarkers), " potential marker genes in CellMatch database for ",
                species, " ", cancer1, " on ", tissue1, "!", sep = ""))
        }
        if (length(cellmarkers) >= 100) {
            message(paste("There are ", length(cellmarkers), " potential marker genes in CellMatch database for ", species,
                " ", cancer1, " on ", tissue1, ".", sep = ""))
        }
    }
    return(marker)
}

.get_clu_pair <- function(meta, cluster) {
    clu_num <- unique(meta$cluster)
    clu_pair <- NULL
    for (i in 1:length(clu_num)) {
        d1 <- data.frame(cluster1 = rep(clu_num[i], (length(clu_num) - 1)), cluster2 = clu_num[-i], stringsAsFactors = FALSE)
        clu_pair <- rbind(clu_pair, d1)
    }
    if (cluster[1] != "All") {
        cluster_match <- cluster %in% clu_num
        cluster_match <- which(cluster_match == FALSE)
        if (length(cluster_match) > 0) {
            stop(paste(cluster[cluster_match], ", not matched with the cell clusters! Please select one or more related clusters.",
                sep = ""))
        }
        clu_pair <- clu_pair[clu_pair$cluster1 %in% cluster, ]
        clu_num1 <- clu_num[clu_num %in% cluster]
        return(list(clu_pair = clu_pair, clu_num = clu_num1))
    } else {
        return(list(clu_pair = clu_pair, clu_num = clu_num))
    }
}

.get_pct_ndata1 <- function(ndata1) {
    ndata1_pct <- apply(ndata1, 1, .percent_cell)
    return(as.numeric(ndata1_pct))
}

.get_pvalue <- function(ndata1, ndata2) {
    clu_marker_pvalue <- rep(1, nrow(ndata1))
    # assigning pct and pvalue
    for (k in 1:nrow(ndata1)) {
        genedata1 <- as.numeric(ndata1[k, ])
        genedata2 <- as.numeric(ndata2[k, ])
        clu_marker_pvalue[k] <- as.numeric(stats::wilcox.test(genedata1, genedata2, alternative = "greater", paired = FALSE)$p.value)
    }
    return(clu_marker_pvalue)
}

.get_marker <- function(ndata, meta, ndata1, clu_pair1, logfc, pvalue, pct_ndata1, cell_min_pct) {
    clu_marker1 <- NULL
    ndata_temp <- ndata[rownames(ndata1), ]
    for (j in 1:nrow(clu_pair1)) {
        clu_marker2 <- as.data.frame(matrix(data = 0, nrow = nrow(ndata_temp), ncol = 6))
        colnames(clu_marker2) <- c("cluster", "gene", "pct", "comp_cluster", "logfc", "pvalue")
        clu_marker2$cluster <- unique(clu_pair1$cluster1)
        clu_marker2$gene <- rownames(ndata_temp)
        clu_marker2$comp_cluster <- clu_pair1$cluster2[j]
        clu_marker2$pct <- pct_ndata1
        # extract normalized data of pair wise clusters
        ndata2 <- ndata_temp[, meta[meta$cluster == clu_pair1$cluster2[j], ]$cell]
        clu_marker2$logfc <- as.numeric(apply(ndata1, 1, mean)) - as.numeric(apply(ndata2, 1, mean))
        clu_marker2$pvalue <- .get_pvalue(ndata1, ndata2)
        # filtering genes with cell_min_pct logfc
        clu_marker2 <- clu_marker2[clu_marker2$pct >= cell_min_pct & clu_marker2$logfc >= logfc & clu_marker2$pvalue < pvalue, ]
        # combine the results for each cluster
        clu_marker1 <- rbind(clu_marker1, clu_marker2)
    }
    return(clu_marker1)
}

.get_marker_scCATCH1 <- function(ndata, meta, cluster, clu_num, clu_pair, cell_min_pct, logfc, pvalue, comp_cluster, verbose) {
    # generating result file
    clu_marker <- NULL
    # calculate the p value and log fold change
    pb <- progress::progress_bar$new(format = "[:bar] Finished::percent Time :elapsedfull", total = length(clu_num), clear = FALSE,
        width = 60, complete = "+", incomplete = "-")
    for (i in 1:length(clu_num)) {
        clu_pair1 <- clu_pair[clu_pair$cluster1 == clu_num[i], ]
        ndata1 <- ndata[, meta[meta$cluster == clu_num[i], ]$cell]
        pct_ndata1 <- .get_pct_ndata1(ndata1)
        ndata1 <- ndata1[which(pct_ndata1 >= cell_min_pct), ]
        pct_ndata1 <- .get_pct_ndata1(ndata1)
        clu_marker1 <- .get_marker(ndata, meta, ndata1, clu_pair1, logfc, pvalue, pct_ndata1, cell_min_pct)
        # generating result file for each cluster
        if (is.null(comp_cluster)) {
            comp_cluster <- length(cluster) - 1
        } else {
            if (comp_cluster > (length(cluster) - 1)) {
                stop("comp_cluster must be less than the length of unique cluster number!")
            }
        }
        if (!is.null(clu_marker1)) {
            if (nrow(clu_marker1) > 0) {
                d1 <- as.data.frame(table(clu_marker1$gene), stringsAsFactors = FALSE)
                d1 <- d1[d1$Freq >= comp_cluster, ]
                if (nrow(d1) > 0) {
                  clu_marker1 <- clu_marker1[clu_marker1$gene %in% d1$Var1, ]
                  clu_marker <- rbind(clu_marker, clu_marker1)
                }
            }
        }
        if (verbose) {
          pb$tick()
        }
    }
    return(clu_marker)
}

.get_marker_scCATCH2 <- function(ndata, meta, cluster, clu_num, clu_pair, cell_min_pct, logfc, pvalue, verbose) {
    # generating result file
    clu_marker <- NULL
    # calculate the p value and log fold change
    pb <- progress::progress_bar$new(format = "[:bar] Finished::percent Time :elapsedfull", total = length(clu_num), clear = FALSE,
        width = 60, complete = "+", incomplete = "-")
    for (i in 1:length(clu_num)) {
        clu_pair1 <- clu_pair[clu_pair$cluster1 == clu_num[i], ]
        ndata1 <- ndata[, meta[meta$cluster == clu_num[i], ]$cell]
        pct_ndata1 <- .get_pct_ndata1(ndata1)
        ndata1 <- ndata1[which(pct_ndata1 >= cell_min_pct), ]
        pct_ndata1 <- .get_pct_ndata1(ndata1)
        ndata_temp <- ndata[rownames(ndata1), ]
        clu_marker2 <- as.data.frame(matrix(data = 0, nrow = nrow(ndata_temp), ncol = 5))
        colnames(clu_marker2) <- c("cluster", "gene", "pct", "logfc", "pvalue")
        clu_marker2$cluster <- unique(clu_pair1$cluster1)
        clu_marker2$gene <- rownames(ndata_temp)
        clu_marker2$pct <- pct_ndata1
        # extract normalized data of pair wise clusters
        ndata2 <- ndata_temp[, meta[meta$cluster %in% clu_pair1$cluster2, ]$cell]
        clu_marker2$logfc <- as.numeric(apply(ndata1, 1, mean)) - as.numeric(apply(ndata2, 1, mean))
        clu_marker2$pvalue <- .get_pvalue(ndata1, ndata2)
        # filtering genes with cell_min_pct logfc
        clu_marker2 <- clu_marker2[clu_marker2$pct >= cell_min_pct & clu_marker2$logfc >= logfc & clu_marker2$pvalue < pvalue, ]
        # combine the results for each cluster
        clu_marker <- rbind(clu_marker, clu_marker2)
        if (verbose) {
            pb$tick()
        }
    }
    Sys.sleep(2)
    return(clu_marker)
}

.get_score <- function(clu_ann) {
    # [1]
    clu_ann_cellname <- as.data.frame(table(clu_ann$celltype), stringsAsFactors = FALSE)
    m <- clu_ann_cellname$Freq + 1
    clu_ann_cellname$Freq <- clu_ann_cellname$Freq/m
    clu_ann_cellname <- clu_ann_cellname[order(clu_ann_cellname$Var1), ]
    # [2]
    clu_ann_article <- unique(clu_ann[, c("pmid", "celltype")])
    clu_ann_article <- as.data.frame(table(clu_ann_article$celltype), stringsAsFactors = FALSE)
    m <- clu_ann_article$Freq + 1
    clu_ann_article$Freq <- clu_ann_article$Freq/m
    clu_ann_article <- clu_ann_article[order(clu_ann_article$Var1), ]
    # scoring and select
    clu_ann_cellname$Freq <- sqrt(clu_ann_article$Freq * clu_ann_cellname$Freq)
    clu_ann_cellname <- clu_ann_cellname[clu_ann_cellname$Freq == max(clu_ann_cellname$Freq), ]
    return(clu_ann_cellname)
}

.get_score_subtype <- function(clu_ann1) {
    clu_ann_subcellname <- as.data.frame(table(clu_ann1$value), stringsAsFactors = FALSE)
    m <- clu_ann_subcellname$Freq + 1
    clu_ann_subcellname$Freq <- clu_ann_subcellname$Freq/m
    clu_ann_subcellname <- clu_ann_subcellname[order(clu_ann_subcellname$Var1), ]
    clu_ann_subarticle <- unique(clu_ann1[, c("pmid", "value")])
    clu_ann_subarticle <- as.data.frame(table(clu_ann_subarticle$value), stringsAsFactors = FALSE)
    m <- clu_ann_subarticle$Freq + 1
    clu_ann_subarticle$Freq <- clu_ann_subarticle$Freq/m
    clu_ann_subarticle <- clu_ann_subarticle[order(clu_ann_subarticle$Var1), ]
    clu_ann_subcellname$Freq <- sqrt(clu_ann_subcellname$Freq * clu_ann_subarticle$Freq)
    clu_ann_subcellname <- clu_ann_subcellname[clu_ann_subcellname$Freq == max(clu_ann_subcellname$Freq), ]
    clu_ann_subcellname <- clu_ann_subcellname[clu_ann_subcellname$Freq > 0.5, ]
    return(clu_ann_subcellname)
}

.get_clu_ann <- function(clu_markergene, marker) {
    clu_ann <- marker[marker$gene %in% clu_markergene, ]
    clu_ann_cellname <- .get_score(clu_ann)
    celltype <- clu_ann_cellname$Var1
    celltype_score <- round(clu_ann_cellname$Freq, digits = 2)
    clu_marker <- unique(clu_ann[clu_ann$celltype %in% clu_ann_cellname$Var1, ]$gene)
    PMID <- unique(clu_ann[clu_ann$celltype %in% clu_ann_cellname$Var1, ]$pmid)
    if (nrow(clu_ann_cellname) == 1) {
        res <- .get_clu_ann_subtype(clu_ann, clu_ann_cellname, celltype, celltype_score, clu_marker, PMID)
    } else {
        res <- list(cellsubtype1 = "NA",cellsubtype2 = "NA",cellsubtype3 = "NA",celltype = celltype,celltype_score = celltype_score,clu_marker = clu_marker,PMID = PMID)
    }
    return(res)
}

.get_clu_ann_subtype <- function(clu_ann, clu_ann_cellname, celltype, celltype_score, clu_marker, PMID) {
    cellsubtype1 <- "NA"
    cellsubtype2 <- "NA"
    cellsubtype3 <- "NA"
    # matching cell subtype to determine the start
    clu_ann1 <- clu_ann[clu_ann$celltype == clu_ann_cellname$Var1, ]
    clu_ann1 <- clu_ann1[, c("gene", "pmid", "subtype3", "subtype2", "subtype1")]
    clu_ann1$row <- 1:nrow(clu_ann1)
    clu_ann1$row <- as.character(clu_ann1$row)
    clu_ann1 <- clu_ann1[, c(6, 1:5)]
    clu_ann1 <- reshape2::melt(data = clu_ann1, id.vars = 1:3, measure.vars = c("subtype3", "subtype2", "subtype1"))
    clu_ann1 <- as.data.frame(clu_ann1)
    clu_ann1$variable <- as.character(clu_ann1$variable)
    clu_ann1 <- clu_ann1[!is.na(clu_ann1$value), ]
    # exist cell subtype
    if (nrow(clu_ann1) > 0) {
        clu_ann_subcellname <- .get_score_subtype(clu_ann1)
        if (nrow(clu_ann_subcellname) == 1) {
            clu_ann_for_det <- clu_ann1[clu_ann1$value == clu_ann_subcellname$Var1, ]
            # max cell subtype label exist in subtype1
            if (unique(clu_ann_for_det$variable == "subtype1")) {
                cellsubtype1 <- unique(clu_ann_for_det$value)
                clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype1, ]$gene)
                PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype1, ]$pmid)
                clu_ann_for_second <- clu_ann1[clu_ann1$variable == "subtype2", ]
                # exist cell second subtype2
                if (nrow(clu_ann_for_second) > 0) {
                  # calculting the second cell subtype2 max score of the same cell type -- short name.
                  clu_second <- clu_ann1[clu_ann1$variable == "subtype2", ]
                  clu_second_name <- as.data.frame(table(clu_second$value), stringsAsFactors = FALSE)
                  m <- clu_second_name$Freq + 1
                  clu_second_name$Freq <- clu_second_name$Freq/m
                  clu_second_name <- clu_second_name[order(clu_second_name$Var1), ]

                  clu_second_article <- unique(clu_second[, c("pmid", "value")])
                  clu_second_article <- as.data.frame(table(clu_second_article$value), stringsAsFactors = FALSE)
                  m <- clu_second_article$Freq + 1
                  clu_second_article$Freq <- clu_second_article$Freq/m
                  clu_second_article <- clu_second_article[order(clu_second_article$Var1), ]
                  clu_second_name$Freq <- sqrt(clu_second_name$Freq * clu_second_article$Freq)
                  # matching cell second subtype2
                  clu_ann_second <- clu_ann[clu_ann$first == cellsubtype1, ]
                  clu_ann_second <- clu_ann_second[!is.na(clu_ann_second$subtype1), ]
                  clu_ann_secondname <- as.data.frame(table(clu_ann_second$subtype2), stringsAsFactors = FALSE)
                  if (nrow(clu_ann_secondname) > 0) {
                    m <- clu_ann_secondname$Freq + 1
                    clu_ann_secondname$Freq <- clu_ann_secondname$Freq/m
                    clu_ann_secondname <- clu_ann_secondname[order(clu_ann_secondname$Var1), ]
                    clu_ann_secondarticle <- unique(clu_ann_second[, c("pmid", "subtype2")])
                    clu_ann_secondarticle <- as.data.frame(table(clu_ann_secondarticle$subtype2), stringsAsFactors = FALSE)
                    m <- clu_ann_secondarticle$Freq + 1
                    clu_ann_secondarticle$Freq <- clu_ann_secondarticle$Freq/m
                    clu_ann_secondarticle <- clu_ann_secondarticle[order(clu_ann_secondarticle$Var1), ]
                    # cell second subtype2 scoring and compare with max second subtype2 score
                    clu_ann_secondname$Freq <- sqrt(clu_ann_secondname$Freq * clu_ann_secondarticle$Freq)
                    clu_ann_secondname <- clu_ann_secondname[clu_ann_secondname$Freq == max(clu_ann_secondname$Freq), ]
                    clu_ann_secondname <- clu_ann_secondname[(clu_ann_secondname$Freq > 0.5) & (clu_ann_secondname$Freq >= max(clu_second_name$Freq)), ]
                    if (nrow(clu_ann_secondname) == 1) {
                      cellsubtype2 <- clu_ann_secondname$Var1
                      clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2, ]$gene)
                      PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2, ]$pmid)
                      clu_ann_for_third <- clu_ann1[clu_ann1$variable == "subtype3", ]
                      # exist cell third subtype3
                      if (nrow(clu_ann_for_third) > 0) {
                        # calculting the third cell subtype3 max score of the same cell type -- short name.
                        clu_third <- clu_ann1[clu_ann1$variable == "subtype3", ]
                        clu_third_name <- as.data.frame(table(clu_third$value), stringsAsFactors = FALSE)
                        m <- clu_third_name$Freq + 1
                        clu_third_name$Freq <- clu_third_name$Freq/m
                        clu_third_name <- clu_third_name[order(clu_third_name$Var1), ]
                        clu_third_article <- unique(clu_third[, c("pmid", "value")])
                        clu_third_article <- as.data.frame(table(clu_third_article$value), stringsAsFactors = FALSE)
                        m <- clu_third_article$Freq + 1
                        clu_third_article$Freq <- clu_third_article$Freq/m
                        clu_third_article <- clu_third_article[order(clu_third_article$Var1), ]
                        clu_third_name$Freq <- sqrt(clu_third_name$Freq * clu_third_article$Freq)
                        # matching cell third subtype3
                        clu_ann_third <- clu_ann_second[clu_ann_second$subtype2 == cellsubtype2, ]
                        clu_ann_third <- clu_ann_third[!is.na(clu_ann_third$subtype2), ]
                        clu_ann_thirdname <- as.data.frame(table(clu_ann_third$subtype3), stringsAsFactors = FALSE)
                        if (nrow(clu_ann_thirdname) > 0) {
                          m <- clu_ann_thirdname$Freq + 1
                          clu_ann_thirdname$Freq <- clu_ann_thirdname$Freq/m
                          clu_ann_thirdname <- clu_ann_thirdname[order(clu_ann_thirdname$Var1), ]
                          clu_ann_thirdarticle <- unique(clu_ann_third[, c("pmid", "subtype3")])
                          clu_ann_thirdarticle <- as.data.frame(table(clu_ann_thirdarticle$subtype3), stringsAsFactors = FALSE)
                          m <- clu_ann_thirdarticle$Freq + 1
                          clu_ann_thirdarticle$Freq <- clu_ann_thirdarticle$Freq/m
                          clu_ann_thirdarticle <- clu_ann_thirdarticle[order(clu_ann_thirdarticle$Var1), ]
                          # cell third subtype3 scoring and compare with max first subtype score
                          clu_ann_thirdname$Freq <- sqrt(clu_ann_thirdname$Freq * clu_ann_thirdarticle$Freq)
                          clu_ann_thirdname <- clu_ann_thirdname[clu_ann_thirdname$Freq == max(clu_ann_thirdname$Freq), ]
                          clu_ann_thirdname <- clu_ann_thirdname[(clu_ann_thirdname$Freq > 0.5) & (clu_ann_thirdname$Freq >= max(clu_third_name$Freq)), ]
                          # Number of (cell third subtype3 with max score) >= 1 -- output -- last
                          if (nrow(clu_ann_thirdname) >= 1) {
                            cellsubtype3 <- clu_ann_thirdname$Var1
                            clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype3, ]$gene)
                            PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype3, ]$pmid)
                          }
                        }
                      }
                    }
                  }
                }
            }
            # max cell subtype label exist in subtype2
            if (unique(clu_ann_for_det$variable == "subtype2")) {
                cellsubtype2 <- unique(clu_ann_for_det$value)
                clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype2, ]$gene)
                PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype2, ]$pmid)
                clu_ann_for_first <- clu_ann1[clu_ann1$variable == "subtype1", ]
                # exist cell first subtype1
                if (nrow(clu_ann_for_first) > 0) {
                  # calculting the cell first subtype1 max score of the same cell type -- short name.
                  clu_first <- clu_ann1[clu_ann1$variable == "subtype1", ]
                  clu_first_name <- as.data.frame(table(clu_first$value), stringsAsFactors = FALSE)
                  m <- clu_first_name$Freq + 1
                  clu_first_name$Freq <- clu_first_name$Freq/m
                  clu_first_name <- clu_first_name[order(clu_first_name$Var1), ]
                  clu_first_article <- unique(clu_first[, c("pmid", "value")])
                  clu_first_article <- as.data.frame(table(clu_first_article$value), stringsAsFactors = FALSE)
                  m <- clu_first_article$Freq + 1
                  clu_first_article$Freq <- clu_first_article$Freq/m
                  clu_first_article <- clu_first_article[order(clu_first_article$Var1), ]
                  clu_first_name$Freq <- sqrt(clu_first_name$Freq * clu_first_article$Freq)
                  # matching cell first subtype1
                  clu_ann_first <- clu_ann[clu_ann$subtype2 == cellsubtype2, ]
                  clu_ann_first <- clu_ann_first[!is.na(clu_ann_first$subtype2), ]
                  clu_ann_firstname <- as.data.frame(table(clu_ann_first$subtype1), stringsAsFactors = FALSE)
                  if (nrow(clu_ann_firstname) > 0) {
                    m <- clu_ann_firstname$Freq + 1
                    clu_ann_firstname$Freq <- clu_ann_firstname$Freq/m
                    clu_ann_firstname <- clu_ann_firstname[order(clu_ann_firstname$Var1), ]
                    clu_ann_firstarticle <- unique(clu_ann_first[, c("pmid", "subtype1")])
                    clu_ann_firstarticle <- as.data.frame(table(clu_ann_firstarticle$subtype1), stringsAsFactors = FALSE)
                    m <- clu_ann_firstarticle$Freq + 1
                    clu_ann_firstarticle$Freq <- clu_ann_firstarticle$Freq/m
                    clu_ann_firstarticle <- clu_ann_firstarticle[order(clu_ann_firstarticle$Var1), ]
                    # cell first subtype1 scoring and compare with max first subtype1 score
                    clu_ann_firstname$Freq <- sqrt(clu_ann_firstname$Freq * clu_ann_firstarticle$Freq)
                    clu_ann_firstname <- clu_ann_firstname[clu_ann_firstname$Freq == max(clu_ann_firstname$Freq), ]
                    clu_ann_firstname <- clu_ann_firstname[(clu_ann_firstname$Freq > 0.5) & (clu_ann_firstname$Freq >= max(clu_ann_firstname$Freq)), ]
                    if (nrow(clu_ann_firstname) == 1) {
                      cellsubtype1 <- clu_ann_firstname$Var1
                      clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype1, ]$gene)
                      PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype1, ]$pmid)
                      clu_ann_for_third <- clu_ann1[clu_ann1$variable == "subtype3", ]
                      # exist cell third subtype3
                      if (nrow(clu_ann_for_third) > 0) {
                        # calculting the third cell subtype3 max score of the same cell type -- short name.
                        clu_third <- clu_ann1[clu_ann1$variable == "subtype3", ]
                        clu_third_name <- as.data.frame(table(clu_third$value), stringsAsFactors = FALSE)
                        m <- clu_third_name$Freq + 1
                        clu_third_name$Freq <- clu_third_name$Freq/m
                        clu_third_name <- clu_third_name[order(clu_third_name$Var1), ]
                        clu_third_article <- unique(clu_third[, c("pmid", "value")])
                        clu_third_article <- as.data.frame(table(clu_third_article$value), stringsAsFactors = FALSE)
                        m <- clu_third_article$Freq + 1
                        clu_third_article$Freq <- clu_third_article$Freq/m
                        clu_third_article <- clu_third_article[order(clu_third_article$Var1), ]
                        clu_third_name$Freq <- sqrt(clu_third_name$Freq * clu_third_article$Freq)
                        # matching cell third subtype3
                        clu_ann_third <- clu_ann_first[clu_ann_first$subtype1 == cellsubtype1, ]
                        clu_ann_third <- clu_ann_third[!is.na(clu_ann_third$subtype1), ]
                        clu_ann_thirdname <- as.data.frame(table(clu_ann_third$subtype3), stringsAsFactors = FALSE)
                        if (nrow(clu_ann_thirdname) > 0) {
                          m <- clu_ann_thirdname$Freq + 1
                          clu_ann_thirdname$Freq <- clu_ann_thirdname$Freq/m
                          clu_ann_thirdname <- clu_ann_thirdname[order(clu_ann_thirdname$Var1), ]
                          clu_ann_thirdarticle <- unique(clu_ann_third[, c("pmid", "subtype3")])
                          clu_ann_thirdarticle <- as.data.frame(table(clu_ann_thirdarticle$subtype3), stringsAsFactors = FALSE)
                          m <- clu_ann_thirdarticle$Freq + 1
                          clu_ann_thirdarticle$Freq <- clu_ann_thirdarticle$Freq/m
                          clu_ann_thirdarticle <- clu_ann_thirdarticle[order(clu_ann_thirdarticle$Var1), ]
                          # cell third subtype3 scoring and compare with max first subtype score
                          clu_ann_thirdname$Freq <- sqrt(clu_ann_thirdname$Freq * clu_ann_thirdarticle$Freq)
                          clu_ann_thirdname <- clu_ann_thirdname[clu_ann_thirdname$Freq == max(clu_ann_thirdname$Freq), ]
                          clu_ann_thirdname <- clu_ann_thirdname[(clu_ann_thirdname$Freq > 0.5) & (clu_ann_thirdname$Freq >=
                            max(clu_third_name$Freq)), ]
                          # Number of (cell third subtype3 with max score) >= 1 -- output -- last
                          if (nrow(clu_ann_thirdname) >= 1) {
                            cellsubtype3 <- clu_ann_thirdname$Var1
                            clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype3, ]$gene)
                            PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype3, ]$pmid)
                          }
                        }
                      }
                    }
                  }
                }
            }
            # max cell subtype label exist in subtype3
            if (unique(clu_ann_for_det$variable == "subtype3")) {
                cellsubtype3 <- unique(clu_ann_for_det$value)
                clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype3, ]$gene)
                PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype3, ]$pmid)
                clu_ann_for_second <- clu_ann1[clu_ann1$variable == "subtype2", ]
                # exist cell second subtype2
                if (nrow(clu_ann_for_second) > 0) {
                  # calculting the second cell subtype2 max score of the same cell type
                  clu_second <- clu_ann1[clu_ann1$variable == "subtype2", ]
                  clu_second_name <- as.data.frame(table(clu_second$value), stringsAsFactors = FALSE)
                  m <- clu_second_name$Freq + 1
                  clu_second_name$Freq <- clu_second_name$Freq/m
                  clu_second_name <- clu_second_name[order(clu_second_name$Var1), ]
                  clu_second_article <- unique(clu_second[, c("pmid", "value")])
                  clu_second_article <- as.data.frame(table(clu_second_article$value), stringsAsFactors = FALSE)
                  m <- clu_second_article$Freq + 1
                  clu_second_article$Freq <- clu_second_article$Freq/m
                  clu_second_article <- clu_second_article[order(clu_second_article$Var1), ]
                  clu_second_name$Freq <- sqrt(clu_second_name$Freq * clu_second_article$Freq)
                  # matching cell second subtype2
                  clu_ann_second <- clu_ann[clu_ann$subtype3 == cellsubtype3, ]
                  clu_ann_second <- clu_ann_second[!is.na(clu_ann_second$subtype3), ]
                  clu_ann_secondname <- as.data.frame(table(clu_ann_second$subtype2), stringsAsFactors = FALSE)
                  if (nrow(clu_ann_secondname) > 0) {
                    m <- clu_ann_secondname$Freq + 1
                    clu_ann_secondname$Freq <- clu_ann_secondname$Freq/m
                    clu_ann_secondname <- clu_ann_secondname[order(clu_ann_secondname$Var1), ]
                    clu_ann_secondarticle <- unique(clu_ann_second[, c("pmid", "subtype2")])
                    clu_ann_secondarticle <- as.data.frame(table(clu_ann_secondarticle$subtype2), stringsAsFactors = FALSE)
                    m <- clu_ann_secondarticle$Freq + 1
                    clu_ann_secondarticle$Freq <- clu_ann_secondarticle$Freq/m
                    clu_ann_secondarticle <- clu_ann_secondarticle[order(clu_ann_secondarticle$Var1), ]
                    # cell second subtype2 scoring and compare with max second subtype2 score
                    clu_ann_secondname$Freq <- sqrt(clu_ann_secondname$Freq * clu_ann_secondarticle$Freq)
                    clu_ann_secondname <- clu_ann_secondname[clu_ann_secondname$Freq == max(clu_ann_secondname$Freq), ]
                    clu_ann_secondname <- clu_ann_secondname[(clu_ann_secondname$Freq > 0.5) & (clu_ann_secondname$Freq >= max(clu_second_name$Freq)), ]
                    if (nrow(clu_ann_secondname) == 1) {
                      cellsubtype2 <- clu_ann_secondname$Var1
                      clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2, ]$gene)
                      PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2, ]$pmid)
                      clu_ann_for_first <- clu_ann1[clu_ann1$variable == "subtype1", ]
                      # exist cell first subtype1
                      if (nrow(clu_ann_for_first) > 0) {
                        # calculting the first cell subtype1 max score of the same cell type
                        clu_first <- clu_ann1[clu_ann1$variable == "subtype1", ]
                        clu_first_name <- as.data.frame(table(clu_first$value), stringsAsFactors = FALSE)
                        m <- clu_first_name$Freq + 1
                        clu_first_name$Freq <- clu_first_name$Freq/m
                        clu_first_name <- clu_first_name[order(clu_first_name$Var1), ]
                        clu_first_article <- unique(clu_first[, c("pmid", "value")])
                        clu_first_article <- as.data.frame(table(clu_first_article$value), stringsAsFactors = FALSE)
                        m <- clu_first_article$Freq + 1
                        clu_first_article$Freq <- clu_first_article$Freq/m
                        clu_first_article <- clu_first_article[order(clu_first_article$Var1), ]
                        clu_first_name$Freq <- sqrt(clu_first_name$Freq * clu_first_article$Freq)
                        # matching cell first subtype1
                        clu_ann_first <- clu_ann_second[clu_ann_second$subtype2 == cellsubtype2, ]
                        clu_ann_first <- clu_ann_first[!is.na(clu_ann_first$subtype2), ]
                        clu_ann_firstname <- as.data.frame(table(clu_ann_first$subtype1), stringsAsFactors = FALSE)
                        if (nrow(clu_ann_firstname) > 0) {
                          m <- clu_ann_firstname$Freq + 1
                          clu_ann_firstname$Freq <- clu_ann_firstname$Freq/m
                          clu_ann_firstname <- clu_ann_firstname[order(clu_ann_firstname$Var1), ]
                          clu_ann_firstarticle <- unique(clu_ann_first[, c("pmid", "subtype1")])
                          clu_ann_firstarticle <- as.data.frame(table(clu_ann_firstarticle$subtype1), stringsAsFactors = FALSE)
                          m <- clu_ann_firstarticle$Freq + 1
                          clu_ann_firstarticle$Freq <- clu_ann_firstarticle$Freq/m
                          clu_ann_firstarticle <- clu_ann_firstarticle[order(clu_ann_firstarticle$Var1), ]
                          # cell first subtype1 scoring and compare with max first subtype1 score
                          clu_ann_firstname$Freq <- sqrt(clu_ann_firstname$Freq * clu_ann_firstarticle$Freq)
                          clu_ann_firstname <- clu_ann_firstname[clu_ann_firstname$Freq == max(clu_ann_firstname$Freq), ]
                          clu_ann_firstname <- clu_ann_firstname[(clu_ann_firstname$Freq > 0.5) & (clu_ann_firstname$Freq >=
                            max(clu_first_name$Freq)), ]
                          # Number of (cell first subtype1 with max score) >= 1 -- output -- last
                          if (nrow(clu_ann_firstname) >= 1) {
                            cellsubtype1 <- clu_ann_firstname$Var1
                            clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype1, ]$gene)
                            PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype1, ]$pmid)
                          }
                        }
                      }
                    }
                  }
                }
            }
        }
    }
    res <- list(cellsubtype1 = cellsubtype1,cellsubtype2 = cellsubtype2,cellsubtype3 = cellsubtype3,celltype = celltype,
                celltype_score = celltype_score,clu_marker = clu_marker,PMID = PMID)
    return(res)
}
