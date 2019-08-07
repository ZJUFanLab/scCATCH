#' Find marker genes for each cluster
#'
#' @description Identify marker genes for each cluster from a Seurat object (>= 3.0.0) after log10 normaliztion and cluser analysis. The marker gene in each cluster is identified according to its expression level compared to it in every other clusters. Only significantly highly expressed one in all pair-wise comparison of the cluster will be selected as a cluster marker gene.
#' @param object Seurat object (>= 3.0.0) after log10 normalization and cluster analysis.
#' @param cell_min_pct Include the gene detected in at least this many cells in each cluster.
#' @param logfc Include the gene with at least this fold change of average gene expression compared to every other clusters.
#' @param pvalue Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters.
#' @return A data.frame containing cluster marker genes and the corresponding expressed cells percentage and average fold change for each cluster.
#' @examples findmarkergenes(mouse_kidney_203_Seurat, 0.25, 0.25, 0.05)
#' @import Seurat stats
#' @export findmarkergenes

findmarkergenes <- function(object, cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05) {
    
    # extract normalized data from Seurat object
    ndata <- as.data.frame(object[["RNA"]]@data)
    
    # extract cluster information of all single cells
    clu_info <- Seurat::Idents(object = object)
    clu_info <- data.frame(cell = names(clu_info), cluster = as.character(clu_info), 
        stringsAsFactors = F)
    clu_info[, 1] <- as.character(clu_info[, 1])
    clu_info[, 2] <- as.character(clu_info[, 2])
    
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    
    clu_num <- c(clu_num1, clu_num2)
    
    # generating pair-wise clusters
    clu_pair <- NULL
    for (i in 1:length(clu_num)) {
        d1 <- data.frame(cluster1 = rep(clu_num[i], (length(clu_num) - 1)), cluster2 = clu_num[-i], 
            stringsAsFactors = F)
        clu_pair <- rbind(clu_pair, d1)
    }
    
    # generating result file
    clu_marker <- NULL
    
    # calculate the p value and log fold change
    for (i in 1:length(clu_num)) {
        
        clu_pair1 <- clu_pair[clu_pair$cluster1 == clu_num[i], ]
        
        # generating result file for each cluster
        clu_marker1 <- NULL
        for (j in 1:nrow(clu_pair1)) {
            clu_marker2 <- as.data.frame(matrix(data = 0, nrow = nrow(ndata), ncol = 6))
            colnames(clu_marker2) <- c("cluster", "gene", "pct", "comp_cluster", 
                "avg_logfc", "pvalue")
            clu_marker2[, 2] <- rownames(ndata)
            
            # extract normalized data of pair wise clusters
            ndata1 <- ndata[, clu_info[clu_info$cluster == clu_pair1$cluster1[j], 
                ]$cell]
            ndata2 <- ndata[, clu_info[clu_info$cluster == clu_pair1$cluster2[j], 
                ]$cell]
            
            clu_marker_pct <- rep(0, nrow(clu_marker2))
            clu_marker_pvalue <- rep(1, nrow(clu_marker2))
            
            # assigning clusters and avg_logfc
            clu_marker2$cluster <- clu_pair1$cluster1[j]
            clu_marker2$comp_cluster <- clu_pair1$cluster2[j]
            clu_marker2$avg_logfc <- rowMeans(ndata1) - rowMeans(ndata2)
            
            # assigning pct and pvalue
            for (k in 1:nrow(clu_marker2)) {
                genedata1 <- as.numeric(ndata1[k, ])
                if (length(genedata1[genedata1 > 0]) >= (length(genedata1) * cell_min_pct)) {
                  clu_marker_pct[k] <- (length(genedata1[genedata1 > 0]))/(length(genedata1))
                  genedata2 <- as.numeric(ndata2[k, ])
                  clu_marker_pvalue[k] <- stats::wilcox.test(genedata1, genedata2, 
                    )$p.value
                }
            }
            clu_marker2$pct <- clu_marker_pct
            clu_marker2$pvalue <- clu_marker_pvalue
            
            # filtering genes with cell_min_pct logfc
            clu_marker2 <- clu_marker2[clu_marker2$pct >= cell_min_pct & clu_marker2$avg_logfc >= 
                logfc & clu_marker2$pvalue < pvalue, ]
            
            # combine the results for each cluster
            clu_marker1 <- rbind(clu_marker1, clu_marker2)
        }
        d1 <- as.data.frame(table(clu_marker1$gene), stringsAsFactors = F)
        d1 <- d1[d1$Freq == (length(clu_num) - 1), ]
        clu_marker1 <- clu_marker1[clu_marker1$gene %in% d1$Var1, ]
        clu_marker_gene <- unique(clu_marker1$gene)
        clu_marker_gene <- clu_marker_gene[order(clu_marker_gene)]
        avg_logfc <- NULL
        for (j in 1:length(clu_marker_gene)) {
            clu_marker_avg_logfc <- clu_marker1[clu_marker1$gene == clu_marker_gene[j], 
                ]
            avg_logfc[j] <- mean(clu_marker_avg_logfc$avg_logfc)
        }
        clu_marker1 <- unique(clu_marker1[, c("cluster", "gene", "pct")])
        clu_marker1 <- clu_marker1[order(clu_marker1$gene), ]
        clu_marker1$avg_logfc <- avg_logfc
        clu_marker <- rbind(clu_marker, clu_marker1)
    }
    rownames(clu_marker) <- 1:nrow(clu_marker)
    return(clu_marker)
}





