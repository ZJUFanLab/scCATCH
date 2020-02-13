#' Find potential marker genes for each cluster
#'
#' @description Identify potential marker genes for each cluster from a Seurat object (>= 3.0.0) after log10 normalization and cluster analysis. The potential marker genes in each cluster are identified according to its expression level compared to it in every other clusters. Only significantly highly expressed one in all pair-wise comparison of the cluster will be selected as a potential marker gene for the cluster. Genes will be revised according to NCBI Gene symbols (updated in Jan. 10, 2020, \url{https://www.ncbi.nlm.nih.gov/gene}) and no matched genes and duplicated genes will be removed. 
#' @param object Seurat object (>= 3.0.0) after log10 normalization and cluster analysis. Please ensure data is log10 normalized data and data has been clustered before running scCATCH pipeline.
#' @param species The specie of cells. The species must be defined. 'Human' or 'Mouse'.
#' @param cluster Select which clusters for potential marker genes identification. e.g. '1', '2', etc. The default is 'All' to find potential makrer genes for each cluster.
#' @param match_CellMatch For large datasets containg > 10,000 cells or > 15 clusters, it is strongly recommended to set match_CellMatch 'TRUE' to match CellMatch database first to include potential marker genes in terms of large system memory it may take.
#' @param cancer If match_CellMatch is set TRUE and the sample is from cancer tissue, then the cancer type may be defined. Select one or more related cancer types in 1.2 of \strong{Details} for human and 2.2 of \strong{Details} for mouse. The dafult is NULL for tissues without cancer.
#' @param tissue If match_CellMatch is set TRUE, then the tissue origin of cells must be defined. Select one or more related tissue types in \strong{Details}. For tissues without cancer, please refer to 1.1 of \strong{Details} for human tissue types and 2.1 of \strong{Details} for mouse tissue types. For tissues with cancer, please refer to 1.2 of \strong{Details} for human tissue types and 2.2 of \strong{Details} for mouse tissue types.
#' @param cell_min_pct Include the gene detected in at least this many cells in each cluster.
#' @param logfc Include the gene with at least this fold change of average gene expression compared to every other clusters.
#' @param pvalue Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters.
#' @return A list include a new data matrix wherein genes are revised by official gene symbols according to NCBI Gene symbols (updated in Jan. 10, 2020, \url{https://www.ncbi.nlm.nih.gov/gene}) and no matched genes and duplicated genes are removed as well as a data.frame containing potential marker genes of each selected cluster and the corresponding expressed cells percentage and average fold change for each cluster.
#' @examples clu_markers <- findmarkergenes(object = mouse_kidney_203_Seurat, species = 'Mouse')
#'
#' @examples clu_markers <- findmarkergenes(object = mouse_kidney_203_Seurat,
#'                                species = 'Mouse',
#'                                cluster = '1',
#'                                match_CellMatch = TRUE,
#'                                tissue = 'Kidney',
#'                                cell_min_pct = 0.1,
#'                                logfc = 0.1,
#'                                pvalue = 0.01)
#'
#' @examples clu_markers <- findmarkergenes(object = mouse_kidney_203_Seurat,
#'                                species = 'Mouse',
#'                                cluster = c('1','2'),
#'                                match_CellMatch = TRUE,
#'                                tissue = c('Kidney','Mesonephros'))
#' @details \strong{1.1} For \strong{Human} tissue, tissue types are listed as follows:
#' @details \strong{Adipose tissue-related}: Abdominal adipose tissue; Adipose tissue; Brown adipose tissue; Fat pad; Subcutaneous adipose tissue; Visceral adipose tissue; White adipose tissue.
#' @details \strong{Bladder-related}: Bladder; Urine.
#' @details \strong{Blood-related}: Blood; Peripheral blood; Plasma; Serum; Umbilical cord blood; Venous blood.
#' @details \strong{Bone-related}: Anterior cruciate ligament; Bone; Bone marrow; Cartilage; Intervertebral disc; Meniscus; Nucleus pulposus; Osteoarthritic cartilage; Periosteum; Skeletal muscle; Spinal cord; Synovial fluid; Synovium.
#' @details \strong{Brain-related}: Brain; Dorsolateral prefrontal cortex; Embryonic brain; Embryonic prefrontal cortex; Fetal brain; Hippocampus; Inferior colliculus; Midbrain.
#' @details \strong{Breast-related}: Breast; Mammary epithelium.
#' @details \strong{Embryo-related}: Embryo; Embryoid body; Embryonic brain; Embryonic prefrontal cortex; Embryonic stem cell; Germ; Primitive streak.
#' @details \strong{Esophagus-related}: Esophagus.
#' @details \strong{Eye-related}: Cornea; Corneal endothelium; Corneal epithelium; Eye; Lacrimal gland; Limbal epithelium; Optic nerve; Retina; Retinal pigment epithelium; Sclerocorneal tissue.
#' @details \strong{Fetus-related}: Amniotic fluid; Amniotic membrane; Fetal brain; Fetal gonad; Fetal kidney; Fetal liver; Fetal pancreas; Placenta; Umbilical cord; Umbilical cord blood; Umbilical vein.
#' @details \strong{Gonad-related}: Corpus luteum; Fetal gonad; Foreskin; Gonad; Ovarian cortex; Ovarian follicle; Ovary; Seminal plasma; Testis.
#' @details \strong{Hair-related}: Chorionic villus; Hair follicle; Scalp.
#' @details \strong{Heart-related}: Heart; Myocardium.
#' @details \strong{Intestine-related}: Colon; Colorectum; Gastrointestinal tract; Gut; Intestinal crypt; Intestine; Jejunum; Large intestine; Small intestinal crypt; Small intestine.
#' @details \strong{Kidney-related}: Adrenal gland; Fetal kidney; Kidney; Renal glomerulus.
#' @details \strong{Liver-related}: Fetal liver; Liver.
#' @details \strong{Lung-related}:Airway epithelium; Alveolus; Bronchoalveolar system; Lung.
#' @details \strong{Lymph-related}: Lymph; Lymph node; Lymphoid tissue.
#' @details \strong{Muscle-related}: Muscle; Skeletal muscle.
#' @details \strong{Nose-related}: Nasal concha; Nasal epithelium; Sinonasal mucosa.
#' @details \strong{Oral cavity-related}: Laryngeal squamous epithelium; Oral mucosa; Salivary gland; Sputum; Submandibular gland; Thyroid; Tonsil; Vocal fold. 
#' @details \strong{Ovary-related}: Corpus luteum; Ovarian cortex; Ovarian follicle; Ovary; Oviduct.
#' @details \strong{Pancreas-related}: Fetal pancreas; Pancreas; Pancreatic acinar tissue; Pancreatic islet.
#' @details \strong{Prostate-related}: Prostate.
#' @details \strong{Skin-related}: Dermis; Skin.
#' @details \strong{Spleen-related}: Spleen; Splenic red pulp.
#' @details \strong{Stomach-related}: Gastric corpus; Gastric epithelium; Gastric gland; Gastrointestinal tract; Pyloric gland; Stomach.
#' @details \strong{Testis-related}: Testis.
#' @details \strong{Tooth-related}: Deciduous tooth; Dental pulp; Gingiva; Molar; Periodontal ligament; Premolar; Tooth.
#' @details \strong{Uterus-related}: Endometrium; Endometrium stroma; Myometrium; Uterus.
#' @details \strong{Vessel-related}: Adventitia; Antecubital vein; Artery; Blood vessel; Umbilical vein.
#' @details \strong{Others}: Ascites; Epithelium; Ligament; Pluripotent stem cell; Thymus; Whartons jelly.
#' 
#' ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @details \strong{1.2} For \strong{Human} tissue about cancer, cancer types and the corresponding tissue types are listed as follows:
#' @details Acute Myelogenous Leukemia: Blood.
#' @details Acute Myeloid Leukemia: Bone marrow.
#' @details Adenoid Cystic Carcinoma: Salivary gland.
#' @details Alveolar Cell Carcinoma: Serum.
#' @details Astrocytoma: Brain.
#' @details B-Cell Lymphoma: Lymphoid tissue.
#' @details Bladder Cancer: Bladder.
#' @details Brain Cancer: Blood vessel; Brain.
#' @details Breast Cancer: Blood: Breast; Mammary gland.
#' @details Cholangiocarcinoma: Liver; Serum.
#' @details Chronic Lymphocytic Leukemia: Blood.
#' @details Chronic Myeloid Leukemia: Blood.
#' @details CNS Primitive Neuroectodermal Tumor: Brain.
#' @details Colon Cancer: Blood; Colon; Serum.
#' @details Colorectal Cancer: Blood; Colon; Colorectum; Gastrointestinal tract; Intestine; Liver; Lung; Venous blood.
#' @details Cutaneous Squamous Cell Carcinoma: Skin.
#' @details Endometrial Cancer: Endometrium.
#' @details Ependymoma: Brain.
#' @details Esophageal Adenocarcinoma: Esophagus.
#' @details Fibroid: Myometrium.
#' @details Follicular Lymphoma: Lymph node.
#' @details Gallbladder Cancer: Gall bladder; Gastrointestinal tract.
#' @details Gastric Cancer: Blood; Peripheral blood; Serum; Stomach.
#' @details Glioblastoma: Blood; Brain.
#' @details Glioma: Blood vessel; Brain.
#' @details Gonadoblastoma: Embryo.
#' @details Head and Neck Cancer: Blood; Brain; Oral cavity.
#' @details Hepatoblastoma: Liver.
#' @details Hepatocellular Cancer: Blood; Bone marrow; Embryo; Liver.
#' @details High-grade glioma: Brain.
#' @details Infantile Hemangiomas: Placenta.
#' @details Intestinal Cancer: Gastrointestinal tract.
#' @details Intracranial Aneurysm: Brain.
#' @details Kaposi's Sarcoma: Lymph node.
#' @details Larynx Cancer: Larynx.
#' @details Leukemia: Bone marrow; Peripheral blood.
#' @details Lipoma: Adipose tissue.
#' @details Liver Cancer: Blood; Liver.
#' @details Lung Adenocarcinoma: Lung.
#' @details Lung Cancer: Blood; Lung.
#' @details Lung Squamous Cell Carcinoma: Lung.
#' @details Lymphoma: Blood; Brain; Kidney; Liver; Lymph; Lymph node.
#' @details Malignant Insulinoma: Pancreas.
#' @details Malignant Mesothelioma: Lung; Pleura.
#' @details Malignant Peripheral Nerve Sheath Tumor: Brain.
#' @details Medulloblastoma: Brain.
#' @details Melanoma: Blood; Peripheral blood; Skin.
#' @details Mucoepidermoid Carcinoma: Salivary gland.
#' @details Multiple Myeloma: Bone marrow; Peripheral blood.
#' @details Myeloma: Bone marrow.
#' @details Natural Killer Cell Lymphoma: Lymph node.
#' @details Nephroblastoma: Kidney.
#' @details Non-Small Cell Lung Cancer: Blood; Lung; Peripheral blood.
#' @details Oesophageal Cancer: Blood.
#' @details Oligodendroglioma: Brain.
#' @details Oral Cancer: Oral cavity.
#' @details Oral Squamous Cell Carcinoma: Oral cavity; Peripheral blood.
#' @details Osteosarcoma: Bone.
#' @details Ovarian Cancer: Ascites; Ovarian cortex; Ovary; Peripheral blood.
#' @details Pancreatic Cancer: Blood vessel; Pancreas.
#' @details Pancreatic Ductal Adenocarcinomas: Pancreas.
#' @details Papillary Thyroid Carcinoma: Thyroid.
#' @details Prostate Cancer: Blood; Peripheral blood; Prostate.
#' @details Renal Cell Carcinoma: Kidney; Serum.
#' @details Renal Clear Cell Carcinoma: Lymph node.
#' @details Retinoblastoma: Eye.
#' @details Salivary Gland Tumor: Parotid gland; Salivary gland.
#' @details Sarcoma: Muscle.
#' @details Small Cell Lung Cancer: Lung.
#' @details Testicular Germ Cell Tumor: Peripheral blood; Testis.
#' @details Thyroid Cancer: Thyroid.
#' @details Tongue Cancer: Tongue.
#' @details Uterine Leiomyoma: Uterus.
#' @details Vascular Tumour: Lymph node.
#' 
#' ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @details \strong{2.1} For \strong{Mouse} tissue, tissue types are listed as follows:
#' @details \strong{Adipose tissue-related}: Adipose tissue; White adipose tissue.
#' @details \strong{Bladder-related}: Bladder.
#' @details \strong{Blood-related}: Blood; Peripheral blood; Serum; Umbilical cord blood.
#' @details \strong{Bone-related}: Bone; Bone marrow; Meniscus; Skeletal muscle; Spinal cord.
#' @details \strong{Brain-related}: Brain; Cerebellum; Fetal brain; Hippocampus; Neural tube.
#' @details \strong{Breast-related}: Mammary epithelium; Mammary gland.
#' @details \strong{Calvaria-related}: Calvaria.
#' @details \strong{Ear-related}: Cochlea; Inner Ear.
#' @details \strong{Embryo-related}: Embryo; Embryoid body; Embryonic heart; Embryonic stem cell.
#' @details \strong{Esophagus-related}: Esophagus.
#' @details \strong{Eye-related}: Corneal epithelium; Eye; Ganglion cell layer of retina; Inner nuclear layer of retina; Lacrimal gland; Retina.
#' @details \strong{Fetus-related}: Fetal brain; Fetal intestine; Fetal liver; Fetal lung; Fetal stomach; Placenta; Umbilical cord; Umbilical cord blood.
#' @details \strong{Gonad-related}: Gonad; Ovary; Testis; Yolk sac.
#' @details \strong{Hair-related}: Hair follicle.
#' @details \strong{Heart-related}: Embryonic heart; Heart; Heart muscle; Neonatal heart.
#' @details \strong{Intestine-related}: Colon; Colon epithelium; Fetal intestine; Gastrointestinal tract; Ileum; Intestinal crypt; Intestine; Mesenteric lymph node; Small intestine.
#' @details \strong{Kidney-related}: Kidney; Mesonephros.
#' @details \strong{Liver-related}: Fetal liver; Liver.
#' @details \strong{Lung-related}: Bronchiole; Fetal lung; Lung; Trachea.
#' @details \strong{Lymph-related}: Lymph node; Lymphoid tissue; Mesenteric lymph node; Peyer patch.
#' @details \strong{Muscle-related}: Heart muscle; Muscle; Neonatal muscle; Skeletal muscle.
#' @details \strong{Neonate-related}: Neonatal calvaria; Neonatal heart; Neonatal muscle; Neonatal pancreas; Neonatal rib; Neonatal skin.
#' @details \strong{Oral cavity-related}: Submandibular gland; Taste bud. 
#' @details \strong{Ovary-related}: Ovary; Yolk sac.
#' @details \strong{Pancreas-related}: Neonatal pancreas; Pancreas; Pancreatic islet.
#' @details \strong{Prostate-related}: Prostate.
#' @details \strong{Skin-related}: Dermis; Epidermis; Neonatal skin; Skin.
#' @details \strong{Spleen-related}: Spleen.
#' @details \strong{Stomach-related}: Fetal stomach; Gastrointestinal tract; Stomach.
#' @details \strong{Testis-related}: Testis.
#' @details \strong{Uterus-related}: Uterus.
#' @details \strong{Vessel-related}: Aorta; Artery; Blood vessel; Carotid artery.
#' @details \strong{Others}: Basilar membrane; Epithelium; Peritoneal cavity; Thymus.
#' 
#' ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @details \strong{2.2} For \strong{Mouse} tissue about cancer, cancer types and the corresponding tissue types are listed as follows:
#' @details Breast Cancer: Lymph node; Breast; Lung.
#' @details Chronic Myeloid Leukemia: Blood.
#' @details Colon Cancer: Colon.
#' @details Colorectal Cancer: Lymph node; Colon; Colorectum.
#' @details Cutaneous Squamous Cell Carcinoma: Skin.
#' @details Hepatocellular Cancer: Blood; Liver.
#' @details Liver Cancer: Liver.
#' @details Lung Cancer: Lung.
#' @details Melanoma: Lung.
#' @details Pancreatic Cancer: Blood.
#' @details Papillary Thyroid Carcinoma: Thyroid.
#' @details Prostate Cancer: Prostate.
#' @details Renal Cell Carcinoma: Kidney.
#' @details Supratentorial Primitive Neuroectodermal Tumor: Brain.
#' @seealso \url{https://github.com/ZJUFanLab/scCATCH}
#' @import Seurat stats
#' @export findmarkergenes

findmarkergenes <- function(object, species = NULL, cluster = "All", match_CellMatch = FALSE, cancer = NULL, 
    tissue = NULL, cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05) {
    # extract normalized data from Seurat object
    ndata <- as.data.frame(object[["RNA"]]@data)
    # extract cluster information of all single cells
    clu_info <- Seurat::Idents(object = object)
    clu_info <- data.frame(cell = names(clu_info), cluster = as.character(clu_info), stringsAsFactors = F)
    clu_info[, 1] <- as.character(clu_info[, 1])
    clu_info[, 2] <- as.character(clu_info[, 2])
    clu_num <- unique(clu_info[, 2])
    clu_num1 <- clu_num[nchar(clu_num) == 1]
    clu_num2 <- clu_num[nchar(clu_num) == 2]
    clu_num1 <- clu_num1[order(clu_num1)]
    clu_num2 <- clu_num2[order(clu_num2)]
    clu_num <- c(clu_num1, clu_num2)
    rm(object)
    cat("Note: the raw data matrix includes", ncol(ndata), "cells and", nrow(ndata), "genes.", "\n")
    Sys.sleep(2)
    clu_num1 <- clu_num
    if (is.null(species)) {
        cat('\n')
        stop("Please define the species! 'Human' or 'Mouse'.")
    }
    geneinfo <- geneinfo[geneinfo$species == species, ]
    # revise gene symbol
    cat('\n')
    cat("---Revising gene symbols according to NCBI Gene symbols (updated in Jan. 10, 2020, https://www.ncbi.nlm.nih.gov/gene) and no matched genes and duplicated genes will be removed.", 
        "\n")
    Sys.sleep(2)
    genename <- rownames(ndata)
    genename1 <- genename[genename %in% geneinfo$Symbol]
    genename2 <- genename[!genename %in% geneinfo$Symbol]
    genename3 <- genename2[genename2 %in% geneinfo$Synonyms]
    genename4 <- rep("NA", length(genename3))
    for (i in 1:length(genename3)) {
        d1 <- geneinfo[geneinfo$Synonyms == genename3[i], ]$Symbol
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
    ndata <- ndata[genedata$raw_name, ]
    rownames(ndata) <- genedata$new_name
    cat('\n')
    cat("Note: the new data matrix includes", ncol(ndata), "cells and", nrow(ndata), "genes.", "\n")
    Sys.sleep(2)
    if (match_CellMatch) {
        ndata3 <- ndata
        # tissue without cancer
        if (is.null(cancer)) {
            CellMatch <- CellMatch[CellMatch$speciesType == species & CellMatch$cancerType == "Normal", 
                ]
            # check tissue
            if (is.null(tissue)) {
                cat('\n')
                stop("Please define the origin tissue of cells! Select one or more related tissue types.")
            }
            tissue_match <- NULL
            # check the tissue
            tissue_match <- tissue %in% CellMatch$tissueType
            tissue_match <- which(tissue_match == "FALSE")
            if (length(tissue_match) > 0) {
                cat('\n')
                stop(paste(tissue[tissue_match], ", not matched with the tissue types in CellMatch database! Please select one or more related tissue types.", 
                  sep = ""))
            }
            CellMatch <- CellMatch[CellMatch$tissueType %in% tissue, ]
            cellmarker_num <- CellMatch$geneSymbol
            tissue1 <- tissue[1]
            if (length(tissue) > 1) {
                for (i in 2:length(tissue)) {
                  tissue1 <- paste(tissue1, tissue[i], sep = ", ")
                }
            }
            if (length(cellmarker_num) < 100) {
                cat('\n')
                cat(paste("Warning: there are only ", length(cellmarker_num), " potential marker genes in CellMatch database for ", 
                  species, " on ", tissue1, "!", sep = ""), "\n")
            }
            if (length(cellmarker_num) >= 100) {
                cat('\n')
                cat(paste("Note: there are ", length(cellmarker_num), " potential marker genes in CellMatch database for ", 
                  species, " on ", tissue1, ".", sep = ""), "\n")
            }
            ndata <- ndata[rownames(ndata) %in% cellmarker_num, ]
        }
        # tissue with cancer
        if (!is.null(cancer)) {
            CellMatch <- CellMatch[CellMatch$speciesType == species, ]
            # check cancer
            cancer_match <- NULL
            # check the tissue
            cancer_match <- cancer %in% CellMatch$cancerType
            cancer_match <- which(cancer_match == "FALSE")
            if (length(cancer_match) > 0) {
                cat('\n')
                stop(paste(cancer[cancer_match], ", not matched with the cancer types in CellMatch database! Please select one or more related cancer types.", 
                  sep = ""))
            }
            CellMatch <- CellMatch[CellMatch$cancerType %in% cancer, ]
            # check tissue
            if (is.null(tissue)) {
                cat('\n')
                stop("Please define the origin tissue of cells! Select one or more related tissue types.")
            }
            tissue_match <- NULL
            # check the tissue
            tissue_match <- tissue %in% CellMatch$tissueType
            tissue_match <- which(tissue_match == "FALSE")
            if (length(tissue_match) > 0) {
                cat('\n')
                stop(paste(tissue[tissue_match], ", not matched with the tissue types in CellMatch database! Please select one or more related tissue types.", 
                  sep = ""))
            }
            CellMatch <- CellMatch[CellMatch$tissueType %in% tissue, ]
            cellmarker_num <- CellMatch$geneSymbol
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
            if (length(cellmarker_num) < 100) {
                cat('\n')
                cat(paste("Warning: there are only ", length(cellmarker_num), " potential marker genes in CellMatch database for ", 
                  species, " ", cancer1, " on ", tissue1, "!", sep = ""), "\n")
            }
            if (length(cellmarker_num) >= 100) {
                cat('\n')
                cat(paste("Note: there are ", length(cellmarker_num), " potential marker genes in CellMatch database for ", 
                  species, " ", cancer1, " on ", tissue1, ".", sep = ""), "\n")
            }
            ndata <- ndata[rownames(ndata) %in% cellmarker_num, ]
        }
    }
    if (nrow(ndata) == 0) {
        cat('\n')
        stop("There is no matched potential marker genes in the matrix! Please try again to find differential expressed genes by setting match_CellMatch as FALSE!")
    }
    # generating pair-wise clusters
    clu_pair <- NULL
    for (i in 1:length(clu_num)) {
        d1 <- data.frame(cluster1 = rep(clu_num[i], (length(clu_num) - 1)), cluster2 = clu_num[-i], 
            stringsAsFactors = F)
        clu_pair <- rbind(clu_pair, d1)
    }
    if (cluster[1] != "All") {
        cluster_match <- cluster %in% clu_num
        cluster_match <- which(cluster_match == "FALSE")
        if (length(cluster_match) > 0) {
            cat('\n')
            stop(paste(cluster[cluster_match], ", not matched with the cell clusters! Please select one or more related clusters.", 
                sep = ""))
        }
        clu_pair <- clu_pair[clu_pair$cluster1 %in% cluster, ]
        clu_num1 <- clu_num[clu_num %in% cluster]
    }
    # generating result file
    clu_marker <- NULL
    Sys.sleep(2)
    cat('\n')
    # calculate the p value and log fold change
    for (i in 1:length(clu_num1)) {
        cat(paste("Finding potential marker genes for cluster", clu_num1[i], sep = " "), "\n")
        clu_pair1 <- clu_pair[clu_pair$cluster1 == clu_num[i], ]
        # generating result file for each cluster
        clu_marker1 <- NULL
        for (j in 1:nrow(clu_pair1)) {
            clu_marker2 <- as.data.frame(matrix(data = 0, nrow = nrow(ndata), ncol = 6))
            colnames(clu_marker2) <- c("cluster", "gene", "pct", "comp_cluster", "avg_logfc", "pvalue")
            clu_marker2[, 2] <- rownames(ndata)
            # extract normalized data of pair wise clusters
            ndata1 <- ndata[, clu_info[clu_info$cluster == clu_pair1$cluster1[j], ]$cell]
            ndata2 <- ndata[, clu_info[clu_info$cluster == clu_pair1$cluster2[j], ]$cell]
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
                  clu_marker_pvalue[k] <- stats::wilcox.test(genedata1, genedata2, )$p.value
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
        if (nrow(clu_marker1) > 0) {
            d1 <- as.data.frame(table(clu_marker1$gene), stringsAsFactors = F)
            d1 <- d1[d1$Freq == (length(clu_num) - 1), ]
            clu_marker1 <- clu_marker1[clu_marker1$gene %in% d1$Var1, ]
            clu_marker_gene <- unique(clu_marker1$gene)
            clu_marker_gene <- clu_marker_gene[order(clu_marker_gene)]
            avg_logfc <- NULL
            for (j in 1:length(clu_marker_gene)) {
                clu_marker_avg_logfc <- clu_marker1[clu_marker1$gene == clu_marker_gene[j], ]
                avg_logfc[j] <- mean(clu_marker_avg_logfc$avg_logfc)
            }
            clu_marker1 <- unique(clu_marker1[, c("cluster", "gene", "pct")])
            clu_marker1 <- clu_marker1[order(clu_marker1$gene), ]
            clu_marker1$avg_logfc <- avg_logfc
            clu_marker <- rbind(clu_marker, clu_marker1)
        }
        cat("***Done***", "\n")
        Sys.sleep(2)
    }
    rm(ndata1, ndata2)
    res <- list()
    if (match_CellMatch) {
        res[[1]] <- ndata3
        rm(ndata3)
    }
    if (!match_CellMatch) {
        res[[1]] <- ndata
        rm(ndata)
    }
    res[[2]] <- "NA"
    names(res) <- c("new_data_matrix", "clu_markers")
    if (is.null(clu_marker)) {
        cat('\n')
        cat("Warning: there is no potential marker genes found in the matrix! You may adjust the cell clusters by applying other clustering algorithm and try again.")
        return(res)
        stop()
    }
    if (nrow(clu_marker) == 0) {
        cat('\n')
        cat("Warning: there is no potential marker genes found in the matrix! You may adjust the cell clusters by applying other clustering algorithm and try again.")
        return(res)
        stop()
    }
    if (nrow(clu_marker) > 0) {
        res[[2]] <- clu_marker
        return(res)
    }
}
