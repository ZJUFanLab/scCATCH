#' Evidence-based score and annotation for each cluster
#'
#' @description Evidence-based score and annotation for each cluster by matching the potential marker genes generated from \code{\link{findmarkergenes}} with known cell marker genes in tissue-specific cell taxonomy reference database (CellMatch).
#' @param object The data.frame containing potential marker genes and the corresponding expressed cells percentage and average fold change for each cluster from the output of \code{\link{findmarkergenes}}.
#' @param species The specie of cells.'Human' or 'Mouse'.
#' @param cancer If the sample is from cancer tissue and you want to match cell marker genes of cancer tissues in CellMatch, then the cancer type may be defined. Select one or more related cancer types in 1.2 of \strong{Details} for human and 2.2 of \strong{Details} for mouse. The dafult is NULL for tissues without cancer.
#' @param tissue The tissue origin of cells. Select one or more related tissue types in \strong{Details}. For tissues without cancer, please refer to 1.1 of \strong{Details} for human tissue types and 2.1 of \strong{Details} for mouse tissue types. For tissues with cancer, please refer to 1.2 of \strong{Details} for human tissue types and 2.2 of \strong{Details} for mouse tissue types.
#' @return A data.frame containing matched cell type for each cluster, related marker genes, evidence-based score and PMID.
#' @examples clu_markers <- findmarkergenes(object = mouse_kidney_203_Seurat, species = 'Mouse')
#' clu_annotation <- scCATCH(object = clu_markers$clu_markers,species = 'Mouse',tissue = 'Kidney')
#' 
#' @details \strong{1.1} For \strong{Human} tissue, tissue types are listed as follows:
#' @details \strong{Adipose tissue-related}: Abdominal adipose tissue; Adipose tissue; Brown adipose tissue; Fat pad; Subcutaneous adipose tissue; Visceral adipose tissue; White adipose tissue.
#' @details \strong{Bladder-related}: Bladder; Urine.
#' @details \strong{Blood-related}: Blood; Peripheral blood; Plasma; Serum; Umbilical cord blood; Venous blood.
#' @details \strong{Bone-related}: Anterior cruciate ligament; Bone; Bone marrow; Cartilage; Intervertebral disc; Meniscus; Nucleus pulposus; Osteoarthritic cartilage; Periosteum; Skeletal muscle; Spinal cord; Synovial fluid; Synovium.
#' @details \strong{Brain-related}: Brain; Dorsolateral prefrontal cortex; Embryonic brain; Embryonic prefrontal cortex; Fetal brain; Hippocampus; Inferior colliculus; Midbrain; Sympathetic ganglion.
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
#' @details \strong{Uterus-related}: Endometrium; Endometrium stroma; Myometrium; Uterus; Vagina.
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
#' @details \strong{Calvaria-related}: Neonatal calvaria.
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
#' \url{https://github.com/ZJUFanLab/scCATCH}
#' @import data.table progress
#' @export scCATCH

scCATCH <- function(object, species = NULL, cancer = NULL, tissue = NULL) {
    if (is.null(species)) {
        stop("Please define the species! 'Human' or 'Mouse'.")
    }
    # cellmarkers matching species tissue
    CellMatch <- CellMatch[CellMatch$speciesType == species, ]
    if (!is.data.frame(object)) {
        if (names(object)[2] != "clu_markers") {
            stop("Please select the clu_clusters as the input!")
        }
        object <- object$clu_markers
    }
    
    # revise gene symbols of object
    object$genesymbol <- object$gene
    # tissue without cancer
    if (is.null(cancer)) {
        CellMatch <- CellMatch[CellMatch$cancerType == "Normal", ]
        # check tissue
        if (is.null(tissue)) {
            stop("Please define the origin tissue of cells! Select one or more related tissue types.")
        }
        tissue_match <- NULL
        # check the tissue
        tissue_match <- tissue %in% CellMatch$tissueType
        tissue_match <- which(tissue_match == "FALSE")
        if (length(tissue_match) > 0) {
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
            cat(paste("Warning: there are only ", length(cellmarker_num), " potential marker genes in CellMatch database for ", 
                species, " on ", tissue1, "!", sep = ""), "\n")
        }
        if (length(cellmarker_num) >= 100) {
            cat(paste("Note: there are ", length(cellmarker_num), " potential marker genes in CellMatch database for ", 
                species, " on ", tissue1, "!", sep = ""), "\n")
        }
        Sys.sleep(2)
        cellmarker_num <- unique(CellMatch$cellName)
        if (length(cellmarker_num) < 10) {
            cat(paste("Warning: there are only ", length(cellmarker_num), " cell types in CellMatch database for ", 
                species, " on ", tissue1, "!", sep = ""), "\n")
        }
        if (length(cellmarker_num) >= 10) {
            cat(paste("Note: there are ", length(cellmarker_num), " cell types in CellMatch database for ", 
                species, " on ", tissue1, "!", sep = ""), "\n")
        }
    }
    # tissue with cancer
    if (!is.null(cancer)) {
        # check cancer
        cancer_match <- NULL
        # check the tissue
        cancer_match <- cancer %in% CellMatch$cancerType
        cancer_match <- which(cancer_match == "FALSE")
        if (length(cancer_match) > 0) {
            stop(paste(cancer[cancer_match], ", not matched with the cancer types in CellMatch database! Please select one or more related cancer types.", 
                sep = ""))
        }
        CellMatch <- CellMatch[CellMatch$cancerType %in% cancer, ]
        # check tissue
        if (is.null(tissue)) {
            stop("Please define the origin tissue of cells! Select one or more related tissue types.")
        }
        tissue_match <- NULL
        # check the tissue
        tissue_match <- tissue %in% CellMatch$tissueType
        tissue_match <- which(tissue_match == "FALSE")
        if (length(tissue_match) > 0) {
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
            cat(paste("Warning: there are only ", length(cellmarker_num), " potential marker genes in CellMatch database for ", 
                species, " ", cancer1, " on ", tissue1, "!", sep = ""), "\n")
        }
        if (length(cellmarker_num) >= 100) {
            cat(paste("Note: there are ", length(cellmarker_num), " potential marker genes in CellMatch database for ", 
                species, " ", cancer1, " on ", tissue1, "!", sep = ""), "\n")
        }
        Sys.sleep(2)
        cellmarker_num <- unique(CellMatch$cellName)
        if (length(cellmarker_num) < 10) {
            cat(paste("Warning: there are only ", length(cellmarker_num), " cell types in CellMatch database for ", 
                species, " on ", tissue1, "!", sep = ""), "\n")
        }
        if (length(cellmarker_num) >= 10) {
            cat(paste("Note: there are ", length(cellmarker_num), " cell types in CellMatch database for ", 
                species, " on ", tissue1, "!", sep = ""), "\n")
        }
    }
    # extract cluster information
    clu.num <- as.character(unique(object$cluster))
    clu.num1 <- clu.num[nchar(clu.num) == 1]
    clu.num2 <- clu.num[nchar(clu.num) == 2]
    clu.num1 <- clu.num1[order(clu.num1)]
    clu.num2 <- clu.num2[order(clu.num2)]
    clu.num <- c(clu.num1, clu.num2)
    # Evidence-based scoring and annotation
    Sys.sleep(2)
    cat("Beginning evidence-based scoring and annotation", "\n")
    Sys.sleep(2)
    pb <- progress_bar$new(format = "[:bar] Finished::percent Remaining::eta", total = length(clu.num), 
        clear = FALSE, width = 60, complete = "+", incomplete = "-")
    
    clu_ann_res <- NULL
    for (i in 1:length(clu.num)) {
        Sys.sleep(1)
        # Filtering cluster marker
        rescluster_marker1 <- object[object$cluster == clu.num[i], ]$genesymbol
        cellsubtype1 <- "NA"
        cellsubtype2 <- "NA"
        cellsubtype3 <- "NA"
        celltype <- "NA"
        celltype_score <- "NA"
        clu_marker <- "NA"
        PMID <- "NA"
        res <- "I"
        # exists marker
        if (length(rescluster_marker1) > 0) {
            # matching CellMarker database
            rescluster_marker2 <- rescluster_marker1[rescluster_marker1 %in% CellMatch$geneSymbol]
            res <- "II"
            # exists matched marker
            if (length(rescluster_marker2) > 0) {
                # matching cell type
                clu_ann <- CellMatch[CellMatch$geneSymbol %in% rescluster_marker2, ]
                clu_ann_cellname <- as.data.frame(table(clu_ann$shortname), stringsAsFactors = F)
                m <- clu_ann_cellname$Freq + 1
                clu_ann_cellname$Freq <- clu_ann_cellname$Freq/m
                clu_ann_cellname <- clu_ann_cellname[order(clu_ann_cellname$Var1), ]
                # conuting related article number
                clu_ann_article <- unique(clu_ann[, c("PMID", "shortname")])
                clu_ann_article <- as.data.frame(table(clu_ann_article$shortname), stringsAsFactors = F)
                m <- clu_ann_article$Freq + 1
                clu_ann_article$Freq <- clu_ann_article$Freq/m
                clu_ann_article <- clu_ann_article[order(clu_ann_article$Var1), ]
                # scoring
                clu_ann_cellname$Freq <- sqrt(clu_ann_article$Freq * clu_ann_cellname$Freq)
                clu_ann_cellname <- clu_ann_cellname[clu_ann_cellname$Freq == max(clu_ann_cellname$Freq), 
                  ]
                # Number of (cell type with max score) > 1 -- output -- possilble
                if (nrow(clu_ann_cellname) > 1) {
                  celltype <- clu_ann_cellname$Var1
                  celltype_score <- round(clu_ann_cellname$Freq, digits = 2)
                  clu_marker <- unique(clu_ann[clu_ann$shortname %in% clu_ann_cellname$Var1, 
                    ]$geneSymbol)
                  PMID <- unique(clu_ann[clu_ann$shortname %in% clu_ann_cellname$Var1, ]$PMID)
                  res <- "III"
                }
                # Number of (cell type with max score) = 1 -- output -- determinate cell type short name
                if (nrow(clu_ann_cellname) == 1) {
                  celltype <- clu_ann_cellname$Var1
                  celltype_score <- round(clu_ann_cellname$Freq, digits = 2)
                  clu_marker <- unique(clu_ann[clu_ann$shortname %in% clu_ann_cellname$Var1, 
                    ]$geneSymbol)
                  PMID <- unique(clu_ann[clu_ann$shortname %in% clu_ann_cellname$Var1, ]$PMID)
                  res <- "IV"
                  # matching cell subtype to determine the start
                  clu_ann1 <- clu_ann[clu_ann$shortname == clu_ann_cellname$Var1, ]
                  clu_ann1 <- clu_ann1[, c("geneSymbol", "PMID", "third", "second", "first")]
                  clu_ann1$row <- 1:nrow(clu_ann1)
                  clu_ann1$row <- as.character(clu_ann1$row)
                  clu_ann1 <- clu_ann1[, c(6, 1:5)]
                  clu_ann1 <- data.table::as.data.table(clu_ann1)
                  clu_ann1 <- melt(data = clu_ann1, id.vars = 1:3, measure.vars = c("third", 
                    "second", "first"))
                  clu_ann1 <- as.data.frame(clu_ann1)
                  clu_ann1$variable <- as.character(clu_ann1$variable)
                  clu_ann1 <- clu_ann1[!is.na(clu_ann1$value), ]
                  
                  # exist cell subtype
                  if (nrow(clu_ann1) > 0) {
                    clu_ann_subcellname <- as.data.frame(table(clu_ann1$value), stringsAsFactors = F)
                    m <- clu_ann_subcellname$Freq + 1
                    clu_ann_subcellname$Freq <- clu_ann_subcellname$Freq/m
                    clu_ann_subcellname <- clu_ann_subcellname[order(clu_ann_subcellname$Var1), 
                      ]
                    
                    clu_ann_subarticle <- unique(clu_ann1[, c("PMID", "value")])
                    
                    clu_ann_subarticle <- as.data.frame(table(clu_ann_subarticle$value), stringsAsFactors = F)
                    m <- clu_ann_subarticle$Freq + 1
                    clu_ann_subarticle$Freq <- clu_ann_subarticle$Freq/m
                    clu_ann_subarticle <- clu_ann_subarticle[order(clu_ann_subarticle$Var1), 
                      ]
                    
                    clu_ann_subcellname$Freq <- sqrt(clu_ann_subcellname$Freq * clu_ann_subarticle$Freq)
                    
                    clu_ann_subcellname <- clu_ann_subcellname[clu_ann_subcellname$Freq == max(clu_ann_subcellname$Freq), 
                      ]
                    clu_ann_subcellname <- clu_ann_subcellname[clu_ann_subcellname$Freq > 0.5, 
                      ]
                    
                    # Number of (cell subtype with max score) == 1 -- determinate cell subtype first/second/third
                    if (nrow(clu_ann_subcellname) == 1) {
                      clu_ann_for_det <- clu_ann1[clu_ann1$value == clu_ann_subcellname$Var1, 
                        ]
                      
                      # max cell subtype label exist in first subtype1
                      if (unique(clu_ann_for_det$variable == "first")) {
                        cellsubtype1 <- unique(clu_ann_for_det$value)
                        clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype1, 
                          ]$geneSymbol)
                        PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype1, 
                          ]$PMID)
                        
                        clu_ann_for_second <- clu_ann1[clu_ann1$variable == "second", ]
                        
                        # exist cell second subtype2
                        if (nrow(clu_ann_for_second) > 0) {
                          
                          # calculting the second cell subtype2 max score of the same cell type -- short name.
                          clu_second <- clu_ann1[clu_ann1$variable == "second", ]
                          clu_second_name <- as.data.frame(table(clu_second$value), stringsAsFactors = F)
                          m <- clu_second_name$Freq + 1
                          clu_second_name$Freq <- clu_second_name$Freq/m
                          clu_second_name <- clu_second_name[order(clu_second_name$Var1), ]
                          
                          clu_second_article <- unique(clu_second[, c("PMID", "value")])
                          clu_second_article <- as.data.frame(table(clu_second_article$value), 
                            stringsAsFactors = F)
                          m <- clu_second_article$Freq + 1
                          clu_second_article$Freq <- clu_second_article$Freq/m
                          clu_second_article <- clu_second_article[order(clu_second_article$Var1), 
                            ]
                          
                          clu_second_name$Freq <- sqrt(clu_second_name$Freq * clu_second_article$Freq)
                          
                          # matching cell second subtype2
                          clu_ann_second <- clu_ann[clu_ann$first == cellsubtype1, ]
                          clu_ann_second <- clu_ann_second[!is.na(clu_ann_second$first), ]
                          
                          clu_ann_secondname <- as.data.frame(table(clu_ann_second$second), stringsAsFactors = F)
                          
                          if (nrow(clu_ann_secondname) > 0) {
                            
                            m <- clu_ann_secondname$Freq + 1
                            clu_ann_secondname$Freq <- clu_ann_secondname$Freq/m
                            clu_ann_secondname <- clu_ann_secondname[order(clu_ann_secondname$Var1), 
                              ]
                            
                            clu_ann_secondarticle <- unique(clu_ann_second[, c("PMID", "second")])
                            clu_ann_secondarticle <- as.data.frame(table(clu_ann_secondarticle$second), 
                              stringsAsFactors = F)
                            m <- clu_ann_secondarticle$Freq + 1
                            clu_ann_secondarticle$Freq <- clu_ann_secondarticle$Freq/m
                            clu_ann_secondarticle <- clu_ann_secondarticle[order(clu_ann_secondarticle$Var1), 
                              ]
                            
                            # cell second subtype2 scoring and compare with max second subtype2 score
                            clu_ann_secondname$Freq <- sqrt(clu_ann_secondname$Freq * clu_ann_secondarticle$Freq)
                            
                            clu_ann_secondname <- clu_ann_secondname[clu_ann_secondname$Freq == 
                              max(clu_ann_secondname$Freq), ]
                            clu_ann_secondname <- clu_ann_secondname[(clu_ann_secondname$Freq > 
                              0.5) & (clu_ann_secondname$Freq >= max(clu_second_name$Freq)), 
                              ]
                            
                            if (nrow(clu_ann_secondname) == 1) {
                              cellsubtype2 <- clu_ann_secondname$Var1
                              clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2, 
                                ]$geneSymbol)
                              PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2, ]$PMID)
                              
                              clu_ann_for_third <- clu_ann1[clu_ann1$variable == "third", ]
                              
                              # exist cell third subtype3
                              if (nrow(clu_ann_for_third) > 0) {
                                
                                # calculting the third cell subtype3 max score of the same cell type -- short name.
                                clu_third <- clu_ann1[clu_ann1$variable == "third", ]
                                clu_third_name <- as.data.frame(table(clu_third$value), stringsAsFactors = F)
                                m <- clu_third_name$Freq + 1
                                clu_third_name$Freq <- clu_third_name$Freq/m
                                clu_third_name <- clu_third_name[order(clu_third_name$Var1), 
                                  ]
                                
                                clu_third_article <- unique(clu_third[, c("PMID", "value")])
                                clu_third_article <- as.data.frame(table(clu_third_article$value), 
                                  stringsAsFactors = F)
                                m <- clu_third_article$Freq + 1
                                clu_third_article$Freq <- clu_third_article$Freq/m
                                clu_third_article <- clu_third_article[order(clu_third_article$Var1), 
                                  ]
                                
                                clu_third_name$Freq <- sqrt(clu_third_name$Freq * clu_third_article$Freq)
                                
                                # matching cell third subtype3
                                clu_ann_third <- clu_ann_second[clu_ann_second$second == cellsubtype2, 
                                  ]
                                clu_ann_third <- clu_ann_third[!is.na(clu_ann_third$second), 
                                  ]
                                
                                clu_ann_thirdname <- as.data.frame(table(clu_ann_third$third), 
                                  stringsAsFactors = F)
                                
                                if (nrow(clu_ann_thirdname) > 0) {
                                  
                                  m <- clu_ann_thirdname$Freq + 1
                                  clu_ann_thirdname$Freq <- clu_ann_thirdname$Freq/m
                                  clu_ann_thirdname <- clu_ann_thirdname[order(clu_ann_thirdname$Var1), 
                                    ]
                                  
                                  clu_ann_thirdarticle <- unique(clu_ann_third[, c("PMID", "third")])
                                  clu_ann_thirdarticle <- as.data.frame(table(clu_ann_thirdarticle$third), 
                                    stringsAsFactors = F)
                                  m <- clu_ann_thirdarticle$Freq + 1
                                  clu_ann_thirdarticle$Freq <- clu_ann_thirdarticle$Freq/m
                                  clu_ann_thirdarticle <- clu_ann_thirdarticle[order(clu_ann_thirdarticle$Var1), 
                                    ]
                                  
                                  # cell third subtype3 scoring and compare with max first subtype score
                                  clu_ann_thirdname$Freq <- sqrt(clu_ann_thirdname$Freq * clu_ann_thirdarticle$Freq)
                                  
                                  clu_ann_thirdname <- clu_ann_thirdname[clu_ann_thirdname$Freq == 
                                    max(clu_ann_thirdname$Freq), ]
                                  clu_ann_thirdname <- clu_ann_thirdname[(clu_ann_thirdname$Freq > 
                                    0.5) & (clu_ann_thirdname$Freq >= max(clu_third_name$Freq)), 
                                    ]
                                  
                                  # Number of (cell third subtype3 with max score) >= 1 -- output -- last
                                  if (nrow(clu_ann_thirdname) >= 1) {
                                    cellsubtype3 <- clu_ann_thirdname$Var1
                                    clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype3, 
                                      ]$geneSymbol)
                                    PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype3, 
                                      ]$PMID)
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      
                      # max cell subtype label exist in second subtype2
                      if (unique(clu_ann_for_det$variable == "second")) {
                        cellsubtype2 <- unique(clu_ann_for_det$value)
                        clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype2, 
                          ]$geneSymbol)
                        PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype2, 
                          ]$PMID)
                        
                        clu_ann_for_first <- clu_ann1[clu_ann1$variable == "first", ]
                        
                        # exist cell first subtype1
                        if (nrow(clu_ann_for_first) > 0) {
                          
                          # calculting the cell first subtype1 max score of the same cell type -- short name.
                          clu_first <- clu_ann1[clu_ann1$variable == "first", ]
                          clu_first_name <- as.data.frame(table(clu_first$value), stringsAsFactors = F)
                          m <- clu_first_name$Freq + 1
                          clu_first_name$Freq <- clu_first_name$Freq/m
                          clu_first_name <- clu_first_name[order(clu_first_name$Var1), ]
                          
                          clu_first_article <- unique(clu_first[, c("PMID", "value")])
                          clu_first_article <- as.data.frame(table(clu_first_article$value), 
                            stringsAsFactors = F)
                          m <- clu_first_article$Freq + 1
                          clu_first_article$Freq <- clu_first_article$Freq/m
                          clu_first_article <- clu_first_article[order(clu_first_article$Var1), 
                            ]
                          
                          clu_first_name$Freq <- sqrt(clu_first_name$Freq * clu_first_article$Freq)
                          
                          # matching cell first subtype1
                          clu_ann_first <- clu_ann[clu_ann$second == cellsubtype2, ]
                          clu_ann_first <- clu_ann_first[!is.na(clu_ann_first$second), ]
                          
                          clu_ann_firstname <- as.data.frame(table(clu_ann_first$first), stringsAsFactors = F)
                          
                          if (nrow(clu_ann_firstname) > 0) {
                            
                            m <- clu_ann_firstname$Freq + 1
                            clu_ann_firstname$Freq <- clu_ann_firstname$Freq/m
                            clu_ann_firstname <- clu_ann_firstname[order(clu_ann_firstname$Var1), 
                              ]
                            
                            clu_ann_firstarticle <- unique(clu_ann_first[, c("PMID", "first")])
                            clu_ann_firstarticle <- as.data.frame(table(clu_ann_firstarticle$first), 
                              stringsAsFactors = F)
                            m <- clu_ann_firstarticle$Freq + 1
                            clu_ann_firstarticle$Freq <- clu_ann_firstarticle$Freq/m
                            clu_ann_firstarticle <- clu_ann_firstarticle[order(clu_ann_firstarticle$Var1), 
                              ]
                            
                            # cell first subtype1 scoring and compare with max first subtype1 score
                            clu_ann_firstname$Freq <- sqrt(clu_ann_firstname$Freq * clu_ann_firstarticle$Freq)
                            
                            clu_ann_firstname <- clu_ann_firstname[clu_ann_firstname$Freq == 
                              max(clu_ann_firstname$Freq), ]
                            clu_ann_firstname <- clu_ann_firstname[(clu_ann_firstname$Freq > 
                              0.5) & (clu_ann_firstname$Freq >= max(clu_ann_firstname$Freq)), 
                              ]
                            
                            if (nrow(clu_ann_firstname) == 1) {
                              cellsubtype1 <- clu_ann_firstname$Var1
                              clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype1, 
                                ]$geneSymbol)
                              PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype1, ]$PMID)
                              
                              clu_ann_for_third <- clu_ann1[clu_ann1$variable == "third", ]
                              
                              # exist cell third subtype3
                              if (nrow(clu_ann_for_third) > 0) {
                                
                                # calculting the third cell subtype3 max score of the same cell type -- short name.
                                clu_third <- clu_ann1[clu_ann1$variable == "third", ]
                                clu_third_name <- as.data.frame(table(clu_third$value), stringsAsFactors = F)
                                m <- clu_third_name$Freq + 1
                                clu_third_name$Freq <- clu_third_name$Freq/m
                                clu_third_name <- clu_third_name[order(clu_third_name$Var1), 
                                  ]
                                
                                clu_third_article <- unique(clu_third[, c("PMID", "value")])
                                clu_third_article <- as.data.frame(table(clu_third_article$value), 
                                  stringsAsFactors = F)
                                m <- clu_third_article$Freq + 1
                                clu_third_article$Freq <- clu_third_article$Freq/m
                                clu_third_article <- clu_third_article[order(clu_third_article$Var1), 
                                  ]
                                
                                clu_third_name$Freq <- sqrt(clu_third_name$Freq * clu_third_article$Freq)
                                
                                # matching cell third subtype3
                                clu_ann_third <- clu_ann_first[clu_ann_first$first == cellsubtype1, 
                                  ]
                                clu_ann_third <- clu_ann_third[!is.na(clu_ann_third$first), ]
                                
                                clu_ann_thirdname <- as.data.frame(table(clu_ann_third$third), 
                                  stringsAsFactors = F)
                                
                                if (nrow(clu_ann_thirdname) > 0) {
                                  
                                  m <- clu_ann_thirdname$Freq + 1
                                  clu_ann_thirdname$Freq <- clu_ann_thirdname$Freq/m
                                  clu_ann_thirdname <- clu_ann_thirdname[order(clu_ann_thirdname$Var1), 
                                    ]
                                  
                                  clu_ann_thirdarticle <- unique(clu_ann_third[, c("PMID", "third")])
                                  clu_ann_thirdarticle <- as.data.frame(table(clu_ann_thirdarticle$third), 
                                    stringsAsFactors = F)
                                  m <- clu_ann_thirdarticle$Freq + 1
                                  clu_ann_thirdarticle$Freq <- clu_ann_thirdarticle$Freq/m
                                  clu_ann_thirdarticle <- clu_ann_thirdarticle[order(clu_ann_thirdarticle$Var1), 
                                    ]
                                  
                                  # cell third subtype3 scoring and compare with max first subtype score
                                  clu_ann_thirdname$Freq <- sqrt(clu_ann_thirdname$Freq * clu_ann_thirdarticle$Freq)
                                  
                                  clu_ann_thirdname <- clu_ann_thirdname[clu_ann_thirdname$Freq == 
                                    max(clu_ann_thirdname$Freq), ]
                                  clu_ann_thirdname <- clu_ann_thirdname[(clu_ann_thirdname$Freq > 
                                    0.5) & (clu_ann_thirdname$Freq >= max(clu_third_name$Freq)), 
                                    ]
                                  
                                  # Number of (cell third subtype3 with max score) >= 1 -- output -- last
                                  if (nrow(clu_ann_thirdname) >= 1) {
                                    cellsubtype3 <- clu_ann_thirdname$Var1
                                    clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype3, 
                                      ]$geneSymbol)
                                    PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype3, 
                                      ]$PMID)
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                      
                      # max cell subtype label exist in first subtype1
                      if (unique(clu_ann_for_det$variable == "third")) {
                        cellsubtype3 <- unique(clu_ann_for_det$value)
                        clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype3, 
                          ]$geneSymbol)
                        PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in% cellsubtype3, 
                          ]$PMID)
                        
                        clu_ann_for_second <- clu_ann1[clu_ann1$variable == "second", ]
                        
                        # exist cell second subtype2
                        if (nrow(clu_ann_for_second) > 0) {
                          
                          # calculting the second cell subtype2 max score of the same cell type -- short name.
                          clu_second <- clu_ann1[clu_ann1$variable == "second", ]
                          clu_second_name <- as.data.frame(table(clu_second$value), stringsAsFactors = F)
                          m <- clu_second_name$Freq + 1
                          clu_second_name$Freq <- clu_second_name$Freq/m
                          clu_second_name <- clu_second_name[order(clu_second_name$Var1), ]
                          
                          clu_second_article <- unique(clu_second[, c("PMID", "value")])
                          clu_second_article <- as.data.frame(table(clu_second_article$value), 
                            stringsAsFactors = F)
                          m <- clu_second_article$Freq + 1
                          clu_second_article$Freq <- clu_second_article$Freq/m
                          clu_second_article <- clu_second_article[order(clu_second_article$Var1), 
                            ]
                          
                          clu_second_name$Freq <- sqrt(clu_second_name$Freq * clu_second_article$Freq)
                          
                          # matching cell second subtype2
                          clu_ann_second <- clu_ann[clu_ann$third == cellsubtype3, ]
                          clu_ann_second <- clu_ann_second[!is.na(clu_ann_second$third), ]
                          
                          clu_ann_secondname <- as.data.frame(table(clu_ann_second$second), stringsAsFactors = F)
                          
                          if (nrow(clu_ann_secondname) > 0) {
                            
                            m <- clu_ann_secondname$Freq + 1
                            clu_ann_secondname$Freq <- clu_ann_secondname$Freq/m
                            clu_ann_secondname <- clu_ann_secondname[order(clu_ann_secondname$Var1), 
                              ]
                            
                            clu_ann_secondarticle <- unique(clu_ann_second[, c("PMID", "second")])
                            clu_ann_secondarticle <- as.data.frame(table(clu_ann_secondarticle$second), 
                              stringsAsFactors = F)
                            m <- clu_ann_secondarticle$Freq + 1
                            clu_ann_secondarticle$Freq <- clu_ann_secondarticle$Freq/m
                            clu_ann_secondarticle <- clu_ann_secondarticle[order(clu_ann_secondarticle$Var1), 
                              ]
                            
                            # cell second subtype2 scoring and compare with max second subtype2 score
                            clu_ann_secondname$Freq <- sqrt(clu_ann_secondname$Freq * clu_ann_secondarticle$Freq)
                            
                            clu_ann_secondname <- clu_ann_secondname[clu_ann_secondname$Freq == 
                              max(clu_ann_secondname$Freq), ]
                            clu_ann_secondname <- clu_ann_secondname[(clu_ann_secondname$Freq > 
                              0.5) & (clu_ann_secondname$Freq >= max(clu_second_name$Freq)), 
                              ]
                            
                            if (nrow(clu_ann_secondname) == 1) {
                              cellsubtype2 <- clu_ann_secondname$Var1
                              clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2, 
                                ]$geneSymbol)
                              PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2, ]$PMID)
                              
                              clu_ann_for_first <- clu_ann1[clu_ann1$variable == "first", ]
                              
                              # exist cell first subtype1
                              if (nrow(clu_ann_for_first) > 0) {
                                
                                # calculting the first cell subtype1 max score of the same cell type -- short name.
                                clu_first <- clu_ann1[clu_ann1$variable == "first", ]
                                clu_first_name <- as.data.frame(table(clu_first$value), stringsAsFactors = F)
                                m <- clu_first_name$Freq + 1
                                clu_first_name$Freq <- clu_first_name$Freq/m
                                clu_first_name <- clu_first_name[order(clu_first_name$Var1), 
                                  ]
                                
                                clu_first_article <- unique(clu_first[, c("PMID", "value")])
                                clu_first_article <- as.data.frame(table(clu_first_article$value), 
                                  stringsAsFactors = F)
                                m <- clu_first_article$Freq + 1
                                clu_first_article$Freq <- clu_first_article$Freq/m
                                clu_first_article <- clu_first_article[order(clu_first_article$Var1), 
                                  ]
                                
                                clu_first_name$Freq <- sqrt(clu_first_name$Freq * clu_first_article$Freq)
                                
                                # matching cell first subtype1
                                clu_ann_first <- clu_ann_second[clu_ann_second$second == cellsubtype2, 
                                  ]
                                clu_ann_first <- clu_ann_first[!is.na(clu_ann_first$second), 
                                  ]
                                
                                clu_ann_firstname <- as.data.frame(table(clu_ann_first$first), 
                                  stringsAsFactors = F)
                                
                                if (nrow(clu_ann_firstname) > 0) {
                                  
                                  m <- clu_ann_firstname$Freq + 1
                                  clu_ann_firstname$Freq <- clu_ann_firstname$Freq/m
                                  clu_ann_firstname <- clu_ann_firstname[order(clu_ann_firstname$Var1), 
                                    ]
                                  
                                  clu_ann_firstarticle <- unique(clu_ann_first[, c("PMID", "first")])
                                  clu_ann_firstarticle <- as.data.frame(table(clu_ann_firstarticle$first), 
                                    stringsAsFactors = F)
                                  m <- clu_ann_firstarticle$Freq + 1
                                  clu_ann_firstarticle$Freq <- clu_ann_firstarticle$Freq/m
                                  clu_ann_firstarticle <- clu_ann_firstarticle[order(clu_ann_firstarticle$Var1), 
                                    ]
                                  
                                  # cell first subtype1 scoring and compare with max first subtype1 score
                                  clu_ann_firstname$Freq <- sqrt(clu_ann_firstname$Freq * clu_ann_firstarticle$Freq)
                                  
                                  clu_ann_firstname <- clu_ann_firstname[clu_ann_firstname$Freq == 
                                    max(clu_ann_firstname$Freq), ]
                                  clu_ann_firstname <- clu_ann_firstname[(clu_ann_firstname$Freq > 
                                    0.5) & (clu_ann_firstname$Freq >= max(clu_first_name$Freq)), 
                                    ]
                                  
                                  # Number of (cell first subtype1 with max score) >= 1 -- output -- last
                                  if (nrow(clu_ann_firstname) >= 1) {
                                    cellsubtype1 <- clu_ann_firstname$Var1
                                    clu_marker <- unique(clu_ann1[clu_ann1$value %in% cellsubtype1, 
                                      ]$geneSymbol)
                                    PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype1, 
                                      ]$PMID)
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
            }
        }
        
        if (length(rescluster_marker1) == 0) {
            rescluster_marker1 <- "NA"
        }
        
        # processing the format of cell type
        d1 <- rescluster_marker1[1]
        if (length(rescluster_marker1) > 1) {
            for (j in 2:length(rescluster_marker1)) {
                d1 <- paste(d1, rescluster_marker1[j], sep = ", ")
            }
        }
        rescluster_marker1 <- d1
        
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
        
        clu_ann <- data.frame(cluster = clu.num[i], cluster_marker = rescluster_marker1, cellsubtype3 = cellsubtype3, 
            cellsubtype2 = cellsubtype2, cellsubtype1 = cellsubtype1, cell_type = celltype, celltype_score = celltype_score, 
            celltype_related_marker = clu_marker, PMID = PMID, stringsAsFactors = F)
        
        clu_ann_res <- rbind(clu_ann_res, clu_ann)
        pb$tick()
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
    cat("---Done---", "\n")
    Sys.sleep(2)
    return(clu_ann_res)
}
