#' Evidence-based score and annotation for each cluster
#'
#' @description Evidence-based score and annotation for each cluster generated from \code{\link{findmarkergenes}} by matching the marker genes with known cell markers in tissue-specific cell taxonomy reference database.
#' @param object The data.frame containing marker genes and the corresponding expressed cells percentage and average fold change for each cluster from the output of \code{\link{findmarkergenes}}.
#' @param species The specie of cells.'Human' or 'Mouse'.
#' @param tissue The tissue origin of cells. Select one or more related tissue types in Details
#' @return A data.frame containing matched cell type for each cluster, related marker genes, evidence-based score and PMID.
#' @examples
#' clu_markers <- findmarkergenes(mouse_kidney_203_Seurat, 0.25, 0.25, 0.05)
#' clu_ann <- scCATCH(clu_markers, 'Mouse', 'Kidney')
#' @details \strong{For human}, tissue include Abdominal adipose tissue, Adipose tissue, Adrenal gland, Adventitia, Airway epithelium, Alveolus, Amniotic fluid, Amniotic membrane, Antecubital vein, Anterior cruciate ligament, Artery, Ascites, Bladder, Blood, Blood vessel, Bone, Bone marrow, Brain, Breast, Bronchoalveolar system, Brown adipose tissue, Cartilage, Chorionic villus, Colon, Colorectum, Cornea, Corneal endothelium, Corneal epithelium, Corpus luteum, Deciduous tooth, Dental pulp, Dermis, Dorsolateral prefrontal cortex, Embryo, Embryoid body, Embryonic brain, Embryonic prefrontal cortex, Embryonic stem cell, Endometrium, Endometrium stroma, Epithelium, Esophagus, Eye, Fat pad, Fetal brain, Fetal gonad, Fetal kidney, Fetal liver, Fetal pancreas, Foreskin, Gastric corpus, Gastric epithelium, Gastric gland, Gastrointestinal tract, Germ, Gingiva, Gonad, Gut, Hair follicle, Heart, Hippocampus, Inferior colliculus, Intervertebral disc, Intestinal crypt, Intestine, Jejunum, Kidney, Lacrimal gland, Large intestine, Laryngeal squamous epithelium, Ligament, Limbal epithelium, Liver, Lung, Lymph, Lymph node, Lymphoid tissue, Mammary epithelium, Meniscus, Midbrain, Molar, Muscle, Myocardium, Myometrium, Nasal concha, Nasal epithelium, Nucleus pulposus, Optic nerve, Oral mucosa, Osteoarthritic cartilage, Ovarian cortex, Ovarian follicle, Ovary, Oviduct, Pancreas, Pancreatic acinar tissue, Pancreatic islet, Periodontal ligament, Periosteum, Peripheral blood, Placenta, Plasma, Pluripotent stem cell, Premolar, Primitive streak, Prostate, Pyloric gland, Renal glomerulus, Retina, Retinal pigment epithelium, Salivary gland, Scalp, Sclerocorneal tissue, Seminal plasma, Serum, Sinonasal mucosa, Skeletal muscle, Skin, Small intestinal crypt, Small intestine, Spinal cord, Spleen, Splenic red pulp, Sputum, Stomach, Subcutaneous adipose tissue, Submandibular gland, Sympathetic ganglion, Synovial fluid, Synovium, Testis, Thymus, Thyroid, Tonsil, Tooth, Umbilical cord, Umbilical cord blood, Umbilical vein, Undefined, Urine, Uterus, Vagina, Venous blood, Visceral adipose tissue, Vocal fold, Whartons jelly, White adipose tissue;
#' \strong{For Mouse}, Adipose tissue, Aorta, Artery, Basilar membrane, Bladder, Blood, Blood vessel, Bone, Bone marrow, Brain, Bronchiole, Carotid artery, Cerebellum, Cochlea, Colon, Colon epithelium, Corneal epithelium, Dermis, Embryo, Embryoid body, Embryonic heart, Embryonic stem cell, Epidermis, Epithelium, Esophagus, Eye, Fetal liver, Ganglion cell layer of retina, Gastrointestinal tract, Gonad, Hair follicle, Heart, Heart muscle, Hippocampus, Ileum, Inner Ear, Inner nuclear layer of retina, Intestinal crypt, Intestine, Kidney, Lacrimal gland, Liver, Lung, Lymph node, Lymphoid tissue, Mammary epithelium, Mammary gland, Meniscus, Mesenteric lymph node, Mesonephros, Muscle, Neural tube, Ovary, Pancreas, Pancreatic islet, Peripheral blood, Peritoneal cavity, Peyer patch, Prostate, Retina, Serum, Skeletal muscle, Skin, Small intestine, Spinal cord, Spleen, Stomach, Submandibular gland, Taste bud, Testis, Thymus, Trachea, Umbilical cord, Umbilical cord blood, Undefined, White adipose tissue, Yolk sac.
#' @import data.table
#' @export scCATCH

scCATCH <- function(object, species, tissue) {

    # cellmarkers matching species tissue
    cell_markers1 <- cellmarkers[cellmarkers$speciesType %in% species & cellmarkers$tissueType %in%
        tissue, ]
    cell_markers1<- cell_markers1[cell_markers1$cellType == 'Normal cell',]
    
    # geneinfo matching species
    humangeneinfo <- geneinfo[geneinfo$species == species, ]

    # object gene standarlization
    object$genesymbol <- object$gene

    for (k in 1:nrow(object)) {
        d1 <- object[k, ]
        if (!d1$genesymbol %in% humangeneinfo$Symbol) {
            d2 <- humangeneinfo[humangeneinfo$Synonyms == d1$genesymbol, ]
            if (nrow(d2) == 0) {
                object[k, "genesymbol"] <- "NA"
            } else {
                if (length(unique(d2$Symbol)) == 1) {
                  object[k, "genesymbol"] <- unique(d2$Symbol)
                }
            }
        }
    }
    object <- object[!object$genesymbol == "NA", ]

    clu_ann_res <- NULL
    clu.num <- as.character(unique(object$cluster))
    clu.num1 <- clu.num[nchar(clu.num) == 1]
    clu.num2 <- clu.num[nchar(clu.num) == 2]

    clu.num1 <- clu.num1[order(clu.num1)]
    clu.num2 <- clu.num2[order(clu.num2)]

    clu.num <- c(clu.num1, clu.num2)

    for (i in 1:length(clu.num)) {
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
            rescluster_marker2 <- rescluster_marker1[rescluster_marker1 %in% cell_markers1$geneSymbol]
            res <- "II"

            # exists matched marker
            if (length(rescluster_marker2) > 0) {

                # matching cell type
                clu_ann <- cell_markers1[cell_markers1$geneSymbol %in% rescluster_marker2,
                  ]

                clu_ann_cellname <- as.data.frame(table(clu_ann$shortname), stringsAsFactors = F)
                m <- clu_ann_cellname$Freq + 1
                clu_ann_cellname$Freq <- clu_ann_cellname$Freq/m
                clu_ann_cellname <- clu_ann_cellname[order(clu_ann_cellname$Var1),
                  ]

                # conuting related article number
                clu_ann_article <- unique(clu_ann[, c("PMID", "shortname")])
                clu_ann_article <- as.data.frame(table(clu_ann_article$shortname),
                  stringsAsFactors = F)
                m <- clu_ann_article$Freq + 1
                clu_ann_article$Freq <- clu_ann_article$Freq/m
                clu_ann_article <- clu_ann_article[order(clu_ann_article$Var1),
                  ]

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
                  PMID <- unique(clu_ann[clu_ann$shortname %in% clu_ann_cellname$Var1,
                    ]$PMID)
                  res <- "III"
                }

                # Number of (cell type with max score) = 1 -- output -- determinate cell type
                # short name
                if (nrow(clu_ann_cellname) == 1) {
                  celltype <- clu_ann_cellname$Var1
                  celltype_score <- round(clu_ann_cellname$Freq, digits = 2)
                  clu_marker <- unique(clu_ann[clu_ann$shortname %in% clu_ann_cellname$Var1,
                    ]$geneSymbol)
                  PMID <- unique(clu_ann[clu_ann$shortname %in% clu_ann_cellname$Var1,
                    ]$PMID)
                  res <- "IV"

                  # matching cell subtype to determine the start
                  clu_ann1 <- clu_ann[clu_ann$shortname == clu_ann_cellname$Var1,
                    ]
                  clu_ann1 <- clu_ann1[, c("geneSymbol", "PMID", "third", "second",
                    "first")]
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
                    clu_ann_subcellname <- as.data.frame(table(clu_ann1$value),
                      stringsAsFactors = F)
                    m <- clu_ann_subcellname$Freq + 1
                    clu_ann_subcellname$Freq <- clu_ann_subcellname$Freq/m
                    clu_ann_subcellname <- clu_ann_subcellname[order(clu_ann_subcellname$Var1),
                      ]

                    clu_ann_subarticle <- unique(clu_ann1[, c("PMID", "value")])

                    clu_ann_subarticle <- as.data.frame(table(clu_ann_subarticle$value),
                      stringsAsFactors = F)
                    m <- clu_ann_subarticle$Freq + 1
                    clu_ann_subarticle$Freq <- clu_ann_subarticle$Freq/m
                    clu_ann_subarticle <- clu_ann_subarticle[order(clu_ann_subarticle$Var1),
                      ]

                    clu_ann_subcellname$Freq <- sqrt(clu_ann_subcellname$Freq *
                      clu_ann_subarticle$Freq)

                    clu_ann_subcellname <- clu_ann_subcellname[clu_ann_subcellname$Freq ==
                      max(clu_ann_subcellname$Freq), ]
                    clu_ann_subcellname <- clu_ann_subcellname[clu_ann_subcellname$Freq >
                      0.5, ]

                    # Number of (cell subtype with max score) == 1 -- determinate cell subtype
                    # first/second/third
                    if (nrow(clu_ann_subcellname) == 1) {
                      clu_ann_for_det <- clu_ann1[clu_ann1$value == clu_ann_subcellname$Var1,
                        ]

                      # max cell subtype label exist in first subtype1
                      if (unique(clu_ann_for_det$variable == "first")) {
                        cellsubtype1 <- unique(clu_ann_for_det$value)
                        clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                          cellsubtype1, ]$geneSymbol)
                        PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                          cellsubtype1, ]$PMID)

                        clu_ann_for_second <- clu_ann1[clu_ann1$variable == "second",
                          ]

                        # exist cell second subtype2
                        if (nrow(clu_ann_for_second) > 0) {

                          # calculting the second cell subtype2 max score of the same cell type -- short
                          # name.
                          clu_second <- clu_ann1[clu_ann1$variable == "second",
                            ]
                          clu_second_name <- as.data.frame(table(clu_second$value),
                            stringsAsFactors = F)
                          m <- clu_second_name$Freq + 1
                          clu_second_name$Freq <- clu_second_name$Freq/m
                          clu_second_name <- clu_second_name[order(clu_second_name$Var1),
                            ]

                          clu_second_article <- unique(clu_second[, c("PMID", "value")])
                          clu_second_article <- as.data.frame(table(clu_second_article$value),
                            stringsAsFactors = F)
                          m <- clu_second_article$Freq + 1
                          clu_second_article$Freq <- clu_second_article$Freq/m
                          clu_second_article <- clu_second_article[order(clu_second_article$Var1),
                            ]

                          clu_second_name$Freq <- sqrt(clu_second_name$Freq * clu_second_article$Freq)

                          # matching cell second subtype2
                          clu_ann_second <- clu_ann[clu_ann$first == cellsubtype1,
                            ]
                          clu_ann_second <- clu_ann_second[!is.na(clu_ann_second$first),
                            ]

                          clu_ann_secondname <- as.data.frame(table(clu_ann_second$second),
                            stringsAsFactors = F)

                          if (nrow(clu_ann_secondname) > 0) {

                            m <- clu_ann_secondname$Freq + 1
                            clu_ann_secondname$Freq <- clu_ann_secondname$Freq/m
                            clu_ann_secondname <- clu_ann_secondname[order(clu_ann_secondname$Var1),
                              ]

                            clu_ann_secondarticle <- unique(clu_ann_second[, c("PMID",
                              "second")])
                            clu_ann_secondarticle <- as.data.frame(table(clu_ann_secondarticle$second),
                              stringsAsFactors = F)
                            m <- clu_ann_secondarticle$Freq + 1
                            clu_ann_secondarticle$Freq <- clu_ann_secondarticle$Freq/m
                            clu_ann_secondarticle <- clu_ann_secondarticle[order(clu_ann_secondarticle$Var1),
                              ]

                            # cell second subtype2 scoring and compare with max second subtype2 score
                            clu_ann_secondname$Freq <- sqrt(clu_ann_secondname$Freq *
                              clu_ann_secondarticle$Freq)

                            clu_ann_secondname <- clu_ann_secondname[clu_ann_secondname$Freq ==
                              max(clu_ann_secondname$Freq), ]
                            clu_ann_secondname <- clu_ann_secondname[(clu_ann_secondname$Freq >
                              0.5) & (clu_ann_secondname$Freq >= max(clu_second_name$Freq)),
                              ]

                            if (nrow(clu_ann_secondname) == 1) {
                              cellsubtype2 <- clu_ann_secondname$Var1
                              clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                                cellsubtype2, ]$geneSymbol)
                              PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2,
                                ]$PMID)

                              clu_ann_for_third <- clu_ann1[clu_ann1$variable ==
                                "third", ]

                              # exist cell third subtype3
                              if (nrow(clu_ann_for_third) > 0) {

                                # calculting the third cell subtype3 max score of the same cell type -- short
                                # name.
                                clu_third <- clu_ann1[clu_ann1$variable == "third",
                                  ]
                                clu_third_name <- as.data.frame(table(clu_third$value),
                                  stringsAsFactors = F)
                                m <- clu_third_name$Freq + 1
                                clu_third_name$Freq <- clu_third_name$Freq/m
                                clu_third_name <- clu_third_name[order(clu_third_name$Var1),
                                  ]

                                clu_third_article <- unique(clu_third[, c("PMID",
                                  "value")])
                                clu_third_article <- as.data.frame(table(clu_third_article$value),
                                  stringsAsFactors = F)
                                m <- clu_third_article$Freq + 1
                                clu_third_article$Freq <- clu_third_article$Freq/m
                                clu_third_article <- clu_third_article[order(clu_third_article$Var1),
                                  ]

                                clu_third_name$Freq <- sqrt(clu_third_name$Freq *
                                  clu_third_article$Freq)

                                # matching cell third subtype3
                                clu_ann_third <- clu_ann_second[clu_ann_second$second ==
                                  cellsubtype2, ]
                                clu_ann_third <- clu_ann_third[!is.na(clu_ann_third$second),
                                  ]

                                clu_ann_thirdname <- as.data.frame(table(clu_ann_third$third),
                                  stringsAsFactors = F)

                                if (nrow(clu_ann_thirdname) > 0) {

                                  m <- clu_ann_thirdname$Freq + 1
                                  clu_ann_thirdname$Freq <- clu_ann_thirdname$Freq/m
                                  clu_ann_thirdname <- clu_ann_thirdname[order(clu_ann_thirdname$Var1),
                                    ]

                                  clu_ann_thirdarticle <- unique(clu_ann_third[,
                                    c("PMID", "third")])
                                  clu_ann_thirdarticle <- as.data.frame(table(clu_ann_thirdarticle$third),
                                    stringsAsFactors = F)
                                  m <- clu_ann_thirdarticle$Freq + 1
                                  clu_ann_thirdarticle$Freq <- clu_ann_thirdarticle$Freq/m
                                  clu_ann_thirdarticle <- clu_ann_thirdarticle[order(clu_ann_thirdarticle$Var1),
                                    ]

                                  # cell third subtype3 scoring and compare with max first subtype score
                                  clu_ann_thirdname$Freq <- sqrt(clu_ann_thirdname$Freq *
                                    clu_ann_thirdarticle$Freq)

                                  clu_ann_thirdname <- clu_ann_thirdname[clu_ann_thirdname$Freq ==
                                    max(clu_ann_thirdname$Freq), ]
                                  clu_ann_thirdname <- clu_ann_thirdname[(clu_ann_thirdname$Freq >
                                    0.5) & (clu_ann_thirdname$Freq >= max(clu_third_name$Freq)),
                                    ]

                                  # Number of (cell third subtype3 with max score) >= 1 -- output -- last
                                  if (nrow(clu_ann_thirdname) >= 1) {
                                    cellsubtype3 <- clu_ann_thirdname$Var1
                                    clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                                      cellsubtype3, ]$geneSymbol)
                                    PMID <- unique(clu_ann1[clu_ann1$value %in%
                                      cellsubtype3, ]$PMID)
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
                        clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                          cellsubtype2, ]$geneSymbol)
                        PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                          cellsubtype2, ]$PMID)

                        clu_ann_for_first <- clu_ann1[clu_ann1$variable == "first",
                          ]

                        # exist cell first subtype1
                        if (nrow(clu_ann_for_first) > 0) {

                          # calculting the cell first subtype1 max score of the same cell type -- short
                          # name.
                          clu_first <- clu_ann1[clu_ann1$variable == "first", ]
                          clu_first_name <- as.data.frame(table(clu_first$value),
                            stringsAsFactors = F)
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
                          clu_ann_first <- clu_ann[clu_ann$second == cellsubtype2,
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

                            clu_ann_firstarticle <- unique(clu_ann_first[, c("PMID",
                              "first")])
                            clu_ann_firstarticle <- as.data.frame(table(clu_ann_firstarticle$first),
                              stringsAsFactors = F)
                            m <- clu_ann_firstarticle$Freq + 1
                            clu_ann_firstarticle$Freq <- clu_ann_firstarticle$Freq/m
                            clu_ann_firstarticle <- clu_ann_firstarticle[order(clu_ann_firstarticle$Var1),
                              ]

                            # cell first subtype1 scoring and compare with max first subtype1 score
                            clu_ann_firstname$Freq <- sqrt(clu_ann_firstname$Freq *
                              clu_ann_firstarticle$Freq)

                            clu_ann_firstname <- clu_ann_firstname[clu_ann_firstname$Freq ==
                              max(clu_ann_firstname$Freq), ]
                            clu_ann_firstname <- clu_ann_firstname[(clu_ann_firstname$Freq >
                              0.5) & (clu_ann_firstname$Freq >= max(clu_ann_firstname$Freq)),
                              ]

                            if (nrow(clu_ann_firstname) == 1) {
                              cellsubtype1 <- clu_ann_firstname$Var1
                              clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                                cellsubtype1, ]$geneSymbol)
                              PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype1,
                                ]$PMID)

                              clu_ann_for_third <- clu_ann1[clu_ann1$variable ==
                                "third", ]

                              # exist cell third subtype3
                              if (nrow(clu_ann_for_third) > 0) {

                                # calculting the third cell subtype3 max score of the same cell type -- short
                                # name.
                                clu_third <- clu_ann1[clu_ann1$variable == "third",
                                  ]
                                clu_third_name <- as.data.frame(table(clu_third$value),
                                  stringsAsFactors = F)
                                m <- clu_third_name$Freq + 1
                                clu_third_name$Freq <- clu_third_name$Freq/m
                                clu_third_name <- clu_third_name[order(clu_third_name$Var1),
                                  ]

                                clu_third_article <- unique(clu_third[, c("PMID",
                                  "value")])
                                clu_third_article <- as.data.frame(table(clu_third_article$value),
                                  stringsAsFactors = F)
                                m <- clu_third_article$Freq + 1
                                clu_third_article$Freq <- clu_third_article$Freq/m
                                clu_third_article <- clu_third_article[order(clu_third_article$Var1),
                                  ]

                                clu_third_name$Freq <- sqrt(clu_third_name$Freq *
                                  clu_third_article$Freq)

                                # matching cell third subtype3
                                clu_ann_third <- clu_ann_first[clu_ann_first$first ==
                                  cellsubtype1, ]
                                clu_ann_third <- clu_ann_third[!is.na(clu_ann_third$first),
                                  ]

                                clu_ann_thirdname <- as.data.frame(table(clu_ann_third$third),
                                  stringsAsFactors = F)

                                if (nrow(clu_ann_thirdname) > 0) {

                                  m <- clu_ann_thirdname$Freq + 1
                                  clu_ann_thirdname$Freq <- clu_ann_thirdname$Freq/m
                                  clu_ann_thirdname <- clu_ann_thirdname[order(clu_ann_thirdname$Var1),
                                    ]

                                  clu_ann_thirdarticle <- unique(clu_ann_third[,
                                    c("PMID", "third")])
                                  clu_ann_thirdarticle <- as.data.frame(table(clu_ann_thirdarticle$third),
                                    stringsAsFactors = F)
                                  m <- clu_ann_thirdarticle$Freq + 1
                                  clu_ann_thirdarticle$Freq <- clu_ann_thirdarticle$Freq/m
                                  clu_ann_thirdarticle <- clu_ann_thirdarticle[order(clu_ann_thirdarticle$Var1),
                                    ]

                                  # cell third subtype3 scoring and compare with max first subtype score
                                  clu_ann_thirdname$Freq <- sqrt(clu_ann_thirdname$Freq *
                                    clu_ann_thirdarticle$Freq)

                                  clu_ann_thirdname <- clu_ann_thirdname[clu_ann_thirdname$Freq ==
                                    max(clu_ann_thirdname$Freq), ]
                                  clu_ann_thirdname <- clu_ann_thirdname[(clu_ann_thirdname$Freq >
                                    0.5) & (clu_ann_thirdname$Freq >= max(clu_third_name$Freq)),
                                    ]

                                  # Number of (cell third subtype3 with max score) >= 1 -- output -- last
                                  if (nrow(clu_ann_thirdname) >= 1) {
                                    cellsubtype3 <- clu_ann_thirdname$Var1
                                    clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                                      cellsubtype3, ]$geneSymbol)
                                    PMID <- unique(clu_ann1[clu_ann1$value %in%
                                      cellsubtype3, ]$PMID)
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
                        clu_marker <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                          cellsubtype3, ]$geneSymbol)
                        PMID <- unique(clu_ann_for_det[clu_ann_for_det$value %in%
                          cellsubtype3, ]$PMID)

                        clu_ann_for_second <- clu_ann1[clu_ann1$variable == "second",
                          ]

                        # exist cell second subtype2
                        if (nrow(clu_ann_for_second) > 0) {

                          # calculting the second cell subtype2 max score of the same cell type -- short
                          # name.
                          clu_second <- clu_ann1[clu_ann1$variable == "second",
                            ]
                          clu_second_name <- as.data.frame(table(clu_second$value),
                            stringsAsFactors = F)
                          m <- clu_second_name$Freq + 1
                          clu_second_name$Freq <- clu_second_name$Freq/m
                          clu_second_name <- clu_second_name[order(clu_second_name$Var1),
                            ]

                          clu_second_article <- unique(clu_second[, c("PMID", "value")])
                          clu_second_article <- as.data.frame(table(clu_second_article$value),
                            stringsAsFactors = F)
                          m <- clu_second_article$Freq + 1
                          clu_second_article$Freq <- clu_second_article$Freq/m
                          clu_second_article <- clu_second_article[order(clu_second_article$Var1),
                            ]

                          clu_second_name$Freq <- sqrt(clu_second_name$Freq * clu_second_article$Freq)

                          # matching cell second subtype2
                          clu_ann_second <- clu_ann[clu_ann$third == cellsubtype3,
                            ]
                          clu_ann_second <- clu_ann_second[!is.na(clu_ann_second$third),
                            ]

                          clu_ann_secondname <- as.data.frame(table(clu_ann_second$second),
                            stringsAsFactors = F)

                          if (nrow(clu_ann_secondname) > 0) {

                            m <- clu_ann_secondname$Freq + 1
                            clu_ann_secondname$Freq <- clu_ann_secondname$Freq/m
                            clu_ann_secondname <- clu_ann_secondname[order(clu_ann_secondname$Var1),
                              ]

                            clu_ann_secondarticle <- unique(clu_ann_second[, c("PMID",
                              "second")])
                            clu_ann_secondarticle <- as.data.frame(table(clu_ann_secondarticle$second),
                              stringsAsFactors = F)
                            m <- clu_ann_secondarticle$Freq + 1
                            clu_ann_secondarticle$Freq <- clu_ann_secondarticle$Freq/m
                            clu_ann_secondarticle <- clu_ann_secondarticle[order(clu_ann_secondarticle$Var1),
                              ]

                            # cell second subtype2 scoring and compare with max second subtype2 score
                            clu_ann_secondname$Freq <- sqrt(clu_ann_secondname$Freq *
                              clu_ann_secondarticle$Freq)

                            clu_ann_secondname <- clu_ann_secondname[clu_ann_secondname$Freq ==
                              max(clu_ann_secondname$Freq), ]
                            clu_ann_secondname <- clu_ann_secondname[(clu_ann_secondname$Freq >
                              0.5) & (clu_ann_secondname$Freq >= max(clu_second_name$Freq)),
                              ]

                            if (nrow(clu_ann_secondname) == 1) {
                              cellsubtype2 <- clu_ann_secondname$Var1
                              clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                                cellsubtype2, ]$geneSymbol)
                              PMID <- unique(clu_ann1[clu_ann1$value %in% cellsubtype2,
                                ]$PMID)

                              clu_ann_for_first <- clu_ann1[clu_ann1$variable ==
                                "first", ]

                              # exist cell first subtype1
                              if (nrow(clu_ann_for_first) > 0) {

                                # calculting the first cell subtype1 max score of the same cell type -- short
                                # name.
                                clu_first <- clu_ann1[clu_ann1$variable == "first",
                                  ]
                                clu_first_name <- as.data.frame(table(clu_first$value),
                                  stringsAsFactors = F)
                                m <- clu_first_name$Freq + 1
                                clu_first_name$Freq <- clu_first_name$Freq/m
                                clu_first_name <- clu_first_name[order(clu_first_name$Var1),
                                  ]

                                clu_first_article <- unique(clu_first[, c("PMID",
                                  "value")])
                                clu_first_article <- as.data.frame(table(clu_first_article$value),
                                  stringsAsFactors = F)
                                m <- clu_first_article$Freq + 1
                                clu_first_article$Freq <- clu_first_article$Freq/m
                                clu_first_article <- clu_first_article[order(clu_first_article$Var1),
                                  ]

                                clu_first_name$Freq <- sqrt(clu_first_name$Freq *
                                  clu_first_article$Freq)

                                # matching cell first subtype1
                                clu_ann_first <- clu_ann_second[clu_ann_second$second ==
                                  cellsubtype2, ]
                                clu_ann_first <- clu_ann_first[!is.na(clu_ann_first$second),
                                  ]

                                clu_ann_firstname <- as.data.frame(table(clu_ann_first$first),
                                  stringsAsFactors = F)

                                if (nrow(clu_ann_firstname) > 0) {

                                  m <- clu_ann_firstname$Freq + 1
                                  clu_ann_firstname$Freq <- clu_ann_firstname$Freq/m
                                  clu_ann_firstname <- clu_ann_firstname[order(clu_ann_firstname$Var1),
                                    ]

                                  clu_ann_firstarticle <- unique(clu_ann_first[,
                                    c("PMID", "first")])
                                  clu_ann_firstarticle <- as.data.frame(table(clu_ann_firstarticle$first),
                                    stringsAsFactors = F)
                                  m <- clu_ann_firstarticle$Freq + 1
                                  clu_ann_firstarticle$Freq <- clu_ann_firstarticle$Freq/m
                                  clu_ann_firstarticle <- clu_ann_firstarticle[order(clu_ann_firstarticle$Var1),
                                    ]

                                  # cell first subtype1 scoring and compare with max first subtype1 score
                                  clu_ann_firstname$Freq <- sqrt(clu_ann_firstname$Freq *
                                    clu_ann_firstarticle$Freq)

                                  clu_ann_firstname <- clu_ann_firstname[clu_ann_firstname$Freq ==
                                    max(clu_ann_firstname$Freq), ]
                                  clu_ann_firstname <- clu_ann_firstname[(clu_ann_firstname$Freq >
                                    0.5) & (clu_ann_firstname$Freq >= max(clu_first_name$Freq)),
                                    ]

                                  # Number of (cell first subtype1 with max score) >= 1 -- output -- last
                                  if (nrow(clu_ann_firstname) >= 1) {
                                    cellsubtype1 <- clu_ann_firstname$Var1
                                    clu_marker <- unique(clu_ann1[clu_ann1$value %in%
                                      cellsubtype1, ]$geneSymbol)
                                    PMID <- unique(clu_ann1[clu_ann1$value %in%
                                      cellsubtype1, ]$PMID)
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

        d1 <- celltype_score[1]
        if (length(celltype_score) > 1) {
            for (j in 2:length(celltype_score)) {
                d1 <- paste(d1, celltype_score[j], sep = ", ")
            }
        }
        celltype_score <- d1

        d1 <- PMID[1]
        if (length(PMID) > 1) {
            for (j in 2:length(PMID)) {
                d1 <- paste(d1, PMID[j], sep = ", ")
            }
        }
        PMID <- d1

        clu_ann <- data.frame(cluster = clu.num[i], cluster_marker = rescluster_marker1,
            cellsubtype3 = cellsubtype3, cellsubtype2 = cellsubtype2, cellsubtype1 = cellsubtype1,
            cell_type = celltype, celltype_score = celltype_score, celltype_related_marker = clu_marker,
            PMID = PMID, stringsAsFactors = F)

        clu_ann_res <- rbind(clu_ann_res, clu_ann)
    }

    for (i in 1:nrow(clu_ann_res)) {
        d1 <- as.character(clu_ann_res[i, c("cellsubtype3", "cellsubtype2", "cellsubtype1",
            "cell_type")])
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
    return(clu_ann_res)
}
