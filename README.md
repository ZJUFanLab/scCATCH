# Updated scCATCH 2.0
[![R >3.6](https://img.shields.io/badge/R-%3E3.6-brightgreen.svg)](https://www.r-project.org/) <a href='#devtools'>![installed with devtools](https://img.shields.io/badge/installed%20with-devtools-blue.svg)</a> [![source package scCATCH__2.0.tar.gz](https://img.shields.io/badge/source%20package-scCATCH__2.0.tar.gz-yellowgreen.svg)](https://github.com/ZJUFanLab/scCATCH_performance_comparison/blob/master/scCATCH_2.0.tar.gz)

### Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data

<img src='https://github.com/ZJUFanLab/scCATCH_performance_comparison/blob/master/Overview.png'>

Recent advance in single-cell RNA sequencing (scRNA-seq) has enabled large-scale transcriptional characterization of thousands of cells in multiple complex tissues, in which accurate cell type identification becomes the prerequisite and vital step for scRNA-seq studies. Currently, the common practice in cell type annotation is to map the highly expressed marker genes with known cell markers manually based on the identified clusters, which requires the priori knowledge and tends to be subjective on the choice of which marker genes to use. Besides, such manual annotation is usually time-consuming.

To address these problems, we introduce a __single cell Cluster-based Annotation Toolkit for Cellular Heterogeneity (scCATCH)__ from cluster marker genes identification to cluster annotation based on evidence-based score by matching the identified potential marker genes with known cell markers in tissue-specific cell taxonomy reference database (CellMatch).

[![download CellMatch](https://img.shields.io/badge/download-CellMatch-orange.svg)](https://github.com/ZJUFanLab/scCATCH_performance_comparison/blob/master/cellmatch.rds)

__CellMatch includes a panel of 353 cell types and related 686 subtypes associated with 184 tissue types, 20,792 cell-specific marker genes and 2,097 references of human and mouse.__

The scCATCH mainly includes two function `findmarkergenes` and `scCATCH` to realize the automatic annotation for each identified cluster. Usage and Examples are detailed below.

# Cite
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.isci.2020.100882-brightgreen.svg)](https://www.sciencedirect.com/science/article/pii/S2589004220300663) [![PMID:32062421](https://img.shields.io/badge/PMID-32062421-blue.svg)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7031312/)

Shao et al., scCATCH:Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data, Volume 23, Issue 3, 27 March 2020. PMID:[32062421](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7031312/)

# <a name='devtools'>News</a>
<font size = 3>1. scCATCH can handle large single-cell transcriptomic dataset containing more than __10,000 cells__ and more than __15 clusters.__</font>

<font size = 3>2. scCATCH can also be used to annotate scRNA-seq data from __tissue with cancer.__</font>

# Install
`devtools::install_github('ZJUFanLab/scCATCH')`

# Usage

`library(scCATCH)`

<font size=5>1. Cluster marker genes identification.</font>

```(r)
clu_markers <- findmarkergenes(object,
                               species = NULL,
                               cluster = 'All',
                               match_CellMatch = FALSE,
                               cancer = NULL,
                               tissue = NULL,
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)
```

Identify potential marker genes for each cluster from a Seurat object (>= 3.0.0) after log10 normalization and cluster analysis. The potential marker genes in each cluster are identified according to its expression level compared to it in every other clusters. Only significantly highly expressed one in all pair-wise comparison of the cluster will be selected as a potential marker gene for the cluster. Genes will be revised according to NCBI Gene symbols (updated in Jan. 10, 2020, [https://www.ncbi.nlm.nih.gov/gene](https://www.ncbi.nlm.nih.gov/gene)) and no matched genes and duplicated genes will be removed.

`object` 
Seurat object (>= 3.0.0) after log10 normalization and cluster analysis. <font color=red>Please ensure data is log10 normalized data and data has been clustered before running scCATCH pipeline.</font>

`species`
The specie of cells. The species must be defined. 'Human' or 'Mouse'.

`cluster`
Select which clusters for potential marker genes identification. e.g. '1', '2', etc. The default is 'All' to find potential makrer genes for each cluster.

`match_CellMatch`
For large datasets containg > 10,000 cells or > 15 clusters, it is strongly recommended to set match_CellMatch 'TRUE' to match CellMatch database first to include potential marker genes in terms of large system memory it may take.

`cancer`
If match_CellMatch is set TRUE and the sample is from cancer tissue, then the cancer type may be defined. Select one or more related cancer types in 3.2 of Details for human and 3.2 of Details for mouse. The dafult is NULL for tissues without cancer.

`tissue`
If match_CellMatch is set TRUE, then the tissue origin of cells must be defined. Select one or more related tissue types in Details. For tissues without cancer, please refer to 3.1.1 of Details for human tissue types and 3.2.1 of Details for mouse tissue types. For tissues with cancer, please refer to 3.1.2 of Details for human tissue types and 3.2.2 of Details for mouse tissue types.

`cell_min_pct`
Include the gene detected in at least this many cells in each cluster. Default is 0.25.

`logfc`
Include the gene with at least this fold change of average gene expression compared to every other clusters. Default is 0.25.

`pvalue`
Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters. Default is 0.05.

<font size=4>Output</font>

`clu_markers` 
A list include a new data matrix wherein genes are revised by official gene symbols according to NCBI Gene symbols (updated in Jan. 10, 2020, [https://www.ncbi.nlm.nih.gov/gene](https://www.ncbi.nlm.nih.gov/gene)) and no matched genes and duplicated genes are removed as well as a data.frame containing potential marker genes of each selected cluster and the corresponding expressed cells percentage and average fold change for each cluster.

<font size=5>2. Cluster annotation</font>

```(r)
clu_ann <- scCATCH(object,
                   species = NULL,
                   cancer = NULL,
                   tissue = NULL)
```

Evidence-based score and annotation for each cluster by matching the potential marker genes generated from `findmarkergenes` with known cell marker genes in tissue-specific cell taxonomy reference database (CellMatch).

`object` 
The data.frame containing marker genes and the corresponding expressed cells percentage and average fold change for each cluster from the output of `findmarkergenes`.

`species`
The species of cells. Select 'Human' or 'Mouse'.

`cancer`
If the sample is from cancer tissue and you want to match cell marker genes of cancer tissues in CellMatch, then the cancer type may be defined. Select one or more related cancer types in 3.1.2 of Details for human and 3.2.2 of Details for mouse. The dafult is NULL for tissues without cancer.

`tissue`
The tissue origin of cells. Select one or more related tissue types in Details. For tissues without cancer, please refer to 3.1.1 of Details for human tissue types and 3.2.1 of Details for mouse tissue types. For tissues with cancer, please refer to 3.1.2 of Details for human tissue types and 3.2.2 of Details for mouse tissue types.

<font size=4>Output</font>

`clu_ann`
A data.frame containing matched cell type for each cluster, related marker genes, evidence-based score and PMID.

<font size=5>3. Details</font>

<font size=3 color=green>3.1.1 For __Human__ tissue, tissue types are listed as follows:</font>

__Adipose tissue-related__: Abdominal adipose tissue; Adipose tissue; Brown adipose tissue; Fat pad; Subcutaneous adipose tissue; Visceral adipose tissue; White adipose tissue.

__Bladder-related__: Bladder; Urine.

__Blood-related__: Blood; Peripheral blood; Plasma; Serum; Umbilical cord blood; Venous blood.

__Bone-related__: Anterior cruciate ligament; Bone; Bone marrow; Cartilage; Intervertebral disc; Meniscus; Nucleus pulposus; Osteoarthritic cartilage; Periosteum; Skeletal muscle; Spinal cord; Synovial fluid; Synovium.

__Brain-related__: Brain; Dorsolateral prefrontal cortex; Embryonic brain; Embryonic prefrontal cortex; Fetal brain; Hippocampus; Inferior colliculus; Midbrain; Sympathetic ganglion.

__Breast-related__: Breast; Mammary epithelium.

__Embryo-related__: Embryo; Embryoid body; Embryonic brain; Embryonic prefrontal cortex; Embryonic stem cell; Germ; Primitive streak.

__Esophagus-related__: Esophagus.

__Eye-related__: Cornea; Corneal endothelium; Corneal epithelium; Eye; Lacrimal gland; Limbal epithelium; Optic nerve; Retina; Retinal pigment epithelium; Sclerocorneal tissue.

__Fetus-related__: Amniotic fluid; Amniotic membrane; Fetal brain; Fetal gonad; Fetal kidney; Fetal liver; Fetal pancreas; Placenta; Umbilical cord; Umbilical cord blood; Umbilical vein.

__Gonad-related__: Corpus luteum; Fetal gonad; Foreskin; Gonad; Ovarian cortex; Ovarian follicle; Ovary; Seminal plasma; Testis.

__Hair-related__: Chorionic villus; Hair follicle; Scalp.

__Heart-related__: Heart; Myocardium.

__Intestine-related__: Colon; Colorectum; Gastrointestinal tract; Gut; Intestinal crypt; Intestine; Jejunum; Large intestine; Small intestinal crypt; Small intestine.

__Kidney-related__: Adrenal gland; Fetal kidney; Kidney; Renal glomerulus.

__Liver-related__: Fetal liver; Liver.

__Lung-related__:Airway epithelium; Alveolus; Bronchoalveolar system; Lung.

__Lymph-related__: Lymph; Lymph node; Lymphoid tissue.

__Muscle-related__: Muscle; Skeletal muscle.

__Nose-related__: Nasal concha; Nasal epithelium; Sinonasal mucosa.

__Oral cavity-related__: Laryngeal squamous epithelium; Oral mucosa; Salivary gland; Sputum; Submandibular gland; Thyroid; Tonsil; Vocal fold. 

__Ovary-related__: Corpus luteum; Ovarian cortex; Ovarian follicle; Ovary; Oviduct.

__Pancreas-related__: Fetal pancreas; Pancreas; Pancreatic acinar tissue; Pancreatic islet.

__Prostate-related__: Prostate.

__Skin-related__: Dermis; Skin.

__Spleen-related__: Spleen; Splenic red pulp.

__Stomach-related__: Gastric corpus; Gastric epithelium; Gastric gland; Gastrointestinal tract; Pyloric gland; Stomach.

__Testis-related__: Testis.

__Tooth-related__: Deciduous tooth; Dental pulp; Gingiva; Molar; Periodontal ligament; Premolar; Tooth.

__Uterus-related__: Endometrium; Endometrium stroma; Myometrium; Uterus; Vagina.

__Vessel-related__: Adventitia; Antecubital vein; Artery; Blood vessel; Umbilical vein.

__Others__: Ascites; Epithelium; Ligament; Pluripotent stem cell; Thymus; Whartons jelly.

<font size=3 color=green>3.1.2 For __Human__ tissue about cancer, cancer types and the corresponding tissue types are listed as follows:</font>

Acute Myelogenous Leukemia: Blood.

Acute Myeloid Leukemia: Bone marrow.

Adenoid Cystic Carcinoma: Salivary gland.

Alveolar Cell Carcinoma: Serum.

Astrocytoma: Brain.

B-Cell Lymphoma: Lymphoid tissue.

Bladder Cancer: Bladder.

Brain Cancer: Blood vessel; Brain.

Breast Cancer: Blood: Breast; Mammary gland.

Cholangiocarcinoma: Liver; Serum.

Chronic Lymphocytic Leukemia: Blood.

Chronic Myeloid Leukemia: Blood.

CNS Primitive Neuroectodermal Tumor: Brain.

Colon Cancer: Blood; Colon; Serum.

Colorectal Cancer: Blood; Colon; Colorectum; Gastrointestinal tract; Intestine; Liver; Lung; Venous blood.

Cutaneous Squamous Cell Carcinoma: Skin.

Endometrial Cancer: Endometrium.

Ependymoma: Brain.

Esophageal Adenocarcinoma: Esophagus.

Fibroid: Myometrium.

Follicular Lymphoma: Lymph node.

Gallbladder Cancer: Gall bladder; Gastrointestinal tract.

Gastric Cancer: Blood; Peripheral blood; Serum; Stomach.

Glioblastoma: Blood; Brain.

Glioma: Blood vessel; Brain.

Gonadoblastoma: Embryo.

Head and Neck Cancer: Blood; Brain; Oral cavity.

Hepatoblastoma: Liver.

Hepatocellular Cancer: Blood; Bone marrow; Embryo; Liver.

High-grade glioma: Brain.

Infantile Hemangiomas: Placenta.

Intestinal Cancer: Gastrointestinal tract.

Intracranial Aneurysm: Brain.

Kaposi's Sarcoma: Lymph node.

Larynx Cancer: Larynx.

Leukemia: Bone marrow; Peripheral blood.

Lipoma: Adipose tissue.

Liver Cancer: Blood; Liver.

Lung Adenocarcinoma: Lung.

Lung Cancer: Blood; Lung.

Lung Squamous Cell Carcinoma: Lung.

Lymphoma: Blood; Brain; Kidney; Liver; Lymph; Lymph node.

Malignant Insulinoma: Pancreas.

Malignant Mesothelioma: Lung; Pleura.

Malignant Peripheral Nerve Sheath Tumor: Brain.

Medulloblastoma: Brain.

Melanoma: Blood; Peripheral blood; Skin.

Mucoepidermoid Carcinoma: Salivary gland.

Multiple Myeloma: Bone marrow; Peripheral blood.

Myeloma: Bone marrow.

Natural Killer Cell Lymphoma: Lymph node.

Nephroblastoma: Kidney.

Non-Small Cell Lung Cancer: Blood; Lung; Peripheral blood.

Oesophageal Cancer: Blood.

Oligodendroglioma: Brain.

Oral Cancer: Oral cavity.

Oral Squamous Cell Carcinoma: Oral cavity; Peripheral blood.

Osteosarcoma: Bone.

Ovarian Cancer: Ascites; Ovarian cortex; Ovary; Peripheral blood.

Pancreatic Cancer: Blood vessel; Pancreas.

Pancreatic Ductal Adenocarcinomas: Pancreas.

Papillary Thyroid Carcinoma: Thyroid.

Prostate Cancer: Blood; Peripheral blood; Prostate.

Renal Cell Carcinoma: Kidney; Serum.

Renal Clear Cell Carcinoma: Lymph node.

Retinoblastoma: Eye.

Salivary Gland Tumor: Parotid gland; Salivary gland.

Sarcoma: Muscle.

Small Cell Lung Cancer: Lung.

Testicular Germ Cell Tumor: Peripheral blood; Testis.

Thyroid Cancer: Thyroid.

Tongue Cancer: Tongue.

Uterine Leiomyoma: Uterus.

Vascular Tumour: Lymph node.

<font size=3 color=green>3.2.1 For __Mouse__ tissue, tissue types are listed as follows:</font>

__Adipose tissue-related__: Adipose tissue; White adipose tissue.

__Bladder-related__: Bladder.

__Blood-related__: Blood; Peripheral blood; Serum; Umbilical cord blood.

__Bone-related__: Bone; Bone marrow; Meniscus; Skeletal muscle; Spinal cord.

__Brain-related__: Brain; Cerebellum; Fetal brain; Hippocampus; Neural tube.

__Breast-related__: Mammary epithelium; Mammary gland.

__Calvaria-related__: Neonatal calvaria.

__Ear-related__: Cochlea; Inner Ear.

__Embryo-related__: Embryo; Embryoid body; Embryonic heart; Embryonic stem cell.

__Esophagus-related__: Esophagus.

__Eye-related__: Corneal epithelium; Eye; Ganglion cell layer of retina; Inner nuclear layer of retina; Lacrimal gland; Retina.

__Fetus-related__: Fetal brain; Fetal intestine; Fetal liver; Fetal lung; Fetal stomach; Placenta; Umbilical cord; Umbilical cord blood.

__Gonad-related__: Gonad; Ovary; Testis; Yolk sac.

__Hair-related__: Hair follicle.

__Heart-related__: Embryonic heart; Heart; Heart muscle; Neonatal heart.

__Intestine-related__: Colon; Colon epithelium; Fetal intestine; Gastrointestinal tract; Ileum; Intestinal crypt; Intestine; Mesenteric lymph node; Small intestine.

__Kidney-related__: Kidney; Mesonephros.

__Liver-related__: Fetal liver; Liver.

__Lung-related__: Bronchiole; Fetal lung; Lung; Trachea.

__Lymph-related__: Lymph node; Lymphoid tissue; Mesenteric lymph node; Peyer patch.

__Muscle-related__: Heart muscle; Muscle; Neonatal muscle; Skeletal muscle.

__Neonate-related__: Neonatal calvaria; Neonatal heart; Neonatal muscle; Neonatal pancreas; Neonatal rib; Neonatal skin.

__Oral cavity-related__: Submandibular gland; Taste bud. 

__Ovary-related__: Ovary; Yolk sac.

__Pancreas-related__: Neonatal pancreas; Pancreas; Pancreatic islet.

__Prostate-related__: Prostate.

__Skin-related__: Dermis; Epidermis; Neonatal skin; Skin.

__Spleen-related__: Spleen.

__Stomach-related__: Fetal stomach; Gastrointestinal tract; Stomach.

__Testis-related__: Testis.

__Uterus-related__: Uterus.

__Vessel-related__: Aorta; Artery; Blood vessel; Carotid artery.

__Others__: Basilar membrane; Epithelium; Peritoneal cavity; Thymus.

<font size=3 color=green>3.2.2 For __Mouse__ tissue about cancer, cancer types and the corresponding tissue types are listed as follows:</font>

Breast Cancer: Lymph node; Breast; Lung.

Chronic Myeloid Leukemia: Blood.

Colon Cancer: Colon.

Colorectal Cancer: Lymph node; Colon; Colorectum.

Cutaneous Squamous Cell Carcinoma: Skin.

Hepatocellular Cancer: Blood; Liver.

Liver Cancer: Liver.

Lung Cancer: Lung.

Melanoma: Lung.

Pancreatic Cancer: Blood.

Papillary Thyroid Carcinoma: Thyroid.

Prostate Cancer: Prostate.

Renal Cell Carcinoma: Kidney.

Supratentorial Primitive Neuroectodermal Tumor: Brain.

# Examples
```(r)
# Step 1: prepare a Seurat object containing log10 normalized single-cell transcriptomic data matrix and the information of cell clusters.
# Note: please define the species for revising gene symbols. Human or Mouse. The default is to find potential marker genes for all clusters with the percentage of expressed cells (≥25%), using WRS test (P<0.05) and a log10 fold change of ≥0.25. These parameters are adjustable for users.

clu_markers <- findmarkergenes(object = mouse_kidney_203_Seurat,
                               species = 'Mouse'
                               cluster = 'All',
                               match_CellMatch = FALSE,
                               cancer = NULL,
                               tissue = NULL,
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)
                               
# Note: for large datasets, please set match_CellMatch as TRUE and provided tissue types. For tissue with cancer, users may provided the cancer types and corresponding tissue types. See Details. 
```

```(r)
# Step 2: evidence-based scoring and annotaion for identified potential marker genes of each cluster generated from findmarkergenes function.

clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Mouse',
                   cancer = NULL,
                   tissue = 'Kidney')

```

```(r)
# Users can also use scCATCH by selecting multiple cluster, cancer types, tissue types as follows:
clu_markers <- findmarkergenes(object = mouse_kidney_203_Seurat,
                               species = 'Mouse'
                               cluster = '1',
                               match_CellMatch = TRUE,
                               cancer = NULL,
                               tissue = 'Kidney',
                               cell_min_pct = 0.1,
                               logfc = 0.1,
                               pvalue = 0.01)
                               
clu_markers <- findmarkergenes(object = mouse_kidney_203_Seurat,
                               species = 'Mouse'
                               cluster = c('1','2'),
                               match_CellMatch = TRUE,
                               cancer = NULL,
                               tissue = c('Kidney','Mesonephros'))
Note: please select the right cancer type and the corresponding tissue type (See Details).
```
# Contributors
scCATCH was developed by Xin Shao. Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn
