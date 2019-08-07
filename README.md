# scCATCH
Recent advance in single-cell RNA sequencing (scRNA-seq) has enabled large-scale transcriptional characterization of thousands of cells in multiple complex tissues, in which accurate cell type identification becomes the prerequisite and vital step for scRNA-seq studies. Currently, the common practice in cell type annotation is to map the highly expressed marker genes with known cell markers manually based on the identified clusters, which requires the priori knowledge and tends to be subjective on the choice of which marker genes to use. Besides, such manual annotation is usually time-consuming. To address these problems, we introduce a single cell Cluster-based auto-Annotation Toolkit for Cellular Heterogeneity (scCATCH) from cluster marker genes identification to cluster annotation based on evidence-based score by matching the marker genes with known cell markers in tissue-specific cell taxonomy reference database. The scCATCH mainly includes two function `findmarkergenes` and `scCATCH` to realize the automated annotation for each identified clusters. 

# Install
`devtools::install_github('ZJUFanLab/scCATCH')`

# Usage

`library(sxtest)`

__1. Cluster marker genes identification.__
`findmarkergenes(object, cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05)`

Identify marker genes for each cluster from a Seurat object (>= 3.0.0) after log10 normaliztion and cluser analysis. The marker gene in each cluster is identified according to its expression level compared to it in every other clusters. Only significantly highly expressed one in all pair-wise comparison of the cluster will be selected as a cluster marker gene.

`object` 
Seurat object (>= 3.0.0) after log10 normalization and cluster analysis.

`cell_min_pct`
Include the gene detected in at least this many cells in each cluster. Default is 0.25.

`logfc`
Include the gene with at least this fold change of average gene expression compared to every other clusters. Default is 0.25.

`pvalue`
Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters. Default is 0.05.

__2. Cluster annotation__
`scCATCH(object, species, tissue)`

Evidence-based score and annotation for each cluster generated from `findmarkergenes` by matching the marker genes with known cell markers in tissue-specific cell taxonomy reference database.

`object` 
The data.frame containing marker genes and the corresponding expressed cells percentage and average fold change for each cluster from the output of `findmarkergenes`.

`species`
The species of cells. Select 'Human' or 'Mouse'.

`tissue`
The tissue origin of cells. Select one or more related tissue types in Details.

__Details__

__For human cells__, `tissue` include Abdominal adipose tissue, Adipose tissue, Adrenal gland, Adventitia, Airway epithelium, Alveolus, Amniotic fluid, Amniotic membrane, Antecubital vein, Anterior cruciate ligament, Artery, Ascites, Bladder, Blood, Blood vessel, Bone, Bone marrow, Brain, Breast, Bronchoalveolar system, Brown adipose tissue, Cartilage, Chorionic villus, Colon, Colorectum, Cornea, Corneal endothelium, Corneal epithelium, Corpus luteum, Deciduous tooth, Dental pulp, Dermis, Dorsolateral prefrontal cortex, Embryo, Embryoid body, Embryonic brain, Embryonic prefrontal cortex, Embryonic stem cell, Endometrium, Endometrium stroma, Epithelium, Esophagus, Eye, Fat pad, Fetal brain, Fetal gonad, Fetal kidney, Fetal liver, Fetal pancreas, Foreskin, Gastric corpus, Gastric epithelium, Gastric gland, Gastrointestinal tract, Germ, Gingiva, Gonad, Gut, Hair follicle, Heart, Hippocampus, Inferior colliculus, Intervertebral disc, Intestinal crypt, Intestine, Jejunum, Kidney, Lacrimal gland, Large intestine, Laryngeal squamous epithelium, Ligament, Limbal epithelium, Liver, Lung, Lymph, Lymph node, Lymphoid tissue, Mammary epithelium, Meniscus, Midbrain, Molar, Muscle, Myocardium, Myometrium, Nasal concha, Nasal epithelium, Nucleus pulposus, Optic nerve, Oral mucosa, Osteoarthritic cartilage, Ovarian cortex, Ovarian follicle, Ovary, Oviduct, Pancreas, Pancreatic acinar tissue, Pancreatic islet, Periodontal ligament, Periosteum, Peripheral blood, Placenta, Plasma, Pluripotent stem cell, Premolar, Primitive streak, Prostate, Pyloric gland, Renal glomerulus, Retina, Retinal pigment epithelium, Salivary gland, Scalp, Sclerocorneal tissue, Seminal plasma, Serum, Sinonasal mucosa, Skeletal muscle, Skin, Small intestinal crypt, Small intestine, Spinal cord, Spleen, Splenic red pulp, Sputum, Stomach, Subcutaneous adipose tissue, Submandibular gland, Sympathetic ganglion, Synovial fluid, Synovium, Testis, Thymus, Thyroid, Tonsil, Tooth, Umbilical cord, Umbilical cord blood, Umbilical vein, Undefined, Urine, Uterus, Vagina, Venous blood, Visceral adipose tissue, Vocal fold, Whartons jelly, White adipose tissue.

__For Mouse cells__, `tissue` include 
Adipose tissue, Aorta, Artery, Basilar membrane, Bladder, Blood, Blood vessel, Bone, Bone marrow, Brain, Bronchiole, Carotid artery, Cerebellum, Cochlea, Colon, Colon epithelium, Corneal epithelium, Dermis, Embryo, Embryoid body, Embryonic heart, Embryonic stem cell, Epidermis, Epithelium, Esophagus, Eye, Fetal liver, Ganglion cell layer of retina, Gastrointestinal tract, Gonad, Hair follicle, Heart, Heart muscle, Hippocampus, Ileum, Inner Ear, Inner nuclear layer of retina, Intestinal crypt, Intestine, Kidney, Lacrimal gland, Liver, Lung, Lymph node, Lymphoid tissue, Mammary epithelium, Mammary gland, Meniscus, Mesenteric lymph node, Mesonephros, Muscle, Neural tube, Ovary, Pancreas, Pancreatic islet, Peripheral blood, Peritoneal cavity, Peyer patch, Prostate, Retina, Serum, Skeletal muscle, Skin, Small intestine, Spinal cord, Spleen, Stomach, Submandibular gland, Taste bud, Testis, Thymus, Trachea, Umbilical cord, Umbilical cord blood, Undefined, White adipose tissue, Yolk sac.

# Examples
```(r)
clu_markers <- findmarkergenes(mouse_kidney_203_Seurat, 0.25, 0.25, 0.05)  
clu_ann <- scCATCH(clu_markers, 'Mouse', 'Kidney')

clu_markers <- findmarkergenes(mouse_kidney_203_Seurat, 0.5, 0.25, 0.01)
clu_ann <- scCATCH(clu_markers, 'Mouse', c('Kidney','Fetal liver'))
```
# Contributors
scCATCH was developed by Xin Shao. Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn