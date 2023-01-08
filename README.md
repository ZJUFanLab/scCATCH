# scCATCH v3.2.1
[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/) <a href='#cran'>![installed with CRAN](https://img.shields.io/badge/installed%20with-CRAN-blue)</a> [![download CellMatch](https://img.shields.io/badge/download-CellMatch-orange.svg)](https://github.com/ZJUFanLab/scCATCH/tree/master/data)

### Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data

<img src='https://github.com/ZJUFanLab/scCATCH_performance_comparison/blob/master/Overview.png'>

Recent advance in single-cell RNA sequencing (scRNA-seq) has enabled large-scale transcriptional characterization of thousands of cells in multiple complex tissues, in which accurate cell type identification becomes the prerequisite and vital step for scRNA-seq studies. Currently, the common practice in cell type annotation is to map the highly expressed marker genes with known cell markers manually based on the identified clusters, which requires the priori knowledge and tends to be subjective on the choice of which marker genes to use. Besides, such manual annotation is usually time-consuming.

To address these problems, we introduce a __single cell Cluster-based Annotation Toolkit for Cellular Heterogeneity (scCATCH)__ from cluster marker genes identification to cluster annotation based on evidence-based score by matching the identified potential marker genes with known cell markers in tissue-specific cell taxonomy reference database (CellMatch).

__CellMatch includes a panel of 353 cell types and related 686 subtypes associated with 184 tissue types, and 2,096 references of human and mouse.__

# Install
```
# install from cran

install.packages("scCATCH")
```
OR
```
# install devtools and install

install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/scCATCH')
```

# Usage
Please refer to the [document](https://raw.githack.com/ZJUFanLab/scCATCH/master/vignettes/scCATCH.pdf) and tutorial [vignette](https://raw.githack.com/ZJUFanLab/scCATCH/master/vignettes/tutorial.html). Available tissues and cancers see the [wiki page](https://github.com/ZJUFanLab/scCATCH/wiki)

# <a name='cran'>News</a>
### v3.2.1
- __Now available on [CRAN](https://CRAN.R-project.org/package=scCATCH)__
- __Allow users to [use custom `cellmatch`](https://github.com/ZJUFanLab/scCATCH/wiki/Use-custom-cellmatch)__
- __Allow users to [select different combination of tissues or cancers for annotation](https://github.com/ZJUFanLab/scCATCH/wiki/Select-different-combination-of-tissues-or-cancers-for-annotation)__
- __Allow users to add more marker genes to `cellmatch` for annotation__
- __Allow users to use markers from different species other than human and mouse__
- __Allow users to use more methods to identify highly expressed genes__
- __Create scCATCH object from Seurat object with the following code:__

`obj <- createscCATCH(data = Seurat_obj[['RNA']]@data, cluster = as.character(Idents(Seurat_obj)))`

# Cite
Please cite us as __Shao et al., scCATCH:Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data, iScience, Volume 23, Issue 3, 27 March 2020. [doi: 10.1016/j.isci.2020.100882](https://www.sciencedirect.com/science/article/pii/S2589004220300663). [PMID:32062421](https://pubmed.ncbi.nlm.nih.gov/32062421/)__
