# scCATCH v3.0
[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/) <a href='#devtools'>![installed with devtools](https://img.shields.io/badge/installed%20with-devtools-blue.svg)</a> 

### Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data

<img src='https://github.com/ZJUFanLab/scCATCH_performance_comparison/blob/master/Overview.png'>

Recent advance in single-cell RNA sequencing (scRNA-seq) has enabled large-scale transcriptional characterization of thousands of cells in multiple complex tissues, in which accurate cell type identification becomes the prerequisite and vital step for scRNA-seq studies. Currently, the common practice in cell type annotation is to map the highly expressed marker genes with known cell markers manually based on the identified clusters, which requires the priori knowledge and tends to be subjective on the choice of which marker genes to use. Besides, such manual annotation is usually time-consuming.

To address these problems, we introduce a __single cell Cluster-based Annotation Toolkit for Cellular Heterogeneity (scCATCH)__ from cluster marker genes identification to cluster annotation based on evidence-based score by matching the identified potential marker genes with known cell markers in tissue-specific cell taxonomy reference database (CellMatch).

[![download CellMatch](https://img.shields.io/badge/download-CellMatch-orange.svg)](https://github.com/ZJUFanLab/scCATCH/tree/master/data)

__CellMatch includes a panel of 353 cell types and related 686 subtypes associated with 184 tissue types, and 2,097 references of human and mouse.__

The scCATCH mainly includes two function `findmarkergene()` and `findcelltype()` to realize the automatic annotation for each identified cluster. Usage and Examples are detailed below.

# Cite
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.isci.2020.100882-brightgreen.svg)](https://www.sciencedirect.com/science/article/pii/S2589004220300663) [![PMID:32062421](https://img.shields.io/badge/PMID-32062421-blue.svg)](https://pubmed.ncbi.nlm.nih.gov/32062421/)

Shao et al., scCATCH:Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data, iScience, Volume 23, Issue 3, 27 March 2020. [doi: 10.1016/j.isci.2020.100882](https://www.sciencedirect.com/science/article/pii/S2589004220300663). [PMID:32062421](https://pubmed.ncbi.nlm.nih.gov/32062421/)

# <a name='devtools'>News</a>
### v3.0
- __`scCATCH` is available on [CRAN](https://CRAN.R-project.org/package=scCATCH)__
- Update Gene symbols in CellMatch according to NCBI Gene symbols (updated in Jan. 2, 2022, https://www.ncbi.nlm.nih.gov/gene).
- __Allow users to use custom `cellmatch`__
- __Allow users to select different combination of tissues or cancers for annotation.__
- __Allow users to add more marker genes to `cellmatch` for annotation.__
- __Allow users to use markers from different species other than human and mouse.__
- __Allow users to use more methods to identify highly expressed genes.__

# Install
```
install.packages("scCATCH")
```

# Usage
Please refer to the [document](https://github.com/ZJUFanLab/scCATCH/blob/master/vignettes/scCATCH.pdf) and tutorial [vignette](https://raw.githack.com/ZJUFanLab/scCATCH/master/vignettes/tutorial.html). Available tissues and cancers see the [wiki page](https://github.com/ZJUFanLab/scCATCH/wiki)

# Issues
[![bug](https://img.shields.io/badge/reported-bug-orange.svg)](https://github.com/ZJUFanLab/scCATCH/issues?q=is%3Aissue+is%3Aclosed) [![error](https://img.shields.io/badge/reported-error-red.svg)](https://github.com/ZJUFanLab/scCATCH/issues?q=is%3Aissue+is%3Aclosed)

Solutions for possilble bugs and errors. Please refer to closed [Issues1](https://github.com/ZJUFanLab/scCATCH/issues?q=is%3Aissue+is%3Aclosed) and [Issues2](https://github.com/ZJUFanLab/scCATCH_performance_comparison/issues?q=is%3Aissue+is%3Aclosed)

# About
scCATCH was developed by Xin Shao. Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn
