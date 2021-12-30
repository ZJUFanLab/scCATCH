# scCATCH v3.0
[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/) <a href='#devtools'>![installed with devtools](https://img.shields.io/badge/installed%20with-devtools-blue.svg)</a> 

### Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data

<img src='https://github.com/ZJUFanLab/scCATCH_performance_comparison/blob/master/Overview.png'>

Recent advance in single-cell RNA sequencing (scRNA-seq) has enabled large-scale transcriptional characterization of thousands of cells in multiple complex tissues, in which accurate cell type identification becomes the prerequisite and vital step for scRNA-seq studies. Currently, the common practice in cell type annotation is to map the highly expressed marker genes with known cell markers manually based on the identified clusters, which requires the priori knowledge and tends to be subjective on the choice of which marker genes to use. Besides, such manual annotation is usually time-consuming.

To address these problems, we introduce a __single cell Cluster-based Annotation Toolkit for Cellular Heterogeneity (scCATCH)__ from cluster marker genes identification to cluster annotation based on evidence-based score by matching the identified potential marker genes with known cell markers in tissue-specific cell taxonomy reference database (CellMatch).

[![download CellMatch](https://img.shields.io/badge/download-CellMatch-orange.svg)](https://github.com/ZJUFanLab/scCATCH/blob/master/cellmatch.rda)

__CellMatch includes a panel of 353 cell types and related 686 subtypes associated with 184 tissue types, and 2,097 references of human and mouse.__

The scCATCH mainly includes two function `findmarkergene()` and `scCATCH()` to realize the automatic annotation for each identified cluster. Usage and Examples are detailed below.

# Cite
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.isci.2020.100882-brightgreen.svg)](https://www.sciencedirect.com/science/article/pii/S2589004220300663) [![PMID:32062421](https://img.shields.io/badge/PMID-32062421-blue.svg)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7031312/)

Shao et al., scCATCH:Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data, iScience, Volume 23, Issue 3, 27 March 2020. [doi: 10.1016/j.isci.2020.100882](https://www.sciencedirect.com/science/article/pii/S2589004220300663). PMID:[32062421](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7031312/)

# <a name='devtools'>News</a>
### v3.0
- Update Gene symbols in CellMatch according to NCBI Gene symbols (updated in June. 30, 2021, https://www.ncbi.nlm.nih.gov/gene).
- Allow users to use custom `cellmatch`
- Allow users to select different combination of tissues or cancers for annotation.
- Allow users to add more marker genes to `cellmatch` for annotation.
- Allow users to use markers from different species other than human and mouse.
- Allow users to use more methods to identify highly expressed genes.

### v2.1
- Update Gene symbols in CellMatch according to NCBI Gene symbols (updated in June 19, 2020, https://www.ncbi.nlm.nih.gov/gene). Unmatched marker genes FLJ42102, LOC101928100, LOC200772 and BC017158 are removed.
- Fix the `Error in intI(j, n = x@Dim[2], dn[[2]], give.dn = FALSE) : invalid character indexing` in `findmarkergenes()` by adding a check of cluster number. Refer to [issue 14](https://github.com/ZJUFanLab/scCATCH/issues/14)
- Fix the `Error in object[object$cluster == clu.num[i], ] : wrong number of dimensions` in `scCATCH()` by adding a check of type of input. Refer to [issue 13](https://github.com/ZJUFanLab/scCATCH/issues/13)
- Add a progress bar for `findmarkergenes()` and `scCATCH()`.
- __scCATCH for R > 4.0.0 can be downloaded in [Release](https://github.com/ZJUFanLab/scCATCH/releases/tag/v2.1) page__. 

### v2.0
- Add `cluster` and `match_CellMatch` parameters to handle large scRNA-seq datasets.
- Add `cancer` parameters to annotate scRNA-seq data from tissue with cancer.

# Install
```
# install devtools and install scCATCH
install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/scCATCH')
```
or

```
# download the source package of scCATCH-3.0.tar.gz and install it
# ensure the right directory for scCATCH-3.0.tar.gz
install.packages(pkgs = 'scCATCH-3.0.tar.gz',repos = NULL, type = "source")
```

# Usage
Please refer to the tutorial [vignette](https://raw.githack.com/ZJUFanLab/scCATCH/master/vignettes/tutorial.html). Available tissues and cancers see the [wiki page](https://github.com/ZJUFanLab/scCATCH/wiki)

# Issues
[![bug](https://img.shields.io/badge/reported-bug-orange.svg)](https://github.com/ZJUFanLab/scCATCH/issues?q=is%3Aissue+is%3Aclosed) [![error](https://img.shields.io/badge/reported-error-red.svg)](https://github.com/ZJUFanLab/scCATCH/issues?q=is%3Aissue+is%3Aclosed)

Solutions for possilble bugs and errors. Please refer to closed [Issues1](https://github.com/ZJUFanLab/scCATCH/issues?q=is%3Aissue+is%3Aclosed) and [Issues2](https://github.com/ZJUFanLab/scCATCH_performance_comparison/issues?q=is%3Aissue+is%3Aclosed)

# About
scCATCH was developed by Xin Shao. Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn
