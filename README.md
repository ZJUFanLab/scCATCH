# scCATCH v2.1
[![R >3.5](https://img.shields.io/badge/R-%3E%3D3.6-brightgreen)](https://www.r-project.org/) <a href='#devtools'>![installed with devtools](https://img.shields.io/badge/installed%20with-devtools-blue.svg)</a> 

### Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data

<img src='https://github.com/ZJUFanLab/scCATCH_performance_comparison/blob/master/Overview.png'>

Recent advance in single-cell RNA sequencing (scRNA-seq) has enabled large-scale transcriptional characterization of thousands of cells in multiple complex tissues, in which accurate cell type identification becomes the prerequisite and vital step for scRNA-seq studies. Currently, the common practice in cell type annotation is to map the highly expressed marker genes with known cell markers manually based on the identified clusters, which requires the priori knowledge and tends to be subjective on the choice of which marker genes to use. Besides, such manual annotation is usually time-consuming.

To address these problems, we introduce a __single cell Cluster-based Annotation Toolkit for Cellular Heterogeneity (scCATCH)__ from cluster marker genes identification to cluster annotation based on evidence-based score by matching the identified potential marker genes with known cell markers in tissue-specific cell taxonomy reference database (CellMatch).

[![download CellMatch](https://img.shields.io/badge/download-CellMatch-orange.svg)](https://github.com/ZJUFanLab/scCATCH_performance_comparison/blob/master/cellmatch.rds)

__CellMatch includes a panel of 353 cell types and related 686 subtypes associated with 184 tissue types, 20,792 cell-specific marker genes and 2,097 references of human and mouse.__

The scCATCH mainly includes two function `findmarkergenes()` and `scCATCH()` to realize the automatic annotation for each identified cluster. Usage and Examples are detailed below.

# Cite
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.isci.2020.100882-brightgreen.svg)](https://www.sciencedirect.com/science/article/pii/S2589004220300663) [![PMID:32062421](https://img.shields.io/badge/PMID-32062421-blue.svg)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7031312/)

Shao et al., scCATCH:Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data, iScience, Volume 23, Issue 3, 27 March 2020. [doi: 10.1016/j.isci.2020.100882](https://www.sciencedirect.com/science/article/pii/S2589004220300663). PMID:[32062421](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7031312/)

# <a name='devtools'>News</a>

### v2.0
- Add `cluster` and `match_CellMatch` parameters to handle large scRNA-seq datasets.
- Add `cancer` parameters to annotate scRNA-seq data from tissue with cancer.

### v2.1
- Update Gene symbols in CellMatch according to NCBI Gene symbols (updated in June 19, 2020, https://www.ncbi.nlm.nih.gov/gene). Unmatched marker genes FLJ42102, LOC101928100, LOC200772 and BC017158 are removed.
- Fix the `Error in intI(j, n = x@Dim[2], dn[[2]], give.dn = FALSE) : invalid character indexing` in `findmarkergenes()` by adding a check of cluster number. Refer to [issue 14](https://github.com/ZJUFanLab/scCATCH/issues/14)
- Fix the `Error in object[object$cluster == clu.num[i], ] : wrong number of dimensions` in `scCATCH()` by adding a check of type of input. Refer to [issue 13](https://github.com/ZJUFanLab/scCATCH/issues/13)
- Add a progress bar for `findmarkergenes()` and `scCATCH()`.

# Install
[![source package scCATCH-2.1.tar.gz](https://img.shields.io/badge/source%20package-scCATCH--2.1.tar.gz-yellowgreen)](https://codeload.github.com/ZJUFanLab/scCATCH/tar.gz/v2.1)

```
# download the source package of scCATCH-2.1.tar.gz and install it
# ensure the right directory for scCATCH-2.1.tar.gz
install.packages(pkgs = 'scCATCH-2.1.tar.gz',repos = NULL, type = "source")
```

or

```
# install devtools and install scCATCH
install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/scCATCH')
```

# Usage

`library(scCATCH)`

### Cluster marker genes identification

```
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

Identify potential marker genes for each cluster from a Seurat object (>= 3.0.0) after the default log1p normalization and cluster analysis. The potential marker genes in each cluster are identified according to its expression level compared to it in every other clusters. Only significantly highly expressed one in all pair-wise comparison of the cluster will be selected as a potential marker gene for the cluster. Genes will be revised according to NCBI Gene symbols (updated in June 19, 2020, [https://www.ncbi.nlm.nih.gov/gene](https://www.ncbi.nlm.nih.gov/gene)) and no matched genes and duplicated genes will be removed.

`object` Seurat object (>= 3.0.0) after the default log1p normalization and cluster analysis. <font color=red>Please ensure data is log1p normalized data and data has been clustered before running scCATCH pipeline.</font>

`species`Species of cells. Species must be defined. `'Human'` or `'Mouse'`.

`cluster`Select which clusters for potential marker genes identification. e.g. '1', '2', etc. Default is `'All'` to find potential marker genes for each cluster.

`match_CellMatch`For large datasets containg > 10,000 cells or > 15 clusters, it is strongly recommended to set match_CellMatch `TRUE` to match CellMatch database first to include potential marker genes in terms of large system memory it may take.

`cancer`If match_CellMatch is set TRUE and the sample is from cancer tissue, then the cancer type may be defined. Select one or more related cancer types detailed in [wiki page](https://github.com/ZJUFanLab/scCATCH/wiki). Default is `NULL`.

`tissue`If match_CellMatch is set TRUE, then the tissue origin of cells must be defined. Select one or more related tissue types detailed in [wiki page](https://github.com/ZJUFanLab/scCATCH/wiki)

`cell_min_pct`Include the gene detected in at least this many cells in each cluster. Default is `0.25`.

`logfc`Include the gene with at least this fold change of average gene expression compared to every other clusters. Default is `0.25`.

`pvalue`Include the significantly highly expressed gene with this cutoff of p value from wilcox test compared to every other clusters. Default is `0.05`.

__Output__

`clu_markers`A list include a new data matrix wherein genes are revised by official gene symbols according to NCBI Gene symbols (updated in June 19, 2020, [https://www.ncbi.nlm.nih.gov/gene](https://www.ncbi.nlm.nih.gov/gene)) and no matched genes and duplicated genes are removed as well as a data.frame containing potential marker genes of each selected cluster and the corresponding expressed cells percentage and average fold change for each cluster.

### Cluster annotation

```
clu_ann <- scCATCH(object,
                   species = NULL,
                   cancer = NULL,
                   tissue = NULL)
```

Evidence-based score and annotation for each cluster by matching the potential marker genes generated from `findmarkergenes()` with known cell marker genes in tissue-specific cell taxonomy reference database (CellMatch).

`object`Data.frame containing marker genes and the corresponding expressed cells percentage and average fold change for each cluster from the output of `findmarkergenes()`.

`species`Species of cells. Select `'Human'` or `'Mouse'`

`cancer`If the sample is from cancer tissue and you want to match cell marker genes of cancer tissues in CellMatch, then the cancer type may be defined. Select one or more related cancer types detailed in [wiki page](https://github.com/ZJUFanLab/scCATCH/wiki)

`tissue`The tissue origin of cells. Select one or more related tissue types in detailed in [wiki page](https://github.com/ZJUFanLab/scCATCH/wiki)

__Output__

`clu_ann`A data.frame containing matched cell type for each cluster, related marker genes, evidence-based score and PMID.

# Examples
```
# Step 1: prepare a Seurat object containing log1p normalized single-cell transcriptomic data matrix and the information of cell clusters.
# Note: please define the species for revising gene symbols. Human or Mouse. The default is to find potential marker genes for all clusters with the percentage of expressed cells (≥25%), using WRS test (P<0.05) and a log1p fold change of ≥0.25. These parameters are adjustable for users.

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

```
# Step 2: evidence-based scoring and annotaion for identified potential marker genes of each cluster generated from findmarkergenes function.

clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Mouse',
                   cancer = NULL,
                   tissue = 'Kidney')

```

```
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
Note: please select the right cancer type and the corresponding tissue type (See wiki page.
```
# Issues
[![bug](https://img.shields.io/badge/reported-bug-orange.svg)](https://github.com/ZJUFanLab/scCATCH/issues?q=is%3Aissue+is%3Aclosed) [![error](https://img.shields.io/badge/reported-error-red.svg)](https://github.com/ZJUFanLab/scCATCH/issues?q=is%3Aissue+is%3Aclosed)

Solutions for possilble bugs and errors. Please refer to closed [Issues1](https://github.com/ZJUFanLab/scCATCH/issues?q=is%3Aissue+is%3Aclosed) and [Issues2](https://github.com/ZJUFanLab/scCATCH_performance_comparison/issues?q=is%3Aissue+is%3Aclosed)

# Contributors
![Xin Shao](https://img.shields.io/badge/Xin-Shao-brightgreen.svg) ![Xiaohui Fan](https://img.shields.io/badge/Xiaohui-Fan-green.svg)

scCATCH was developed by Xin Shao. Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn
