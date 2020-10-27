# Single-Cell Tumor Immune Atlas project: Patient Stratification

These are the scripts necessary for the creation of the patients cell-proportion dataset, for the hierarchical k-means clustering of the patients and to evaluate the clustering with a random forest model trained on the clustering output.

## Package versions
Required packages for the code in this folder and versions used:

* [Tidyverse 1.3.0](https://cran.r-project.org/web/packages/tidyverse/vignettes/paper.html)
* [Seurat 3.2.2](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419305598%3Fshowall%3Dtrue)
* [MClust 5.4.6](https://doi.org/10.32614/RJ-2016-021)
* [dendextend 1.14](10.1093/bioinformatics/btv428)
* [magrittr 1.5](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html)
* [ggsci 2.9](https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html)
* [pals 1.6](https://kwstat.github.io/pals/)
* [RColorBrewer 1.1](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
* [patchwork 1.1.0](https://github.com/thomasp85/patchwork)
* [clValid 0.5](http://dx.doi.org/10.18637/jss.v025.i04)
* [factoextra 1.0.7](https://cran.r-project.org/package=factoextra)
* [NbClust 3.0](http://dx.doi.org/10.18637/jss.v061.i06)
* [GGally 2.0.0](https://ggobi.github.io/ggally/)
* [Rtsne 0.15](https://github.com/jkrijthe/Rtsne)
* [ComplexHeatmap 2.4.3](http://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
* [circlize 0.4.10](https://doi.org/10.1093/bioinformatics/btu393)
* [kohonen 3.0.10](https://cran.r-project.org/web/packages/kohonen/A)
* [caTools 1.18](https://cran.r-project.org/web/packages/caTools/index.html)
* [rpart 4.1](https://cran.r-project.org/web/packages/rpart/)
* [rpart.plot 3.0.9](https://cran.r-project.org/web/packages/rpart.plot/)
* [randomForest 4.6](https://cran.r-project.org/web/packages/randomForest/)
* [rattle 5.4.0](https://cran.r-project.org/web/packages/rattle/)
* [caret 6.0](http://topepo.github.io/caret/index.html)
* [scales 1.1.1]()
