# Projection of different datasets on the atlas

To demonstrate the power of our atlas, we predicted the cell types on datasets from different cancer types and varying experimental designs. We make use of FindTransferAnchors utility provided by Seurat [(Stuart et al, 2019)](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) for transfering cell types on TICA to the query objects. The following scripts shows the necesary commands to perform this and create the plots on the article.


## Dependencies

* [R 3.6.0](https://cran.r-project.org/)
* [Seurat 3.2.0](https://cran.r-project.org/web/packages/Seurat/index.html)
* [tidyverse 1.3.0](https://cran.r-project.org/web/packages/tidyverse/index.html)
* [magritts 1.5.0](https://cloud.r-project.org/package=magrittr)
* [pals 1.6.0](https://kwstat.github.io/pals/)
* [flextable 0.5.10](https://davidgohel.github.io/flextable/)
* [ComplexHeatmap 2.4.3](https://github.com/jokergoo/ComplexHeatmap)
* [patchwork 1.0.0](https://patchwork.data-imaginist.com/)
* [matchSCore2 0.1.0](https://github.com/elimereu/matchSCore2)
* [RColorBrewer 1.1.2](https://cloud.r-project.org/package=RColorBrewer)


## Data

* TICA Seurat object: download it as specified in the publication (Zenodo link)
* Query samples: two human uveal melanoma cancers, one human ovarian cancer, one human uveal melanoma liver metastasis, one human non-small cell lung cancer brain metastasis (including TCR) and two mice colorectal cancers (one full and one only T cells and TCR) 
