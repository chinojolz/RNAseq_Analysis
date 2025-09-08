# RNAseq_Analysis
Differential gene expression analysis of RNA-seq data using base R, including volcano plots and functional insights into top regulated genes


### Project Overview
This project explores an RNA-seq dataset comparing diseased cell lines and diseased cell lines treated with compound X. The analysis involves differential expression, visualization with a volcano plot, and functional annotation of top regulated genes.


### Task
- Generate a volcano plot.
- Determine the upregulated genes (Genes with Log2FC > 1 and pvalue < 0.01)
- Determine the upregulated genes (Genes with Log2FC < -1 and pvalue < 0.01)
- What are the functions of the top 5 upregulated genes and top 5 downregulated genes. (Use genecards)


### Datasource
he dataset contains an experiment between a diseased cell line and diseased cell lines treated with compound X. The difference in expression change between the two health status is computed as Fold change to log 2 (Log2FC) and the significance of each is computed in p-value. Access [Dataset](https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt) here.

