# RNAseq_Analysis
Differential gene expression analysis of RNA-seq data using R, including volcano plots and functional insights into top regulated genes


### Project Overview
This project explores an RNA-seq dataset comparing diseased cell lines and diseased cell lines treated with compound X. The analysis involves differential expression, visualization with a volcano plot, and functional annotation of top regulated genes.


### Task
- Generate a volcano plot.
- Determine the upregulated genes (Genes with Log2FC > 1 and pvalue < 0.01)
- Determine the downregulated genes (Genes with Log2FC < -1 and pvalue < 0.01)
- What are the functions of the top 5 upregulated genes and top 5 downregulated genes. (Use genecards)


### Datasource
The dataset contains an experiment between a diseased cell line and diseased cell lines treated with compound X. The difference in expression change between the two health status is computed as Fold change to log 2 (Log2FC) and the significance of each is computed in p-value. Access [Dataset](https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt) here.


### Methods

#### Data Import
```r
link_to_rnaseq <- "https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt"
rna_seq <- read.table(file = link_to_rnaseq, header = TRUE)
names(rna_seq)
row(rna_seq)
head(rna_seq)
```

#### Volcano Plot
```r
rna_seq$negLogP <- -log10(rna_seq$pvalue)
plot(rna_seq$log2FoldChange, rna_seq$negLogP,
     main = "Volcano Plot of RNA-seq Data",
     xlab = "log2 Fold Change",
     ylab = "-log10(p-value)",
     pch = 20, col = "black")
abline(v = c(-1, 1), col = "red", lty = 2)
abline(h = -log10(0.01), col = "blue", lty = 2)
```

<img width="704" height="432" alt="Screenshot 2025-09-08 015144" src="https://github.com/user-attachments/assets/881b2870-0155-42c0-a477-9e88e2a2c66c" />


#### Gene Classification
```r
rna_seq$diffexpressed <- 'NO'
rna_seq$diffexpressed[rna_seq$log2FoldChange > 1 & rna_seq$pvalue < 0.01] <- 'UP'
rna_seq$diffexpressed[rna_seq$log2FoldChange < -1 & rna_seq$pvalue < 0.01] <- 'DOWN'
head(rna_seq)
```

<img width="1265" height="240" alt="Screenshot 2025-09-08 015253" src="https://github.com/user-attachments/assets/5c2a5cab-0e34-47da-8d3d-0393a9e90b31" />


#### Volcano Plot showing Upregulated and Downregulated Genes
```r
plot(rna_seq$log2FoldChange, rna_seq$negLogP,
     main = "Volcano Plot with Highlighted Genes",
     xlab = "log2 Fold Change",
     ylab = "-log10(p-value)",
     pch = 20,
     col = ifelse(rna_seq$diffexpressed == "UP", "red",
            ifelse(rna_seq$diffexpressed == "DOWN", "blue", "grey")))
abline(v = c(-1, 1), col = "grey", lty = 2)
abline(h = -log10(0.01), col = "grey", lty = 2)
```

<img width="702" height="437" alt="Screenshot 2025-09-08 015338" src="https://github.com/user-attachments/assets/01d65c30-cf79-454e-b70a-5b53cebec1e1" />

The volcano plot shows the distribution of genes based on their log2 fold change (x-axis) and statistical significance (-log10 p-value, y-axis). Genes on the right side (red dots) represent upregulated genes in the treated diseased cells (compound X vs untreated); Genes on the left side (blue dots) represent downregulated genes after treatment; Grey dots represent genes with no significant differential expression.

 **Interpretation**: Compound X treatment induces both upregulation and downregulation of multiple genes, suggesting it influences disease-related molecular pathways.


#### Top 5 Upregulated Genes
```r
up_reg <- rna_seq %>%
  filter(diffexpressed == "UP") %>%
  arrange(desc(log2FoldChange)) %>%
  head(5) %>%
   select(Gene, log2FoldChange, pvalue) 
print(up_reg)
```

<img width="1265" height="216" alt="Screenshot 2025-09-08 015405" src="https://github.com/user-attachments/assets/d125019d-41f2-4610-9d8e-f60a8c0896d6" />



#### Top 5 Downregulated Genes
```r
down_genes <- rna_seq %>%
  filter(diffexpressed == "DOWN") %>%
  arrange(log2FoldChange) %>%
  head(5) %>%
  select(Gene, log2FoldChange, pvalue)

print(down_genes)
```

<img width="1267" height="221" alt="Screenshot 2025-09-08 015450" src="https://github.com/user-attachments/assets/d3630212-3fbb-4a18-a8a6-aff9f40dbff8" />

### Summary of Analysis
The analysis revealed distinct sets of genes upregulated and downregulated upon Compound X treatment. Upregulated genes (e.g., DTHD1, EMILIN2) suggest enhanced apoptosis and matrix remodeling, while downregulated genes (e.g., TBX5, IFITM1) indicate suppression of immune-related and developmental transcriptional programs. These findings highlight potential molecular mechanisms by which Compound X exerts its therapeutic effects.
### Reference
Genecard: [See here](https://www.genecards.org/) 

Task: [HackBio](https://course.thehackbio.com/classroom/2)
