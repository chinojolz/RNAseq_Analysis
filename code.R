link_to_rnaseq <- "https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt"
rna_seq <- read.table(file = link_to_rnaseq, header = TRUE)
names(rna_seq)
row(rna_seq)
head(rna_seq)
rna_seq$negLogP <- -log10(rna_seq$pvalue)
plot(rna_seq$log2FoldChange, rna_seq$negLogP,
     main = "Volcano Plot of RNA-seq Data",
     xlab = "log2 Fold Change",
     ylab = "-log10(p-value)",
     pch = 20, col = "black")
abline(v = c(-1, 1), col = "red", lty = 2)
abline(h = -log10(0.01), col = "blue", lty = 2)
rna_seq$diffexpressed <- 'NO'
rna_seq$diffexpressed[rna_seq$log2FoldChange > 1 & rna_seq$pvalue < 0.01] <- 'UP'
rna_seq$diffexpressed[rna_seq$log2FoldChange < -1 & rna_seq$pvalue < 0.01] <- 'DOWN'
head(rna_seq)
plot(rna_seq$log2FoldChange, rna_seq$negLogP,
     main = "Volcano Plot with Highlighted Genes",
     xlab = "log2 Fold Change",
     ylab = "-log10(p-value)",
     pch = 20,
     col = ifelse(rna_seq$diffexpressed == "UP", "red",
            ifelse(rna_seq$diffexpressed == "DOWN", "blue", "grey")))
abline(v = c(-1, 1), col = "grey", lty = 2)
abline(h = -log10(0.01), col = "grey", lty = 2)
up_reg <- rna_seq %>%
  filter(diffexpressed == "UP") %>%
  arrange(desc(log2FoldChange)) %>%
  head(5) %>%
   select(Gene, log2FoldChange, pvalue) 
print(up_reg)
down_genes <- rna_seq %>%
  filter(diffexpressed == "DOWN") %>%
  arrange(log2FoldChange) %>%
  head(5) %>%
  select(Gene, log2FoldChange, pvalue)

print(down_genes)
