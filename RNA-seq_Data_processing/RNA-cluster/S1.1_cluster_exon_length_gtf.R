# ！/usr/bin/R

# The exon length of each reads was obtained from the reference genome

library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/home/songjia/reference/human/Homo_sapiens.GRCh38.98.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
n=t(as.data.frame(exons_gene_lens))
write.table(n,"cluster_exon_length_gtf.txt", sep = "\t", row.names = TRUE)

