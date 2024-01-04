library("ggplot2")
library("PieGlyph")
library("ggpmisc")
library("ggpubr")

celltypecsvf <- "/home/disk/linxin/lxr/SRR25213403/ccRCC-2-L.csv" ## 细胞类型解卷积结果，记得替换

bc2posf <- "/home/disk/linxr/231207_clist/result/bc2pos/decoder_bc2pos.txt" ## spot空间坐标文件，06/08/10都用这个就行
## 23的样本是10x的数据，用/home/disk/linxr/231207_clist/result/bc2pos/10x_bc2pos.txt

umif <- "/home/disk/linxr/231207_clist/result/ccRCC_p1_filter_umi_25/spot_uminum.txt" ## 每个spot的umi数量
## 23的样本用/home/disk/linxr/231207_clist/result/cSCC_p4_filter_umi_25/spot_uminum.txt
## 06的样本用/home/disk/linxr/231207_clist/result/ccRCC_p2_filter_umi_25/spot_uminum.txt
## 08的样本用/home/disk/linxr/231207_clist/result/chRCC_p1_filter_umi_25/spot_uminum.txt
## 10的样本用/home/disk/linxr/231207_clist/result/ccRCC_p3_filter_umi_25/spot_uminum.txt

ctdf <- read.table(celltypecsvf, h = T, row.names = 1, sep = ",")
ctdf <- t(ctdf)
cts <- colnames(ctdf)
bc2posdf <- read.table(bc2posf, h = F, row.names = 1)
colnames(bc2posdf) <- c("x", "y")
umidf <- read.table(umif, h = F, row.names = 1)
colnames(umidf) <- c("umi")

merdf <- merge(ctdf, bc2posdf, by = "row.names", all = T)
rownames(merdf) <- merdf$Row.names
merdf <- merdf[, -1]
merdf <- merge(merdf, umidf, by = "row.names", all = T)
merdf[is.na(merdf)] <- 0

## 要选取空间的x/y坐标的上下区间
x1 <- 33
x2 <- 42
y1 <- 13
y2 <- 37

subdf <- subset(merdf, x >= x1 & x <= x2 & y >= y1 & y <= y2)
subdf$x <- factor(subdf$x, levels = c(x1:x2))
p1 <- ggplot(subdf, aes(x = x, y = y, color = umi)) +
     geom_point(size = 5) +
     scale_color_gradient(low = "gray", high = "red") +  
     labs(title = paste0("ccRCC_p1 MT UMI count"), x = "Coord X", y = "Coord Y", color = "UMI count") + ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
     theme_classic() +
     theme(plot.title = element_text(hjust = 0.5))
pdf(file=paste0("../result/ccRCC_p1_x", x1, "_", x2, "_y", y1, "_", y2, "_celltype_spot_umi_heatmap.pdf"), height = 5, width = 3.2) ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
print(p1)
dev.off()

highmt <- ifelse(subdf$umi >= 100, "High", "Low") ## 
subdf$highmt <- highmt
p1 <- ggplot(subdf, aes(x = x, y = y, color = highmt)) +
     geom_point(size = 5) +  
     labs(title = paste0("ccRCC_p1 high MT UMI"), x = "Coord X", y = "Coord Y", color = "MT UMI") + ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
     theme_classic() +
     theme(plot.title = element_text(hjust = 0.5))
pdf(file=paste0("../result/ccRCC_p1_x", x1, "_", x2, "_y", y1, "_", y2, "_celltype_spot_umi_highlow.pdf"), height = 5, width = 3.2) ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
print(p1)
dev.off()

p1 <- ggplot(data = subdf, aes(x = x, y = y))+
  geom_pie_glyph(slices = cts)+
  labs(title = paste0("ccRCC_p1 MT UMI celltype"), x = "Coord X", y = "Coord Y") + ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
  theme_classic()
pdf(paste0("../result/ccRCC_p1_x", x1, "_", x2, "_y", y1, "_", y2, "_MT_celltype_pie.pdf"), height = 5.6, width = 5.3) ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
print(p1)
dev.off()

ssubdf <- subset(subdf, highmt == "High")
p1 <- ggplot(data = ssubdf, aes(x = x, y = y))+
  geom_pie_glyph(slices = cts)+
  labs(title = paste0("ccRCC_p1 high MT UMI celltype"), x = "Coord X", y = "Coord Y") + ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
  theme_classic()
pdf(paste0("../result/ccRCC_p1_x", x1, "_", x2, "_y", y1, "_", y2, "_high_MT_celltype_pie.pdf"), height = 4.8, width = 4.8) ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
print(p1)
dev.off()

ssubdf <- subset(subdf, highmt == "Low")
p1 <- ggplot(data = ssubdf, aes(x = x, y = y))+
  geom_pie_glyph(slices = cts)+
  labs(title = paste0("ccRCC_p1 low MT UMI celltype"), x = "Coord X", y = "Coord Y") + ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
  theme_classic()
pdf(paste0("../result/ccRCC_p1_x", x1, "_", x2, "_y", y1, "_", y2, "_low_MT_celltype_pie.pdf"), height = 5.6, width = 5.3) ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
print(p1)
dev.off()

for(eachct in cts){
    ssubdf <- subdf[, c("x", "y", eachct, "umi", "highmt")]
    colnames(ssubdf) <- c("x", "y", "ct", "umi", "highmt")
    ssubdf$x <- factor(ssubdf$x, levels = c(x1:x2))
   p1 <- ggplot(ssubdf, aes(x = x, y = y, color = ct)) +
     geom_point(size = 5) +
     scale_color_gradient(low = "gray", high = "red") +  
     labs(title = paste0("ccRCC_p1 ", eachct), x = "Coord X", y = "Coord Y", color = "Proportion") + ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
     theme_classic() +
     theme(plot.title = element_text(hjust = 0.5))
pdf(file=paste0("../result/ccRCC_p1_x", x1, "_", x2, "_y", y1, "_", y2, "_", eachct,"_spot_umi_heatmap.pdf"), height = 5, width = 3.2) ## ## 这里的ccRCC_p1替换成对应的ccRCC_p2/ccRCC_p3/chRCC_p1/cSCC_p4
print(p1)
dev.off() 

 p2 <- ggplot(ssubdf, aes(x = umi, y = ct)) +
          geom_point() +
          geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
          #stat_poly_eq(formula = y ~ x,
          #aes(label = paste(..eq.label..,
          #                  ..rr.label..,
          #                  sep = "~~~")),
          #parse = TRUE) +
          stat_correlation(method = "pearson",
                           aes(label = paste(after_stat(cor.label),
                                             after_stat(p.value.label),
                                             after_stat(n.label),
                                             sep = '*";"*')),
                           parse = TRUE) +
          labs(title = paste0("ccRCC_p1 ", eachct, " vs MT UMI"), x = "MT UMI count", y = "Proportion of cell type", color = "MT UMI") +
          theme_classic() +
          theme(plot.title = element_text(hjust = 0.5))
pdf(file=paste0("../result/ccRCC_p1_x", x1, "_", x2, "_y", y1, "_", y2, "_", eachct,"_spot_umi_vs_", eachct, "_corr_onegroup.pdf"), height = 5, width = 5)
print(p2)
dev.off() 

}

library(reshape2) 
subdf1 = subdf[,1:32] ## 取决于有多少种细胞类型
subdf1$highmt = subdf$highmt
df_long = melt(subdf1)
head(df_long)


p1 <- ggplot(df_long, aes(x = variable, y = value, fill = highmt)) +
  geom_violin(scale = "width") +
  stat_compare_means(method = "t.test", label = "p.signif") +
  xlab("Cell type") +
  ylab("Proportion") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(angle = 45,vjust=1,hjust = 1))
pdf(file=paste0("../result/ccRCC_p1_x", x1, "_", x2, "_y", y1, "_", y2, "_", eachct,"_spot_umi_celltype_vlnplot.pdf"),height = 6,width = 14)
print(p1)
dev.off()