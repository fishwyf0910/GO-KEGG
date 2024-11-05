#conda activate r4.3

#BiocManager::install("org.Dr.eg.db")
library(ggplot2)
library(stringr)
library(org.Dr.eg.db)
library(topGO)
library(enrichplot)
library(clusterProfiler)

#导入数据
df <- read.csv("genelist.txt-2",header=F)

#转换数据类型
df02 <- as.character(df$V1)

#symbol to entrez ID
DEG.gene_symbol = df02
DEG.entrez_id = mapIds(x=org.Dr.eg.db,keys=DEG.gene_symbol,keytype="SYMBOL",column = "ENTREZID")

gene = bitr(DEG.gene_symbol, fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Dr.eg.db")

#GO富集分析
ego <- enrichGO(gene=gene$ENTREZID,
        keyType="ENTREZID",#输入基因类型
        OrgDb=org.Dr.eg.db,#导入背景基因
        ont = "all",#GO的类型：BP,CC,MF
        pAdjustMethod="BH",#FDM矫正
        pvalueCutoff=0.05,#p值过滤值
        readable=TRUE)

#GO注释柱形图
pdf(file="GO_bar.pdf",width=10,height=9)
barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#GO注释气泡图
pdf(file="GO_bubble.pdf",width = 10,height = 8)
dotplot(ego,showCategory = 10,split="ONTOLOGY",orderBy = "GeneRatio") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
#导出GO结果
write.csv(ego,"GO.output.csv")

#KEGG富集
KEGG=enrichKEGG(gene = gene$ENTREZID,
                organism = "dre",#斑马鱼
                keyType = "kegg",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
#保存结果
write.csv(KEGG, file = "KEGG_output.csv", row.names = FALSE)

#barplot_kegg_top15图绘制--柱状图
barplot(KEGG,
        x = "GeneRatio",
        color = "p.adjust",
        size = NULL,
        split = NULL,
        font.size = 8,
        showCategory = 15,
        title = "KEGG_enrichment")
ggsave(file = "barplot_kegg_top15.pdf", width = 8, height = 7)

#dotplot_kegg_top15图绘制---气泡图
dotplot(
  KEGG,
  x = "GeneRatio",
  color = "p.adjust",
  showCategory = 15,
  size = NULL,
  split = NULL,
  font.size = 8,
  title = "",
  orderBy = "x",
  label_format = 30)
ggsave(file = "dotplot_kegg_top15.pdf", width = 8, height = 7)
