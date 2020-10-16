library(ggplot2)
library(edgeR)
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(topGO)
font_import()

theme_science <- function (base_size = 12, base_family = "Arial Black") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(DESeq2)


WT_KO_RNA <- read.table("WT_KO_RNA.txt")

row.names(WT_KO_RNA) <- WT_KO_RNA$gene

WT_KO_RNA$gene <- NULL

WT_KO_RNA_matrix <- as.matrix(WT_KO_RNA)

head(WT_KO_RNA_matrix)

coldata <- DataFrame(condition=factor(c("1","1","2","2")), 
                     row.names=as.character(colnames(WT_KO_RNA_matrix)))
coldata

ddsHTSeq <- DESeqDataSetFromMatrix(countData= WT_KO_RNA_matrix, 
                                   colData = coldata, 
                                   design = ~ condition )

ddsHTSeq

dds <- DESeq(ddsHTSeq)

res <- results(dds)

res

RNA_deseq_summary_Wk <- as.data.frame(summary(res))

sum(res$padj <0.01, na.rm=TRUE)

res_wk <- as.data.frame(res)

res_wk$gene <- row.names(res_wk)

res_wk$Sig <- ifelse(res_wk$padj<0.05, "<0.05", ">0.05")
write.csv(res_wk, "RNA_deseq_WK.csv")

res_wk_up <- subset(res_wk, log2FoldChange > 1 & Sig == "<0.05")
write.csv(res_wk_up, "RNA_sig_up_deseq_WK.csv")

res_wk_down <- subset(res_wk, log2FoldChange < -1 & Sig == "<0.05")
write.csv(res_wk_down, "RNA_sig_down_deseq_WK.csv")

ggplot(res_wk, aes(log2FoldChange, -log10(padj), colour=Sig)) + geom_point(alpha=0.2) +
  theme_science() + scale_colour_manual(values=cbPalette) + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="Green") +
  annotate('rect', xmin = -0.95, xmax = 0.995, ymin = -0.1, ymax = Inf,alpha = .75, fill="white")
ggsave("volcano_wt-ko_RNA_deseq2.tiff", height=4, width=4, unit="in")

Norm_counts <- as.data.frame(counts(dds, normalized=TRUE))

#plot replicate to replicate comparison
labels = data.frame(x = 1, y = 5, label = format(round(cor(Norm_counts$WTshSCR1, Norm_counts$WTshSCR2), 5)))
ggplot(Norm_counts, aes(log10(WTshSCR1), log10(WTshSCR2))) + geom_point(alpha = 0.1) + theme_science() +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("replicates_WTshScr_RNA_deseq2.tiff", height=4, width=4, unit="in")

#determine correlation between replicates
cor.test(Norm_counts$WTshSCR1, Norm_counts$WTshSCR2)

#plot replicate to replicate comparison
labels = data.frame(x = 1, y = 5, label = format(round(cor(Norm_counts$KOshSCR1, Norm_counts$KOshSCR2), 5)))
ggplot(Norm_counts, aes(log10(KOshSCR1), log10(KOshSCR2))) + geom_point(alpha = 0.1) + theme_science() +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("replicates_KOshSCR_RNA_deseq2.tiff", height=4, width=4, unit="in")

#determine correlation between replicates
cor.test(Norm_counts$KOshSCR1, Norm_counts$KOshSCR2)


WT_WT_RNA <- read.table("WT_WT_RNA.txt")

row.names(WT_WT_RNA) <- WT_WT_RNA$gene

WT_WT_RNA <- WT_WT_RNA[!(row.names(WT_WT_RNA) == "Cdkn2a"),]

WT_WT_RNA$gene <- NULL

WT_WT_RNA_matrix <- as.matrix(WT_WT_RNA)

head(WT_WT_RNA_matrix)

coldata <- DataFrame(condition=factor(c("1","1","2","2")), 
                     row.names=as.character(colnames(WT_WT_RNA_matrix)))
coldata

ddsHTSeq <- DESeqDataSetFromMatrix(countData= WT_WT_RNA_matrix, 
                                   colData = coldata, 
                                   design = ~ condition )
ddsHTSeq

dds <- DESeq(ddsHTSeq)

res <- results(dds)

res

RNA_deseq_summary_WW <- as.data.frame(summary(res))

sum(res$padj <0.01, na.rm=TRUE)

res_ww <- as.data.frame(res)

res_ww$gene <- row.names(res_ww)

res_ww$Sig <- ifelse(res_ww$padj<0.05, "<0.05", ">0.05")
write.csv(res_ww, "RNA_deseq_WW.csv")

res_ww_up <- subset(res_ww, log2FoldChange > 1 & Sig == "<0.05")
write.csv(res_ww_up, "RNA_sig_up_deseq_WW.csv")

res_ww_down <- subset(res_ww, log2FoldChange < -1 & Sig == "<0.05")
write.csv(res_ww_down, "RNA_sig_down_deseq_WW.csv")

ggplot(res_ww, aes(log2FoldChange, -log10(padj), colour=Sig)) + geom_point(alpha=0.2) +
  theme_science() + scale_colour_manual(values=cbPalette) + xlab("Log2 Fold Change") + ylab("-log10 (FDR)") + geom_vline(xintercept = -1, linetype="dashed")+ geom_vline(xintercept = 1, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="Green") +
  annotate('rect', xmin = -0.95, xmax = 0.995, ymin = -0.1, ymax = Inf,alpha = .75, fill="white")
ggsave("volcano_wt-wt_RNA_deseq2.tiff", height=4, width=4, unit="in")

Norm_counts_ww <- as.data.frame(counts(dds, normalized=TRUE))

#plot replicate to replicate comparison
labels = data.frame(x = 1, y = 5, label = format(round(cor(Norm_counts_ww$WTshSCR1, Norm_counts_ww$WTshSCR2), 5)))
ggplot(Norm_counts_ww, aes(log10(WTshARF1), log10(WTshARF2))) + geom_point(alpha = 0.1) + theme_science() +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("replicates_WTshARF_RNA_deseq2.tiff", height=4, width=4, unit="in")

#determine correlation between replicates
cor.test(Norm_counts_ww$WTshSCR1, Norm_counts_ww$WTshSCR2)


summary_deseq2 <- merge(res_wk, res_ww, by = "gene")

cor.test(summary_deseq2$log2FoldChange.x, summary_deseq2$log2FoldChange.y)

labels = data.frame(x = 2.5, y = -5, label = 0.61)

ggplot(summary_deseq2, aes(log2FoldChange.x, log2FoldChange.y)) + geom_point(alpha=0.1) + theme_science() + geom_abline() +
  xlab("Log2 Fold Change WT-KO") + ylab("Log2 Fold Change WTshSCR-WTshARF") +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("deg_corr_RNA_deseq2.tiff", width=4, height=4, unit="in")

go_wk_up <- read.delim("RNA_go_wk_up.txt")
go_wk_down <- read.delim("RNA_go_wk_down.txt")

go_wk_up$GO.biological.process.complete <- gsub("\\(GO:.......\\)", "", go_wk_up$GO.biological.process.complete)
go_wk_down$GO.biological.process.complete <- gsub("\\(GO:.......\\)", "", go_wk_down$GO.biological.process.complete)

go_wk_up <- subset(go_wk_up, !GO.biological.process.complete == "Unclassified (UNCLASSIFIED)")
go_wk_down <- subset(go_wk_down, !GO.biological.process.complete == "Unclassified (UNCLASSIFIED)")

ggplot(subset(go_wk_up, upload_1..P.value. < 3E-15), aes(reorder(GO.biological.process.complete, -upload_1..P.value.), -log10(upload_1..P.value.))) + 
  geom_bar(stat="identity", fill = "#E69F00") + coord_flip() + xlab("") + ylab("-log10(Adjusted P.Value)") + theme_science()
ggsave("go_wk_up.tiff", height = 5, width = 7, units = "in")

ggplot(subset(go_wk_down, upload_1..P.value. < 6E-5), aes(reorder(GO.biological.process.complete, -upload_1..P.value.), -log10(upload_1..P.value.))) + 
  geom_bar(stat="identity", fill = "#56B4E9") + coord_flip() + xlab("") + ylab("-log10(Adjusted P.Value)") + theme_science()
ggsave("go_wk_down.tiff", height = 5, width = 7, units = "in")




