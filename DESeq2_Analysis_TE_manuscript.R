library(ggplot2)
library(edgeR)
library(dplyr)
library(ggplot2)
library(extrafont)
library(gridExtra)
library(RGraphics)
library(ggridges)
font_import()
loadfonts()

loadfonts(dev="win")
windowsFonts()

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

setwd("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/")

#load HTSeq counts for WT-shSCR replicates
WTshSCR1_ribo <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_WTshScr1_ribo.txt")
colnames(WTshSCR1_ribo) <- c("gene", "WTshSCR1")
WTshSCR2_ribo <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_WTshScr2_ribo.txt")
colnames(WTshSCR2_ribo) <- c("gene", "WTshSCR2")

#merge counts from replicates
WTshScr_rep_ribo <- merge(WTshSCR1_ribo, WTshSCR2_ribo, by = "gene")

WTshScr_rep_ribo <- subset(WTshScr_rep_ribo, !gene %in% c("no_feature", "alignment_not_unique", "ambiguous"))


#load HTSeq counts for WT-shARF replicates
WTshARF1_ribo <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_WTshARF1_ribo.txt")
colnames(WTshARF1_ribo) <- c("gene", "WTshARF1")
WTshARF2_ribo <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_WTshARF2_ribo.txt")
colnames(WTshARF2_ribo) <- c("gene", "WTshARF2")

#merge counts from replicates
WTshARF_rep_ribo <- merge(WTshARF1_ribo, WTshARF2_ribo, by = "gene")

WTshARF_rep_ribo <- subset(WTshARF_rep_ribo, !gene %in% c("no_feature", "alignment_not_unique", "ambiguous"))

#load HTSeq counts for KO-shSCR replicates
KOshSCR1_ribo <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_KOshSCR1_ribo.txt")
colnames(KOshSCR1_ribo) <- c("gene", "KOshSCR1")
KOshSCR2_ribo <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_KOshSCR2_ribo.txt")
colnames(KOshSCR2_ribo) <- c("gene", "KOshSCR2")

#merge counts from replicates
KOshScr_rep_ribo <- merge(KOshSCR1_ribo, KOshSCR2_ribo, by = "gene")

KOshScr_rep_ribo <- subset(KOshScr_rep_ribo, !gene %in% c("no_feature", "alignment_not_unique", "ambiguous"))


#merge WTshSCR and KOshSCR counts, write table and read in
wk_ribo <- merge(WTshScr_rep_ribo, KOshScr_rep_ribo, by = "gene")

write.table(wk_ribo, "WT_KO_ribo.txt")

WT_KO_ribo <- read.table("WT_KO_ribo.txt")

#load HTSeq counts for WT-shSCR replicates
WTshSCR1 <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_WTshSCR1_RNA.txt")
colnames(WTshSCR1) <- c("gene", "WTshSCR1")
WTshSCR2 <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_WTshSCR2_RNA.txt")
colnames(WTshSCR2) <- c("gene", "WTshSCR2")

#merge counts from replicates
WTshScr_rep <- merge(WTshSCR1, WTshSCR2, by = "gene")

WTshScr_rep <- subset(WTshScr_rep, !gene %in% c("no_feature", "alignment_not_unique", "ambiguous"))


#load HTSeq counts for WT-shARF replicates
WTshARF1 <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_WTshARF1_RNA.txt")
colnames(WTshARF1) <- c("gene", "WTshARF1")
WTshARF2 <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_WTshARF2_RNA.txt")
colnames(WTshARF2) <- c("gene", "WTshARF2")

#merge counts from replicates
WTshARF_rep <- merge(WTshARF1, WTshARF2, by = "gene")

WTshARF_rep <- subset(WTshARF_rep, !gene %in% c("no_feature", "alignment_not_unique", "ambiguous"))

#load HTSeq counts for KO-shSCR replicates
KOshSCR1 <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_KOshSCR1_RNA.txt")
colnames(KOshSCR1) <- c("gene", "KOshSCR1")
KOshSCR2 <- read.table("/Users/cottr/Box Sync/Results/Ribo_RNA_Seq_HiSeq/htseq_results_KOshSCR2_RNA.txt")
colnames(KOshSCR2) <- c("gene", "KOshSCR2")

#merge counts from replicates
KOshScr_rep <- merge(KOshSCR1, KOshSCR2, by = "gene")

KOshScr_rep <- subset(KOshScr_rep, !gene %in% c("no_feature", "alignment_not_unique", "ambiguous"))


#merge WTshSCR and KOshSCR counts, write table and read in
wk <- merge(WTshScr_rep, KOshScr_rep, by = "gene")

write.table(wk, "WT_KO_RNA.txt")

WT_KO_rna <- read.table("WT_KO_RNA.txt")


#merge WT_KO_rna and WT_KO_ribo

WT_KO_TE <- merge(WT_KO_rna, WT_KO_ribo, by = "gene")

row.names(WT_KO_TE) <- WT_KO_TE$gene

WT_KO_TE$gene <- NULL

WT_KO_TE_matrix <- as.matrix(WT_KO_TE)

head(WT_KO_TE_matrix)

coldata <- DataFrame(assay=factor(c("mRNA", "mRNA", "mRNA", "mRNA", "Ribo", "Ribo", "Ribo", "Ribo")), 
                     condition=factor(c("1","1","2","2","1","1","2","2")), 
                     row.names=as.character(colnames(WT_KO_TE_matrix)))
coldata

dds_wt_ko <- DESeqDataSetFromMatrix(countData= WT_KO_TE_matrix, 
                              colData = coldata, 
                              design = ~ assay + condition + assay:condition)

dds_wt_ko

dds_wt_ko <- DESeq(dds_wt_ko, test="LRT", reduced = ~ assay + condition)


Norm_counts <- as.data.frame(counts(dds_wt_ko, normalized=TRUE))

colnames(Norm_counts) <- c("WTshSCR1_RNA", "WTshSCR2_RNA", "KOshSCR1_RNA", "KOshSCR2_RNA", "WTshSCR1_Ribo", "WTshSCR2_Ribo", "KOshSCR1_Ribo", "KOshSCR2_Ribo")

Norm_counts$WTshSCR_RNA <- (Norm_counts$WTshSCR1_RNA + Norm_counts$WTshSCR2_RNA)/2
Norm_counts$WTshSCR_Ribo <- (Norm_counts$WTshSCR1_Ribo + Norm_counts$WTshSCR2_Ribo)/2
Norm_counts$KOshSCR_RNA <- (Norm_counts$KOshSCR1_RNA + Norm_counts$KOshSCR2_RNA)/2
Norm_counts$KOshSCR_Ribo <- (Norm_counts$KOshSCR1_Ribo + Norm_counts$KOshSCR2_Ribo)/2

#plot replicate to replicate comparison
labels = data.frame(x = 1, y = 5, label = format(round(cor(Norm_counts$WTshSCR1_RNA, Norm_counts$WTshSCR2_RNA), 5)))
ggplot(Norm_counts, aes(log10(WTshSCR1_RNA), log10(WTshSCR2_RNA))) + geom_point(alpha = 0.1) + theme_science() +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("replicates_WTshScr_RNA.tiff", height=4, width=4, unit="in")

#plot replicate to replicate comparison
labels = data.frame(x = 1, y = 5, label = format(round(cor(Norm_counts$WTshSCR1_Ribo, Norm_counts$WTshSCR2_Ribo), 5)))
ggplot(Norm_counts, aes(log10(WTshSCR1_Ribo), log10(WTshSCR2_Ribo))) + geom_point(alpha = 0.1) + theme_science() +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("replicates_WTshScr_Ribo.tiff", height=4, width=4, unit="in")

#plot replicate to replicate comparison
labels = data.frame(x = 1, y = 5, label = format(round(cor(Norm_counts$KOshSCR1_RNA, Norm_counts$KOshSCR2_RNA), 5)))
ggplot(Norm_counts, aes(log10(KOshSCR1_RNA), log10(KOshSCR2_RNA))) + geom_point(alpha = 0.1) + theme_science() +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("replicates_KOshScr_RNA.tiff", height=4, width=4, unit="in")

#plot replicate to replicate comparison
labels = data.frame(x = 1, y = 5, label = format(round(cor(Norm_counts$KOshSCR1_Ribo, Norm_counts$KOshSCR2_Ribo), 5)))
ggplot(Norm_counts, aes(log10(KOshSCR1_Ribo), log10(KOshSCR2_Ribo))) + geom_point(alpha = 0.1) + theme_science() +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("replicates_KOshScr_Ribo.tiff", height=4, width=4, unit="in")

ggplot(Norm_counts, aes(log2(WTshSCR_RNA), log2(WTshSCR_Ribo))) + geom_point()
ggplot(Norm_counts, aes(log2(KOshSCR_RNA), log2(KOshSCR_Ribo))) + geom_point()


res <- results(dds_wt_ko)

res

sum(res$padj <0.01, na.rm=TRUE)

res_wk <- as.data.frame(res)

res_wk$gene <- row.names(res_wk)


res_wk$Sig <- ifelse(res_wk$padj<0.05, "<0.05", ">0.05")

ggplot(subset(res_wk, !is.na(padj)), aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha=0.5) +
  theme_science() + scale_colour_manual("FDR", values=cbPalette) + xlab("Log2 Fold Change") + ylab("-log10 FDR") + 
  geom_vline(xintercept = -0.5, linetype="dashed")+ geom_vline(xintercept = 0.5, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed") + scale_x_continuous(limits = c(-3,3)) + theme(legend.position="bottom")
ggsave("volcano_wt-ko_TE_deseq2.tiff", height=4, width=4, unit="in")

ggplot(res_wk, aes(log2FoldChange)) + geom_histogram() + theme_science()

TOP_genes <-read.delim("known_top.txt")

TOP_genes$TOP <- TRUE

res_wk <- merge(res_wk, TOP_genes, by.x = "gene", by.y = "Gene.name", all = TRUE)

res_wk$TOP[is.na(res_wk$TOP)] <- FALSE

write.csv(res_wk, "res_wk.csv")

res_wk_up <- subset(res_wk, log2FoldChange > 0 & Sig == "<0.05")

write.csv(res_wk_up, "TE_up_WK.csv")

res_wk_up_top <- subset(res_wk_up, TOP == TRUE)

res_wk_top <- subset(res_wk, TOP == TRUE)

res_wk_top_up <- subset(res_wk_top, log2FoldChange > 0 & Sig == "<0.05")

phyper(q =length(row.names(res_wk_top_up)), 
       m = length(row.names(res_wk_up)), 
       n = length(row.names(res_wk))-length(row.names(res_wk_top)), 
       k =length(row.names(res_wk_top)), log.p = TRUE)


ggplot(res_wk, aes(TOP, log2FoldChange, fill = TOP)) + geom_violin() + theme_science() + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = cbPalette, guide = FALSE) + 
  scale_x_discrete(breaks=c(FALSE, TRUE),labels=c("Non-5'TOP", "5'TOP"))+
  xlab("") + ylab("Log2 Fold Change")
ggsave("violin.tiff", height = 4, width = 3, unit = "in")

median(na.omit(res_wk$log2FoldChange[res_wk$TOP == TRUE]))
median(na.omit(res_wk$log2FoldChange[res_wk$TOP == FALSE]))

t.test(res_wk$log2FoldChange[res_wk$TOP == TRUE], res_wk$log2FoldChange[res_wk$TOP == FALSE])

#GO terms from  PANTHER Overrepresentation Tes, GO Ontology database Released 2019-12-09, http://geneontology.org/
go_wk_up <- read.delim("TE_go_wk_up.txt")

go_wk_up$GO.biological.process.complete <- gsub("\\(GO:.......\\)", "", go_wk_up$GO.biological.process.complete)

go_wk_up <- subset(go_wk_up, !GO.biological.process.complete == "Unclassified (UNCLASSIFIED)")

ggplot(go_wk_up, aes(reorder(GO.biological.process.complete, -upload_1..FDR.), -log10(upload_1..FDR.))) + 
  geom_bar(stat="identity",fill = "#E69F00") + coord_flip() + xlab("") + ylab("-log10 FDR") + theme_science()
ggsave("go_wk_up_te.tiff", height = 7, width = 7, units = "in")


torin <- read.delim("NIHMS367466-supplement-3.txt")

res_wk_torin <- merge(res_wk, torin, by.x = "gene", by.y = "Gene.symbol")

ggplot(res_wk_torin, aes(log2FoldChange, log2_FC_TE, colour = TOP)) + geom_point(alpha = 0.2) + theme_science() +
  labs(y = "log2 Fold Change TE - Torin", x = "log2 Fold Change TE - Arf-KO") + 
  geom_smooth(method = "lm") + scale_x_continuous(limits = c(-2.5,2.5)) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
ggsave("cor_torin_wk.tiff")

cor.test(res_wk_torin$log2FoldChange[res_wk_torin$TOP == TRUE], res_wk_torin$log2_FC_TE[res_wk_torin$TOP == TRUE])



ww_ribo <- merge(WTshScr_rep_ribo, WTshARF_rep_ribo, by = "gene")
write.table(ww_ribo, "WT_WT_ribo.txt")

WT_WT_ribo <- read.table("WT_WT_ribo.txt")

ww <- merge(WTshScr_rep, WTshARF_rep, by = "gene")
write.table(ww, "WT_WT_RNA.txt")

WT_WT_rna <- read.table("WT_WT_RNA.txt")

WT_WT_TE <- merge(WT_WT_rna, WT_WT_ribo, by = "gene")

row.names(WT_WT_TE) <- WT_WT_TE$gene

WT_WT_TE$gene <- NULL

WT_WT_TE_matrix <- as.matrix(WT_WT_TE)

head(WT_WT_TE_matrix)

coldata <- DataFrame(assay=factor(c("mRNA", "mRNA", "mRNA", "mRNA", "Ribo", "Ribo", "Ribo", "Ribo")), 
                     condition=factor(c("1","1","2","2","1","1","2","2")), 
                     row.names=as.character(colnames(WT_WT_TE_matrix)))
coldata

dds_wt_wt <- DESeqDataSetFromMatrix(countData= WT_WT_TE_matrix, 
                                   colData = coldata, 
                                   design = ~ assay + condition + assay:condition)

dds_wt_wt

dds_wt_wt <- DESeq(dds_wt_wt, test="LRT", reduced = ~ assay + condition)

res_ww <- results(dds_wt_wt)

res_ww

sum(res_ww$padj <0.01, na.rm=TRUE)

res_ww <- as.data.frame(res_ww)

res_ww$gene <- row.names(res_ww)


res_ww$Sig <- ifelse(res_ww$padj<0.05, "<0.05", ">0.05")

write.csv(res_ww, "res_ww.csv")


Norm_counts_ww <- as.data.frame(counts(dds_wt_wt, normalized=TRUE))

colnames(Norm_counts_ww) <- c("WTshSCR1_RNA", "WTshSCR2_RNA", "WTshARF1_RNA", "WTshARF2_RNA", "WTshSCR1_Ribo", "WTshSCR2_Ribo", "WTshARF1_Ribo", "WTshARF2_Ribo")

Norm_counts_ww$WTshSCR_RNA <- (Norm_counts_ww$WTshSCR1_RNA + Norm_counts_ww$WTshSCR2_RNA)/2
Norm_counts_ww$WTshSCR_Ribo <- (Norm_counts_ww$WTshSCR1_Ribo + Norm_counts_ww$WTshSCR2_Ribo)/2
Norm_counts_ww$WTshARF_RNA <- (Norm_counts_ww$WTshARF1_RNA + Norm_counts_ww$WTshARF2_RNA)/2
Norm_counts_ww$WTshARF_Ribo <- (Norm_counts_ww$WTshARF1_Ribo + Norm_counts_ww$WTshARF2_Ribo)/2

#plot replicate to replicate comparison
labels = data.frame(x = 1, y = 5, label = format(round(cor(Norm_counts_ww$WTshARF1_RNA, Norm_counts_ww$WTshARF2_RNA), 5)))
ggplot(Norm_counts_ww, aes(log10(WTshARF1_RNA), log10(WTshARF2_RNA))) + geom_point(alpha = 0.1) + theme_science() +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("replicates_WTshARF_RNA.tiff", height=4, width=4, unit="in")

#plot replicate to replicate comparison
labels = data.frame(x = 1, y = 5, label = format(round(cor(Norm_counts_ww$WTshARF1_Ribo, Norm_counts_ww$WTshARF2_Ribo), 5)))
ggplot(Norm_counts_ww, aes(log10(WTshARF1_Ribo), log10(WTshARF2_Ribo))) + geom_point(alpha = 0.1) + theme_science() +
  geom_text(data = labels, aes(x = x, y = y, label = label))
ggsave("replicates_WTshARF_Ribo.tiff", height=4, width=4, unit="in")



ggplot(subset(res_ww, !is.na(padj)), aes(log2FoldChange, -log10(padj), colour = Sig)) + geom_point(alpha=0.5) +
  theme_science() + scale_colour_manual("FDR", values=cbPalette) + xlab("Log2 Fold Change") + ylab("-log10 FDR") + 
  geom_vline(xintercept = -0.5, linetype="dashed")+ geom_vline(xintercept = 0.5, linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed") + scale_x_continuous(limits = c(-3,3)) + theme(legend.position="bottom")
ggsave("volcano_wt-wt_TE_deseq2.tiff", height=4, width=4, unit="in")

ggplot(res_ww, aes(log2FoldChange)) + geom_histogram() + theme_science()

wk_ww <- merge(res_wk, res_ww, by = "gene")

cor.test(wk_ww$log2FoldChange.x, wk_ww$log2FoldChange.y)

ggplot(wk_ww, aes(log2FoldChange.x, log2FoldChange.y)) + geom_point(alpha = 0.2) + geom_smooth(method = "lm", colour = "#56B4E9") + theme_science() +
  labs(x = "Log2 Fold Change TE - KO", y = "Log2 Fold Change TE - Knockdown") + annotate("text", x = -5, y = c(6, 5.3), label = c("r = 0.534", "p-value < 2.2E-16"))
ggsave("ww_wk_te.tiff", height = 4.5, width = 4.5, units = "in")


Examples_wk <- subset(res_wk, gene %in% c("Pabpc1", "Rpl22", "Rpl23a", "Tpt1", "Eef2"))


Examples_wk$FC <- 2^Examples_wk$log2FoldChange
Examples_wk$FC_se <- 2^Examples_wk$lfcSE

ggplot(Examples_wk, aes(gene, FC)) + geom_col() + theme_science()



Examples_ww <- subset(res_ww, gene %in% c("Pabpc1", "Rpl22", "Rpl23a", "Tpt1", "Eef2"))


Examples_ww$FC <- 2^Examples_ww$log2FoldChange

Examples_ww$FC_se <- 2^Examples_ww$lfcSE

ggplot(Examples_ww, aes(gene, FC)) + geom_col() + theme_science()

Examples_ww$Type <- "ARF Knockdown"
Examples_wk$Type <- "Arf Knockout"

Examples_wk$TOP <- NULL


Examples <- rbind(Examples_wk, Examples_ww)

Examples$gene <- toupper(Examples$gene)

ggplot(Examples, aes(gene, FC, fill = Type)) + geom_col(position = "dodge", colour = "Black") + theme_science() + 
  scale_fill_manual(values = c("white", "grey")) + labs(x = "", y = "Fold Change TE") + scale_y_continuous(limits = c(0,2.5), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.position = "bottom", legend.title = element_blank())
ggsave("FC_TE_Examples.tiff", height = 3, width = 4, units = "in")  
