library(tximport)
library(rhdf5)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Input results folder name
dir = "p53_mir34a_wt_genes_results" #Change to your ExperimentID

wd = file.path("C:", "cygwin64", "home", "marti", "RNA-seq", "P53_mir34a_wt", dir)
setwd(wd)

# Annotation databases
endb <- EnsDb.Mmusculus.v79
hendb <- EnsDb.Hsapiens.v86

df <- transcripts(endb, return.type = "data.frame")
h_df <- transcripts(hendb, return.type = "data.frame")

tx2gene <- data.frame(df$tx_id, df$gene_id)
htx2gene <- data.frame(h_df$tx_id, h_df$gene_id)

tx2gene_rmsk <- read.table("tx2gene_rmsk.txt", header = TRUE)

df1 <- genes(endb, return.type = "data.frame")
h_df1 <- genes(hendb, return.type = "data.frame")

eng2gs <- data.frame(df1$gene_id, df1$gene_name)
heng2gs <- data.frame(h_df1$gene_id, h_df1$gene_name)

colnames(eng2gs) <- c("GENEID", "Gene_Symbol")
colnames(heng2gs) <- c("GENEID", "Gene_Symbol")


# results direcrotry samples.txt is a table of all processed SRR_IDs.
samples <- read.table("paul_samples.txt", header = TRUE)
files <- file.path("C:", "cygwin64", "home", "marti", "RNA-seq", "P53_mir34a_wt", dir, samples$Run, "abundance.h5")
names(files) <- samples$Run

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
abundances <- txi.kallisto$abundance

heatmap(abundances[rowSds(abundances) > 10,], scale = "row", Colv = NA, distfun = function(x) as.dist((1-cor(t(x)))))

# a_heng2gs <- heng2gs[heng2gs$GENEID %in% rownames(abundances),]
# abundances <- abundances[a_heng2gs$GENEID, ]
# rownames(abundances) <- a_heng2gs$Gene_Symbol


# enter in replicate information
column_data <- data.frame(cell.type = as.factor(samples$details))
rownames(column_data) <- colnames(txi.kallisto$counts)


# DEseq
dds <- DESeqDataSetFromTximport(txi.kallisto, column_data, ~cell.type)
rld <- rlogTransformation(dds)  #used for PCA
res <- DESeq(dds)

#PCA
ntop = 500
Pvars <- rowVars(assay(rld))
select <- order(Pvars, decreasing = T)[seq_len(min(ntop,
                                                   length(Pvars)))]

PCA <- prcomp(t(assay(rld)[select,]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
dataGG = data.frame(PC1 = PCA$x[,1], PC2=PCA$x[,2],
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4],
                    condition = colData(rld)$cell.type)
qplot(PC1,PC2, data=dataGG, color = condition,
      main = "PC1 vs PC2, top variable repeats", size = I(2))+theme_classic()



#Change to samples you want to compare
res.vis <- results(res, c("cell.type", "mir-34a", "wt"))
vis_eng2gs <- eng2gs[eng2gs$GENEID %in% rownames(res.vis),]
res.vis <- res.vis[vis_eng2gs$GENEID, ]
rownames(res.vis) <- vis_eng2gs$Gene_Symbol



# Plotting
label.genes <- c() #What genes/repeats do you want to see?
#res.vis <- res.vis[abs(res.vis$log2FoldChange) < 15,]
res.vis <- res.vis[log2(res.vis$baseMean) > 1, ]
with(res.vis, qplot(main = "34a_KO_vs_WT", log2(baseMean), log2FoldChange, label=ifelse(log2(baseMean) > 5 & abs(log2FoldChange) > 2 & pvalue < 0.01, as.character(rownames(res.vis)), "") ,color = ifelse(abs(log2FoldChange) > 1 & pvalue < 0.05, "significant", 
                                                                   "not-significant"))) + scale_color_manual(name = "", values = c("darkgrey", "red"))+ geom_text(color = "black", size = 4)+theme_classic() 

