library(edgeR)
library(scran)
library(scater)


# Pre-processing and normalisation of Tonsils dataset (GSE70580)-----

scRNAseq_path <- "/stornext/General/data/academic/lab_davis/Belz/ILC2_Eosinophil/"
load(paste0(scRNAseq_path,"GSE70580_RAW_counts_allCells.RData"))
test_cell <- read.table(paste0(scRNAseq_path,"GSM1810510_T74_P1_A9_ILC1_expression.txt"))

y_Ilc2[1:6,1:6]
library(SingleCellExperiment)
library(scran)
library(scater)

sce <- SingleCellExperiment(assays = list(counts = y_Ilc2))
rowData(sce) <- data.frame(test_cell[,1:2])


mito <- grep("^MT-", rowData(sce)$V1)
percent.mito <- colSums(counts(sce)[mito, ]) / colSums(counts(sce))
nGenes <- colSums(counts(sce) != 0)




keep1 <- rowSums(counts(sce) > 0) >= ncol(counts(sce))*0.05
keep2 <- !is.na(rowData(sce)$V1)
keep3 <- !duplicated(rowData(sce)$V1)
keep4 <-  percent.mito < 0.05
keep <- keep1 & keep2 & keep3 
table(keep)

sce <- sce[keep1 & keep2 & keep3, keep4]

clst <- quickCluster(sce, method="igraph", min.mean=0.5)
table(clst)
sce <- computeSumFactors(sce, cluster=clst)
summary(sizeFactors(sce))


sce <- normalize(sce)
cpm(sce) <- calculateCPM(sce, use_size_factors = TRUE)


fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
dec <- decomposeVar(fit=fit)
top.dec.hs <- dec[order(dec$bio, decreasing=TRUE), ]
top.dec.hs <- as.data.frame(top.dec.hs)
#head(top.dec, n=10L)

top <- 5000
hvg_hs_sce <- rownames(top.dec.hs)[1:top]



## Pre-process and normalise scRNAseq WT Skin ILC2 mouse (GSE117568) -------

WT_skin_ILC2 <- "GSE117568_RAW"
annot_path <- "../../genomes/DATA/Mus_musculus.gene_info"

dge_skin_ILC2 <- read10X(mtx = "GSM3303972_Adult_skin_matrix.mtx.gz", 
                  genes = "GSM3303972_Adult_skin_genes.tsv.gz",
                                    barcodes = "GSM3303972_Adult_skin_barcodes.tsv.gz",
                                                      path = WT_skin_ILC2,
                                                                        DGEList = TRUE)

ann <- alias2SymbolUsingNCBI(dge_skin_ILC2$genes$Symbol, 
                             required.columns=c("GeneID","Symbol","description"), 
                                                          gene.info.file= annot_path)
dge_skin_ILC2$genes <- cbind(dge_skin_ILC2$genes,
                             Official=ann$Symbol, GeneID=ann$GeneID,
                                                          description = ann$description)
head(dge_skin_ILC2$genes)


# convert mouse to human --------
musGenes <- dge_skin_ILC2$genes$Symbol
 
 MGI_homologs <- read.table("HMD_HumanPhenotype.rpt", sep = "\t", header = FALSE)

 dge_skin_ILC2$genes$humanHomolog <- MGI_homologs$V1[match(dge_skin_ILC2$genes$Symbol,MGI_homologs$V5)]


# proceed with normalisation then subset, replacing the missing genes with 0s and apply the classifier

mito <- grep("mitochondrial", dge_skin_ILC2$genes$description)
percent.mito <- colSums(dge_skin_ILC2$counts[mito, ]) / dge_skin_ILC2$samples$lib.size
nGenes <- colSums(dge_skin_ILC2$counts != 0)
dge_skin_ILC2$samples <- cbind(dge_skin_ILC2$samples, percent.mito=percent.mito, nGenes=nGenes)
head(dge_skin_ILC2$samples, n=10)
plot(dge_skin_ILC2$samples[,c("lib.size","nGenes")], pch=16, cex=0.7)


o <- order(rowSums(dge_skin_ILC2$counts), decreasing=TRUE)
dge_skin_ILC2 <- dge_skin_ILC2[o, ]
keep1 <- rowSums(dge_skin_ILC2$counts > 0) >= ncol(dge_skin_ILC2)*0.01
keep2 <- !is.na(dge_skin_ILC2$genes$Official)
keep3 <- !duplicated(dge_skin_ILC2$genes$Official)
keep4 <- dge_skin_ILC2$samples$percent.mito < 0.08
keep <- keep1 & keep2 & keep3 
table(keep)

dge_skin_ILC2 <- dge_skin_ILC2[keep1 & keep2 & keep3, keep4]




sce2 <- SingleCellExperiment(list(counts=dge_skin_ILC2$counts))
sce2

clst <- quickCluster(sce2, method="igraph", min.mean=0.5)
table(clst)
sce2 <- computeSumFactors(sce2, cluster=clst)
summary(sizeFactors(sce2))
libSize <- dge_skin_ILC2$samples$lib.size

plot(libSize/1e3, sizeFactors(sce2), log="xy", pch=16, cex=0.7,
xlab="Library size (thousands)", ylab="Size factor")

sce2 <- normalize(sce2)
cpm(sce2) <- calculateCPM(sce2, use_size_factors = TRUE)
rowData(sce2) <- dge_skin_ILC2$genes



fit <- trendVar(sce2, use.spikes=FALSE, loess.args=list(span=0.05))
dec <- decomposeVar(fit=fit)
top.dec.mm <- dec[order(dec$bio, decreasing=TRUE), ]
top.dec.mm <- as.data.frame(top.dec.mm)
#head(top.dec, n=10L)

top <- 5000
hvg_mm_sce2 <- rownames(top.dec.mm)[1:top]

# feature selection ----

shared_hvg <- intersect(as.character(rowData(sce)$V1[match(hvg_hs_sce, rowData(sce)$V2)]),
                 as.character(rowData(sce2)$humanHomolog[match(hvg_mm_sce2, rownames(sce2))]))


#filter highly variable genes based on the list of 700 genes
nanoString_genes <- xlsx::read.xlsx("Nanostring_PanCancer_Immune_Profiling_Panel_Gene_List.xlsx", sheetIndex = 2,
                                    rowIndex = 2:732)



shared_hvg <- shared_hvg[shared_hvg %in% nanoString_genes$Gene.Name]


# compute log2CPMs and ensure non of the entries are missing -----
table(shared_hvg %in% rowData(sce)$V1)

# shared_hvg <- shared_hvg[shared_hvg %in% rowData(sce)$V1] #

tonsils_data <- log2(cpm(sce)[match(shared_hvg, rowData(sce)$V1),] + 2)
table(complete.cases(tonsils_data))


wt_skin_ILC2 <- log2(cpm(sce2)[match(shared_hvg, rowData(sce2)$humanHomolog),] + 2)
table(complete.cases(wt_skin_ILC2))


rownames(tonsils_data) <- shared_hvg
rownames(wt_skin_ILC2) <- shared_hvg

y <- gsub("(.*)_(.*)_(.*)_(.*)", "\\4", colnames(tonsils_data))

# need to trun this into a binary classification problem
y <- as.factor(ifelse(y == "ILC2", "ILC2", "other"))
y <- relevel(y, ref = "other")




## Augment number of observations from ILC2 class to mitigate class imbalance ----

wt_skin_l2 <- wt_skin_ILC2/sqrt(colSums(wt_skin_ILC2^2))
tonsils_l2 <- tonsils_data[,y %in% "ILC2"]/sqrt(colSums(tonsils_data[,y %in% "ILC2"]^2))


library(FNN)
nn <- get.knnx(data =  t(tonsils_l2), 
               query = t(wt_skin_l2), k = 3) 


# confirm distances are actually cosine similarities
# select cells ----

keep <- rowSums(nn$nn.dist > 1.1) > 2
table(keep)



x <- cbind(tonsils_data, wt_skin_ILC2[,keep])



# augment y by cells chosen from wt skin mouse
y <- as.factor(c(ifelse(y == "ILC2", "ILC2", "other"), rep("ILC2", sum(keep))))
y <- relevel(y, ref = "other")
table(y)



### Training -----

library(xgboost)
set.seed(2329)

validation <- sample(1:ncol(x), size = floor(0.2*ncol(x)))


dtrain <- xgb.DMatrix(data = t(x[,-validation]), label= as.numeric(y[-validation]) - 1) # 1 if ILC2, 0 if other
dtest <- xgb.DMatrix(data = t(x[,validation]), label= as.numeric(y[validation]) - 1)


watchlist <- list(train=dtrain, test=dtest)


bst <- xgb.train(data=dtrain,
                 max.depth=5,
                 eta=0.05, 
                 nthread = 2, nrounds=9, watchlist=watchlist,
                 objective = "binary:logistic")



