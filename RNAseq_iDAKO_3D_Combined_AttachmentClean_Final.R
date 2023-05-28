# ================
#  Load libraries
# ================

library("DESeq2")
library("ggplot2")
library("EnhancedVolcano")
library("clusterProfiler")
library(pheatmap)
library(grid)
library(VennDiagram)
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("topGO")
library(tibble)
library(EnhancedVolcano)
library(dplyr)
library(openxlsx)
library(readxl)

# Set Convert (used to convert Gene names to ID numbers)
Convert <- as.data.frame(org.Mm.egSYMBOL)

# Set path to data
PATH = "C:/Users/RuneInglev/OneDrive - RSP Systems A S/Skrivebord/R stuff/"



# ========================
#  iDAKO Data Preparation 
# ========================

rawdata.iDAKO <- read.delim(paste0(PATH,"iDAKO_RNAseq.txt"))
rownames(rawdata.iDAKO) <- rawdata.iDAKO$Symbol

# Prepare the names for the different replicates
coldata.iDAKO <- matrix(nrow=22, ncol=3) 
colnames(coldata.iDAKO) <- c("Condition", "Replicate", "Names")

condition.iDAKO <- c("gWA_iDAKO", "gWA_Vehicle", "iWA_iDAKO", "iWA_Vehicle")

coldata.iDAKO[,1] <- c(rep(condition.iDAKO[1],7),rep(condition.iDAKO[2],4),rep(condition.iDAKO[3],7),rep(condition.iDAKO[4],4))
coldata.iDAKO[1:7,2] <- 1:7
coldata.iDAKO[8:11,2] <- 1:4
coldata.iDAKO[12:18,2] <- 1:7
coldata.iDAKO[19:22,2] <- 1:4
coldata.iDAKO[,3] <- paste0(coldata.iDAKO[,1],"_", coldata.iDAKO[,2])

# Extract the countdata from the raw dataframe 
countdata.iDAKO <- rawdata.iDAKO[,9:30] 

# Give names to the columns
colnames(countdata.iDAKO) <- coldata.iDAKO[,3] 
dim(countdata.iDAKO)

#Sorting out everything which in none of the conditions exceeds 20 reads 
countdata.iDAKO <- as.matrix(countdata.iDAKO[apply(countdata.iDAKO,1,max)>20,]) 
dim(countdata.iDAKO)

# Remove bad samples (found from PCA)
countdata.iDAKO <- countdata.iDAKO[, -c(2, 14)]
coldata.iDAKO <- coldata.iDAKO[-c(2,14),]



# ===============================
#  3D Model Data Preparation
# ===============================

rawdata.3D <- read.delim(paste0(PATH,"iWAT-3D-LipExp.txt"))
rownames(rawdata.3D) <- rawdata.3D$Symbol

# Prepare names for the replicates
coldata.3D <- matrix(nrow=12, ncol=3)
colnames(coldata.3D) <- c("Condition", "Replicate", "Names")

condition.3D <- c("DMSO_lean", "DMSO_obese", "Lip_lean", "Lip_obese")

coldata.3D[,1] <- c(rep(condition.3D[1],3),rep(condition.3D[2],3),rep(condition.3D[3],3),rep(condition.3D[4],3))
coldata.3D[,2] <- rep(c("a","b","c","a","b","c","a","b","c","a","b","c"))
coldata.3D[,3] <- paste0(coldata.3D[,1],"_", coldata.3D[,2])

# Extract the countdata from the raw dataframe
countdata.3D <- rawdata.3D[,9:20] 

# Naming the columns. Using information from previous coldata file (above)
colnames(countdata.3D) <- coldata.3D[,3] 
dim(countdata.3D)

#Sorting out everything which in none of the conditions exceeds 100 reads 
countdata.3D <- as.matrix(countdata.3D[apply(countdata.3D,1,max)>100,]) 
dim(countdata.3D)



# ================================
# Combined dataset preparation
# ================================

# First, since the same genes are not necesarily in both datasets so we must find the genes that are in common
genes_common <- intersect(rownames(countdata.3D),rownames(countdata.iDAKO))

# Grab the genes (rows) from each count dataset and 
countdata.combined <- cbind(countdata.3D[genes_common,], countdata.iDAKO[genes_common,])

# Create the combined column data
coldata.combined <- rbind(coldata.3D,coldata.iDAKO)

# Combined conditions
condition.combined <- rbind(condition.3D,condition.iDAKO)



# ===========================
#  DESeq2 Analysis on iDAKO
# ===========================

# Create the model
DDS.iDAKO <- DESeqDataSetFromMatrix(countdata.iDAKO, coldata.iDAKO, design = ~Condition)
DDS.iDAKO <- DESeq(DDS.iDAKO)
rlD.iDAKO <- rlog(DDS.iDAKO)
rlD.iDAKO.df <- assay(rlD.iDAKO)

colnames(rlD.iDAKO.df) <- paste0("rlog_",coldata.iDAKO[,3])

#Making normalized counts.iDAKO to be used later (are nice for plotting different genes)
counts.iDAKO.norm <- counts(DDS.iDAKO, normalized = TRUE) 

# Calculate the means for each of the 4 Conditions
for(i in c(1,2,3,4)){
  tmp_rlog <- rowMeans(rlD.iDAKO.df[,coldata.iDAKO[,1] == condition.iDAKO[i]])
  tmp_norm <- rowMeans(counts.iDAKO.norm[,coldata.iDAKO[,1] == condition.iDAKO[i]])
  
  rlD.iDAKO.df <- cbind(rlD.iDAKO.df,tmp_rlog)
  counts.iDAKO.norm <- cbind(counts.iDAKO.norm,tmp_norm)
  
  colnames(rlD.iDAKO.df)[ncol(rlD.iDAKO.df)] <- paste0("rlog_average_", condition.iDAKO[i])
  colnames(counts.iDAKO.norm)[ncol(counts.iDAKO.norm)] <- paste0("counts.average_", condition.iDAKO[i])
}

# Combine the normalized counts and rlog values into one matrix
data.iDAKO <- cbind(counts.iDAKO.norm, rlD.iDAKO.df)

#Making contrasts, and adding them to data.1 dataframe. The data.1 dataframe will end up containing all my information
tmp <- as.data.frame(results(DDS.iDAKO, contrast=c("Condition", "gWA_iDAKO", "gWA_Vehicle")))
colnames(tmp) <- paste0("gWA_iDAKO_vs_gWA-Vehicle", colnames(tmp))
data.iDAKO <- cbind(data.iDAKO, tmp)
tmp <- as.data.frame(results(DDS.iDAKO, contrast=c("Condition", "iWA_iDAKO", "iWA_Vehicle")))
colnames(tmp) <- paste0("iWA_iDAKO_vs_iWA_Vehicle", colnames(tmp))
data.iDAKO <- cbind(data.iDAKO, tmp)
tmp <- as.data.frame(results(DDS.iDAKO, contrast=c("Condition", "gWA_Vehicle", "iWA_Vehicle")))
colnames(tmp) <- paste0("gWA_Vehicle_vs_iWA_Vehicle", colnames(tmp))
data.iDAKO <- cbind(data.iDAKO, tmp)

#Fjern NAs:
data.iDAKO[is.na(data.iDAKO)] <- 0

# Used by volcano plots
results.iDAKO.iWA <- results(DDS.iDAKO,
                               contrast = c("Condition", "iWA_iDAKO", "iWA_Vehicle"))
results.iDAKO.iWA <- lfcShrink(DDS.iDAKO,
                                 contrast = c("Condition", "iWA_iDAKO", "iWA_Vehicle"), res=results.iDAKO.iWA, type = "normal")




# =============================
#  DESeq2 Analysis on 3D model
# =============================

DDS.3D <- DESeqDataSetFromMatrix(countdata.3D, coldata.3D, design = ~Condition)
DDS.3D <- DESeq(DDS.3D)
rlD.3D <- rlog(DDS.3D)
rlD.3D.df <- assay(rlD.3D)

colnames(rlD.3D.df) <- paste0("rlog_",coldata.3D[,3])

#Making normalized counts to be used later (are nice for plotting different genes)
counts.3D.norm <- counts(DDS.3D, normalized = TRUE) 

#Udregner en average rlog-v?rdi for hver condition
for(i in c(1,2,3,4)){
  tmp_rlog <- rowMeans(rlD.3D.df[,coldata.3D[,1] == condition.3D[i]])
  tmp_norm <- rowMeans(counts.3D.norm[,coldata.3D[,1] == condition.3D[i]])
  
  rlD.3D.df <- cbind(rlD.3D.df, tmp_rlog)
  counts.3D.norm <- cbind(counts.3D.norm, tmp_norm)
  
  colnames(rlD.3D.df)[ncol(rlD.3D.df)] <- paste0("rlog_average_", condition.3D[i])
  colnames(counts.3D.norm)[ncol(counts.3D.norm)] <- paste0("counts.average_", condition.3D[i])
}

#Danner en stor matrix af de forskellige matrices
data.3D <- cbind(counts.3D.norm, rlD.3D.df)

#Making contrasts, and adding them to data.1 dataframe. The data.1 dataframe will end up containing all my information
tmp <- as.data.frame(results(DDS.3D, contrast=c("Condition", "Lip_lean", "DMSO_lean")))
colnames(tmp) <- paste0("Lip_lean_vs_DMSO_lean", colnames(tmp))
data.3D <- cbind(data.3D, tmp)
tmp <- as.data.frame(results(DDS.3D, contrast=c("Condition", "DMSO_obese", "DMSO_lean")))
colnames(tmp) <- paste0("DMSO_obese_vs_DMSO_lean", colnames(tmp))
data.3D <- cbind(data.3D, tmp)
tmp <- as.data.frame(results(DDS.3D, contrast=c("Condition", "Lip_obese", "DMSO_obese")))
colnames(tmp) <- paste0("Lip_obese_vs_DMSO_obese", colnames(tmp))
data.3D <- cbind(data.3D, tmp)

#Fjern NAs:
data.3D[is.na(data.3D)] <- 0

# Used by volcano plots
results.3D <- results(DDS.3D,
                        contrast = c("Condition", "Lip_lean", "DMSO_lean"))
results.3D <- lfcShrink(DDS.3D,
                          contrast = c("Condition", "Lip_lean", "DMSO_lean"), res=results.3D, type = "normal")




# ===================================
#  Combined dataset DESeq analysis
# ===================================

DDS.combined <- DESeqDataSetFromMatrix(countdata.combined, coldata.combined, design = ~Condition)
DDS.combined <- DESeq(DDS.combined)

rlD.combined <- rlog(DDS.combined)
rlD.combined.df <- assay(rlD.combined)

#Angiver nu at det er rlog-v?rdier
colnames(rlD.combined.df) <- paste0("rlog_",coldata.combined[,3])

#Making normalized counts.iDAKO to be used later (are nice for plotting different genes)
counts.combined.norm <- counts(DDS.combined, normalized = TRUE) 

#Udregner en average rlog-v?rdi for hver condition
for (i in c(1,2,3,4,5,6,7,8)){
  tmp_rlog <- rowMeans(rlD.combined.df[, coldata.combined[,1] == condition.combined[i]])
  tmp_norm <- rowMeans(counts.combined.norm[, coldata.combined[,1] == condition.combined[i]])
  
  rlD.combined.df <- cbind(rlD.combined.df, tmp_rlog)
  counts.combined.norm <- cbind(counts.combined.norm, tmp_norm)
  
  colnames(rlD.combined.df)[ncol(rlD.combined.df)] <- paste0("rlog_average_", condition.combined[i])
  colnames(counts.combined.norm)[ncol(counts.combined.norm)] <- paste0("counts.average_", condition.combined[i])
  
}

#Danner en stor matrix af de forskellige matrices
data.combined <- cbind(counts.combined.norm, rlD.combined.df)

#Making contrasts, and adding them to data.1 dataframe. The data.1 dataframe will end up containing all my information
tmp <- as.data.frame(results(DDS.combined, contrast=c("Condition", "iWA_iDAKO", "iWA_Vehicle")))
colnames(tmp) <- paste0("iWA_iDAKO_vs_iWA_Vehicle", colnames(tmp))
data.combined <- cbind(data.combined, tmp)
tmp <- as.data.frame(results(DDS.combined, contrast=c("Condition", "Lip_lean", "DMSO_lean")))
colnames(tmp) <- paste0("Lip_lean_vs_DMSO_lean", colnames(tmp))
data.combined <- cbind(data.combined, tmp)
tmp <- as.data.frame(results(DDS.combined, contrast=c("Condition", "DMSO_obese", "DMSO_lean")))
colnames(tmp) <- paste0("DMSO_obese_vs_DMSO_lean", colnames(tmp))
data.combined <- cbind(data.combined, tmp)
tmp <- as.data.frame(results(DDS.combined, contrast=c("Condition", "Lip_obese", "DMSO_obese")))
colnames(tmp) <- paste0("Lip_obese_vs_DMSO_obese", colnames(tmp))
data.combined <- cbind(data.combined, tmp)

data.combined[is.na(data.combined)] <- 0



# =======================================
#  INDUCED AND REPRESSED GENES ANALYSIS
# =======================================

#  3D Model genes
LipMaintGenes.3D <- data.3D[ data.3D$Lip_lean_vs_DMSO_leanlog2FoldChange < -0.58 & data.3D$Lip_lean_vs_DMSO_leanpadj <0.05,] %>%
  arrange((Lip_lean_vs_DMSO_leanlog2FoldChange))

LipRepGenes.3D <- data.3D[ data.3D$Lip_lean_vs_DMSO_leanlog2FoldChange > 0.58 & data.3D$Lip_lean_vs_DMSO_leanpadj <0.05,] %>%
  arrange(desc(Lip_lean_vs_DMSO_leanlog2FoldChange))

LipConstGenes.3D <- data.3D[data.3D$Lip_lean_vs_DMSO_leanpadj > 0.97,]

LipMaintGenes.3D$symbol <- rownames(LipMaintGenes.3D)
LipRepGenes.3D$symbol <- rownames(LipRepGenes.3D)
LipConstGenes.3D$symbol <- rownames(LipConstGenes.3D)


#  iDAKO genes iWA
LipMaintGenes.iWA <- data.iDAKO[ data.iDAKO$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange < -0.5 & data.iDAKO$iWA_iDAKO_vs_iWA_Vehiclepadj <0.05,] %>%
  arrange((iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange))

LipRepGenes.iWA <- data.iDAKO[ data.iDAKO$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange > 0.5 & data.iDAKO$iWA_iDAKO_vs_iWA_Vehiclepadj <0.05,] %>%
  arrange(desc(iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange))

LipConstGenes.iWA <- data.iDAKO[data.iDAKO$iWA_iDAKO_vs_iWA_Vehiclepadj > 0.97,] #Here, DE == 5% FDR and at least two fold change (log2FC > 1)

LipMaintGenes.iWA$symbol <- rownames(LipMaintGenes.iWA)
LipRepGenes.iWA$symbol <- rownames(LipRepGenes.iWA)
LipConstGenes.iWA$symbol <- rownames(LipConstGenes.iWA)


#  Combined Genes IWA
LipMaintGenes.comb.iWA <- data.combined[ data.combined$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange < -0.5 & data.combined$iWA_iDAKO_vs_iWA_Vehiclepadj <0.05,] %>%
  arrange((iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange))

LipRepGenes.comb.iWA <- data.combined[ data.combined$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange > 0.5 & data.combined$iWA_iDAKO_vs_iWA_Vehiclepadj <0.05,] %>%
  arrange(desc(iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange))

LipConstGenes.comb.iWA <- data.combined[data.combined$iWA_iDAKO_vs_iWA_Vehiclepadj > 0.97,] #Here, DE == 5% FDR and at least two fold change (log2FC > 1)

LipMaintGenes.comb.iWA$symbol <- rownames(LipMaintGenes.comb.iWA)
LipRepGenes.comb.iWA$symbol <- rownames(LipRepGenes.comb.iWA)
LipConstGenes.comb.iWA$symbol <- rownames(LipConstGenes.comb.iWA)


# Combined genes 3D
LipMaintGenes.comb.3D <- data.combined[ data.combined$Lip_lean_vs_DMSO_leanlog2FoldChange < -0.5 & data.combined$Lip_lean_vs_DMSO_leanpadj <0.05,] %>%
  arrange((Lip_lean_vs_DMSO_leanlog2FoldChange))

LipRepGenes.comb.3D <- data.combined[ data.combined$Lip_lean_vs_DMSO_leanlog2FoldChange > 0.5 & data.combined$Lip_lean_vs_DMSO_leanpadj <0.05,] %>%
  arrange(desc(Lip_lean_vs_DMSO_leanlog2FoldChange))

LipConstGenes.comb.3D <- data.combined[data.combined$Lip_lean_vs_DMSO_leanpadj > 0.97,]

LipMaintGenes.comb.3D$symbol <- rownames(LipMaintGenes.comb.3D)
LipRepGenes.comb.3D$symbol <- rownames(LipRepGenes.comb.3D)
LipConstGenes.comb.3D$symbol <- rownames(LipConstGenes.comb.3D)


# Genes that are significant and maintained from both sets (combined DESeq2 Analysis)
Significant.comb <- data.combined[data.combined$Lip_lean_vs_DMSO_leanpadj < 0.05 & data.combined$iWA_iDAKO_vs_iWA_Vehiclepadj < 0.05,]

SignificantMaint.comb.iWA <- Significant.comb[ Significant.comb$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange < -0.5,]
SignificantMaint.comb.3D <- Significant.comb[ Significant.comb$Lip_lean_vs_DMSO_leanlog2FoldChange < -0.5,]

total <- union(SignificantMaint.comb.3D,SignificantMaint.comb.iWA)


# Significant genes in both datasets (from each of the individual DESeq Analyses)
data.3D.both <- data.3D[genes_common,]
data.iDAKO.both <- data.iDAKO[genes_common,]

Significant.both.3D <- data.3D.both[data.3D.both$Lip_lean_vs_DMSO_leanpadj < 0.05 & data.iDAKO.both$iWA_iDAKO_vs_iWA_Vehiclepadj < 0.05,]
Significant.both.iWA <- data.iDAKO.both[data.3D.both$Lip_lean_vs_DMSO_leanpadj < 0.05 & data.iDAKO.both$iWA_iDAKO_vs_iWA_Vehiclepadj < 0.05,]

SignificantMaint.both.3D <- Significant.both.3D[Significant.both.3D$Lip_lean_vs_DMSO_leanlog2FoldChange < -0.5,]
SignificantMaint.both.iWA <- Significant.both.iWA[Significant.both.iWA$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange < -0.5,]

SignificantMaint.both <- union(rownames(SignificantMaint.both.iWA),rownames(SignificantMaint.both.3D))


# ====================================================================
#  Plotting 2FoldChange of intersection genes (from combined dataset)
# ======================================================================

maint_inter <- intersect(rownames(LipMaintGenes.comb.3D),rownames(LipMaintGenes.comb.iWA))
rep_inter <- intersect(rownames(LipRepGenes.comb.3D),rownames(LipRepGenes.comb.iWA))

# Subset the data for the two groups
group1_data <- data.combined[maint_inter,]
group2_data <- data.combined[rep_inter, ]

# Create a combined data frame for plotting
combined_data <- rbind(data.frame(group = "Maintained", group1_data),
                       data.frame(group = "Repressed", group2_data))

# Plot the scatterplot-
ggplot(combined_data, aes(x = iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange, y = Lip_lean_vs_DMSO_leanlog2FoldChange, color = group)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  xlab("Variable X") +
  ylab("Variable Y") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  scale_color_manual(values = c("darkred", "darkgreen")) +
  theme_minimal()



# =======================
#  KEGG ANalysis Data
# =======================

LipMaintKEGG.3D <- as.data.frame(enrichKEGG(Convert[ Convert$symbol %in% LipMaintGenes.3D$symbol,1], organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1))
LipMaintKEGG.3D <- LipMaintKEGG.3D[ LipMaintKEGG.3D$qvalue <= 0.01,]

LipRepKEGG.3D <- as.data.frame(enrichKEGG(Convert[ Convert$symbol %in% LipRepGenes.3D$symbol,1], organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1))
LipRepKEGG.3D <- LipRepKEGG.3D[ LipRepKEGG.3D$qvalue <= 0.01,]

write.xlsx(as.data.frame(LipMaintKEGG.3D), file = paste0(PATH,"LipMaint3D_Kegg.xlsx"), sheetName = "3D")
write.xlsx(as.data.frame(LipRepKEGG.3D), file = paste0(PATH,"LipRep3D_Kegg.xlsx"), sheetName = "3D")

LipMaintKEGG.iWA <- as.data.frame(enrichKEGG(Convert[ Convert$symbol %in% LipMaintGenes.iWA$symbol,1], organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1))
LipMaintKEGG.iWA <- LipMaintKEGG.iWA[ LipMaintKEGG.iWA$qvalue <= 0.01,]

LipRepKEGG.iWA <- as.data.frame(enrichKEGG(Convert[ Convert$symbol %in% LipRepGenes.iWA$symbol,1], organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1))
LipRepKEGG.iWA <- LipRepKEGG.iWA[ LipRepKEGG.iWA$qvalue <= 0.01,]

write.xlsx(as.data.frame(LipMaintKEGG.iWA), file = paste0(PATH,"LipMaintiWA_Kegg.xlsx"), sheetName = "iWAT")
write.xlsx(as.data.frame(LipRepKEGG.iWA), file = paste0(PATH,"LipRepiWA_Kegg.xlsx"), sheetName = "iWAT")




# ===================
#  GETTING TOP GENES
# ===================

topgenes.iWAT <- as.numeric(LipMaintGenes.iWA$counts.average_iWA_iDAKO)
print(topgenes.iWAT)

genenames <- rownames(LipMaintGenes.iWA)
LipMaintGenes.iWA$counts.iDAKO_average_iWA_iDAKO <- as.numeric(as.character(LipMaintGenes.iWA$counts.average_iWA_iDAKO))
Lol <- LipMaintGenes.iWA$counts.average_iWA_iDAKO

numeric_Maint <- as.numeric(as.character(Lol))
print(numeric_Maint)
print(LipMaintGenes.iWA$counts.average_iWA_iDAKO)

genenames_counts_iDAKO_iWAT <- cbind(genenames,numeric_Maint)
dataframe_genenames_counts_iDAKO_iWAT <- as.data.frame(genenames_counts_iDAKO_iWAT)
topgenes.iWAT.iDAKO.arranged <- dataframe_genenames_counts_iDAKO_iWAT[order(dataframe_genenames_counts_iDAKO_iWAT$numeric_Maint),]
print(topgenes.iWAT.iDAKO.arranged)

topgenes_arranged_by_numeric_Maint <- topgenes.iWAT.iDAKO.arranged

# Convert the numeric_Maint column to numeric
topgenes.iWAT.iDAKO.arranged$numeric_Maint <- as.numeric(as.character(topgenes.iWAT.iDAKO.arranged$numeric_Maint))

# Sort the data frame based on the numeric_Maint column
topgenes_arranged_by_numeric_Maint <- topgenes.iWAT.iDAKO.arranged[order(topgenes.iWAT.iDAKO.arranged$numeric_Maint), ]

topgenes.iWAT.iDAKO.arranged$numeric_Maint <- as.numeric(as.character(topgenes.iWAT.iDAKO.arranged$numeric_Maint))

# Sort the data frame based on the numeric_Maint column
topgenes_arranged_by_numeric_Maint <- topgenes.iWAT.iDAKO.arranged[order(topgenes.iWAT.iDAKO.arranged[["numeric_Maint"]]), ]



# ====================
# PCA Plots
# ====================

# PCA FOR IDAKO
PCAData.iDAKO <- plotPCA(rlD.iDAKO, intgroup=c("Condition", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(PCAData.iDAKO, "percentVar"))
ggplot(PCAData.iDAKO, aes(PC1, PC2, color=Condition, shape=Replicate)) +
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

# PCA for 3D
PCAData.3D <- plotPCA(rlD.3D, intgroup=c("Condition", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(PCAData.3D, "percentVar"))
ggplot(PCAData.3D, aes(PC1, PC2, color=Condition, shape=Replicate)) +
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

# PCA for Combined
PCAData.combined <- plotPCA(rlD.combined, intgroup=c("Condition", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(PCAData.combined, "percentVar"))
ggplot(PCAData.combined, aes(PC1, PC2, color=Condition)) +
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))



# =============
#  Plotting MA-plot
# ==============

plot(log10(data.3D$Lip_lean_vs_DMSO_leanbaseMean),
     data.3D$Lip_lean_vs_DMSO_leanlog2FoldChange,
     main="",
     xlab = "Log10(Mean of normalised counts)",
     ylab = "",
     title(ylab = "Log2FoldChange", line = 2.5, cex.lab = 1.5),
     col="Darkgrey",
     pch = 20,
     ylim = c(-3,3),
     xlim = c(0.5,6),
     las = 1,
     cex.lab = 1.5) + 
  points(log10(LipRepGenes.3D$Lip_lean_vs_DMSO_leanbaseMean), LipRepGenes.3D$Lip_lean_vs_DMSO_leanlog2FoldChange, col = "Darkgreen", pch = 20, ylim = c(-10,10)) +
  text(log10(LipRepGenes.3D$Lip_lean_vs_DMSO_leanbaseMean[1:15]), LipRepGenes.3D$Lip_lean_vs_DMSO_leanlog2FoldChange[1:15], labels = rownames(LipRepGenes.3D)[1:15], pch = 15, adj = c(-0.1,-1), cex = 0.6, ylim = c(-10,10)) +
  points(log10(LipMaintGenes.3D$Lip_lean_vs_DMSO_leanbaseMean), LipMaintGenes.3D$Lip_lean_vs_DMSO_leanlog2FoldChange, col = "Darkred", pch = 20, ylim = c(-10,10)) +
  text(log10(LipMaintGenes.3D$Lip_lean_vs_DMSO_leanbaseMean[1:15]), LipMaintGenes.3D$Lip_lean_vs_DMSO_leanlog2FoldChange[1:15], labels = rownames(LipMaintGenes.3D)[1:15], pch = 15, adj = c(0.7,2), cex = 0.6, ylim = c(-10,10)) +
  abline(0.58,0,col="black") + abline(-0.58,0,col="black") + text(5.5, -2, paste0(format(nrow(LipMaintGenes.3D), big.mark = " ")," genes"), col="Darkred", cex = 1.5) +
  text(5.5, 2, paste0(format(nrow(LipRepGenes.3D), big.mark = " ")," genes"), col="Darkgreen", cex = 1.5)

plot(log10(data.iDAKO$iWA_iDAKO_vs_iWA_VehiclebaseMean),
     data.iDAKO$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange,
     main="",
     xlab = "Log10(Mean of normalised counts)",
     ylab = "",
     col="Darkgrey",
     pch = 20,
     ylim = c(-5,5),
     xlim = c(0.5,7),
     las = 1,
     cex.lab = 1.5) +
  title(ylab = "Log2FoldChange", line = 2.5, cex.lab = 1.5) +
  points(log10(LipRepGenes.iWA$iWA_iDAKO_vs_iWA_VehiclebaseMean), LipRepGenes.iWA$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange, col = "Darkgreen", pch = 20, ylim = c(-10,10)) +
  text(log10(LipRepGenes.iWA$iWA_iDAKO_vs_iWA_VehiclebaseMean[1:15]), LipRepGenes.iWA$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange[1:15], labels = rownames(LipRepGenes.iWA)[1:15], pch = 15, adj = c(-0.1,-1), cex = 0.6, ylim = c(-10,10)) +
  points(log10(LipMaintGenes.iWA$iWA_iDAKO_vs_iWA_VehiclebaseMean), LipMaintGenes.iWA$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange, col = "Darkred", pch = 20, ylim = c(-10,10)) +
  text(log10(LipMaintGenes.iWA$iWA_iDAKO_vs_iWA_VehiclebaseMean[1:15]), LipMaintGenes.iWA$iWA_iDAKO_vs_iWA_Vehiclelog2FoldChange[1:15], labels = rownames(LipMaintGenes.iWA)[1:15], pch = 15, adj = c(0.7,2), cex = 0.6, ylim = c(-10,10)) +
  abline(0.58,0,col="black") + abline(-0.58,0,col="black") + text(5.5, -2, paste0(format(nrow(LipMaintGenes.iWA), big.mark = " ")," genes"), col="Darkred", cex = 1.5) +
  text(5.5, 2, paste0(format(nrow(LipRepGenes.iWA), big.mark = " ")," genes"), col="Darkgreen", cex = 1.5)



# ===================
#  Volcano Plots
# ===================

# 3D model

VolcanoPlot.3D <- EnhancedVolcano(results.3D,
                               lab = rownames(results.3D),
                               x = "log2FoldChange",
                               y = "pvalue",
                               
                               gridlines.major = FALSE,
                               gridlines.minor = FALSE,
                               ylim = c(0,-log10(10e-49)),
                               xlim = c(-5,5),
                               captionLabSize = 10,
                               
                               
                               xlab = bquote(~Log[2]~ 'FC (Lip/Vehicle)'),
                               FCcutoff = 0.5,
                               title = "",
                               subtitle = "",
                               cutoffLineType = 'dashed',
                               cutoffLineCol = 'black',
                               cutoffLineWidth = 0.5,
                               
                               legendLabels=c('NS','Log2 FC','Adjusted p-value',
                                              'Adjusted p-value & Log2 FC'),
                               legendPosition = 'bottom',
                               legendLabSize = 10,
                               legendIconSize = 3.0)

plot(VolcanoPlot.3D)

# ====
#  iDAKO
# ====

VolcanoPlot.iWA <- EnhancedVolcano(results.iDAKO.iWA,
                               lab = rownames(results.iDAKO.iWA),
                               x = "log2FoldChange",
                               y = "pvalue",
                               
                               gridlines.major = FALSE,
                               gridlines.minor = FALSE,
                               ylim = c(0,-log10(10e-49)),
                               xlim = c(-5,5),
                               captionLabSize = 10,
                               
                               labSize = 3,
                               
                               xlab = bquote(~Log[2]~ 'FC (iDAKO/Vehicle)'),
                               FCcutoff = 0.58,
                               title = "",
                               subtitle = "",
                               cutoffLineType = 'dashed',
                               cutoffLineCol = 'black',
                               cutoffLineWidth = 0.5,
                               
                               legendLabels=c('NS','Log2 FC','Adjusted p-value',
                                              'Adjusted p-value & Log2 FC'),
                               legendPosition = 'bottom',
                               legendLabSize = 10,
                               legendIconSize = 3.0)

plot(VolcanoPlot.iWA)




# ====================
#  BOXPLOTS
# ====================

#  3D model

wilcox.test(LipRepGenes.3D$rlog_average_DMSO_lean, LipMaintGenes.3D$rlog_average_DMSO_lean, alternative = "two.sided")
wilcox.test(LipRepGenes.3D$rlog_average_DMSO_lean, LipConstGenes.3D$rlog_average_DMSO_lean, alternative = "two.sided")
wilcox.test(LipConstGenes.3D$rlog_average_DMSO_lean, LipMaintGenes.3D$rlog_average_DMSO_lean,alternative = "two.sided")

# Define the data and groups
rlogRepr.3D <- LipRepGenes.3D$rlog_average_DMSO_lean
rlogConst.3D <- LipConstGenes.3D$rlog_average_DMSO_lean
rlogMaint.3D <- LipMaintGenes.3D$rlog_average_DMSO_lean
GroupData.3D <- c(rlogMaint.3D, rlogRepr.3D, rlogConst.3D)

# Create a factor with the desired order of levels
groups.3D <- c(rep("Maintained", length(rlogMaint.3D)),
            rep("Repressed", length(rlogRepr.3D)),
            rep("Constant", length(rlogConst.3D)))

groups.3D <- factor(groups.3D, levels = c("Maintained", "Repressed", "Constant"))
            
# Define the colors for each group
colors <- c("darkred", "darkgreen", "grey")

# Plot the boxplot

par(cex.axis = 1.5) # Change the number 1.5 to the desired size factor for group labels

boxplot(GroupData.3D ~ groups.3D, col = colors, main = "", xlab = "", ylab = "", cex = 0)

# If you want to change the size of the x and y titles as well, set the size factor for cex.lab
par(cex.lab = 1.5) # Change the number 1.5 to the desired size factor for axis titles

# Now, add the x and y titles with the modified size:
title(xlab = "", ylab = "rlog")




# ====
#  iDAKO iWAT
# =====

# Define the data and groups
rlogRepr.iWA <- LipRepGenes.iWA$rlog_average_iWA_Vehicle
rlogConst.iWA <- LipConstGenes.iWA$rlog_average_iWA_Vehicle
rlogMaint.iWA <- LipMaintGenes.iWA$rlog_average_iWA_Vehicle

GroupData.iWA <- c(rlogMaint.iWA, rlogRepr.iWA, rlogConst.iWA)
groups.iWA <- c(rep("Maintained", length(rlogMaint.iWA)),
            rep("Repressed", length(rlogRepr.iWA)),
            rep("Constant", length(rlogConst.iWA)))

groups.iWA <- factor(groups.iWA, levels = c("Maintained", "Repressed", "Constant"))

# Define the colors for each group
colors <- c("darkred", "darkgreen", "grey")

# Plot the boxplot
boxplot(GroupData.iWA ~ groups.iWA, col = colors, main = "Expression level of gene group (iWA)", xlab = "", ylab = "rlog")



# ==============================
#  KEGG Analysis Heat map iDAKO
# ==============================

KEGG_analysis.iDAKO <- read_excel(paste0(PATH,"/KeggAnalysis_iDAKO.xlsx"))
KEGG.iDAKO <- as.data.frame(KEGG_analysis.iDAKO)
rownames(KEGG.iDAKO) <- KEGG.iDAKO$Pathway 

#The first column of our excel file containing the name of the pathways has the header "KEGG_repressed". Here we turn those names into the rownames 
KEGG.iDAKO <- KEGG.iDAKO[,2:3] #Column 2:3 in this case contain the -log(qValues) 

#Saturating it a bit. So all genes which a scaled rlog value greater than 2.5 should just be viewed a 2.5.
KEGG.iDAKO[KEGG.iDAKO> 10] <- 10

#Everything lower than -2 should just be called -2
KEGG.iDAKO[KEGG.iDAKO< 0] <- 0

#The first column of our excel file containing the name of the pathways has the header "KEGG_repressed". Here we turn those names into the rownames 
KEGG.iDAKO.maint <- KEGG.iDAKO[1:6,1, drop = FALSE] #Column 2:3 in this case contain the -log(qValues) 

#The first column of our excel file containing the name of the pathways has the header "KEGG_repressed". Here we turn those names into the rownames 
KEGG.iDAKO.rep <- KEGG.iDAKO[7:10,2, drop = FALSE] #Column 2:3 inpheatmap(KEGG, scale="none", cluster_rows = FALSE, cluster_cols = FALSE,  gaps_row = NULL, fontsize_row = 15, face = "boldpheatmap(KEGG, scale="none", cluster_rows = FALSE, cluster_cols = FALSE,  gaps_row = NULL, fontsize_row = 15, face = "bold


#Choosing the colors. 
my_palette <- colorRampPalette(c("white", "darkred"))(1090)
col<- my_palette
pheatmap(KEGG.iDAKO, scale="none", cluster_rows = FALSE, cluster_cols = FALSE,  gaps_row = NULL, color = col)

#Choosing the colors. 
light_red <- rgb(1, 0.8, 0.8)
my_palette <- colorRampPalette(c(light_red,"darkred"))(1090)
col<- my_palette
pheatmap(KEGG.iDAKO.maint, scale="none", cluster_rows = FALSE, cluster_cols = FALSE,  gaps_row = NULL, fontsize_row = 15, face = "bold", width = 15, color = col, border_color = "black")

#Choosing the colors.
my_palette <- colorRampPalette(c( "lightgreen","darkgreen"))(1090)
col<- my_palette
pheatmap(KEGG.iDAKO.rep, scale="none", cluster_rows = FALSE, cluster_cols = FALSE,  gaps_row = NULL, fontsize_row = 20, width = 5, color = col, border_color = "grey")




# ============================
#  KEGG Analysis Heat map 3D
# ============================

KEGG_analysis.3D <- read_excel(paste0(PATH,"KeggAnalysis_3D.xlsx"))
KEGG.3D <- as.data.frame(KEGG_analysis.3D)
rownames(KEGG.3D) <- KEGG.3D$Pathway 

#The first column of our excel file containing the name of the pathways has the header "KEGG_repressed". Here we turn those names into the rownames 
KEGG.3D <- KEGG.3D[,2:3] #Column 2:3 in this case contain the -log(qValues) 

#The first column of our excel file containing the name of the pathways has the header "KEGG_repressed". Here we turn those names into the rownames 
KEGG.3D.maint <- KEGG.3D[1:8,1, drop = FALSE] #Column 2:3 in this case contain the -log(qValues) 

#The first column of our excel file containing the name of the pathways has the header "KEGG_repressed". Here we turn those names into the rownames 
KEGG.3D.rep <- KEGG.3D[9:12,2, drop = FALSE] #Column 2:3 in this case contain the -log(qValues) 

#Saturating it a bit. So all genes wich a scaled rlog value greater than 2.5 should just be viewed a 2.5.
KEGG.3D.rep[KEGG.3D.rep > 5] <- 10

#Everything lower than -2 should just be called -2
KEGG.3D.rep[KEGG.3D.rep< 3] <- 0

#Choosing the colors.
my_palette <- colorRampPalette(c("white", "dark red"))(1090)
col<- my_palette
pheatmap(KEGG.3D, scale="none", cluster_rows = FALSE, cluster_cols = FALSE,  gaps_row = NULL, color = col)

#Choosing the colors.
my_palette <- colorRampPalette(c("white", "indianred", "darkred"))(50)
col<- my_palette
pheatmap(KEGG.3D.maint, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = NULL, fontsize_row = 15, width = 10, color = col)

#Choosing the colors. 
my_palette <- colorRampPalette(c("lightgreen", "darkgreen"))(1090)
col<- my_palette
pheatmap(KEGG.3D.rep, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = NULL, fontsize_row = 15, color = col)




# ====================
# KEGG DOTPLOT ANalysis
# ====================

# KEGG DOT 3D
KeggDotMaint.3D <- LipMaintGenes.3D$symbol # Replace with your actual gene IDs
KeggDotRepr.3D <- LipRepGenes.3D$symbol   # Replace with your actual gene IDs

KeggDot3D <- list(
  mapIds(org.Mm.eg.db, KeggDotMaint.3D, 'ENTREZID', 'SYMBOL'),
  mapIds(org.Mm.eg.db, KeggDotRepr.3D, 'ENTREZID', 'SYMBOL')
)

names(KeggDot3D) <- c("basal lipolysis maintained", "basal lipolysis repressed")

KeggDot3D <- compareCluster(KeggDot3D, fun="enrichKEGG",organism="mmu", pvalueCutoff=0.1)
KeggDot3D@compareClusterResult$Description <- sub("Mus musculus \\(house mouse\\)", "", KeggDot3D@compareClusterResult$Description)

# KEGG DOT iWAT
KeggDotMaint.iWA <- rownames(LipMaintGenes.iWA) # Replace with your actual gene IDs
KeggDotRepr.iWA <- rownames(LipRepGenes.iWA)   # Replace with your actual gene IDs

KeggDotiWA <- list(
  mapIds(org.Mm.eg.db, KeggDotMaint.iWA, 'ENTREZID', 'SYMBOL'),
  mapIds(org.Mm.eg.db, KeggDotRepr.iWA, 'ENTREZID', 'SYMBOL')
)

names(KeggDotiWA) <- c("basal lipolysis maintained", "basal lipolysis repressed")

KeggDotiWA <- compareCluster(KeggDotiWA, fun="enrichKEGG",organism="mmu", pvalueCutoff=0.1)
KeggDotiWA@compareClusterResult$Description <- sub("Mus musculus \\(house mouse\\)", "", KeggDotiWA@compareClusterResult$Description)

# KEGG DOT Combined
KeggDotMaint.3D <- rownames(LipMaintGenes.3D) # Replace with your actual gene IDs
KeggDotMaint.iWA <- rownames(LipMaintGenes.iWA)   # Replace with your actual gene IDs

KeggDotCombinedMaintained <- list(
  mapIds(org.Mm.eg.db, KeggDotMaint.3D, 'ENTREZID', 'SYMBOL'),
  mapIds(org.Mm.eg.db, KeggDotMaint.iWA, 'ENTREZID', 'SYMBOL')
)

names(KeggDotCombinedMaintained) <- c("3D Maintained", "iWAT Maintained")

KeggDotCombinedMaintained <- compareCluster(KeggDotCombinedMaintained, fun="enrichKEGG",organism="mmu", pvalueCutoff=0.1)
KeggDotCombinedMaintained@compareClusterResult$Description <- sub("Mus musculus \\(house mouse\\)", "", KeggDotCombinedMaintained@compareClusterResult$Description)

# Kegg PLow Combined
union(rownames(SignificantMaint.both.3D),rownames(SignificantMaint.both.iWA))

#KeggDotCombPlow <- union(rownames(LipPLow.comb.3D),rownames(LipPLow.comb.iWA))
KeggDotCombPlow <- union(rownames(SignificantMaint.both.3D),rownames(SignificantMaint.both.iWA))
KeggDotCombPlow <- list(
  mapIds(org.Mm.eg.db, KeggDotCombPlow, 'ENTREZID', 'SYMBOL')
)

names(KeggDotCombPlow) <- c("Maintained")

KeggDotCombPlow <- compareCluster(KeggDotCombPlow, fun="enrichKEGG",organism="mmu", pvalueCutoff=0.1)
KeggDotCombPlow@compareClusterResult$Description <- sub("Mus musculus \\(house mouse\\)", "", KeggDotCombPlow@compareClusterResult$Description)



# =========================
# KEgg DOtplot Plots
# =========================

# PLOTTING OF KEGG 3D Maintained (Color adjusted)
png(file=paste0(PATH,"img/keggdot_maintained.png"), width = 1100, height = 1100) ## to automatically save the file (when ended by dev.off() later), you can leave this out.

maintained_pvals <- KeggDot3D@compareClusterResult[KeggDot3D@compareClusterResult$Cluster == "basal lipolysis maintained", "p.adjust"]
min10_pvals <- sort(maintained_pvals)[1:10]
mx <- max(min10_pvals)

dotplot(KeggDot3D,
        x = "Cluster", color = "p.adjust", showCategory = 10, split = NULL, font.size = 20,
        title = "KEGG pathways - basal lipolysis maintained vs repressed",
        by = "p.adjust", size = "count", includeAll = TRUE, label_format = 100) +
  scale_color_gradient(low = "darkred", high = "pink", limits = c(0, mx))+
  theme(axis.text.x = element_text(size = 14),  # Adjust the size of x-axis text
        axis.text.y = element_text(size = 14))  # Adjust the size of y-axis text
dev.off()


# PLOTTING OF KEGG 3D Repressed (Color adjusted)
png(file=paste0(PATH,"img/keggdot3D_repressed.png"), width = 1100, height = 1100) ## to automatically save the file (when ended by dev.off() later), you can leave this out.

repressed_pvals <- KeggDot3D@compareClusterResult[KeggDot3D@compareClusterResult$Cluster == "basal lipolysis repressed", "p.adjust"]
min10_pvals <- sort(repressed_pvals)[1:10]
mx <- max(min10_pvals)

dotplot(KeggDot3D,
        x = "Cluster", color = "p.adjust", showCategory = 10, split = NULL, font.size = 20,
        title = "KEGG pathways - basal lipolysis maintained vs repressed",
        by = "p.adjust", size = "count", includeAll = TRUE, label_format = 100) +
  scale_color_gradient(low = "darkgreen", high = "lightgreen", limits = c(0, mx))+
  theme(axis.text.x = element_text(size = 14),  # Adjust the size of x-axis text
        axis.text.y = element_text(size = 14))  # Adjust the size of y-axis text
dev.off()


# PLOTTING OF KEGG iDAKO Maintained (color adjusted)
png(file=paste0(PATH,"img/keggdotiDAKO_maintained.png"), width = 1100, height = 900) ## to automatically save the file (when ended by dev.off() later), you can leave this out.

maintained_pvals <- KeggDotiWA@compareClusterResult[KeggDotiWA@compareClusterResult$Cluster == "basal lipolysis maintained", "p.adjust"]
min10_pvals <- sort(maintained_pvals)[1:10]
mx <- max(min10_pvals)

dotplot(KeggDotiWA,
        x = "Cluster", color = "p.adjust", showCategory = 15, split = NULL, font.size = 12,
        title = "KEGG pathways - basal lipolysis maintained vs repressed",
        by = "p.adjust", size = "count", includeAll = TRUE, label_format = 100) +
  scale_color_gradient(low = "darkred", high = "pink", limits = c(0, mx)) +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size of x-axis text
        axis.text.y = element_text(size = 14))  # Adjust the size of y-axis text

dev.off()


# PLOTTING OF KEGG iDAKO Repressed (Color adjusted)
png(file=paste0(PATH,"img/keggdotiDAKO_repressed.png"), width = 1100, height = 900) ## to automatically save the file (when ended by dev.off() later), you can leave this out.

repressed_pvals <- KeggDotiWA@compareClusterResult[KeggDotiWA@compareClusterResult$Cluster == "basal lipolysis repressed", "p.adjust"]
min10_pvals <- sort(repressed_pvals)[1:10]
mx <- max(min10_pvals)

dotplot(KeggDotiWA,
        x = "Cluster", color = "p.adjust", showCategory = 10, split = NULL, font.size = 12,
        title = "KEGG pathways - basal lipolysis maintained vs repressed",
        by = "p.adjust", size = "count", includeAll = TRUE, label_format = 100) +
  scale_color_gradient(low = "darkgreen", high = "lightgreen", limits = c(0, mx)) +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size of x-axis text
        axis.text.y = element_text(size = 14))  # Adjust the size of y-axis text

dev.off()

# PLOTTING OF KEGG COMBINED (Just from both sets)
png(file=paste0(PATH,"img/keggdot_comb_maintained.png"), width = 1100, height = 1100) ## to automatically save the file (when ended by dev.off() later), you can leave this out.

maintained_pvals <- KeggDotCombinedMaintained@compareClusterResult[KeggDotCombinedMaintained@compareClusterResult$Cluster == "basal lipolysis maintained", "p.adjust"]
min10_pvals <- sort(maintained_pvals)[1:10]
mx <- max(min10_pvals)

KeggDotCombinedMaintained@compareClusterResult$logval <- -log(KeggDotCombinedMaintained@compareClusterResult$p.adjust)

dotplot(KeggDotCombinedMaintained,
        x = "Cluster", color = "logval", showCategory = 10, split = NULL, font.size = 20,
        title = "KEGG pathways - 3D Maintained vs iDAKO Maintained",
        by = "p.adjust", size = "count", includeAll = TRUE, label_format = 100) +
  scale_color_gradient(low = "pink", high = "darkred")+
  theme(axis.text.x = element_text(size = 14),  # Adjust the size of x-axis text
        axis.text.y = element_text(size = 14))  # Adjust the size of y-axis text
dev.off()


# PLOTTING OF KEGG COMBINED (Selected based off padjust)
png(file=paste0(PATH,"img/keggdot_comb_maintained_plow.png"), width = 1100, height = 1100) ## to automatically save the file (when ended by dev.off() later), you can leave this out.

maintained_pvals <- KeggDotCombPlow@compareClusterResult[KeggDotCombPlow@compareClusterResult$Cluster == "basal lipolysis maintained", "p.adjust"]
min10_pvals <- sort(maintained_pvals)[1:10]
mx <- max(min10_pvals)

KeggDotCombPlow@compareClusterResult$logval <- -log(KeggDotCombPlow@compareClusterResult$p.adjust)

dotplot(KeggDotCombPlow,
        x = "Cluster", color = "logval", showCategory = 10, split = NULL, font.size = 20,
        title = "KEGG pathways - Maintained",
        by = "p.adjust", size = "count", includeAll = TRUE, label_format = 100) +
  scale_color_gradient(low = "pink", high = "darkred")+
  theme(axis.text.x = element_text(size = 14),  # Adjust the size of x-axis text
        axis.text.y = element_text(size = 14))  # Adjust the size of y-axis text
dev.off()




# ======================================
#  BAR PLOT OF ADIPOCYTE SPECIFIC GENES
# ======================================

genes_adi = read_excel(paste0(PATH,"singlenuclei_genelist.xlsx"))

genes_adi_maint.3D = intersect(genes_adi$gene, rownames(LipMaintGenes.3D))
genes_adi_maint.iWA = intersect(genes_adi$gene, rownames(LipMaintGenes.iWA))

genes_adi_rep.3D = intersect(genes_adi$gene, rownames(LipRepGenes.3D))
genes_adi_rep.iWA = intersect(genes_adi$gene, rownames(LipRepGenes.iWA))

genes_adi_const.3D = intersect(genes_adi$gene, rownames(LipConstGenes.3D))
genes_adi_const.iWA = intersect(genes_adi$gene, rownames(LipConstGenes.iWA))

write.xlsx(as.data.frame(genes_adi_maint.3D),file=paste0(PATH,"maintained_adi_3D.xlsx"))
write.xlsx(as.data.frame(genes_adi_maint.iWA),file=paste0(PATH,"maintained_adi_iWA.xlsx"))

length(genes_adi$gene)
length(genes_adi_maint.3D)

r_maint.3D = length(genes_adi_maint.3D)/nrow(genes_adi)
r_maint.iWA = length(genes_adi_maint.iWA)/nrow(genes_adi)

r_rep.3D = length(genes_adi_rep.3D)/nrow(genes_adi)
r_rep.iWA = length(genes_adi_rep.iWA)/nrow(genes_adi)

r_const.3D = length(genes_adi_const.3D)/nrow(genes_adi)
r_const.iWA = length(genes_adi_const.iWA)/nrow(genes_adi)

group <- c("Maintained", "Repressed", "Constant")
ratio <- c(r_maint.3D, r_rep.3D, r_const.3D)
adipocyte_specific.3D <- data.frame(group, ratio)
adipocyte_specific.3D$group <- factor(adipocyte_specific.3D$group, levels = c("Maintained", "Repressed", "Constant"))

group <- c("Maintained", "Repressed", "Constant")
ratio <- c(r_maint.iWA, r_rep.iWA, r_const.iWA)
adipocyte_specific.iWA <- data.frame(group, ratio)
adipocyte_specific.iWA$group <- factor(adipocyte_specific.iWA$group, levels = c("Maintained", "Repressed", "Constant"))

# Define the two vectors (replace with your actual data)
iWAT_3D <- r_maint.3D
iWAT_iDAKO <- r_maint.iWA

# Combine the vectors into a dataframe
adipocyte_specific.combined <- data.frame(
  Value = c(iWAT_3D, iWAT_iDAKO),
  Group = c(rep("iWAT 3D", length(iWAT_3D)), rep("iWAT iDAKO", length(iWAT_iDAKO))),
  Category = factor(rep(1:length(iWAT_3D), 2))
)

# Create the bar plot
ggplot(adipocyte_specific.3D, aes(x = group, y = ratio, fill = group)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "", y = "Ratio of overlap \n(n of of observations/ total number of adipocyte specific genes)") +
  ggtitle("") +
  scale_fill_manual(values = c("Maintained" = "darkred", "Repressed" = "darkgreen", "Constant" = "grey")) +
  theme(
    panel.background = element_rect(fill = "transparent", color = "black"),
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

# Create the bar plot
ggplot(adipocyte_specific.iWA, aes(x = group, y = ratio, fill = group)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "", y = "Ratio of overlap \n(n of of observations/ total number of adipocyte specific genes)") +
  ggtitle("Ratio of Adipocyte Specific Genes (iWAs)") +
  scale_fill_manual(values = c("Maintained" = "darkred", "Repressed" = "darkgreen", "Constant" = "grey")) +
  theme(panel.background = element_rect(fill = "transparent", color = "black")
  )    


# Create the bar plot
ggplot(adipocyte_specific.combined, aes(x = Group, y = Value, fill = Group)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "", y = "Ratio of overlap \n(n of of observations/ total number of adipocyte specific genes)") +
  ggtitle("Ratio of Adipocyte Specific Genes (3D cultured iWAs)") +
  scale_fill_manual(values = c("Maintained" = "darkred", "Repressed" = "darkgreen", "Constant" = "grey")) +
  theme(panel.background = element_rect(fill = "transparent", color = "black")
 )    




# =====================
#  STACKED BAR PLOT (DEGs)
# =====================

group <- c("3D", "3D", "iDAKO", "iDAKO")
category <- c("Repressed", "Maintained", "Repressed", "Maintained")
genes_numbers <- c(nrow(LipRepGenes.3D), nrow(LipMaintGenes.3D), nrow(LipRepGenes.iWA), nrow(LipMaintGenes.iWA))
gene_numbers <- data.frame(group, category, genes_numbers)

# Create the stacked bar plot
options(repr.plot.width = 20, repr.plot.height = 12)

# Calculate the total number of values in each group and category
df_summary <- gene_numbers %>%
  group_by(group, category) %>%
  summarise(total = sum(genes_numbers), .groups = 'drop')

# Create the ggplot2 bar plot
ggplot(gene_numbers, aes(fill = category, y = genes_numbers, x = group)) +
  geom_bar(position = "stack", stat = "identity", width = 0.8) +
  labs(x = "", y = "DEGs", fill = "Category") +
  scale_fill_manual(values = c("Repressed" = "darkgreen", "Maintained" = "darkred")) +
  # Add the total number of values as text on the bars for "Repressed" category
  geom_text(data = df_summary %>% filter(category == "Repressed"),
            aes(label = total, y = total / 2, x = group),
            position = position_stack(vjust = 1), size = 8, fontface= "bold", color = "black") +
  # Add the total number of values as text on the bars for "Maintained" category
  geom_text(data = df_summary %>% filter(category == "Maintained"),
            aes(label = total, y = total / 2, x = group),
            position = position_stack(vjust = 3), size = 8, fontface = "bold", color = "black") +
  theme(
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18, face = "bold"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.margin = margin(t = 0, r = 5, b = 0, l = 5),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),
    panel.border = element_rect(fill = NA, color = "black")
  )


length(intersect(rownames(LipMaintGenes.3D),rownames(LipMaintGenes.iWA)))/length(union(rownames(LipMaintGenes.3D),rownames(LipMaintGenes.iWA)))




# ========================
#  VENN DIAGRAM DATA
# ========================

set.3D_maint <- union(rownames(LipMaintGenes.3D),rownames(LipRepGenes.3D))
set.iWA_maint <- rownames(LipMaintGenes.iWA)

intersect(set.3D_maint,set.iWA_maint)


set.3D_maint <- rownames(LipMaintGenes.3D)
set.iWA_maint <- rownames(LipMaintGenes.iWA)
set_maint_3d_iwa <- list(Set1 = set.3D_maint, Set2 = set.iWA_maint)

set.3D_rep <- rownames(LipRepGenes.3D)
set.iWA_rep <- rownames(LipRepGenes.iWA)
intersect_repr_3d_iwa <- list(Set1 = set.3D_rep, Set2 = set.iWA_rep)




# ========================
#  VENN DIAGRAM iDAKO vs iWAT 3D
# ========================

# Create the Venn diagram

venn.plot <- venn.diagram(
  x = set_maint_3d_iwa,
  category.names = c("iWAT 3D", "iWAT iDAKO"),
  fill = c("darkred", "darkred"),  # Change colors
  filename = NULL,
  cat.cex = 2,  # Adjust the category label size if necessary
  cex = 2,  # Adjust the numbers size if necessary
  cat.default.pos = "text",
  cat.pos = c(-5, 5),
  cat.dist = c(0.25, 0.25),
  cat.fontface = "bold",
)

# Plot the Venn diagram
grid.newpage()
grid.draw(venn.plot)

# Print genes found in the overlap
overlap_combined <- intersect(set_maint_3d_iwa[[1]], set_maint_3d_iwa[[2]])
print(overlap_combined)

intersect(genes_adi$gene,overlap_combined)

# Venn diagram repressed overlap 
venn.plot <- venn.diagram(
  x = intersect_repr_3d_iwa,
  category.names = c("iWAT 3D", "iWAT iDAKO"),
  fill = c("darkgreen", "darkgreen"),  # Change colors
  filename = NULL,
  cat.cex = 1,  # Adjust the category label size if necessary
  cex = 1.5,  # Adjust the numbers size if necessary
  cat.default.pos = "text",
  cat.pos = c(-5, -10),
  cat.dist = c(0.22, 0.22),
  cat.fontface = "bold",
)

# Plot the Venn diagram
grid.newpage()
grid.draw(venn.plot)

# Print overlap
overlap_rep<- intersect(intersect_repr_3d_iwa[[1]], intersect_repr_3d_iwa[[2]])
print(overlap_rep)

intersect(genes_adi$gene,overlap_rep)
