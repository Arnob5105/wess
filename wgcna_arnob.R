rm(list = ls())
setwd("F:/WES/CGGA/")

#BiocManager::install("WGCNA")

library(WGCNA)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)
options(stringsAsFactors = FALSE);


library(preprocessCore)

## loading data

data <- read.csv('F:/WES/CGGA/gbm_only.csv',row.names = 1,header = TRUE, sep = ",")
which(duplicated(data$gene_name) == TRUE)
#data <- normalize.quantiles(as.matrix(data))
dev.new(width=3+ncol(data)/6, height=5)
par(mar=c(7,4,2,1))
#title <- paste ("GSE50161", "/", annotation(data), sep ="")
boxplot(data, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
dev.off()

#Good sample genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

## Dtecting outliers by clustering

htree <- hclust(dist(t(data)), method = "average")
plot(htree)

pdf(file = "cluster.pdf", width = 40, height = 9);
par(cex = 1.3);
par(mar = c(0,4,2,0))
plot(htree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()


# Detting outliers by PCA

pca <- prcomp(t(data))
pca.dat <- pca $x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digit =2)
pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2))+
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0("PC1: ", pca.var.percent[1],'%'),
       y = paste0("PC2: ", pca.var.percent[2], '%'))


# pheno data

coldata <- read.csv('F:/WES/CGGA/clinical_new.csv',row.names = 1,header = TRUE, sep = ",")
all(row.names(coldata)) %in% colnames(data)
all(row.names(coldata)) == colnames(data)

# DEseq2 matrix
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~1)

dds75 <- dds[rowSums(counts(dds)> 15) > 50,]
nrow(dds75)
write.csv(rownames(dds75@assays@data@listData[["counts"]]), 'lowcount_new.csv')

# normalization
dds_norm <- vst(dds75)
dds_count <- assay (dds_norm) 



#### CEMItools
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install( "CEMiTool")
BiocManager::install("clusterProfiler")
install.packages("CEMiTool")
install.packages("tidyverse", type="binary")
options(BioC_mirror = "http://bioconductor.org")

library("ggpmisc")
library("CEMiTool")

cem <- cemitool(data)
main <- as.data.frame(cem@selected_genes)
write.csv(main, 'cem_combind_HD.csv')



## soft power

power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(t(data),
                         powerVector = power,
                         networkType = "unsigned",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.85, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)

## soft power

powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(t(dds_count), powerVector = powers, networkType = "signed", verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
#pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.8;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.83,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

power=sft$powerEstimate #4
soft_power <- 9
# Option 1: automatic
temp_cor <- cor
cor <- WGCNA::cor
net = blockwiseModules(t(data), power = soft_power,
                       maxBlockSize = 6500,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.20,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)

cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- net$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(net$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(net$dendrograms[[1]], cbind(net$unmergedColors, net$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# option 2 (step by step) 
power <- 12
power = power
adjacency = adjacency(t(dds_count), power = power)
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM

# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
#pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(t(dds_count), colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Merge close modules
MEDissThres=0.10
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(t(dds_count), dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
#pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()


# clinical data
colData <- coldata

colData = read.csv(file = "meta.csv",header = T, row.names = 1, sep = ',',)

colData$class <- factor(colData$class, levels = c('control','Case'))

trait <- binarizeCategoricalColumns(colData$class,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)
colnames(trait) <- 'condition'

traits <- cbind(colData[-c(1)],trait)
# Define numbers of genes and samples
nSamples <- nrow(t(dds_count))
nGenes <- ncol(t(dds_count))


module.trait.corr <- cor(mergedMEs, coldata, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

heatmap.data <- merge(mergedMEs, coldata, by = 'row.names')
rownames(heatmap.data) <- heatmap.data[,1]
heatmap.data <- heatmap.data[,-1 ] 


#heatmap.data <- heatmap.data %>% 
  #column_to_rownames(var = 'Row.names')




CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[11],
             y = names(heatmap.data)[1:10],
             col = c("blue1", "skyblue", "white", "pink", "red"))


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(module.trait.corr, 2), "\n(",
                    signif(module.trait.corr.pvals, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.corr)
pdf("module_trait_glioblastoma.pdf", width = 10, height = 15)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.corr,
               xLabels = colnames(coldata),
               yLabels = colnames(mergedMEs),
               ySymbols = colnames(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


module.gene.mapping <- as.data.frame(mergedColors)
module.gene.mapping <- cbind(dds_count[-c(1)],module.gene.mapping)
gene <-module.gene.mapping %>% 
      filter(`mergedColors` == 'MElightgreen') %>% 
      rownames()
#write.csv(gene, 'green.csv')


cancer = as.data.frame(coldata);
names(cancer) = "glioblastoma"


modNames = substring(names(mergedMEs), 3)

#MET = orderMEs(cbind(MEs, Verru))

geneModuleMembership = as.data.frame(cor(t(dds_count), mergedMEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(t(dds_count), cancer, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(coldata), sep="");
names(GSPvalue) = paste("p.GS.", names(coldata), sep="");


module = "MElightgreen"
# Rename to moduleColors
moduleColors = mergedColors
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for glioblastoma",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "purple")
abline(h=0.7,col = "black")
abline(v=0.9,col = "black")
dev.off()

library(dplyr)
#mm_threshold <- 0.85
#gs_threshold <- 0.85
module <-as.data.frame( geneModuleMembership[moduleGenes, column])
significance <- as.data.frame(geneTraitSignificance[moduleGenes, 1])
colnames(significance) <- 'sig'
colnames(module) <- 'MMgreen'


binding <- cbind(gene,module,significance)
rownames(binding) <- gene

filter<-  binding %>% filter(  MMgreen > 0.9 & sig < -0.7)

print(rownames(filter))
view(filter)
# Filter genes based on module membership and gene significance

# Print the filtered genes
write.csv(rownames(filter), 'cyan_gs_mm_gene.csv')
