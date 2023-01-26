# This R script performs one-step individual blockwise WGCNA module
# 
# This R script is called by submit_rWGCNA_individual.sh, which helps
#  set up the R environment, set a working directory, and pass
#  on essential user-defined arguments (see Assign arguments to variables section).
#  Additional hardwired parameters can be set by user (see Define hardwired parameters).
#
# Input:
# 5 arguments:
#	1) path/to/working_directory
#	2) project_name
#	3) input_dataset.txt, log2-transformed gene expression data set with sample_name columns and gene rows
#	4) trait_data.txt, trait data with $Samples column corresponding to column names of input_dataset.txt
#	5) maximum block size (integer), preferably slightly larger than the number of genes in your full dataset, if computational resources allow; may break plotting otherwise
#
# Output:
# Plots:
#   SampleClustering.pdf #basic QC to see how samples cluster
#   SoftPowerThresh_signed.pdf #soft power thresholds for all re-sampled datasets, if calculated
#   ConsensusDendrogram.pdf #module dendrogram
# Objects/Environment:
#   datExpr.RData # dataset
#   Individual-dataInput.RData #multiple objects: datExpr, nGenes, Traits, nSamples, setLabels, shortLabels, exprSize
#   Individual-NetworkConstruction.RData # multiple objects: MEs, moduleLabels, moduleColors, sampleTree
#   MEs.Rdata #module eigengenes calculated on datExpr
#   Final_image.RData # R image containing all objects and environment for additional analysis
# If stdout is saved in job submission (recommended), output from this R script
#  will be saved for reference 
#
# Sys.time() is called at various steps to assist in any debugging
#
# Finish setting up environnment by specifying libPaths and loading libraries
.libPaths(c("/PHShome/ry077/R/x86_64-pc-linux-gnu-library/4.0","/PHShome/hcw11/R/x86_64-pc-linux-gnu-library/4.0","/apps/source/R/R-4.0.2-withX11/lib64/R/library"))

# Load necessary libraries
library(BiocParallel)
library(WGCNA);
Sys.time()

#2) create a directory for output, as well as Plots and Jobs folders
#mkdir working/directory
#cd working/directory
#mkdir Plots Jobs

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);###
allowWGCNAThreads() 
options(warn=1)

# Assign arguments to variables
args<-commandArgs(TRUE)
workingDir <- args[1]
proj_name <- args[2]
inData <- args[3] 
traitData <- args[4] #trait data
maxBlockSize <-as.integer(args[5]) #preferably as large as the number of genes in your full dataset, if computational resources allow

# Define hardwired parameters:
randomSeed<-12345 #set seed for reproducibility
netType<-"signed" #set network type to signed
scaleP = 0.95 #reference percentile for scaling multiExpr #default is 0.95
minModuleSize = 30; #smallest module size

# Print arguments/parameters to stdout for future reference
print(paste0("Seed: ", randomSeed))
print(paste0("Working directory: ", workingDir))
print(paste0("Project name: ", proj_name))
print(paste0("Input data: ", inData))
print(paste0("Trait data: ", traitData))
print(paste0("Maximum block size: ", maxBlockSize))
print(paste0("Network type: ", netType))
print(paste0("Minimum module size: ", minModuleSize))
print(paste0("Reference percentile for scaling multiExpr: ", scaleP)) 

# Set working directory and seed
setwd(workingDir)
set.seed(randomSeed)

# Load/transpose log2-transformed input data
cleany<-read.table(inData)
datExpr<-t(cleany)

names(datExpr) = rownames(cleany); #gene names
rownames(datExpr) = names(cleany); #sample names

##################
# Begin analysis #
##################

# Quality control
print("Performing quality control...")

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

#If the last statement returns TRUE, all genes and samples have passed the cuts.
#If it returns FALSE, the following code removes the offending samples and genes:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Save plot of sample clustering;
sampleTree = hclust(dist(datExpr), method = "average");
pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
print("...quality control done, see Plots/SampleClustering.pdf")

# Load clinical trait data, match the samples for which they were measured to the expression samples
print("Retrieving trait data")
meta<-read.table("/data/talkowski/hwick/projects/tutorials/WGCNA/WGCNA_from_laptop/metaData_mRNA.tab", sep="\t", header=T) #use all wildtypes
allTraits<-read.table(traitData, sep="\t", header=T)
setSamples = rownames(datExpr);
traitRows = match(setSamples, allTraits$Samples);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();

print("Trait data gathered for datExpr; Saving data...")
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
save(datExpr, datTraits, nGenes, nSamples,
     file = "Individual-dataInput.RData")

print("Soft power thresholding...")
Sys.time()

# Choose a set of soft-thresholding powers
powers  =  c(c(2:30))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, 
                        powerVector = powers, 
                        networkType = netType,
                        verbose = 5)
# Plot the results:

pdf("Plots/powerThresh_signed.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
dataDF <- sft$fitIndices
power <- dataDF[dataDF$SFT.R.sq > 0.8,]$Power[1]
print(paste0("Using soft power threshold: ", power))

Sys.time()
softPower<-power

########################################
# One-step blockwise module detection: #
########################################

print("Beginning blockwise consensus module detection")
Sys.time()
net = blockwiseModules(datExpr, 
                       maxBlockSize = maxBlockSize,
                       randomSeed = randomSeed,
                       power = softPower,
                       minModuleSize = minModuleSize,
                       networkType = netType,
                       TOMType = netType,
                       minKMEtoStay = 0.2, #consistent with consensus analysis
                       mergeCutHeight = 0.25,
                       nThreads = 0,
                       numericLabels = TRUE, #consider changing to get color labels...
                       pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       verbose = 3)

Sys.time()
print("Blockwise modules done. Extracting modules....")

# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
#Save MEs
saveRDS(MEs, file="MEs.RData")
geneTree = net$dendrograms[[1]];

# Plot the dendrogram and the module colors underneath
pdf(file = "Plots/IndividualDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

table(moduleColors)
paste0("There are ",length(unique(moduleColors))," modules")

print("Modules extracted. See Plots/IndividualDendrogram-auto.pdf. Saving final objects and environment...")

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Individual-NetworkConstruction-auto.RData")
save.image("Final_image.RData")
print("Success! Analysis complete!")
Sys.time()
