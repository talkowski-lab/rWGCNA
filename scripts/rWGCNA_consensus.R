# This R script performs one-step consensus blockwise WGCNA module
#  detection on a user-defined number of data sets resampled without 
#  replacement from a larger user-supplied dataset, at a sample size of 0.66 of
#  the larger data set.
# 
# This R script is called by submit_rWGCNA_consensus.sh, which helps
#  set up the R environment, set a working directory, and pass
#  on essential user-defined arguments (see Assign arguments to variables section).
#  Additional hardwired parameters can be set by user (see Define hardwired parameters).
#
# by Heather Wick hwick@broadinstitute.org heather.c.wick@gmail.com
#
# Input:
# 7 arguments:
#	1) path/to/working_directory
#	2) project_name
#	3) number representing number of times data set is to be resampled (integer)
#	4) input_dataset.txt, log2-transformed gene expression data set with sample_name columns and gene rows
#	5) trait_data.txt, trait data with $Samples column corresponding to column names of input_dataset.txt
#	6) maximum block size (integer), preferably slightly larger than the number of genes in your full dataset, if computational resources allow; may break plotting otherwise
#	7) boolean, use default soft power threshold? If T, default value is 12; otherwise will calculate soft power threshold
#
# Output:
# Plots:
#   SampleClustering.pdf #basic QC to see how samples cluster
#   SoftPowerThresh_signed.pdf #soft power thresholds for all re-sampled datasets, if calculated
#   ConsensusDendrogram.pdf #module dendrogram
# Objects/Environment:
#   multiExpr.RData # object of re-sampled data sets
#   Consensus-dataInput.RData #multiple objects: multiExpr, nGenes, Traits, nSamples, setLabels, shortLabels, exprSize
#   Consensus-NetworkConstruction.RData # multiple objects: consMEs, moduleLabels, moduleColors, consTree
#   consMEsOnDatExpr.Rdata #module eigengenes calculated on datExpr
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

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);###
allowWGCNAThreads() 
options(warn=1)

# Assign arguments to variables
args<-commandArgs(TRUE)
workingDir <- args[1]
proj_name <- args[2]
nRepTOM <- as.integer(args[3]) #num times data set is sampled
inData <- args[4]
traitData <- args[5] #trait data
maxBlockSize <-args[6] #preferably as large as the number of genes in your full dataset, if computational resources allow
defaultPower <-args[7] #bool. T=use default soft power threshold, which is 12; F=calculate soft power threshold on individual sample iterations and use median 

# Define hardwired parameters:
replace<-F #do not replace samples within bootstrap sets
fraction<-0.66 #fraction of full data set used in each sampling iteration; this is also the fraction that must be in agreement to stay in consensus module instead of being reassigned to grey
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
print(paste0("Number of sampling iterations: ", nRepTOM))
print(paste0("Fraction of full data set used in each sample iteration: ", fraction))
print(paste0("Sample with replacement?: ", replace))
print(paste0("Default soft power threshold used?: ", defaultPower))
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

# Define functions:
#  bootstrap function from Perslab rWGCNA seurat project, https://github.com/perslab/wgcna-toolbox, see author info embedded below:
bootstrap <- function(datExpr,
                      nPermutations,
                      replace,
                      fraction,
                      randomSeed)
  # @Usage: Resample a dataset.
  # @args:
  #       datExpr: Dataset with samples in the rows, is coerced to matrix
  #       nPermutations: runs
  #       replace: sample with replacement?
  #       fraction: sample what fraction of the total each iteration?
  #       randomSeed: initial random seed
  # @return:
  #       result: list of resampled datasets in matrix format. If replace = T, first item is the unpermuted dataset.
  # @author: Jonatan Thompson jjt3f2188@gmail.com
  # @date: 180222

{
  startRunIndex = 1
  endRunIndex = if (replace==T) nPermutations+1 else nPermutations 
  result = vector("list", length= if (replace==T) nPermutations+1 else nPermutations);
  nSamples = nrow(datExpr);
  nFeatures = ncol(datExpr);
  
  try(
    for (run in startRunIndex:endRunIndex)
      
    {
      
      #set.seed(randomSeed + 2*run + 1);
      
      if (run == startRunIndex & replace == T) {
        useSamples = c(1:nSamples) # The first run just returns the full dataset
      } else if (run>startRunIndex | replace == F) 
      {  useSamples = sample(nSamples, as.integer(nSamples * fraction), replace = replace)
      } 
      
      
      samExpr = as.matrix(datExpr[useSamples, ]);
      
      result[[run]]$data <- samExpr
    })
  
  return(result)
}

##################
# Begin analysis #
##################

# Create multi-expression set using bootstrap function
multiExpr <- bootstrap(datExpr=datExpr, 
                       nPermutations = nRepTOM,
                       replace = replace,
                       fraction = fraction,
                       randomSeed = randomSeed)

str(multiExpr)
if (checkSets(multiExpr)$structureOK) {
  print("multiExpr created succesffully; saving...")
} else {
  print("multiExpr not created")
}
Sys.time()

# Save multiExpr
saveRDS(multiExpr, "multiExpr.RData")

# Quality control
print("Performing quality control...")
#get number of sets
nSets = checkSets(multiExpr)$nSets
nSets
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
# In this case, they are bootstraps so simple names will suffice
setLabels<-""
shortLabels<-""
for (i in 1:nSets) {
  setLabels[i]<- paste0("bs",i)
  shortLabels[i]<- paste0("bs",i)
}

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)
exprSize

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

# If the last statement returns TRUE, all genes and samples have passed the cuts.
# If it returns FALSE, the following code removes the offending samples and genes:
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;


# Visualize samples
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

# Save plot of sample clustering;
#note: include more parameters in name?
pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();
print("...quality control done, see Plots/SampleClustering.pdf")
print("Retrieving trait data for multiExpr")

# Load clinical trait data, match the samples for which they were measured to the expression samples
# We will save a multiMeta object that has the traits for each Expr set within the multiExpr

#allTraits<-read.table("traitdata_for_github.txt", sep="\t", header=T) 
allTraits<-read.table(traitData, sep="\t", header=T)

# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$Samples);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage()

print("Trait data gathered for multiExpr; Saving data...")

save(multiExpr, nGenes, Traits, nSamples, setLabels, shortLabels, exprSize, 
     file = "Consensus-dataInput.RData")

print("Soft power thresholding...")
Sys.time()
if (defaultPower) {
  power <-12
  #power <-18
  print(paste0("Using default soft power threshold: ", power))
} else {
  print("Calculating power thresholds...")
  powers  =  c(c(2:30))
  # Initialize a list to hold the results of scale-free analysis
  powerTables = vector(mode = "list", length = nSets);
  # Call the network topology analysis function for each set in turn
  for (set in 1:nSets)
    powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, networkType = netType,
                                                       verbose = 2)[[2]]);
  # Plot the results:
  colors = palette(rainbow(nSets))  
  # Will plot these columns of the returned scale free analysis tables
  plotCols = c(2,5,6,7)
  colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
               "Max connectivity");
  # Get the minima and maxima of the plotted points
  ylim = matrix(NA, nrow = 2, ncol = 4);
  for (set in 1:nSets)
  {
    for (col in 1:length(plotCols))
    {
      ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
      ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    }
  }
  
  pdf("Plots/SoftPowerThresh_signed.pdf")
  par(mfcol = c(2,2));
  par(mar = c(4.2, 4.2 , 2.2, 0.5))
  cex1 = 0.7;
  for (col in 1:length(plotCols)) for (set in 1:nSets)
  {
    if (set==1)
    {
      plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
           xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
           main = colNames[col]);
      addGrid();
    }
    if (col==1)
    {
      text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
           labels=powers,cex=cex1,col=colors[set]);
    } else
      text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
           labels=powers,cex=cex1,col=colors[set]);
    if (col==1)
    {
      legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
    } else
      legend("topright", legend = setLabels, col = colors, pch = 20) ;
  }
  dev.off()
  print("Soft power thresholds for multiExpr calculated; see Plots/SoftPowerThresh_signed.pdf")
  powercon<-""
  for( i in 1:length(powerTables)) {
    powerTables[[i]]$data$SFT.R.sq
    dataDF1 <- powerTables[[i]]$data
    #check for NA values (SFT.R.sq does not reach 0.8), set to max
    if(is.na(dataDF1[dataDF1$SFT.R.sq > 0.8,]$Power[1])){
      powercon[i] <-31
    } else {
      powercon[i] <- dataDF1[dataDF1$SFT.R.sq > 0.8,]$Power[1]
    }
  }
  power<-floor(median(as.integer(powercon)))
  print(paste0("Using median soft power threshold: ", power))
  collectGarbage();
}
Sys.time()
softPower = power

########################################
# One-step blockwise module detection: #
########################################

print("Beginning blockwise consensus module detection")
Sys.time()

# Note: if your data may have outliers, consider setting corType to biweight
# instead of the default pearson

net = blockwiseConsensusModules(
  multiExpr, 
  maxBlockSize = maxBlockSize,
  randomSeed = randomSeed,
  power = softPower, 
  minModuleSize = minModuleSize, 
  networkType = netType,
  TOMType = netType,
  saveIndividualTOMs = F,
  calibrationQuantile = scaleP, 
  deepSplit = 2,
  pamRespectsDendro = FALSE, 
  mergeCutHeight = 0.25, #default is 0.15
  nThreads = 0, #0 is the default; if it can be determined dynamically (should be able to be?) it will use that many threads; otherwise will use 2
  numericLabels = TRUE,
  minKMEtoStay = 0.2, #0 turns off min KMEtoStay; default is 0.2
  saveTOMs = F, #TOMS take up a lot of space
  verbose = 5)

Sys.time()
print("Blockwise consensus modules done. Extracting modules....")
consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)

table(moduleColors)
paste0("There are ",length(unique(moduleColors))," modules")
consTree = net$dendrograms[[1]];

# Plot the gene dendrogram and the corresponding module colors:
pdf(file = "Plots/ConsensusDendrogram.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()

print("Modules extracted. See Plots/ConsensusDendogram.pdf. Saving final objects and environment...")

# Calculate module eigengenes on whole data set
consMEsOnDatExpr<-moduleEigengenes(datExpr, colors=moduleColors)
saveRDS(consMEsOnDatExpr, file="consMEsOnDatExpr.RData")

save(consMEs, moduleLabels, moduleColors, consTree, file = "Consensus-NetworkConstruction.RData")
save.image("Final_image.RData")
print("Success! Analysis complete!")
Sys.time()
