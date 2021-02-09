library(DESeq2)
library(edgeR)
library(DEFormats)
library(RNOmni)
library(recount3)

### This script accesses Geuvadis gene expression data through recount3
### and normalizes it based on the following protocol from GTEx:
### - TMM normalization for all samples
### - Filter to keep genes in which at least 20% of samples must have
###   a TPM value > 0.1 and at least 20% of samples must have at least 6 reads
### - Rank normalization to fit the data into a normal distribution and make it comparable to other datasets

############################################################################


# load Geuvadis data as a ranged summarized experiment object from recount3
human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  project == "ERP001942" & project_type == "data_sources"
)
rse_gene_ERP001942 <- create_rse(proj_info) # creates rse 


### TMM Normalization 

names(colData(rse_gene_ERP001942)) = "group"
# creates a DGEList object with $counts, $genes, and $samples
dgeFullData = DGEList(rse_gene_ERP001942)
# normalizes $counts by TMM
TMMFullData <- calcNormFactors(dgeFullData, method="TMM")
# creates a matrix out of the TMM normalized counts to be cleaned
TMMCounts <- as.matrix (TMMFullData$counts)
# cleans the matrix based on the following criteria: ≥6 reads in at least 20% of samples
countsCleaned <- TMMCounts[rowSums(TMMCounts >= 6) > (ncol(TMMCounts)* .2),]


### TPM Calculation

#creates a function for calculating TPM values
calc_tpm <- function(x, gene.length) {
  x <- as.matrix(x)
  len.norm.lib.size <- colSums(x / gene.length)
  return((t(t(x) / len.norm.lib.size) * 1e06) / gene.length)
}
# creates a matrix with calculated TPM values for each sample from TMM normalized counts and $genes$width
rawTPMVals <- calc_tpm(TMMFullData, gene.length = TMMFullData$genes$width)
# cleans the matrix based on the following criteria:>0.1 TPM in at least 20% of samples
cleanedTPMVals <- rawTPMVals[rowSums(rawTPMVals > 0.1) > (ncol(rawTPMVals)* .2),]
# converts cleaned TMM normalized counts matrix to a dataframe for intersection
cleanCountsDf <- as.data.frame(countsCleaned)
# converts cleaned TPM matrix to dataframe for intersection
cleanTPMDf <- as.data.frame(cleanedTPMVals)


### Intersection between cleaned TMM normalized counts and cleaned TPM values

# intersects the two dataframes by row names and creates new dataframe with genes that have
# >0.1 TPM in at least 20% of samples and ≥6 reads in at least 20% of samples
cleanedIntersection <- cbind( cleanCountsDf[ intersect(rownames(cleanCountsDf), rownames(cleanTPMDf)), ])


### Rank Normalization/ Inverse Normal Transform

# turns intersected cleaned data frame into a matrix
intersectedMatrix <- data.matrix(cleanedIntersection)
# applies rank normalization/ inverse normal transform for each gene, which is each row of the matrix
rankNormMatrix <- t(apply(intersectedMatrix, 1, RankNorm))


### Reassigning phenotype colnames to match genotype sample IDs

# reads txt file with key and converts to a dataframe
myLookupDf <- read.csv(file="filereport_read_run_PRJEB3366_tsv.txt", sep="\t", header=TRUE)
# extracts genotype sample ID names
geuvNames <- myLookupDf[c(2)]
# removes "GEUV:"
actualGeuvNames <- gsub("^.*?:","",geuvNames$sample_alias)
# turns dataframe to matrix to create a space for column name
actualGeuvNames <- data.matrix(actualGeuvNames)
# manually adds right column name
colnames(actualGeuvNames) <- c("sample_alias")
# replaces genotype sample IDs back into dataframe with phenotype IDs
myLookupDf$sample_alias <- actualGeuvNames
# creates phenotype data sample IDs that match genotype data
updatedPhenotypeDataCols <-  myLookupDf$sample_alias[match(colnames(rankNormMatrix),myLookupDf$run_accession)]
# reassigns updated colnames
colnames(rankNormMatrix) <- updatedPhenotypeDataCols


### Cleaning phenotype data to have the same samples as genotype data

# names of samples with genotyped SVs
genotypeSamples <- as.matrix(read.table("samples.txt"))
# returns the samples that are in common
intersectedSamples <- intersect(colnames(rankNormMatrix), genotypeSamples)
# write common samples to a file, for subsetting from SV VCF
write.table(intersectedSamples, file = "intersectedSamples.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# subsets the phenotype data to only include in common samples
cleanedPhenotypeData <- subset(rankNormMatrix, select = intersectedSamples)
# returns the samples that are in the genotype data that are not in phenotype data
uniqueGenotypeSamples <- setdiff(genotypeSamples, intersectedSamples)
# returns the samples that are in the phenotype data that are not in genotype data
uniquePhenotypeSamples <- setdiff(colnames(rankNormMatrix), intersectedSamples)


### Reformatting Geuvadis phenotype data for fastQTL

# extracts subsetted gene IDs
selectedGenes <- as.matrix(rownames(cleanedPhenotypeData))
# extracts gene metadata from original DGEList object
geneInfo <- as.matrix(dgeFullData$genes)
# subsets metadata for start position, end position, and gene ID (same as rowname)
addPhenotypeRows <- subset (geneInfo, select = c("start", "end", "gene_id"))
# removes genes that are not in normalized and cleaned count data
cleanedAddPhenotypeRows <- subset(addPhenotypeRows, rownames(addPhenotypeRows) %in% selectedGenes)
# adds metadata columns to count matrix
phenotypeInput <- cbind(cleanedAddPhenotypeRows, cleanedPhenotypeData)
# write data to file
write.table(phenotypeInput, file = "phenotypeData.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)