suppressWarnings(suppressMessages(library(infercnv, quietly = TRUE)))
suppressWarnings(suppressMessages(library(tidyverse, quietly = TRUE)))
suppressWarnings(suppressMessages(library(Seurat, quietly = TRUE)))
suppressWarnings(suppressMessages(library(phylogram, quietly = TRUE)))
suppressWarnings(suppressMessages(library(ape, quietly = TRUE)))
suppressWarnings(suppressMessages(library(hdf5r, quietly = TRUE)))
library(SpatialInferCNV, quietly = TRUE)

# Get the command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3) {
  InputAnnotationFile <- args[1]
  InputCountFile <- args[2]
  LibraryName <- args[3]
} else {
  stop("Please provide the path of InputAnnotationFile, InputCountFile, and LibraryName")
}

GeneOrder <- "/diskmnt/Projects/Users/efang/SpatialInferCNV/siCNV_GeneOrderFile.tsv"
outputDir <- sprintf("/diskmnt/Projects/Users/efang/SpatialInferCNV/InferCNVRunsOutputValidation/InferCNVrun_outputs-%s", LibraryName)
joinedCountsName <- sprintf("%s-Joined_Counts.tsv", LibraryName)
annotationsName <- sprintf("%s-FinalAnnotationsForExport.tsv", LibraryName)


histologicalAnnotations <- ImportHistologicalAnnotations("H10X", InputAnnotationFile)
dataCounts <- ImportCountData("H10X", InputCountFile)

# Merge annotation with counts
joinedCountsAnnotations <- MergingCountAndAnnotationData("H10X", histologicalAnnotations, dataCounts)
joinedCountsAnnotations <- joinedCountsAnnotations %>% column_to_rownames(var = "Genes")

# remove the counts data (not needed)
rm(dataCounts)

H10FinalAnnotationsForExport <- FinalAnnotations(histologicalAnnotations, joinedCountsAnnotations)

# write the data frames locally
write.table(joinedCountsAnnotations, joinedCountsName, sep = "\t")
write.table(H10FinalAnnotationsForExport, annotationsName, 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

# Create the inferCNV object from the 3 files
H10_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix=joinedCountsName, 
                                             gene_order_file=GeneOrder,
                                             annotations_file=annotationsName,
                                             delim="\t",
                                             ref_group_names=c("normal"),) 
# Run inferCNV
H10_infCNV = infercnv::run(H10_infCNV, 
                           cutoff=0.1, #0.1 means 10X data
                           out_dir=outputDir, 
                           cluster_by_groups=T,
                           HMM=T, 
                           denoise=TRUE) #denoising applies noise reduction for the plot 
file.remove(joinedCountsName)
file.remove(annotationsName)
