library(dplyr)  # data filtration
library(survival)  # time-to-event analysis (used for analyzing the probablity of a patient dying given the time they've lived)
library(ggplot2)
library(biomaRt)  # gene id conversions & other bio data
#library(rCausalMGM)
library(stringr)
library(DESeq2)  # gene expression analysis
library(vsn)  # variance, stabilization, and normalization --> data normalization


# finds genes that represent clusters & subclusters based on correlation patterns
# parameters respectively: correlation matrix, single linkage cuttoff, & average linkage cutoff
    # cutoff- the max distance data points can be from each other to be partitioned into a gene cluster
    # linkage - measuring the distance between clusters/data points
hubSummarize <- function(cor.mat, single.cutoff, average.cutoff) {
    groups <- list()
    
    # hierarchical clustering genes using single linkage --> nearest neighbors grouping
    # (1-abs(cor.mat)) transforms the correlation matrix that measures similarity into dissimilarity matrix
    # as.dist coverts the dissimilarity matrix into a matrix of distances btwn data
    hc.single <- hclust(as.dist(1-abs(cor.mat)), method='single')

    clusts <- cutree(hc.single, h=single.cutoff)  # assigns genes to clusters based on the cutoff value
    nclusts <- length(unique(clusts))  #number of unique clusters

    keep.genes <- c()  # names of the genes kept based on the logic below

    # Iterate over each cluster
    for (idx in seq_len(nclusts)) {
      # Get the genes belonging to the current cluster
      gene.group <- names(clusts[clusts == idx])
      
      # If the cluster contains only one gene
      if (length(gene.group) == 1) {
        # Add the gene to the list of kept genes
        keep.genes <- c(keep.genes, gene.group)
      } else {
        # Calculate the average linkage hierarchical clustering of the genes in the cluster
        hc.ave <- hclust(as.dist(1 - abs(cor.mat[gene.group, gene.group])), method = 'average')
        
        # Assign genes to subclusters based on the average linkage clustering
        clusts.group <- cutree(hc.ave, h = average.cutoff)
        nclusts.group <- length(unique(clusts.group))
        
        # Print the index of the current cluster
        print(idx)
        
        # Iterate over each subcluster within the current cluster
        for (jdx in seq_len(nclusts.group)) {
          # Get the genes belonging to the current subcluster
          gene.group.group <- names(clusts.group[clusts.group == jdx])
          
          # If the subcluster contains more than one gene
          if (length(gene.group.group) > 1) {
            # Calculate the correlation matrix for the genes in the subcluster
            group.cor.mat <- cor.mat[gene.group.group, gene.group.group]
            
            # Calculate the average correlations for the genes in the subcluster
            ave.cors <- (colSums(abs(group.cor.mat)) - 1) / (length(gene.group.group) - 1)
            
            # If the maximum correlation (excluding self-correlation) is below a threshold
            if (max(group.cor.mat[upper.tri(group.cor.mat)]) < (1 - single.cutoff)) {
              # Add the genes in the subcluster to the list of kept genes
              keep.genes <- c(keep.genes, gene.group.group)
            } else {
              # Select the gene with the highest average correlation
              max.cor.gene <- names(which.max(ave.cors))
              
              # Print the average correlations and the gene with the highest correlation
              print(ave.cors)
              print(max.cor.gene)
              
              # Add the gene with the highest correlation to the list of kept genes
              keep.genes <- c(keep.genes, max.cor.gene)
              
              # Store the subcluster genes under the selected gene in a list
              groups[[max.cor.gene]] <- gene.group.group
            }
          } else {
            # Add the gene in the subcluster to the list of kept genes
            keep.genes <- c(keep.genes, gene.group.group)
          }
        }
      }
    }
    
    # Return the list of kept genes and the groups of genes
    return(list(keep = keep.genes, groups = groups))
}


#### LOAD AND PREPARE RNA-SEQ DATA


# Read the gene expression counts from a TSV file
counts <- read.csv('COPD-Gene-Expression-Counts/2021_counts_raw.tsv', sep = '\t', row.names = 1)

# Display the dimensions of the counts matrix
dim(counts)

# Set the configuration to ignore SSL certificate verification
httr::set_config(httr::config(ssl_verifypeer = 0L))

# Use the "ensembl" dataset from the "hsapiens_gene_ensembl" BioMart
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# Retrieve information about genes based on Ensembl IDs
genes <- row.names(counts)
gene.list <- getBM(filters = "ensembl_gene_id_version",
                   attributes = c("ensembl_gene_id_version", "ensembl_gene_id",
                                  "hgnc_symbol", "description", "chromosome_name",
                                  "gene_biotype", "source"),
                   values = genes, mart = mart)

# Display the first few rows and the dimensions of the gene list
head(gene.list)
dim(gene.list)

# Create a table of gene biotypes
table(gene.list$gene_biotype)

# Filter the gene list to include only protein coding, miRNA, and lncRNA genes,
# while excluding ribosomal/mitochondrial genes and genes on the Y chromosome
gene.list <- gene.list %>%
  filter(gene_biotype %in% c('protein_coding', 'miRNA', 'lncRNA')) %>%
  filter(hgnc_symbol != "") %>%
  filter(!(grepl('MT[-]', hgnc_symbol) |
             grepl('^RPS.', hgnc_symbol) |
             grepl('^RPL.', hgnc_symbol) |
             grepl('^MRPS.', hgnc_symbol) |
             grepl('^MRPL.', hgnc_symbol))) %>%
  filter(chromosome_name != "Y")

# Display the dimensions of the filtered gene list
dim(gene.list)

# Create a table of Ensembl gene IDs with their counts
ensembl.counts <- table(gene.list$ensembl_gene_id_version)

# Find duplicate Ensembl gene IDs
dup.ensembl <- names(ensembl.counts[ensembl.counts > 1])

# Filter the gene list to keep only one row per duplicated Ensembl gene ID
gene.list <- gene.list %>% group_by(ensembl_gene_id_version) %>% filter(row_number() == 1)

# Filter the counts matrix to keep only rows corresponding to the filtered gene list
counts <- counts[gene.list$ensembl_gene_id_version,]

# Display the dimensions of the updated counts matrix
dim(counts)

# Create a table of gene counts based on HGNC symbols
gene.counts <- table(gene.list$hgnc_symbol)

# Find duplicate HGNC symbols
dup.genes <- names(gene.counts[gene.counts > 1])

# Filter the gene list to keep only one row per duplicated HGNC symbol
dup.gene.list <- gene.list %>% filter(hgnc_symbol %in% dup.genes)

# Order the duplicated gene list by HGNC symbol
dup.gene.list <- dup.gene.list[order(dup.gene.list$hgnc_symbol),]

# Merge counts for duplicated genes and update the counts matrix
for (gene in dup.genes) {
  dup.ensembl.ids <- dup.gene.list$ensembl_gene_id_version[dup.gene.list$hgnc_symbol == gene]
  counts[dup.ensembl.ids[1],] <- colSums(counts[dup.ensembl.ids,])
  counts <- counts[-which(rownames(counts) == dup.ensembl.ids[-1]),]
}

# Filter the gene list to keep only rows that match the updated counts matrix
gene.list <- gene.list %>% filter(ensembl_gene_id_version %in% rownames(counts))

# Identify the Ensembl gene IDs that are not present in the updated counts matrix
rownames(counts)[!rownames(counts) == gene.list$ensembl_gene_id_version]
gene.list$ensembl_gene_id_version[!rownames(counts) == gene.list$ensembl_gene_id_version]

# Rename the rows in the counts matrix with HGNC symbols
rownames(counts) <- gene.list$hgnc_symbol

# Display the first few rows and the dimensions of the updated counts matrix
head(counts)
dim(counts)

# Write the counts matrix to a CSV file
write.csv(counts, 'COPD-data/raw_counts.csv')

# Read the phenotype information from a TSV file
exprPheno <- read.csv('COPDgene-RNAseq/RNAseq_f4_masterpheno.tsv', sep = '\t')

# Select the unique sample IDs
ids <- unique(exprPheno$sid)

# Remove empty sample IDs
ids <- ids[which(ids != "")]

# Filter the phenotype information to keep only the samples with selected IDs
exprPheno <- exprPheno %>% filter(sid %in% ids)

# Rename the rows in the phenotype information with sample IDs
rownames(exprPheno) <- paste0('X', exprPheno$sid)

# Select the sample IDs that exist in both counts and phenotype matrices
samp.ids <- intersect(colnames(counts), rownames(exprPheno))

# Select the sample IDs that have non-missing values for specific phenotype columns
samp.ids <- rownames(na.omit(exprPheno[samp.ids, c('neutrophl_pct', 'lymphcyt_pct',
                                                   'monocyt_pct', 'eosinphl_pct')]))

# Scale the selected phenotype columns
exprPheno <- exprPheno %>% mutate_at(c('neutrophl_pct', 'lymphcyt_pct',
                                       'monocyt_pct', 'eosinphl_pct'), scale)

# Create a DESeq dataset from the filtered counts and scaled phenotype information
dds <- DESeqDataSetFromMatrix(round(counts[, samp.ids]),
                              exprPheno[samp.ids, ],
                              design = ~neutrophl_pct + lymphcyt_pct +
                                monocyt_pct + eosinphl_pct)

# Display the dimensions of the DESeq dataset
dim(dds)

# Compute quantiles of the row means of the DESeq dataset
quantile(rowMeans(counts(dds)), seq(0, 1, 0.05))

# Keep only the rows with a row mean greater than 0.5
keep <- rowMeans(counts(dds)) > 0.5
sum(keep)
dds <- dds[keep,]

# Display the dimensions of the updated DESeq dataset
dim(dds)

# Run DESeq normalization and differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Compute quantiles of the row means of the normalized DESeq dataset
quantile(rowMeans(counts(dds, normalize = TRUE)), seq(0, 1, 0.05))

# Transform the DESeq dataset using the variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Generate a mean versus standard deviation plot
meanSdPlot(assay(vsd))

# Save the PCA plots before and after ComBat normalization to a PDF file
pdf('pca_uncorrected.pdf', width = 12, height = 8)
plotPCA(vsd, 'Lc.Batch')
plotPCA(vsd, 'CBC_QC')
plotPCA(vsd, 'neutrophl_pct')
plotPCA(vsd, 'lymphcyt_pct')
plotPCA(vsd, 'monocyt_pct')
plotPCA(vsd, 'eosinphl_pct')
dev.off()

# Load the sva package for ComBat normalization
library(sva)

# Perform ComBat normalization on the counts matrix using batch and covariate information
combat.counts <- ComBat_seq(counts(dds),
                            batch = colData(dds)$Lc.Batch,
                            covar_mod = model.matrix(~neutrophl_pct + lymphcyt_pct +
                                                       monocyt_pct + eosinphl_pct,
                                                     colData(dds)))
write.csv(combat.counts, "COPD-data/combat_counts_batch_normalization.csv")

# Create a new DESeq dataset with the normalized counts and the original design formula
dds <- DESeqDataSetFromMatrix(combat.counts,
                              colData(dds),
                              design = ~neutrophl_pct + lymphcyt_pct +
                                monocyt_pct + eosinphl_pct)

write.csv(assay(dds), "COPD-data/DESeq_normalize_counts_data.csv")

# Display the dimensions of the updated DESeq dataset
dim(dds)

# Compute quantiles of the row means of the DESeq dataset after ComBat normalization
quantile(rowMeans(counts(dds)), seq(0, 1, 0.05))

# Run DESeq normalization and differential expression analysis on the normalized data
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Compute quantiles of the row means of the normalized DESeq dataset after ComBat normalization
quantile(rowMeans(counts(dds, normalize = TRUE)), seq(0, 1, 0.05))

# Transform the DESeq dataset after ComBat normalization using the VST
vsd <- vst(dds, blind = FALSE)

# Generate a mean versus standard deviation plot
meanSdPlot(assay(vsd))

# Save the PCA plots after ComBat normalization to a PDF file
pdf('pca_combat.pdf', width = 12, height = 8)
plotPCA(vsd, 'Lc.Batch')
plotPCA(vsd, 'CBC_QC')
plotPCA(vsd, 'neutrophl_pct')
plotPCA(vsd, 'lymphcyt_pct')
plotPCA(vsd, 'monocyt_pct')
plotPCA(vsd, 'eosinphl_pct')
dev.off()

# Write the normalized counts matrix to a CSV file
write.csv(counts(dds), 'COPD-data/combat_counts.csv')

# Write the normalized counts matrix with normalization factors to a CSV file
write.csv(counts(dds, normalize = TRUE), 'COPD-data/combat_norm_counts.csv')


# Finding the sd of each gene
genes_sd <- apply(assay(vsd), 1, function(x) sd(x))

# Sorting the genes based on their sd from most variant to least
sorted_genes <- names(genes_sd)[order(genes_sd, decreasing = TRUE)] # Order gives indices of genes  
top_500_variative_genes <- sorted_genes[1:500]

# Correlation matrix of the top 500 most variative genes
top_500_assay <- assay(vsd)[top_500_variative_genes, ]
top_500_cor_mat <- cor(t(top_500_assay), method = "spearman") # Flipping rows & columns

heatmap(top_500_cor_mat, main = "Spearman Correlation Heatmap")

library(pheatmap)
pheatmap(top_500_cor_mat, breaks=seq(-1,1,length.out=101))

# Set of genes to be used in the model
model_genes <- hubSummarize(top_500_cor_mat, 0.15, 0.2)

pheatmap(top_500_cor_mat[model_genes$keep,model_genes$keep], breaks=seq(-1,1,length.out=101))

median(top_500_cor_mat['BCL2L1', model_genes$groups$BCL2L1])

# Getting the count data for the representative genes being kept
representative_genes_counts <- assay(vsd)[model_genes$keep, ]

# Tagging the representative genes with the identifier "_grp"
for (i in 1:length(rownames(representative_genes_counts))) {
  for (j in 1:length(names(model_genes$groups))) {
    if (rownames(representative_genes_counts)[i] == names(model_genes$groups)[j]) {
      rownames(representative_genes_counts)[i] <- paste(names(model_genes$groups)[j],"_grp", sep = "")
    }
  }
}

write.csv(representative_genes_counts, "COPD-data/representative_gene_counts.csv")

# Sorting through the normalized ComBat corrected counts to only include genes in lm22
library(readxl)
lm22_data <- read_excel("LM22.xlsx")

# Cutting the unnecessary text to get list of lm22 genes
lm22_genes <- lm22_data[-c(1:12), ]
lm22_genes <- lm22_genes[, 3]
lm22_genes <- lm22_genes[-c(1:1), ]
lm22_genes <- as.vector(lm22_genes)$...3

# Getting only the genes in lm22 from the gene counts with batch effects removed to put into CIBERSORT
combat_counts_norm = read.csv("COPD-data/combat_norm_counts.csv", sep = ",", header=TRUE, row.names=1)

in_lm22 <- c()

for (gene in rownames(combat_counts_norm)) {
  if (gene %in% lm22_genes) {
    in_lm22 <- c(in_lm22, gene)
  }
}

print(in_lm22)

# CIBERSORTx deconvolution

# lm22 genes with their count data
combat_counts_norm_cibersort <- combat_counts_norm[rownames(combat_counts_norm) %in% lm22_genes, ]

# Creating the table inputted into CIBERSORTx (must add tab in first row before first character)
write.table(combat_counts_norm_cibersort, "COPD-data/combat_counts_normalized_cibersort.txt", sep="\t", quote=FALSE)


# InstaPrism deconvolution

# # Cleaning the lm22_data
# lm22_prism <-lm22_data
# lm22_prism <- lm22_data[-c(1:10), -c(1:2)]
# 
# # Naming the NA cell types
# for (i in 2:length(colnames(lm22_prism))) {
#   
#   if (is.na(lm22_prism[1, i])) {
#     lm22_prism[1, i] <- lm22_prism[1, i-1]
#   }
# }
# lm22_prism <- lm22_prism[-c(2), ]
# lm22_prism[1,1] <- NA

# write.csv(lm22_prism, "COPD-data/cleaned_lm22.csv")

# Cleaned lm22 from loop alterations & additions above
lm22_cleaned <- read.csv("COPD-data/cleaned_lm22.csv", sep = ",")
lm22_cleaned <- lm22_cleaned[, -c(1)]

# Getting cell states & types
lm22_cell_types <- lm22_cleaned[1, 2:length(lm22_cleaned[1,])]
lm22_cell_states <- lm22_cleaned[2, 2:length(lm22_cleaned[1,])]

# Removing the cell types row and making the appropriate row & column the headers
lm22_cleaned <- lm22_cleaned[-c(1), ]

colnames(lm22_cleaned) <- as.character(lm22_cleaned[1, ])
lm22_cleaned <- lm22_cleaned[-1, ] # Remove the first row after making it the column names

rownames(lm22_cleaned) <- as.character(lm22_cleaned[, 1]) # Remove the first column after making it the row names 
lm22_cleaned <- lm22_cleaned[, -1]

# Converting character matricies into numeric matricies
lm22_cleaned <- lm22_cleaned %>% mutate_all(as.numeric)
combat_counts_norm <- combat_counts_norm %>% mutate_all(as.numeric)


# Running InstaPrism

library(InstaPrism)

InstaPrism.res = InstaPrism(input_type = 'raw', sc_Expr = as.matrix(lm22_cleaned), 
                            bulk_Expr = as.matrix(combat_counts_norm),
                            cell.type.labels = lm22_cell_types, cell.state.labels = lm22_cell_states,
                            n.iter=50, convergence.plot = TRUE)

# InstaPrism cell state output
InstaPrism_cs <- t(InstaPrism.res@Post.ini.cs@theta)

# InstaPrism cell state output
InstaPrism_ct <- t(InstaPrism.res@Post.ini.ct@theta)

# Loading the clinical data
clin <- read.csv('COPDGene_Data-P1P2P3.2021Aug/COPDGene_P1P2P3_SM_NS_Long_Mar20-2.txt', sep='\t')

# Filtering clinical data for phase 2 patients only
clin.p2 <- clin %>% filter(Phase_study==2)

rownames(clin.p2) <- clin.p2$sid

# Loading the CIBERSORT output
cibersort_output <- read.csv('COPD-Data/cibersort_ouput.csv', sep=",")

# Taking out the leading "X" in the patient ids
cibersort_output[, 1] <- gsub("^X", "", cibersort_output[, 1])
rownames(cibersort_output) <- cibersort_output[, 1]
rownames(InstaPrism_cs) <- gsub("^X", "", rownames(InstaPrism_cs))

# Common patients in CIBERSORT & clinical data
common_patients <- intersect(clin.p2$sid, cibersort_output$Mixture)

# Plotting Clinical vs CIBERSORT data
par(mar = c(4, 4, 2, 1))
plot(clin.p2[common_patients,'neutrophl'], cibersort_output[common_patients,'Neutrophils'], main = "Clinical vs CIBERSORT: neutrophl",
     xlab = "clinical", ylab = "cibersort")

plot(clin.p2[common_patients,'monocyt'], cibersort_output[common_patients,'Monocytes'], main = "Clinical vs CIBERSORT: monocyt",
     xlab = "clinical", ylab = "cibersort")

plot(clin.p2[common_patients,'eosinphl'], cibersort_output[common_patients,'Eosinophils'], main = "Clinical vs CIBERSORT: eosinphl",
     xlab = "clinical", ylab = "cibersort")

# Clin & CIBERSORT
# Only keeping rows that have non-NA values in both matricies
keep_rows_neutrophl <- complete.cases(clin.p2[common_patients,'neutrophl'], cibersort_output[common_patients,'Neutrophils'])
keep_rows_monocyt <- complete.cases(clin.p2[common_patients,'monocyt'], cibersort_output[common_patients,'Monocytes'])
keep_rows_eosinphl <- complete.cases(clin.p2[common_patients,'eosinphl'], cibersort_output[common_patients,'Eosinophils'])

pearson_correlation_neutrophl <- cor(clin.p2[common_patients,'neutrophl'][keep_rows_neutrophl], cibersort_output[common_patients,'Neutrophils'][keep_rows_neutrophl], method = "pearson")
pearson_correlation_monocyt <- cor(clin.p2[common_patients,'monocyt'][keep_rows_monocyt], cibersort_output[common_patients,'Monocytes'][keep_rows_monocyt], method = "pearson")
pearson_correlation_eosinphl <- cor(clin.p2[common_patients,'eosinphl'][keep_rows_eosinphl], cibersort_output[common_patients,'Eosinophils'][keep_rows_eosinphl], method = "pearson")

print(paste("Pearson correlations clinical vs CIBERSORT: neutrophl = ", pearson_correlation_neutrophl, ", monocyt = ", pearson_correlation_monocyt, ", eosinphl = ", pearson_correlation_eosinphl))


# Common InstaPrism & CIBERSORT and InstaPrism & clinical
common_patients_ciber <- intersect(rownames(InstaPrism_cs), cibersort_output$Mixture)
common_patients_clin <- intersect(clin.p2$sid, rownames(InstaPrism_cs))

# Plotting Clinical vs InstaPrism data
plot(clin.p2[common_patients_clin, 'neutrophl_pct'], InstaPrism_cs[common_patients_clin, "Neutrophils"], main = "Clinical vs InstaPrism: neutrophl",
     xlab = "clinical", ylab = "InstaPrism")

plot(clin.p2[common_patients_clin, 'monocyt_pct'], InstaPrism_cs[common_patients_clin, "Monocytes"], main = "Clinical vs InstaPrism: monocyt",
     xlab = "clinical", ylab = "InstaPrism")

plot(clin.p2[common_patients_clin, 'eosinphl_pct'], InstaPrism_cs[common_patients_clin, "Eosinophils"], main = "Clinical vs InstaPrism: eosinphl",
     xlab = "clinical", ylab = "InstaPrism")


# Converting the absolute numbers in CIBERSORT output into percentages for InstaPrism comparison
cibersort_output_pct <- cibersort_output

for (row in 1:nrow(cibersort_output_pct)) {
  for (cell_state in 2:(ncol(cibersort_output_pct) - 1)) {
    cibersort_output_pct[row, cell_state] <- cibersort_output_pct[row, cell_state] / cibersort_output_pct[row, "Absolute.score..sig.score."]
  }
}

# Plotting InstaPrism vs CIBERSORT
plot(InstaPrism_cs[common_patients_ciber, "Neutrophils"], cibersort_output_pct[common_patients_ciber, "Neutrophils"], main = "InstaPrism vs CIBERSORT: neutrophl",
     xlab = "InstaPrism", ylab = "cibersort")

plot(InstaPrism_cs[common_patients_ciber, "Monocytes"], cibersort_output_pct[common_patients_ciber, "Monocytes"], main = "InstaPrism vs CIBERSORT: monocyt",
     xlab = "InstaPrism", ylab = "cibersort")

plot(InstaPrism_cs[common_patients_ciber, "Eosinophils"], cibersort_output_pct[common_patients_ciber, "Eosinophils"], main = "InstaPrism vs CIBERSORT: eosinphl",
     xlab = "InstaPrism", ylab = "cibersort")


# Clin & InstaPrism
keep_rows_neutrophl <- complete.cases(clin.p2[common_patients_clin, 'neutrophl_pct'], InstaPrism_cs[common_patients_clin, "Neutrophils"])
keep_rows_monocyt <- complete.cases(clin.p2[common_patients_clin, 'monocyt_pct'], InstaPrism_cs[common_patients_clin, "Monocytes"])
keep_rows_eosinphl <- complete.cases(clin.p2[common_patients_clin, 'eosinphl_pct'], InstaPrism_cs[common_patients_clin, "Eosinophils"])

pearson_correlation_neutrophl <- cor(clin.p2[common_patients,'neutrophl_pct'][keep_rows_neutrophl], InstaPrism_cs[common_patients_clin, "Neutrophils"][keep_rows_neutrophl], method = "pearson")
pearson_correlation_monocyt <- cor(clin.p2[common_patients,'monocyt_pct'][keep_rows_monocyt], InstaPrism_cs[common_patients_clin, "Monocytes"][keep_rows_monocyt], method = "pearson")
pearson_correlation_eosinphl <- cor(clin.p2[common_patients,'eosinphl_pct'][keep_rows_eosinphl], InstaPrism_cs[common_patients_clin, "Eosinophils"][keep_rows_eosinphl], method = "pearson")

print(paste("Pearson correlations clinical vs InstaPrism: neutrophl = ", pearson_correlation_neutrophl, ", monocyt = ", pearson_correlation_monocyt, ", eosinphl = ", pearson_correlation_eosinphl))

# InstaPrism & CIBERSORT
keep_rows_neutrophl <- complete.cases(InstaPrism_cs[common_patients_ciber, "Neutrophils"], cibersort_output[common_patients_ciber,'Neutrophils'])
keep_rows_monocyt <- complete.cases(InstaPrism_cs[common_patients_ciber, "Monocytes"], cibersort_output[common_patients_ciber,'Neutrophils'])
keep_rows_eosinphl <- complete.cases(InstaPrism_cs[common_patients_ciber, "Eosinophils"], cibersort_output[common_patients_ciber,'Neutrophils'])

pearson_correlation_neutrophl <- cor(InstaPrism_cs[common_patients_ciber, "Neutrophils"][keep_rows_neutrophl], cibersort_output[common_patients_ciber,'Neutrophils'][keep_rows_neutrophl], method = "pearson")
pearson_correlation_monocyt <- cor(InstaPrism_cs[common_patients_ciber, "Monocytes"][keep_rows_monocyt], cibersort_output[common_patients_ciber,'Monocytes'][keep_rows_monocyt], method = "pearson")
pearson_correlation_eosinphl <- cor(InstaPrism_cs[common_patients_ciber, "Eosinophils"][keep_rows_eosinphl], cibersort_output[common_patients_ciber,'Eosinophils'][keep_rows_eosinphl], method = "pearson")

print(paste("Pearson correlations InstaPrism vs CIBERSORT: neutrophl = ", pearson_correlation_neutrophl, ", monocyt = ", pearson_correlation_monocyt, ", eosinphl = ", pearson_correlation_eosinphl))


# Looking at the graphs to determine what cell types should be removed
InstaPrism_cs[InstaPrism_cs<0.0001] <- 0

library(reshape2)

plot(InstaPrism_cs[common_patients,'Plasma cells'], cibersort_output_pct[common_patients,'Plasma.cells'])

plotdf <- melt(InstaPrism_cs)

colMeans(InstaPrism_cs<1e-4)[colMeans(InstaPrism_cs<1e-4)>0.5]

colMeans(InstaPrism_cs)


colMeans(cibersort_output_pct==0)

colMeans(clin.p2[common_patients,grep('_pct', colnames(clin.p2))]==0)

head(plotdf)

library(ggplot2)

ggplot(plotdf, aes(x=1, y=value, fill=Var2)) +
  geom_violin() + 
  facet_wrap(vars(Var2), scales='free', ncol=3) +
  theme_classic()

plotdf <- melt(cibersort_output_pct)

head(plotdf)

library(ggplot2)

ggplot(plotdf, aes(x=1, y=value, fill=variable)) +
  geom_violin() + 
  facet_wrap(vars(variable), scales='free', ncol=3) +
  theme_classic()


# Removing cell types that are measured in the clinical data & cell types that are barely detected from InstaPrism
# Columns needed to be removed: colMeans(InstaPrism_cs<1e-4)[colMeans(InstaPrism_cs<1e-4)>0.5]
InstaPrism_cs_for_model <- InstaPrism_cs

columns_deleted <- c('T cells follicular helper', 'T cells gamma delta', 'Dendritic cells resting','Mast cells resting',
                     'Mast cells activated', 'Eosinophils', 'Monocytes', 'Neutrophils')
InstaPrism_cs_for_model <- InstaPrism_cs_for_model[, !(colnames(InstaPrism_cs_for_model) %in% columns_deleted)]

write.csv(InstaPrism_cs_for_model, "COPD-data/InstaPrism_deconvolution_for_model.csv")

#### LOAD AND PREPARE PHASE 2 CLINICAL DATA

library(readxl)

clin.dict <- read_excel('/Users/neil/Documents/Hillman/Hillman23/COPDGene/COPDGene_Data-P1P2P3.2021Aug/COPDGene_P1P2P3_Visitlevel_DataDict_Mar20_rev_16Aug21.xlsx', 1) %>% as.data.frame

rownames(clin.dict) <- clin.dict$VariableName

head(clin.dict)

factor.vars <- c(setdiff(rownames(clin.dict)[clin.dict$CodedVariable=='Y'],
                       c('HealthStatus', 'SchoolCompleted', 'ATS_ERS',
                         'DrugCostCovered', 'LungDiseaseInformed')),
                 c('Internet', 'Insurance', 'CVD'))
numeric.vars <- c(setdiff(rownames(clin.dict)[clin.dict$CodedVariable=='N'],
                          c('Internet', 'Insurance')),
                  c('HealthStatus', 'SchoolCompleted', 'ATS_ERS', 
                    'DrugCostCovered', 'LungDiseaseInformed',
                    'neutrophl.to.lymphcyt.ratio', 'ADI_NATRANK'))


clin <- read.csv('COPDGene_Data-P1P2P3.2021Aug/COPDGene_P1P2P3_SM_NS_Long_Mar20-2.txt', sep='\t')


clin <- clin %>% filter(smoking_status != 0)

cvd.feats <- c('CongestHeartFail', 'CoronaryArtery', 'HeartAttack',
               'PeriphVascular', 'Stroke', 'TIA', 'CABG', 'Angioplasty')

clin$CVD <- as.integer(apply(clin[,cvd.feats], 1,
                             function(x) {
                                 ifelse(all(is.na(x)) |
                                        (any(is.na(x)) & all(x==0, na.rm=T)),
                                        NA,
                                        any(x==1, na.rm=T))
                             }))


clin.p2 <- clin %>% filter(Phase_study==2)

rownames(clin.p2) <- paste0('X', clin.p2$sid)


clin.p2.feats <- c('gender', 'race', 'smoking_status', 'age_visit',
                   'days_since_baseline',
                   'BMI', 'sysBP', 'diasBP', 'HR', 'Resting_SaO2',
                   'finalgold_visit', 
                   'distwalked', 'Diabetes', 'HighBloodPres', 'CVD',
                   'HighCholest', 'OsteoArth', 'Osteoporosis', 
                   'HaveCough', 'HavePhlegm', 'Chronic_Bronchitis',
                   'AwakeByCough', 'AwakeByShrtBrth',
                   'MMRCDyspneaScor', 'Asthma', 'Pneumonia', 
                   'SleepApnea', 
                   'FEV1pp_post', 'FVCpp_post', 'FEV1_FVC_post', 'PEF_post',
                   'BDR', 'ATS_PackYears', 'Duration_Smoking', 
                   'Severe_Exacerbations', 'Exacerbation_Frequency',
                   'neutrophl_pct', 'lymphcyt_pct',
                   'monocyt_pct', 'eosinphl_pct', 'basophl_pct',
                   'wbc','hemoglobin', 'MCV', 'MCHC',  
                   'Platelets', 'Income', 'InternetAccess',
                   'LungDiseaseInformed','SchoolCompleted',
                   'pctEmph_Thirona', 
                   'Pi10_Thirona', 'AWT_seg_Thirona')

colSums(is.na(clin.p2[colnames(counts),clin.p2.feats]))

clin.p2 <- na.omit(clin.p2[colnames(counts),clin.p2.feats])

table(clin.p2[colnames(counts),'finalgold_visit'])

mortality <- read.csv('All-Cause-Mortality-Patient-Data/COPDGene_VitalStatus_SM_NS_Oct22.csv', row.name=1)

rownames(mortality) <- paste0('X', rownames(mortality))

head(mortality)

mortality <- mortality[rownames(clin.p2),]
mortality$days_followed <- mortality$days_followed - clin.p2[rownames(mortality),'days_since_baseline']
mortality <- mortality[!is.na(mortality$days_followed) & mortality$days_followed>0,]

dim(mortality)
sum(mortality[,'vital_status'])

mortality$overall <- Surv(mortality[,'days_followed'], mortality[,'vital_status'])

##########
clin.p2 <- na.omit(cbind(clin.p2[,-which(colnames(clin.p2)=='days_since_baseline')], overall=mortality[rownames(clin.p2),'overall']))
rownames(clin.p2) <- gsub("^X", "", rownames(clin.p2))

# Reading in matricies for combination
# Representative genes
rep_gene_counts <- t(read.csv("COPD-data/representative_gene_counts.csv", sep=","))
colnames(rep_gene_counts) <- rep_gene_counts[1, ]
rep_gene_counts <- rep_gene_counts[-c(1), ]
rownames(rep_gene_counts) <- gsub("^X", "", rownames(rep_gene_counts))

# InstaPrism
InstaPrism_data<- read.csv("COPD-data/InstaPrism_deconvolution_for_model.csv", sep=",")
rownames(InstaPrism_data) <- InstaPrism_data[, 1]
InstaPrism_data <- InstaPrism_data[, -c(1)]

InstaPrism_data[InstaPrism_data < 1e-4] <- 0 # Making everything below 1e-4, 0

dim(clin.p2)
dim(rep_gene_counts)
dim(InstaPrism_data)

common_ids <- intersect(rownames(clin.p2), rownames(rep_gene_counts))
common_ids <- intersect(common_ids, rownames(InstaPrism_data))

# COPD model data that contains the clinical measured metrics, representative genes of gene clusters, and InstaPrism deconvolution for cell type prediction
copd_model_data <- cbind(clin.p2[common_ids, ], rep_gene_counts[common_ids, ], InstaPrism_data[common_ids, ])
write.csv(copd_model_data %>% mutate_at('overall', as.matrix), "COPD-data/copd_model_data.csv")

# Train & Test Data for All GOLD Stages
set.seed(14)
fold_id <- sample(rep(1:10, length.out=nrow(copd_model_data)))

for (i in 1:10) {
  test_data <- copd_model_data[fold_id==i,]
  test_data <- test_data[, colnames(test_data) != 'finalgold_visit'] # Taking gold stages out of model
  write.csv(test_data %>% mutate_at('overall', as.matrix), paste("Cross-Validation-Matricies/All-GOLD-stages/test/", "all_copd_model_data_test_", i, ".csv", sep = ""))
  
  train_data <- copd_model_data[fold_id!=i,]
  train_data <- train_data[, colnames(train_data) != 'finalgold_visit'] # Taking gold stages out of model
  write.csv(train_data %>% mutate_at('overall', as.matrix), paste("Cross-Validation-Matricies/All-GOLD-stages/train/", "all_copd_model_data_train_", i, ".csv", sep = ""))
}

# Train & Test Data for GOLD Stages 2-4

copd_model_data_2_4 <- copd_model_data
copd_model_data_2_4 <- copd_model_data_2_4[copd_model_data_2_4$finalgold_visit %in% c(2,3,4), ]
write.csv(copd_model_data_2_4, "COPD-data/2_4_copd_model_data.csv")

dim(copd_model_data_2_4)

# Adjusted
copd_model_data_2_4 <- read.csv("/Users/neil/Documents/Hillman/Hillman23/COPDGene/COPD-data-new/2_4_copd_model_data_adjusted_v2.csv", sep=",", row.names = 1)

set.seed(24)
fold_id <- sample(rep(1:10, length.out=nrow(copd_model_data_2_4)))

for (i in 1:10) {
  test_data <- copd_model_data_2_4[fold_id==i,]
  test_data <- test_data[, colnames(test_data) != 'finalgold_visit'] # Taking gold stages out of model
  write.csv(test_data, paste("/Users/neil/Documents/Hillman/Hillman23/COPDGene/Cross-Validation-Matricies-New/adjusted/test/", "2_4_copd_model_data_test_", i, ".csv", sep = ""))
  
  train_data <- copd_model_data_2_4[fold_id!=i,]
  train_data <- train_data[, colnames(train_data) != 'finalgold_visit'] # Taking gold stages out of model
  write.csv(train_data, paste("/Users/neil/Documents/Hillman/Hillman23/COPDGene/Cross-Validation-Matricies-New/adjusted/train/", "2_4_copd_model_data_train_", i, ".csv", sep = ""))
}

# Unadjusted
copd_model_data_2_4 <- read.csv("/Users/neil/Documents/Hillman/Hillman23/COPDGene/COPD-data-new/2_4_copd_model_data_unadjusted_v2.csv", sep=",", row.names = 1)

set.seed(24)
fold_id <- sample(rep(1:10, length.out=nrow(copd_model_data_2_4)))

for (i in 1:10) {
  test_data <- copd_model_data_2_4[fold_id==i,]
  test_data <- test_data[, colnames(test_data) != 'finalgold_visit'] # Taking gold stages out of model
  write.csv(test_data, paste("/Users/neil/Documents/Hillman/Hillman23/COPDGene/Cross-Validation-Matricies-New/unadjusted/test/", "2_4_copd_model_data_test_", i, ".csv", sep = ""))
  
  train_data <- copd_model_data_2_4[fold_id!=i,]
  train_data <- train_data[, colnames(train_data) != 'finalgold_visit'] # Taking gold stages out of model
  write.csv(train_data, paste("/Users/neil/Documents/Hillman/Hillman23/COPDGene/Cross-Validation-Matricies-New/unadjusted/train/", "2_4_copd_model_data_train_", i, ".csv", sep = ""))
}


library(rCausalMGM)

copd_model_data_2_4 <- copd_model_data_2_4 %>% mutate_if(colnames(copd_model_data_2_4) %in% factor.vars, factor)

ig.path <- mgmPath(copd_model_data_2_4[,!colnames(copd_model_data_2_4) %in% c('overall','finalgold_visit')],
                   rank=TRUE,
                   verbose=TRUE)

plot(ig.path)

plot(ig.path$graph.bic, 'overall', list(fontsize=54))
