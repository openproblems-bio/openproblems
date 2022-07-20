# Standardized CCC method runs via LIANA
# @param sce SingleCellExperiment
# @param method character name for method to be called via LIANA
# @return a tibble with liana results

# Libraries
library(SingleCellExperiment)
library(tidyverse)
library(liana)

# Parameters
expr_prop <- 0.1
idents_col <- "label"

# Convert data to an R-friendly sparse format
sce@assays@data$counts <- as(as.matrix(sce@assays@data$counts), "sparseMatrix")
sce@assays@data$logcounts <- as(as.matrix(sce@assays@data$logcounts), "sparseMatrix")


# Obtain LIANA's Consensus resource
op_resource <- select_resource("Consensus")[[1]]

# Check if the target organism is human
if(sce@metadata$target_organism!=9606){
    # Generate orthologous resource
    op_resource <- generate_homologs(op_resource = op_resource,
                                     target_organism = sce@metadata$target_organism)
}

# Run LIANA
liana_res <- liana_wrap(sce,
                        resource = 'custom',
                        external_resource = op_resource,
                        method=method,
                        expr_prop = expr_prop,
                        idents_col = idents_col,
                        base = 2.718282 # Relevant only for logfc
                        )

# Return
liana_res
