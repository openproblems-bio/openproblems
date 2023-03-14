# Standardized CCC method runs via LIANA
# @param sce SingleCellExperiment
# @param min_expression_prop minimum proportion of cells per cell type
# expressing the ligands and receptors for a given interaction
# @param idents_col cell identity column name
# @param test logical; if a test run
# @return a tibble with liana results

# Libraries
library(SingleCellExperiment, quietly = TRUE)
library(tibble, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(liana, quietly = TRUE)

# Parameters
nperms <- if (test) 2 else 1000

# Convert data to an R-friendly sparse format
sce@assays@data$counts <-
  as(as.matrix(sce@assays@data$counts), "sparseMatrix")
sce@assays@data$logcounts <-
  as(as.matrix(sce@assays@data$logcounts), "sparseMatrix")



# Run LIANA
liana_res <- liana_wrap(sce,
  resource = "custom",
  external_resource = op_resource,
  expr_prop = min_expression_prop,
  idents_col = idents_col,
  base = 2.718282, # Relevant only for logfc
  permutation.params = list(nperms = nperms), # Relevant only for CPDB
  ...
)

# Aggregate if a run /w multiple methods
if (!is.tibble(liana_res)) {
  liana_res <- liana_res %>%
    liana_aggregate(aggregate_how = aggregate_how) %>%
    # inverse distribution
    mutate(aggregate_rank = 1 - aggregate_rank)
}

# Return (Keep Complexes [not subunits] for Consistency)
liana_res %>%
  dplyr::select(-any_of(c("ligand", "receptor"))) %>%
  dplyr::rename(ligand = ligand.complex, receptor = receptor.complex)
