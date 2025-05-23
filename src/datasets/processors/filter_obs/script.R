requireNamespace("anndata", quietly = TRUE)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)

## VIASH START
#par <- list(
#  input = "resources/neurips-2023-raw/sc_counts_reannotated_with_counts.h5ad",
#  output = "resources/datasets/neurips-2023-data/sc_counts_cleaned.h5ad"
#)
## VIASH END

# Load data
input <- anndata::read_h5ad(par$input)

cat("Input:\n")
print(input)

# set up obs_filt
obs_filt <- rep(TRUE, nrow(input$obs))

# Alvocidib only T cells in only 2 donors, remove
obs_filt <- obs_filt & input$obs$sm_name != "Alvocidib"
# BMS-387032 - one donor with only T cells, two other consistent, but only 2 cell types - leave the 2 cell types in, remove donor 2 with only T cells
obs_filt <- obs_filt & !(input$obs$sm_name == "BMS-387032" & input$obs$donor_id == "Donor 2")
# BMS-387032 remove myeloid cells and B cells
obs_filt <- obs_filt & !(input$obs$sm_name == "BMS-387032" & input$obs$cell_type %in% c("B cells", "Myeloid cells"))
# CGP 60474 has only T cells left, remove
obs_filt <- obs_filt & input$obs$sm_name != "CGP 60474"
# Canertinib - the variation of Myeloid cell proportions is very large, skip Myeloid
obs_filt <- obs_filt & !(input$obs$sm_name == "Canertinib" & input$obs$cell_type == "Myeloid cells")
# Foretinib - large variation in Myeloid cell proportions (some in T cells), skip Myeloid.
obs_filt <- obs_filt & !(input$obs$sm_name == "Foretinib" & input$obs$cell_type == "Myeloid cells")
# Ganetespib (STA-9090) - donor 2 has no Myeloid and small NK cells proportions. Skip Myeloid, remove donor 2
obs_filt <- obs_filt & !(input$obs$sm_name == "Ganetespib (STA-9090)" & input$obs$donor_id == "Donor 2")
# IN1451 - donor 2 has no NK or B, remove Donor 2
obs_filt <- obs_filt & !(input$obs$sm_name == "IN1451" & input$obs$donor_id == "Donor 2")
# Navitoclax - donor 3 doesn't have B cells and has different T and Myeloid proportions, remove donor 3
obs_filt <- obs_filt & !(input$obs$sm_name == "Navitoclax" & input$obs$donor_id == "Donor 3")
# PF-04691502 remove Myeloid (only present in donor 3)
obs_filt <- obs_filt & !(input$obs$sm_name == "PF-04691502" & input$obs$cell_type == "Myeloid cells")
# Proscillaridin A;Proscillaridin-A remove Myeloid, since the variation is very high (4x)
obs_filt <- obs_filt & !(input$obs$sm_name == "Proscillaridin A;Proscillaridin-A" & input$obs$cell_type == "Myeloid cells")
# R428 - skip NK due to high variation (close to 3x)
obs_filt <- obs_filt & !(input$obs$sm_name == "R428" & input$obs$cell_type == "NK cells")
# UNII-BXU45ZH6LI - remove due to large variation across all cell types and missing cell types
obs_filt <- obs_filt & input$obs$sm_name != "UNII-BXU45ZH6LI"

# filter genes
output <- input[obs_filt, ]

cat("Output:\n")
print(output)

cat("Writing output to ", par$output, "\n")
zzz <- output$write_h5ad(par$output, compression = "gzip")
