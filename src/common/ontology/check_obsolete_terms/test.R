library(assertthat)

## VIASH START
meta <- list(
  executable = "target/docker/common/ontology/check_obsolete_terms",
  resources_dir = "resources_test/common/"
)
## VIASH END

input_file <- paste0(meta$resources_dir, "/cellxgene_census/dataset.h5ad")
ontology_file <- paste0(meta$resources_dir, "/cellxgene_census/cl.obo")
temp_file <- tempfile(fileext = ".h5ad")
temp2_file <- tempfile(fileext = ".h5ad")

# add obsolete terms to the dataset
input <- anndata::read_h5ad(input_file)
input$obs$cell_type_ontology_term_id <- as.character(input$obs$cell_type_ontology_term_id)
input$obs$cell_type_ontology_term_id[1:3] <- "CL:0000375" # obsolete, replaced by 'CL:0007010'
input$obs$cell_type_ontology_term_id[4:6] <- "CL:0000399" # obsolete, removed
input$obs$cell_type_ontology_term_id[7:9] <- "CL:0007011" # not obsolete
zzz <- input$write_h5ad(temp_file)

# run component
zzz <- processx::run(
  meta$executable,
  c(
    "--input", temp_file,
    "--struct", "obs",
    "--input_term", "cell_type_ontology_term_id",
    "--ontology", ontology_file,
    "--output", temp2_file,
    "--output_term", "cell_type_ontology_term_id_new",
    "--output_name", "cell_type_new",
    "--output_obsolete", "cell_type_obsolete_new"
  ),
  echo = TRUE
)

# check output
output <- anndata::read_h5ad(temp2_file)

print(output$obs[1:10, , drop = FALSE])

assert_that(
  all(output$obs$cell_type_ontology_term_id_new[1:3] == "CL:0007010"),
  all(is.na(output$obs$cell_type_ontology_term_id_new[4:6])),
  all(output$obs$cell_type_ontology_term_id_new[7:9] == "CL:0007011"),
  all(output$obs$cell_type_new[1:3] == "preosteoblast"),
  all(is.na(output$obs$cell_type_new[4:6])),
  all(output$obs$cell_type_new[7:9] == "enteric neuron"),
  all(!output$obs$cell_type_obsolete_new[1:3]),
  all(output$obs$cell_type_obsolete_new[4:6]),
  all(!output$obs$cell_type_obsolete_new[7:9])
)
