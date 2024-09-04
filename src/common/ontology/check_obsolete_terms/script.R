library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(ontologyIndex, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "resources_test/common/cellxgene_census/dataset.h5ad",
  ontology = "resources_test/common/cellxgene_census/cl.obo",
  input_term = "cell_type_ontology_term_id",
  struct = "obs",
  output = "output.h5ad",
  output_term = "cell_type_ontology_term_id",
  output_name = "cell_type",
  output_obsolete = "cell_type_ontology_obsolete",
  obsolete_as_na = TRUE
)
## VIASH END

cat("Read ontology\n")
ont <- ontologyIndex::get_ontology(
  par$ontology,
  extract_tags = "everything"
)
ont_tib <- ont %>%
  as.data.frame %>%
  select(id, name, obsolete, replaced_by) %>%
  as_tibble

cat("Read anndata\n")
adata <- anndata::read_h5ad(par$input, backed = "r")

cat("Find terms\n")
term_ids <- adata[[par$struct]][[par$input_term]]

unique_term_ids <- as.character(unique(term_ids))

cat("Look for obsolete or replaced terms\n")
ont_map <- ont_tib %>%
  slice(match(unique_term_ids, id)) %>%
  transmute(
    orig_id = id,
    id = case_when(
      !obsolete ~ id,
      replaced_by != "" ~ replaced_by,
      rep(par$obsolete_as_na, length(id)) ~ rep(NA_character_, length(id)),
      TRUE ~ id
    )
  ) %>%
  left_join(ont_tib %>% select(id, name, obsolete), by = "id")

cat("Store new columns in data structure\n")
new_data <- ont_map %>% slice(match(term_ids, orig_id))
adata[[par$struct]][[par$output_term]] <- new_data$id
adata[[par$struct]][[par$output_name]] <- new_data$name
adata[[par$struct]][[par$output_obsolete]] <- ifelse(
  !is.na(new_data$obsolete),
  new_data$obsolete,
  TRUE
)

cat("Write to file\n")
anndata::write_h5ad(adata, par$output)
