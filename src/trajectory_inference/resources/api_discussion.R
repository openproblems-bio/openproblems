library(tidyverse)

# example from: https://github.com/rcannood/phdthesis/blob/master/ch3_dynbenchmark/fig/snote1fig_1.pdf

cell_ids <- c("a", "b", "c", "d", "e")
milestone_ids <- c("W", "X", "Y", "Z")

milestone_network <- tribble(
  ~from, ~to, ~length, ~directed,
  "W", "X", 2, FALSE,
  "X", "Y", 3, FALSE,
  "X", "Z", 4, FALSE
)

milestone_percentages <- tribble(
  ~cell_id, ~milestone_id, ~percentage,
  "a", "W", 0.9,
  "a", "X", 0.1,
  "b", "W", 0.2,
  "b", "X", 0.8,
  "c", "X", 0.8,
  "c", "Z", 0.2,
  "d", "X", 0.2,
  "d", "Y", 0.7,
  "d", "Z", 0.1,
  "e", "X", 0.3,
  "e", "Y", 0.2,
  "e", "Z", 0.5
)

divergence_regions <- tribble(
  ~divergence_id, ~milestone_id, ~is_start,
  "XYZ", "X", TRUE,
  "XYZ", "Y", FALSE,
  "XYZ", "Z", FALSE
)

library(dynwrap)
library(dynplot)

convert_milestone_percentages_to_progressions(
  cell_ids = cell_ids,
  milestone_ids = milestone_ids,
  milestone_network = milestone_network,
  milestone_percentages = milestone_percentages
)

traj <-
  wrap_data(cell_ids = cell_ids) %>%
  add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    divergence_regions = divergence_regions
  )

plot_graph(traj)

traj$progressions









