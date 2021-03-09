+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "result"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 20  # Order that this section will appear.

title = "Multimodal Data Integration"
subtitle = "*Realigning multimodal measurements of the same cell*"

[design]
  # Choose how many columns the section has. Valid values: 1 or 2.
  columns = "1"


  # Text color (true=light or false=dark).
  text_color_light = false

[design.spacing]
  # Customize the section spacing. Order is top, right, bottom, left.
  padding = ["20px", "0", "20px", "0"]

[advanced]
 # Custom CSS.
 css_style = ""

 # CSS class.
 css_class = ""

+++
## The task

Several recently described technologies allow for simultaneous measurement of different aspects of cell state. For example, [sci-CAR](https://doi.org/10.1126/science.aau0730) jointly profiles RNA expression and chromatin accessibility on the same cell and [CITE-seq](https://doi.org/10.1038/nmeth.4380) measures surface protein abundance and RNA expression from each cell. However, these joint profiling methods have several tradeoffs compared to unimodal measurements.

Joint methods can be more expensive or lower throughput or more noisy than measuring a single modality at a time. Therefore it is useful to develop methods that are capable of integrating measurements of the same biological system but obtained using different technologies.

Here the goal is to learn a latent space where observations from the same cell acquired using different modalities. A perfect result has each of the paired observations sharing the same coordinates in the latent space.

## The metrics
Metrics for multimodal data integration aim to characterize how well the aligned datasets correspond to the ground truth.

* **kNN AUC**: Let $f(i) ∈ F$ be the scRNA-seq measurement of cell $i$, and $g(i) ∈ G$ be the scATAC- seq measurement of cell $i$. kNN-AUC calculates the average percentage overlap of neighborhoods of $f(i)$ in $F$ with neighborhoods of $g(i)$ in $G$. Higher is better.
* **MSE**: Mean squared error (MSE) is the average distance between each pair of matched observations of the same cell in the learned latent space. Lower is better.


## The results

{{< task_tables task="multimodal_data_integration" >}}
