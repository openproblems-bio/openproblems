+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 20  # Order that this section will appear.

title = "Results"
subtitle = ""

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

|    | metric   |      value | method     | task                        | dataset             |
|---:|:---------|-----------:|:-----------|:----------------------------|:--------------------|
|  0 | accuracy | 0.505859   | dummy      | Label Projection            | dummy               |
|  0 | accuracy | 1          | cheat      | Label Projection            | dummy               |
|  1 | dummy    | 0          | cheat      | Label Projection            | dummy               |
|  1 | dummy    | 0          | dummy      | Label Projection            | dummy               |
|  0 | knn_auc  | 0.00298303 | procrustes | Multimodal Data Integration | scicar_cell_lines   |
|  0 | knn_auc  | 0.0587961  | cheat      | Multimodal Data Integration | scicar_cell_lines   |
|  0 | knn_auc  | 0.062768   | mnn        | Multimodal Data Integration | scicar_cell_lines   |
|  1 | mse      | 0          | cheat      | Multimodal Data Integration | scicar_cell_lines   |
|  1 | mse      | 0.833141   | procrustes | Multimodal Data Integration | scicar_cell_lines   |
|  1 | mse      | 1.05657    | mnn        | Multimodal Data Integration | scicar_cell_lines   |
|  0 | knn_auc  | 0.00119093 | mnn        | Multimodal Data Integration | scicar_mouse_kidney |
|  0 | knn_auc  | 0.0223434  | procrustes | Multimodal Data Integration | scicar_mouse_kidney |
|  0 | knn_auc  | 0.0424163  | cheat      | Multimodal Data Integration | scicar_mouse_kidney |
|  1 | mse      | 0          | cheat      | Multimodal Data Integration | scicar_mouse_kidney |
|  1 | mse      | 0.939055   | procrustes | Multimodal Data Integration | scicar_mouse_kidney |
|  1 | mse      | 0.980936   | mnn        | Multimodal Data Integration | scicar_mouse_kidney |
