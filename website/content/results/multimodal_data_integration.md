+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 20  # Order that this section will appear.

title = "Multimodal Data Integration"
subtitle = "Realigning multimodal measurements of the same cell using various technologies"

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
|    | metric   |    value | method                               | task                        |   memory_mb |   memory_leaked_mb |   runtime_s | dataset             |
|---:|:---------|---------:|:-------------------------------------|:----------------------------|------------:|-------------------:|------------:|:--------------------|
|  0 | knn_auc  | 0.417506 | mnn_log_cpm                          | Multimodal Data Integration |    35.4883  |         3.46875    |   7.68056   | scicar_cell_lines   |
|  0 | knn_auc  | 0.521865 | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |     7.18359 |         1.53125    |   0.216146  | scicar_cell_lines   |
|  0 | knn_auc  | 0.52371  | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     5.72656 |         0.015625   |   0.214223  | scicar_cell_lines   |
|  0 | knn_auc  | 0.627982 | mnn_log_scran_pooling                | Multimodal Data Integration |    10.2695  |         0.332031   |   3.88503   | scicar_cell_lines   |
|  0 | knn_auc  | 0.769706 | procrustes                           | Multimodal Data Integration |     5.87891 |         0.015625   |   0.0487141 | scicar_cell_lines   |
|  1 | mse      | 0.442196 | procrustes                           | Multimodal Data Integration |     5.87891 |         0.015625   |   0.0487141 | scicar_cell_lines   |
|  1 | mse      | 0.974549 | mnn_log_scran_pooling                | Multimodal Data Integration |    10.2695  |         0.332031   |   3.88503   | scicar_cell_lines   |
|  1 | mse      | 0.984945 | mnn_log_cpm                          | Multimodal Data Integration |    35.4883  |         3.46875    |   7.68056   | scicar_cell_lines   |
|  1 | mse      | 0.998209 | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |     7.18359 |         1.53125    |   0.216146  | scicar_cell_lines   |
|  1 | mse      | 1.01319  | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     5.72656 |         0.015625   |   0.214223  | scicar_cell_lines   |
|  0 | knn_auc  | 0.41964  | mnn_log_cpm                          | Multimodal Data Integration |    14.1328  |         0.414062   |   3.90117   | scicar_mouse_kidney |
|  0 | knn_auc  | 0.549842 | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |     6.13672 |         0.132812   |   0.250555  | scicar_mouse_kidney |
|  0 | knn_auc  | 0.552441 | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     5.85156 |         0.00390625 |   0.2584    | scicar_mouse_kidney |
|  0 | knn_auc  | 0.59586  | mnn_log_scran_pooling                | Multimodal Data Integration |    14.1914  |         0.0703125  |   3.90184   | scicar_mouse_kidney |
|  0 | knn_auc  | 0.846586 | procrustes                           | Multimodal Data Integration |     6.05859 |        -0.00390625 |   0.0877079 | scicar_mouse_kidney |
|  1 | mse      | 0.185888 | procrustes                           | Multimodal Data Integration |     6.05859 |        -0.00390625 |   0.0877079 | scicar_mouse_kidney |
|  1 | mse      | 1.00367  | mnn_log_cpm                          | Multimodal Data Integration |    14.1328  |         0.414062   |   3.90117   | scicar_mouse_kidney |
|  1 | mse      | 1.004    | harmonic_alignment_sqrt_cpm          | Multimodal Data Integration |     5.85156 |         0.00390625 |   0.2584    | scicar_mouse_kidney |
|  1 | mse      | 1.00711  | harmonic_alignment_log_scran_pooling | Multimodal Data Integration |     6.13672 |         0.132812   |   0.250555  | scicar_mouse_kidney |
|  1 | mse      | 1.03358  | mnn_log_scran_pooling                | Multimodal Data Integration |    14.1914  |         0.0703125  |   3.90184   | scicar_mouse_kidney |