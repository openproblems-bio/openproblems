+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 20  # Order that this section will appear.

title = ""
subtitle = ""

[design]
  # Choose how many columns the section has. Valid values: 1 or 2.
  columns = "1"

[design.background]
  # Apply a background color, gradient, or image.
  #   Uncomment (by removing `#`) an option to apply it.
  #   Choose a light or dark text color by setting `text_color_light`.
  #   Any HTML color name or Hex value is valid.

  # Background color.
  # color = "navy"

  # Background gradient.
  # gradient_start = "DeepSkyBlue"
  # gradient_end = "SkyBlue"

  # Background image.


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
# Benchmarking formalized challenges in single-cell analysis

The Open Problems effort is inspired by the progress in machine learning driven by researchers aiming to maximize the performance of their algorithms against standardized, well-defined computational tasks (e.g. [Papers with Code State-of-the-art Leaderboards](https://paperswithcode.com/sota)). Several [Grand Challenges](https://doi.org/10.1186/s13059-020-1926-6) have been identified for single cell analysis, but these challenges require formalization before method developers can attempt to solve them. Our goal is to formalize challenges such as these and create a living community-driven state-of-the-art benchmarking platform to facilitate development of single-cell methods.

We identify several traits that allow these benchmarks to serve as Open Problems that drive innovation:
* Tasks are formally defined with a clear mathematical interpretation
* Easily accessible gold-standard datasets are publicly available in a ready-to-go standardized format
* One or more quantitative metrics are defined for each task to judge success
* State-of-the-art methods are ranked in a continuously updated leaderboard

This project is sponsored by the [Chan Zuckerberg Initiative](https://chanzuckerberg.com/science/), and we are looking for broad input from the scientific community.

If you'd like to learn more, click here: <a href="/contributing"><button type="button" class="btn btn-primary btn-lg">**Get Involved**</button></a>
