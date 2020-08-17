+++
# A Demo section created with the Blank widget.
# Any elements can be added in the body: https://sourcethemes.com/academic/docs/writing-markdown-latex/
# Add more sections by duplicating this file and customizing to your requirements.

widget = "blank"  # See https://sourcethemes.com/academic/docs/page-builder/
headless = true  # This file represents a page section.
active = true  # Activate this widget? true/false
weight = 20  # Order that this section will appear.

title = "Single Cell Open Problems"
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
### Integrating the Human Cell Atlas

Thanks to recent advances in single cell technologies, we now have the opportunity to profile the state of millions of cells an unprecedented resolution. Using off-the-shelf reagents and protocols, individual laboratories can measure global gene expression, DNA accessibility, and high dimensional proteomic data. These innovations promise to reveal new levels of biological heterogeneity across tissues, developmental stages, and organisms. Is is the goal of [the Human Cell Atlas](https://www.humancellatlas.org/) to push these boundaries and characterize the complete cellular diversity in the trillions of cells in the human body.

Assembling a Human Cell Atlas will require development of [new computational methods](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1926-6) developed by groups spanning computational and biological sciences. To facilitate these collaborations, we are aggregating formalized open problems with standardized datasets and test metrics. This approach is inspired by similar challenges in computer science such as the [ImageNet Challenge](http://www.image-net.org/challenges/LSVRC/) or the [Netflix Prize](https://www.netflixprize.com/) that drive innovation in algorithm development for computationally important tasks. However, with a few notable [exceptions](https://predictioncenter.org/), these challenges do not exist for biological problems.

###  Our goal
We want to make it easier for computer scientists, machine learning experts, and computational biologists to identify important problems worth solving, to access standardized datasets, and to benchmark novel methods against the current state of the art using a common set of test metrics.

### Who can get involved
Participating in the development of the Single Cell Open Problems is open to all labs with an emphasis on those participating in the [Chan Zuckerberg Initiative](https://chanzuckerberg.com/)'s [Seed Networks for the Human Cell Atlas](https://chanzuckerberg.com/science/programs-resources/single-cell-biology/seednetworks/).

### How to get involved

We have broken down the development of the Single Cell Open Problems into **tasks**. A task is a specific quantifiable problems that addresses an open problem in integrating the Human Cell Atlas. An example of a task is [Multimodal Data Integration](results/#multimodal_data_integration) in which the goal is to take a set of multimodal measurements of the same set of cells (e.g. joint single cell RNA and ATAC profiling) and identify which pairs of measurements were acquired from the same cell without using cell barcodes.

Each task is composed of three aspects:

* **Datasets** - a group of high quality curated datasets matched to a task
* **Metrics** - a set of quantitative measures that are used to rank methods
* **Methods** - algorithms contributed by the community to perform the task

At this stage, we are looking for two kinds of contributions:
* **Task developers** - these individuals will help develop tasks for labs to solve
* **Method developers** - these groups will contribute new methods that will get displayed on our leaderboards

We have more information on each of these parts of the Open Problems project in our documentation. This website is backed by an open source [GitHub repository](https://github.com/singlecellopenproblems/SingleCellOpenProblems).
