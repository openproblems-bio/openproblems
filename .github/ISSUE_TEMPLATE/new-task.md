---
name: New task
about: Start creating a new benchmarking task in OpenProblems
title: "[New task] ..."
labels: New task
assignees: ''

---

**Task motivation**
Explain the motivation behind your proposed task. Describe the biological or computational problem you aim to address and why it’s important. Discuss the current state of research in this area and any gaps or challenges that your task could help address. This section should convince readers of the significance and relevance of your task.

**Task description**
Provide a clear and concise description of your task, detailing the specific problem it aims to solve. Outline the input data types, the expected output, and any assumptions or constraints. Be sure to explain any terminology or concepts that are essential for understanding the task.

**Proposed ground-truth in datasets**
Describe the datasets you plan to use for your task. OpenProblems offers a standard set of datasets (See [“Common datasets”](https://openproblems.bio/documentation/reference/openproblems-v2/src-datasets.html)) which you can peruse through.

Explain how these datasets will provide the ground-truth for evaluating the methods implemented in your task. If possible, include references or links to the datasets to facilitate reproducibility.

**Initial set of methods to implement**
List the initial set of methods you plan to implement for your task. Briefly describe each method’s core ideas and algorithms, and explain why you think they are suitable for your task.

Consider including both established and cutting-edge methods to provide a comprehensive benchmarking of the state-of-the-art.

**Proposed control methods**
Outline the control methods you propose for your task. These methods serve as a starting point to test the relative accuracy of new methods in the task and as quality control for the defined metrics.

Include both positive controls, which are methods with known outcomes resulting in the best possible metric values, and negative controls, which are simple, naive, or random methods that do not rely on sophisticated techniques or domain knowledge. Explain the rationale for your chosen controls.

**Proposed Metrics**
Describe the metrics you propose for evaluating the performance of methods in your task. Explain the rationale for selecting these metrics and how they will accurately assess the methods’ success in addressing the task’s challenges. Consider including multiple metrics to capture different aspects of method performance.
