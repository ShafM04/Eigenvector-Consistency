-----

# Implementing and Assessing the Inverse Problem in 3D

This repository contains the code and documentation for the MEng thesis project by **Shafayat Mustafa**, supervised by **Dr. Benjamin Cameron** at the University of Southampton.

The project focuses on implementing and validating a 3D inverse problem methodology to determine a material's stress field from a known strain field, without requiring prior knowledge of its material properties.

-----

## Table of Contents

  * Overview
  * The Core Problem: Forward vs. Inverse
  * Project Methodology
      * 1.  Data Generation (Forward Problem)
      * 2.  Data Processing: Eigenvector Field Correction
      * 3.  Solving for Stress (Inverse Problem)
  * Key Results & Findings
  * Repository Contents
  * How to Use the Code
      * Prerequisites
      * Input Data Format
      * Execution
      * Output
  * Future Work

-----

## Overview

In engineering analysis, the Finite Element Method (FEM) is traditionally used to solve the "forward problem": given known material properties, geometry, and loads, we compute the resulting stress and strain. However, this approach is challenging when material properties are unknown, heterogeneous, or difficult to measure, which is common in complex composites or biological tissues.

This project tackles the **"inverse problem"**: using a measured strain field as the input, we compute the corresponding stress field by solving a system of Partial Differential Equations (PDEs). This is particularly valuable for characterizing novel materials and for structural health monitoring where only surface deformations might be available.

*Figure 1: Visual representation of the inverse problem. Measured strain data is used to compute the stress field.*

This study implements the full inverse problem workflow in 3D, from processing raw strain data to solving for and validating the final stress tensor field.

## The Core Problem: Forward vs. Inverse

| Feature | Forward Problem (Traditional FEM) | Inverse Problem (This Project) |
| :--- | :--- | :--- |
| **Known Inputs** | Material Properties (E, ν), Geometry, Loads | Measured Strain Field, Boundary Conditions |
| **Unknown Outputs** | Stress, Strain, Displacement | Stress Field, Material Properties (Implicit) |
| **Governing Principle** | Solve for system response based on known physics. | Infer internal characteristics from observed response. |
| **Challenge** | Requires accurate material property data, which can be a major limitation. | Requires high-quality, full-field strain data (e.g., from DVC) and robust algorithms to handle noisy or ambiguous vector data. |

## Project Methodology

The workflow is divided into three main stages:

### 1\. Data Generation (Forward Problem)

To create a clean, noise-free dataset for validation, an initial forward problem simulation was run in COMSOL Multiphysics.

  * A simple cuboid and a more complex geometry with curved surfaces were modeled.
  * Defined loads and boundary conditions were applied to induce a stress-strain state.
  * The resulting full-field strain tensor data ($\\epsilon\_{xx}, \\epsilon\_{xy}, \\epsilon\_{xz}, \\epsilon\_{yy}, \\epsilon\_{yz}, \\epsilon\_{zz}$) at each mesh node was exported to a text file. This simulated data serves as the input for the inverse problem, acting as a stand-in for experimental data from methods like Digital Volume Correlation (DVC).

### 2\. Data Processing: Eigenvector Field Correction

This is the core of the Python script (`Complexbonecode.py`) in this repository. The exported strain field must be processed to extract the principal strain directions (eigenvectors), which are essential for the inverse solver.

A significant challenge is that the numerical computation of eigenvectors has an inherent sign ambiguity (a vector can point in either direction along the same axis). A field of randomly-flipped eigenvectors is discontinuous and cannot be used. The script corrects this using a custom **Eigenvector Consistency Algorithm**:

1.  **Eigenvalue Decomposition:** For each point in the 3D grid, the 3x3 strain tensor matrix is decomposed to find its eigenvalues (principal strain magnitudes) and eigenvectors (principal strain directions $q\_1, q\_2, q\_3$).
2.  **Consistency Propagation:** A **Breadth-First Search (BFS)** algorithm systematically traverses the 3D grid, starting from the origin.
3.  **Orientation Correction:** For each new point, its eigenvectors are compared to the already-corrected eigenvectors of its neighbor. The orientation that yields the largest dot product (i.e., is most aligned) is chosen. This ensures a smooth, continuous vector field.
4.  **Orthogonality:** The right-hand rule is enforced using the cross-product ($q\_3 = q\_1 \\times q\_2$) to maintain a consistent, orthogonal coordinate system at every point.
5.  **Output:** The corrected, continuous eigenvector field is saved to a text file for use in the next stage.

*Figure 2: Visualization of the corrected eigenvector field after processing, ensuring all vectors are consistently oriented.*

### 3\. Solving for Stress (Inverse Problem)

The corrected eigenvector field is imported back into COMSOL to solve for the stress tensor.

  * A "Wave Form PDE" interface is used to solve the equilibrium equation, which simplifies to $\\nabla\\cdot\\sigma = 0$ for a static problem.
  * The stress tensor $\\sigma$ is cleverly related to a new unknown field variable $u$ and the known eigenvector field $q$ via the expression:
    $$\sigma=\sum_{i=1}^{3}(u\cdot q_{i})q_{i}q_{i}^{T}$$
    This formulation allows the PDE system to be solved for $u$ using the known principal directions from the Python script.
  * The final stress tensor components are calculated from the solved field $u$ and compared against the "ground truth" results from the initial forward problem for validation.

## Key Results & Findings

  * **Success in Simple Geometry:** For a simple cuboid, the inverse problem accurately reproduced the stress fields for components governed by the applied boundary conditions ($\\sigma\_{xx}, \\sigma\_{xy}, \\sigma\_{yy}$) with a percentage uncertainty of less than 0.5%.

    *Figure 3: Comparison of the $\\sigma\_{xx}$ stress component from the (a) forward problem and (b) inverse problem, showing excellent agreement.*

  * **Limitation in Unconstrained Directions:** The method struggled to accurately calculate the stress component in the unconstrained z-direction ($\\sigma\_{zz}$). This highlights that the inverse problem's accuracy is highly dependent on having sufficient boundary conditions to constrain the solution in all directions.

  * **Challenges with Complexity:** For the more complex geometry, inaccuracies were more pronounced, particularly near curved surfaces and where boundary effects converge. However, some components still showed strong correlation, proving the method's potential.

  * **Computational Efficiency:** The Python-based data processing stage exhibits **linear time complexity**, meaning its runtime scales predictably with the number of nodes in the mesh. This is a positive indicator for scaling the method to larger, more detailed models.

    *Figure 4: The eigenvector processing algorithm shows a linear relationship between runtime and the number of nodes.*

## Repository Contents

  * `README.md`: This project overview.
  * `Complexbonecode.py`: The Python script for processing strain data and enforcing eigenvector consistency.
  * `FEEG3003 IP 33509069 report (1).pdf`: The full MEng research report this project is based on.
  * **Sample Data (Not included):** The script requires an input file (e.g., `3D_strain_output_v1.txt`).

## How to Use the Code

The `Complexbonecode.py` script performs the eigenvector consistency correction on a 3D strain field.

### Prerequisites

You will need Python 3 and the following libraries:

  * NumPy
  * Matplotlib

You can install them using pip:

```bash
pip install numpy matplotlib
```

### Input Data Format

The script requires a space-delimited text file as input. This file should not contain headers, but each line must represent a node in the mesh with 9 columns:
`x y z e_xx e_xy e_xz e_yy e_yz e_zz`

For example:

```
0.0 0.0 0.0 0.0012 0.0001 -0.0002 0.0011 0.0003 -0.0021
0.1 0.0 0.0 0.0013 0.0002 -0.0003 0.0012 0.0004 -0.0022
...
```

You will need to update the filename inside the script:

```python
# In Complexbonecode.py
with open('YOUR_DATA_FILE.txt', 'r') as f:
    # ...
```

### Execution

To run the script, navigate to the repository directory in your terminal and execute:

```bash
python Complexbonecode.py
```

### Output

The script will produce two outputs:

1.  **A text file `eigenvectors_bfs_output1.1.txt`**: This file contains the corrected eigenvector components for each node, ready to be imported into a PDE solver like COMSOL. The format is: `X Y Z Q00 Q01 Q02 Q10 Q11 Q12 Q20 Q21 Q22`.
2.  **A Matplotlib 3D plot**: This visualizes a subsample of the corrected eigenvectors, allowing for a quick visual check of their consistency.

## Future Work

Based on the project's findings, future research could focus on:

  * **Improving `σzz` Accuracy:** Investigate methods to better constrain the problem in the z-direction, a tentially by incorporating additional physics or boundary conditions.
  * **Enhancing Robustness for Complex Geometries:** Refine the PDE formulation or boundary condition handling to improve accuracy for models with significant geometric complexity.
  * **Performance Optimization:** While the current algorithm is efficient, explore machine learning approaches to accelerate the simulation and data processing pipeline even further.
