# 3D Inverse Problem for Stress Analysis

This repository contains the Python code for the dissertation titled "Implementing and assessing the validity of the inverse problem in 3D". The code implements an inverse problem approach to determine the stress field in a 3D structure from a given strain field.

The methodology involves two main stages:
1.  **Strain Data Processing**: This stage uses Python to process the input strain data, compute eigenvectors, and ensure their consistency across the structure. An algorithm, specifically a breadth-first search (BFS), is used to handle eigenvector alignment, even in complex geometries with disconnected regions.
2.  **Stress Equation Formation**: This stage utilizes COMSOL Multiphysics to solve a set of Partial Differential Equations (PDEs) that relate the processed eigenvector field to the unknown stress field.

The project validates this approach by comparing the stress tensor plots and uncertainty values from the inverse problem with those from a traditional forward Finite Element Method (FEM) simulation.

## **Installation**

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/your-username/your-repo-name.git](https://github.com/your-username/your-repo-name.git)
    ```
2.  **Install the required Python packages:**
    ```bash
    pip install numpy matplotlib
    ```
3.  **Software Requirements:**
    * A working installation of COMSOL Multiphysics is required to run the PDE simulations for both the forward and inverse problems.

## **Usage**

To run the data processing part of the project, navigate to the `src` directory and execute the `main.py` script:

```bash
python src/main.py
