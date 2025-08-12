# Finite Element Method Applied to 1-Dimensional PDE using Python

## Overview
This project implements the Finite Element Method (FEM) in Python to solve a 1-dimensional partial differential equation with mixed Dirichlet and Neumann boundary conditions.  
We derive the weak form of the PDE, define Lagrange basis functions (quadratic and cubic), assemble the global stiffness matrix and force vector, apply boundary conditions, and compare results to the analytical solution.

The repository includes:
- A step-by-step derivation of the weak form.
- FEM implementation for multiple polynomial orders.
- Visualization of FEM vs. analytical solutions.
- Error analysis as a function of mesh refinement.

---

## Repository Structure
- **`source`**  
  Python scripts for:
  - Weak form discretization  
  - Basis function generation  
  - Matrix assembly  
  - Solving and plotting results  

- **`data`**  
  - Analytical vs. numerical solution data
  - Error metrics for different element sizes
  - Generated plots and result files  

- **`documentation`**  
  - PDF explaining theory and derivations
  - Methodology notes
  - References

---

## Problem Definition
We solve:

$$
-\frac{d}{d} \left( -4 \frac{du}{dx} \right)+8 \frac{du}{dx} +5u+sin x+8cos x, \ \pi<x<\pi
$$


Subject to:
- Neumann BC: $$\frac{du}{dx} = -20, \quad x = -\pi $$
- Dirichlet BC: u = 2, x = $$\pi$$

With analytical solution:

$$
u(x)=2 e^{-x} \left( e^{\pi} \sin \frac{x}{2}+\frac{2(3-e^{2\pi})}{e^{\pi}} \cos \frac{x}{2} \right)+\sin x
$$


---

## Features
- Weak form derivation using integration by parts.
- Support for **quadratic (M=4)** and **cubic (M=3)** basis functions.
- Numerical integration with **Simpson's rule**.
- Global stiffness matrix and force vector assembly.
- Boundary condition enforcement.
- L² error computation.
- Comparison plots for varying element counts.

---

## Results
- **Quadratic basis functions** yield lower error for coarse meshes compared to cubic.
- Increasing the number of elements smooths the approximation and matches the analytical solution more closely.

---

## Getting Started

### **Prerequisites**
Install Python 3.8+ with the following packages:
```bash
pip install numpy sympy scipy matplotlib pandas

