# Truly Cell-Free Vector Nodal Meshless Method (VNMM-MLPG 2D) for $H(\mathbf{curl})$ Electromagnetics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-20%20passed-brightgreen.svg)](tests/)

Official companion codebase and reproducibility package for the research article:

> **"A Truly Cell-Free Vector Nodal Meshless Method via Meshless Local Petrov–Galerkin (VNMM-MLPG) for 2D Electromagnetics"**  
> *Target Journal:* International Journal for Numerical Methods in Engineering (Wiley)

---

## 📖 Overview

The **Vector Nodal Meshless Method via Meshless Local Petrov–Galerkin (VNMM-MLPG)** is a truly cell-free spatial discretization paradigm for electromagnetic vector fields governed by Maxwell's curl-curl equations in $H(\mathbf{curl})$. By combining directional vector point interpolation with weak forms integrated over autonomous circular subdomains $\Omega_s^I$ centered at physical nodes, VNMM-MLPG completely eliminates background integration cells, intra-cell stencil switching, and geometric meshing bottlenecks.

### Key Methodological Pillars:
1. **Truly Cell-Free Integration (Autonomous Disks):** Weak-form integrals are evaluated directly over circular subdomains $\Omega_s^I$ of radius $r_s^I = \alpha_s h_I$ ($\alpha_s = 0.80$), completely decoupling spatial integration from background grids.
2. **Local Galerkin (MLPG6) with Closed-Form Operators:** Because fixed complete linear stencils ($\mathcal{P}^1$) yield constant curl and divergence per subdomain, Stokes' and Gauss' theorems imply that weight test functions vanishing on $\partial\Omega_s^I$ (MLPG1) produce identically zero stiffness, necessitating a Local Galerkin (MLPG6) approach. All local matrices are evaluated in exact closed analytical form.
3. **Overlapping Symmetric Subdomain Galerkin (SSG) Assembly:** Eliminates the structural asymmetry ($K_{Ij} \ne K_{Ji}$) and discrete loss of coercivity of classical row-by-row MLPG/LPIM by assembling additive local subdomain energies, guaranteeing strictly symmetric positive-semidefinite global operators ($\mathbf{K} = \mathbf{K}^T \ge 0, \mathbf{M} = \mathbf{M}^T > 0$) without requiring lens-intersection numerical cubatures.
4. **VMS-Inspired Orthogonal Energy Stabilization:** Resolves the localized rank-2 deficiency of single-disk integration by projecting the two harmonic zero-energy modes (pure shear and hyperbolic) onto directional degrees of freedom, restoring full rank 4 and eliminating mode softening via quadratic subgrid scaling ($\alpha_{\mathrm{stab}} \propto h_{\mathrm{nom}}^2$).
5. **Superlinear Div-Curl Regularization:** Shifts gradient nullspace modes into the high-frequency spectrum via superlinear divergence penalty scaling ($s_{\mathrm{div}} \propto h_{\mathrm{nom}}^{-1.5}$), achieving spectral purity without volumeric locking.

---

## 📁 Repository Structure

```text
VNMM-MLPG2d-IJNME/
├── README.md               # Complete reproducibility documentation (this file)
├── requirements.txt        # Minimal Python dependencies (numpy, scipy, matplotlib)
├── LICENSE                 # MIT Open-Source License
├── run_all.py              # Master CLI runner to reproduce any or all results
│
├── src/                    # Core Python library
│   ├── cavity_mesh_2d.py           # Point cloud and cavity mesh generator (Cartesian & perturbed jitter)
│   ├── mlpg_subdomains.py          # Local subdomains, weight functions, and analytical truncation
│   ├── mlpg_assembler_2d.py        # Exact closed-form local stiffness, mass, and stabilization matrices
│   ├── support_nodes_vnmm_2d_6_p1.py # Scale-invariant 6-node stencil selector via KDTree
│   ├── shape_functions_vnmm_2d_6_p1.py # Complete linear P1 vector shape functions
│   ├── vnmm_assembler_2d.py        # Cell-based Gauss integration (for Stage 6 comparative benchmark)
│   ├── eigen_solver_cavity_2d.py   # Shift-and-invert Lanczos eigensolver for reference Gauss
│   ├── gauss_quadrature_2d.py      # Gauss-Legendre quadrature rules
│   └── fem_edge_2d.py              # Canonical Whitney edge finite element solver
│
├── scripts/                # Self-contained reproduction scripts mapped directly to paper
│   ├── generate_mlpg_figures.py          # Generates all publication-quality vector figures (.pdf & .png)
│   ├── mlpg_patch_test_2d.py             # Stage 2: Galerkin Patch Test across 7 affine vector modes
│   ├── mlpg_source_problem_2d.py         # Stage 3: MMS Dirichlet & Mixed Neumann convergence
│   ├── parametric_analysis_alpha_s_mlpg_2d.py # Stage 4: Subdomain radius factor alpha_s sensitivity
│   ├── compute_alpha_jitter_study.py     # Stage 4: Jitter sensitivity study for alpha_s
│   ├── mlpg_cavity_eigenvalues_2d.py     # Stage 5: Cavity eigenmodes, s_div scaling & c_stab sweep
│   ├── compute_eabe_benchmark.py         # Stage 6: Benchmark runtime and NNZ (VNMM-MLPG vs Gauss)
│   ├── compute_jitter_benchmark.py       # Stage 6: Cutoff convergence under 0% and 25% nodal jitter
│   ├── plot_mms_error_distribution_49x49.py # Auxiliary MMS error distribution analysis
│   └── *.json                            # Precomputed benchmark datasets for instant reproduction
│
├── tests/                  # Automated verification test suite (20 tests)
│   ├── test_mlpg_subdomains_2d.py        # Tests coverage, boundary truncation, and geometric moments
│   ├── test_mlpg_assembler_2d.py         # Tests analytical vs numerical quadrature equivalence (< 10^-13)
│   ├── test_mlpg_patch_test_2d.py        # Tests affine exactness to machine precision (< 10^-13)
│   ├── test_mlpg_source_problem_2d.py    # Tests monotonic MMS convergence (L2 quadratic, broken H(curl) linear)
│   └── test_mlpg_cavity_eigenvalues_2d.py# Tests cavity spectrum, positive-definiteness, and stabilization
│
└── figures/                # Output directory where generated publication figures are saved (.pdf & .png)
```

---

## 🚀 Quick Start

### 1. Requirements & Environment Setup

The code requires Python 3.8+ and standard scientific packages. Install dependencies via:

```bash
git clone https://github.com/rencmbr/VNMM-MLPG2d-IJNME.git
cd VNMM-MLPG2d-IJNME
pip install -r requirements.txt
```

### 2. Verify the Installation (Automated Test Suite)

Run the automated test suite verifying machine-precision affine patch tests, analytical moment closures, MMS convergence, and cavity spectral assembly:

```bash
python3 run_all.py --tests
```
*Expected output:* `Ran 20 tests in ~2.4s - OK`.

---

## 📊 Complete Reproduction Guide

You can reproduce any specific figure or validation stage from the article using `run_all.py`:

### Paper Figures

| Paper Item | Description | Direct Command | Output Artifact |
| :--- | :--- | :--- | :--- |
| **Figure 1** | Schematic of MLPG nodal cloud, autonomous integration subdomains, and 6-node stencils | `python3 run_all.py --fig 1` | `figures/fig_mlpg_subdomains_schematic.pdf` |
| **Figure 2** | Method of Manufactured Solutions (MMS) asymptotic refinement convergence (Dirichlet & Neumann) | `python3 run_all.py --fig 2` | `figures/fig_mlpg_mms_convergence.pdf` |
| **Figure 3** | Parametric sensitivity of subdomain radius scaling factor $\alpha_s$ (coverage vs containment) | `python3 run_all.py --fig 3` | `figures/fig_mlpg_alpha_s_sensitivity.pdf` |
| **Figure 4** | Resonant cavity cutoff wavenumber convergence under $25\%$ stochastic nodal jitter | `python3 run_all.py --fig 4` | `figures/fig_mlpg_cavity_jitter.pdf` |
| **Figure 5** | Orthogonal energy stabilization parametric sweep ($c_{\mathrm{stab}}$) and quadratic rate recovery | `python3 run_all.py --fig 5` | `figures/fig_mlpg_cavity_stabilization.pdf` |
| **Figure 6** | Computational benchmark comparing VNMM-MLPG against cell-based VNMM-Gauss ($0\%$ and $25\%$ jitter) | `python3 run_all.py --fig 6` | `figures/fig_mlpg_benchmark_comparison.pdf` |

To generate **all 6 figures** in both vector PDF and high-resolution PNG formats at once:
```bash
python3 run_all.py --figures
```

### Validation Stages & Tables

| Paper Stage | Corresponding Tables | Description | Direct Command |
| :--- | :--- | :--- | :--- |
| **Stage 2** | Sec 3.1 Table | Galerkin patch test across 7 affine vector modes to machine precision ($\sim 10^{-14}$) | `python3 run_all.py --stage 2` |
| **Stage 3** | Tables 2, 3, 4 | MMS spatial refinement ($L^2$ field and broken $H(\mathbf{curl})$ semi-norm) & Unsymmetric vs Symmetric comparison | `python3 run_all.py --stage 3` |
| **Stage 4** | Table 5 | Subdomain scaling parameter sweep ($\alpha_s \in [0.707, 1.15]$) under $0\%$ and $25\%$ jitter | `python3 run_all.py --stage 4` |
| **Stage 5** | Tables 6, 7, 8 | Cavity eigenmodes, superlinear $s_{\mathrm{div}}$ penalty scaling, and VMS orthogonal stabilization sweep | `python3 run_all.py --stage 5` |
| **Stage 6** | Tables 9, 10 | Computational benchmark: matrix assembly speedup ($2.0\times\text{--}3.7\times$), operator sparsity ($-20\%\text{--}24\%$), and asymptotic rates | `python3 run_all.py --stage 6` |

---

## 🔬 Core Numerical Results

* **Machine-Precision Linear Consistency:** Passes all 7 affine patch tests (translations, dilatation, shear, curl, divergence, and general affine field) with errors on the order of $\sim 10^{-14}\text{--}10^{-15}$.
* **Asymptotic Refinement Rates:** Achieves monotonic quadratic convergence $\mathcal{O}(h_{\mathrm{nom}}^{1.93})$ in $L^2(\vec{E})$ and linear convergence $\mathcal{O}(h_{\mathrm{nom}}^{1.05})$ in the broken $H(\mathbf{curl})$ semi-norm under both pure Dirichlet and mixed Dirichlet--Neumann boundary conditions.
* **Algebraic Symmetry & Coercivity:** Eliminates the unphysical complex eigenvalues and divergence of unsymmetric row-by-row MLPG, guaranteeing $\mathbf{K} = \mathbf{K}^T \ge 0$, $\dim \ker(\mathbf{K}) = 2$ (rigid translations), and strict positive-definiteness $\mathbf{K}_{\mathrm{red}} > 0$ after PEC boundary reduction.
* **Mode Softening Elimination:** Restores near-quadratic spectral convergence ($p = 2.15\text{--}2.39$ on regular grids and $p = 1.62$ under $25\%$ jitter) via closed-form VMS orthogonal energy stabilization ($c_{\mathrm{stab}} = 70$), driving finest-grid cavity errors down to $1.29\%$.
* **Computational Speedup:** Accelerates global matrix assembly by $2.0\times\text{--}3.7\times$ over cell-based $3\times 3$ Gauss quadrature while producing $20\%\text{--}24\%$ sparser stiffness operators.

---

## 📜 Citation

If you find this codebase or methodology useful in your research, please cite our article:

```bibtex
@article{mesquita2026vnmm_mlpg,
  author  = {Mesquita, Renato C. and Ortiz, Willian A. and Gonçalves, Cássio F. and Nicomedes, Wellerson P.},
  title   = {A Truly Cell-Free Vector Nodal Meshless Method via Meshless Local Petrov--Galerkin ({VNMM-MLPG}) for {2D} Electromagnetic Problems},
  journal = {International Journal for Numerical Methods in Engineering},
  year    = {2026},
  note    = {Under Review}
}
```

---

## 📄 License

This project is licensed under the terms of the [MIT License](LICENSE).
