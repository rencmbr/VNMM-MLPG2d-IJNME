# Truly Cell-Free Vector Nodal Meshless Method (VNMM-MLPG 2D) for Electromagnetic Eigenvalue Problems

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-20%20passed-brightgreen.svg)](tests/)

Official companion codebase and reproducibility package for the research article:

> **"A Truly Cell-Free Vector Nodal Meshless Method via Meshless Local Petrov--Galerkin (VNMM-MLPG) for 2D Electromagnetic Eigenvalue Problems"**  
> *Authors:* Renato C. Mesquita, Williams L. de Nicomedes, Bárbara M. F. Gonçalves, and Luilly A. G. Ortiz  
> *Target Journal:* International Journal for Numerical Methods in Engineering (Wiley)

---

## 📖 Overview

The **Vector Nodal Meshless Method via Meshless Local Petrov--Galerkin (VNMM-MLPG)** is a truly cell-free spatial discretization formulated for source-free two-dimensional electromagnetic eigenvalue problems governed by Maxwell's curl-curl equations in $H(\mathbf{curl})$. By combining directional vector point interpolation ($\mathcal{P}^1$) with weak forms integrated over autonomous circular subdomains $\Omega_s^I$ centered at physical nodes, VNMM-MLPG completely eliminates background integration cells, intra-cell stencil switching, and geometric meshing bottlenecks.

### Key Methodological Pillars:
1. **Truly Cell-Free Integration (Autonomous Disks):** Weak-form integrals are evaluated directly over circular subdomains $\Omega_s^I$ of radius $r_s^I = \alpha_s h_I$ ($\alpha_s = 0.80$), completely decoupling spatial integration from background grids.
2. **Local Galerkin (MLPG6) with Closed-Form Operators:** Because fixed complete linear stencils ($\mathcal{P}^1$) yield constant curl and divergence per subdomain, Stokes' and Gauss' theorems imply that weight test functions vanishing on $\partial\Omega_s^I$ (MLPG1) produce identically zero stiffness, necessitating a Local Galerkin (MLPG6) approach. All local curl, divergence, mass, and stabilization matrices---on both full interior disks and boundary-truncated circular segments---are evaluated in exact closed analytical form.
3. **Overlapping Symmetric Subdomain Galerkin (SSG) Assembly:** Eliminates the structural asymmetry ($K_{Ij} \ne K_{Ji}$) and loss of stability of classical row-by-row MLPG/LPIM by assembling additive local subdomain energies, guaranteeing real spectra and symmetric positive-semidefinite global operators ($\mathbf{K} = \mathbf{K}^T \ge 0$, $\mathbf{M} = \mathbf{M}^T > 0$, strictly $\mathbf{K}_{\mathrm{red}} > 0$ on fixed connected PEC clouds) without requiring lens-intersection numerical cubatures, while its continuous multiplicity-weighted Rayleigh quotient homogenizes at second order on regular rectangular lattices.
4. **VMS-Inspired Linear Orthogonal Energy Stabilization:** Resolves the localized rank-2 deficiency of single-disk $\mathcal{P}^1$ integration by projecting the two curl- and divergence-free harmonic linear gradient modes (cross-derivative and hyperbolic quadrupole) onto directional degrees of freedom, restoring rank 4 on the gradient subspace. Across extended $15\times 15$ to $81\times 81$ refinements, calibrated linear scaling ($\alpha_{\mathrm{stab}} = c_{\mathrm{stab}}\,(h_{\mathrm{char}}/L)$, baseline $\gamma = 1, c_{\mathrm{stab}} = 4.0$) sustains monotonic cutoff decay without the fine-mesh rank re-softening exhibited by quadratic decay ($\gamma = 2$) above $49\times 49$.
5. **Superlinear Div-Curl Regularization:** Shifts gradient nullspace modes outside the operational 10-mode spectrum via resolution-scaled divergence penalization ($s_{\mathrm{div}} = c_{\mathrm{div}}(h_{\mathrm{char}}/L)^{-1.5}$, $\beta = 1.5$, anchored at $s_{\mathrm{div}}^{(0)} = 15.0$ on $21\times 21$), preserving physical cavity resonances without volumetric locking.

---

## 📁 Repository Structure

```text
VNMM-MLPG2d-IJNME/
├── README.md                                   # Complete reproducibility documentation (this file)
├── requirements.txt                            # Minimal Python dependencies (numpy, scipy, matplotlib)
├── LICENSE                                     # MIT Open-Source License
├── run_all.py                                  # Master CLI runner to reproduce tests, figures, and tables
│
├── src/                                        # Core Python library
│   ├── __init__.py
│   ├── cavity_mesh_2d.py                       # Point cloud and cavity mesh generator (regular & jittered)
│   ├── mlpg_subdomains.py                      # Autonomous disks, analytical truncation, and geometric moments
│   ├── mlpg_assembler_2d.py                    # Closed-form local curl, div, mass, and stabilization operators
│   ├── support_nodes_vnmm_2d_6_p1.py           # Scale-invariant 6-node stencil selector via KDTree
│   ├── shape_functions_vnmm_2d_6_p1.py         # Complete linear P1 vector shape functions
│   ├── vnmm_assembler_2d.py                    # Cell-based Gauss integration (Stage 6 comparative baseline)
│   ├── eigen_solver_cavity_2d.py               # Reference cavity eigensolver utilities & analytical modes
│   ├── gauss_quadrature_2d.py                  # Gauss-Legendre quadrature rules
│   └── fem_edge_2d.py                          # Canonical Whitney edge finite element solver
│
├── scripts/                                    # Self-contained reproduction scripts mapped directly to paper
│   ├── generate_mlpg_figures.py                # Generates all 7 publication figures (.pdf & .png)
│   ├── mlpg_patch_test_2d.py                   # Stage 2: Galerkin Patch Test across 7 affine modes (Table 2)
│   ├── mlpg_source_problem_2d.py               # Stage 3: Auxiliary MMS convergence (Tables 3-4) & M(x) study
│   ├── parametric_analysis_alpha_s_mlpg_2d.py  # Stage 4: Subdomain radius alpha_s sensitivity (Fig. 3)
│   ├── compute_alpha_jitter_study.py           # Stage 4/5: alpha_s x jitter x refinement study
│   ├── mlpg_cavity_eigenvalues_2d.py           # Stage 5: 21x21 s_div sweep & baseline spectrum (Tables 5-7)
│   ├── study_fine_mesh_beta_cstab_2d.py        # Stage 5: 8-grid (15..81) coupled (beta, c_stab) & M(x) sweep
│   ├── study_alpha_gamma_and_beta2_2d.py       # Stage 5: Linear vs quadratic gamma & beta in {1.5,2.0,2.5} sweep
│   ├── reproduce_81x81_spectrum_and_runtime.py # Stage 5/6: 81x81 Table 10 spectrum & 49x49 end-to-end runtime
│   ├── compute_eabe_benchmark.py               # Stage 6: Assembly runtime & NNZ vs VNMM-Gauss (Table 11)
│   ├── compute_jitter_benchmark.py             # Stage 5/6: 0% and 25% jitter cavity convergence (15..49)
│   ├── plot_mms_error_distribution_49x49.py    # Auxiliary boundary vs interior MMS error breakdown
│   └── *.json                                  # Precomputed benchmark & parameter-sweep datasets
│
├── tests/                                      # Automated verification test suite (20 tests)
│   ├── test_mlpg_subdomains_2d.py              # Coverage, boundary truncation, and exact geometric moments
│   ├── test_mlpg_assembler_2d.py               # Analytical vs high-order quadrature equivalence (< 10^-13)
│   ├── test_mlpg_patch_test_2d.py              # Affine exactness to machine precision (< 10^-13)
│   ├── test_mlpg_source_problem_2d.py          # Monotonic MMS convergence (L2 quadratic, broken H(curl) linear)
│   └── test_mlpg_cavity_eigenvalues_2d.py      # Cavity spectrum, positive-definiteness, and stabilization
│
└── figures/                                    # Generated publication figures (Figs. 1-7, .pdf & .png)
```

---

## 🚀 Quick Start

### 1. Requirements & Environment Setup

The code requires Python 3.8+ and standard scientific packages (`numpy`, `scipy`, `matplotlib`):

```bash
git clone https://github.com/rencmbr/VNMM-MLPG2d-IJNME.git
cd VNMM-MLPG2d-IJNME
pip install -r requirements.txt
```

### 2. Verify the Installation (Automated Test Suite)

Run the automated test suite verifying analytical moment closures, quadrature equivalence, machine-precision affine patch tests, MMS convergence, and cavity spectral assembly:

```bash
python3 run_all.py --tests
```
*Expected output:* `Ran 20 tests ... OK`.

---

## 📊 Complete Reproduction Guide

Use `run_all.py` to inspect available targets or reproduce any manuscript figure, validation stage, or table:

```bash
python3 run_all.py --list
```

### Paper Figures (Figs. 1–7)

| Paper Item | Description | Direct Command | Output Artifact |
| :--- | :--- | :--- | :--- |
| **Figure 1** | Schematic of MLPG nodal cloud, autonomous integration subdomains, and 6-node stencils | `python3 run_all.py --fig 1` | `figures/fig_mlpg_subdomains_schematic.pdf` |
| **Figure 2** | Method of Manufactured Solutions (MMS) refinement convergence (Dirichlet & Mixed Neumann) | `python3 run_all.py --fig 2` | `figures/fig_mlpg_mms_convergence.pdf` |
| **Figure 3** | Parametric sensitivity of subdomain radius factor $\alpha_s$ (coverage vs. containment) | `python3 run_all.py --fig 3` | `figures/fig_mlpg_alpha_s_sensitivity.pdf` |
| **Figure 4** | Visual comparison of analytical vs. global VNMM-MLPG modal fields ($TE_{10}, TE_{11}, TE_{21}$) on $21\times 21$ and $81\times 81$ regular grids | `python3 run_all.py --fig 4` | `figures/fig_mlpg_modal_fields.pdf` |
| **Figure 5** | Unstabilized cavity cutoff wavenumber error under $0\%$ and $25\%$ nodal jitter ($15\times 15$ to $49\times 49$) | `python3 run_all.py --fig 5` | `figures/fig_mlpg_cavity_jitter.pdf` |
| **Figure 6** | Orthogonal energy stabilization convergence across 8 grids ($15^2$ to $81^2$): unstabilized vs. $\gamma=2$ vs. $\gamma=1$ | `python3 run_all.py --fig 6` | `figures/fig_mlpg_cavity_stabilization.pdf` |
| **Figure 7** | Computational benchmark comparing linear VNMM-MLPG ($15^2\text{--}81^2$) vs. VNMM-Gauss ($15^2\text{--}49^2$) under $0\%$ and $25\%$ jitter | `python3 run_all.py --fig 7` | `figures/fig_mlpg_benchmark_comparison.pdf` |

To regenerate **all 7 figures** in both vector PDF and 300-DPI PNG formats:
```bash
python3 run_all.py --figures
```

### Validation Stages & Tables

| Paper Stage | Corresponding Tables | Description | Direct Command |
| :--- | :--- | :--- | :--- |
| **Stage 2** | Table 2 | Polynomial reproduction: Galerkin patch test across 7 affine vector modes to machine precision ($\sim 10^{-14}\text{--}10^{-15}$) | `python3 run_all.py --stage 2` |
| **Stage 3** | Tables 3, 4 | Auxiliary local-equilibrium MMS verification (Dirichlet & Mixed Neumann; SSG vs. unsymmetric row assembly) + Sec. 3.5.7 $M(\mathbf{x})$ check | `python3 run_all.py --stage 3` |
| **Stage 4** | Fig. 3 | Subdomain geometry, coverage, and containment bounds ($\alpha_s = 0.80$ within $[0.70, 1.20]$) under $0\%$ and $25\%$ jitter | `python3 run_all.py --stage 4` |
| **Stage 5 (Baseline)** | Tables 5–7 | $21\times 21$ $s_{\mathrm{div}}$ sweep, 10-mode spectrum at $s_{\mathrm{div}} = 15.0$, and fixed vs. scaled $s_{\mathrm{div}}$ refinement ($15^2\text{--}49^2$) | `python3 run_all.py --stage 5` |
| **Stage 5 (Fine Sweeps)** | Tables 8, 9 | 8-grid ($15^2\text{--}81^2$) coupled $(\beta, \gamma, c_{\mathrm{stab}})$ parameter sweeps & $M(\mathbf{x})$ Rayleigh shifts (regenerates precomputed JSONs) | `python3 run_all.py --stage 5f` |
| **Stage 5/6 ($81^2$ & $49^2$)** | Table 10 & Sec. 3.6 | Full 10-mode $81\times 81$ spectrum ($0\%$ & $25\%$ jitter) and paired $49\times 49$ end-to-end runtime & $\mathrm{cond}_2$ vs. VNMM-Gauss | `python3 run_all.py --stage 5s` |
| **Stage 6 (Benchmark)** | Tables 11, 12 | Matrix assembly speedup ($2.0\times\text{--}3.7\times$), stiffness NNZ sparsity ($-20\%\text{--}24\%$), and cutoff accuracy vs. VNMM-Gauss | `python3 run_all.py --stage 6` |

---

## 🔬 Core Numerical Results

* **Fine-Mesh Spectral Cavity Convergence ($15\times 15$ to $81\times 81$):** Pairing linear orthogonal stabilization ($\gamma = 1, c_{\mathrm{stab}} = 4.0$, giving $\alpha_{\mathrm{stab}} = 0.050$ at $81\times 81$) with $\beta = 1.5$ ($s_{\mathrm{div}} = 120.00$ at $81\times 81$) eliminates spurious modes across all 240 tabulated cavity runs and all 654 parameter-grid evaluations, achieving monotonic cutoff wavenumber decay down to **$0.78\%$** on regular grids ($p_{15\text{--}81} = 1.72$, accelerating to $p_{36\text{--}81} = 2.48$) and **$0.79\%$** under $25\%$ jitter ($p_{15\text{--}81} = 1.61$, $p_{36\text{--}81} = 1.62$), with discrete divergence ratio $\max R_{\mathrm{div}} \le 0.06\%$ at $81\times 81$.
* **Modal Field Reconstruction:** On $21\times 21$ and $81\times 81$ regular grids ($0\%$ jitter), multiplicity-weighted global reconstructions $\vec{E}_{\mathrm{global}}$ for $TE_{10}$, $TE_{11}$, and $TE_{21}$ cleanly recover analytical standing-wave lobes, nodal lines, tangential PEC boundary behavior, and vector field orientations without checkerboard artifacts (Fig. 4).
* **Spatial Multiplicity Homogenization:** On regular rectangular lattices, the continuous multiplicity-weighted Rayleigh quotient shift decays at rate $\mathcal{O}(h_{\mathrm{nom}}^{1.89})$ (down to $0.024\%$ on $81\times 81$), whereas $94.4\%$ coverage under $25\%$ jitter breaks periodic cell cancellation and leaves a slowly decaying geometric floor ($0.63\%\text{--}1.33\%$ across $49^2\text{--}81^2$, explaining the $0.79\%$ fine-mesh plateau).
* **Algebraic Symmetry & Fixed-Mesh Definiteness:** SSG assembly avoids the complex spectra and divergence of unsymmetric row-by-row MLPG, ensuring $\mathbf{K} = \mathbf{K}^T \ge 0$, $\mathbf{M} = \mathbf{M}^T > 0$, $\dim\ker(\mathbf{K}) = 2$, and strict positive-definiteness $\mathbf{K}_{\mathrm{red}} > 0$ after PEC elimination---verified across all 16 benchmark clouds ($15^2\text{--}81^2$, $0\%$ & $25\%$ jitter) with full-rank spanning connectivity ($\sigma_2 \ge 1.0 \times 10^{-3}$, MST bottleneck $\ge 0.395$).
* **Auxiliary Local Verification:** Passes all 7 affine patch tests to machine precision ($\sim 10^{-14}\text{--}10^{-15}$) and attains $\mathcal{O}(h_{\mathrm{nom}}^{1.93})$ broken $L^2$ and $\mathcal{O}(h_{\mathrm{nom}}^{1.05})$ broken $H(\mathbf{curl})$ rates under exact local weak equilibrium.
* **Computational Efficiency:** Closed-form single-disk assembly cuts matrix assembly runtime by $2.0\times\text{--}3.7\times$ across $15^2\text{--}49^2$, achieves **$3.35\times\text{--}3.94\times$ end-to-end (assembly + Shift-and-Invert Arnoldi eigensolve) speedup** at $49\times 49$, and reduces stiffness matrix non-zeros by $20\%\text{--}24\%$ relative to calibrated $3\times 3$ cell-based VNMM-Gauss.

---

## 📜 Citation

If you find this codebase or methodology useful in your research, please cite our article:

```bibtex
@article{mesquita2026vnmm_mlpg,
  author  = {Mesquita, Renato C. and de Nicomedes, Williams L. and Gon{\c{c}}alves, B{\'a}rbara M. F. and Ortiz, Luilly A. G.},
  title   = {A Truly Cell-Free Vector Nodal Meshless Method via Meshless Local Petrov--Galerkin ({VNMM-MLPG}) for {2D} Electromagnetic Eigenvalue Problems},
  journal = {International Journal for Numerical Methods in Engineering},
  year    = {2026},
  note    = {Under Review}
}
```

---

## 📄 License

This project is licensed under the terms of the [MIT License](LICENSE).
