#!/usr/bin/env python3
"""
reproduce_81x81_spectrum_and_runtime.py
=======================================
Reproduces two key IJNME manuscript validations for the calibrated linear
Orthogonal Energy Stabilization baseline (beta = 1.5, gamma = 1, c_stab = 4.0):

1. Table 10 (Section 3.5.6): Full 10-mode physical PEC cavity spectrum at 81x81
   (6,241 active DoFs, h_nom = 0.0393 m, s_div = 120.00, alpha_stab = 0.050,
   alpha_s = 0.80, seed = 42) under both 0% and 25% nodal jitter.
2. Stage 6 End-to-End Wall-Clock Benchmark (Section 3.6): Assembly, sparse
   Shift-and-Invert Arnoldi eigensolve (k = 25, sigma = -0.05), conditioning
   cond_2(K_red), and total speedup at 49x49 vs. calibrated VNMM-Gauss (s_div=6.0).
"""

import os
import sys
import time
import numpy as np
from scipy.sparse.linalg import eigsh

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(CURRENT_DIR)
for p in [ROOT_DIR, os.path.join(ROOT_DIR, "src"), CURRENT_DIR]:
    if p not in sys.path:
        sys.path.insert(0, p)

from src.cavity_mesh_2d import generate_cavity_mesh
from src.mlpg_subdomains import build_mlpg_subdomains_2d
from src.mlpg_assembler_2d import assemble_mlpg_matrices_2d
from src.vnmm_assembler_2d import assemble_vnmm_matrices_2d
from src.eigen_solver_cavity_2d import apply_pec_boundary_conditions
from scripts.study_fine_mesh_beta_cstab_2d import compute_sparse_cond2

MODE_LIST = [
    ("TE10", 1, 0),
    ("TE01", 0, 1),
    ("TE11", 1, 1),
    ("TE20", 2, 0),
    ("TE02", 0, 2),
    ("TE21", 2, 1),
    ("TE12", 1, 2),
    ("TE22", 2, 2),
    ("TE30", 3, 0),
    ("TE03", 0, 3),
]
LX = np.pi
SPUR_THRESHOLD = 50.0


def solve_mlpg_spectrum(N: int, jitter_frac: float, alpha_s: float = 0.80, seed: int = 42):
    """Assembles and solves the calibrated linear VNMM-MLPG cavity spectrum."""
    s_div = float(15.0 * ((N - 1) / 20.0) ** 1.5)
    alpha_stab = float(4.0 / (N - 1))

    t0 = time.perf_counter()
    coords, vectors, is_bnd = generate_cavity_mesh(
        Nx=N, Ny=N, Lx=LX, Ly=LX, interior_type="random", jitter_frac=jitter_frac, seed=seed
    )
    subs = build_mlpg_subdomains_2d(
        coords=coords, vectors=vectors, is_boundary=is_bnd, alpha_s=alpha_s
    )
    K_tot, M_mat, _ = assemble_mlpg_matrices_2d(
        subdomains=subs, s_div=s_div, alpha_stab=alpha_stab, assembly_mode="subdomain_galerkin"
    )
    K_red, M_red, _ = apply_pec_boundary_conditions(K_tot, M_mat, is_bnd)
    t_asm = time.perf_counter() - t0

    # Assemble curl-only reference to evaluate R_div exactly
    K_curl, _, _ = assemble_mlpg_matrices_2d(
        subdomains=subs, s_div=0.0, alpha_stab=0.0, assembly_mode="subdomain_galerkin"
    )
    Kc_red, _, _ = apply_pec_boundary_conditions(K_curl, M_mat, is_bnd)
    K_stab_only, _, _ = apply_pec_boundary_conditions(
        assemble_mlpg_matrices_2d(
            subdomains=subs, s_div=0.0, alpha_stab=alpha_stab, assembly_mode="subdomain_galerkin"
        )[0],
        M_mat,
        is_bnd,
    )
    Kd_red = (K_red - K_stab_only) / s_div

    t1 = time.perf_counter()
    vals, vecs = eigsh(K_red, k=25, M=M_red, sigma=-0.05, which="LM", tol=1e-10)
    order = np.argsort(vals)
    vals = np.real(vals[order])
    vecs = np.real(vecs[:, order])
    mask = vals > 0.05
    vals = vals[mask]
    vecs = vecs[:, mask]
    t_sol = time.perf_counter() - t1

    phys = []
    spur_in_band = 0
    for i in range(len(vals)):
        v = vecs[:, i]
        ec = float(v @ (Kc_red @ v))
        ed = float(v @ (Kd_red @ v))
        r_div = (ed / max(ec, 1e-12)) * 100.0
        if i < 10 and r_div > SPUR_THRESHOLD:
            spur_in_band += 1
        if r_div <= SPUR_THRESHOLD and len(phys) < 10:
            phys.append((float(np.sqrt(vals[i])), float(r_div)))

    return {
        "modes": phys,
        "spur_in_band": spur_in_band,
        "s_div": s_div,
        "alpha_stab": alpha_stab,
        "nnz": int(K_red.nnz),
        "dofs": int(K_red.shape[0]),
        "t_asm": t_asm,
        "t_sol": t_sol,
        "t_tot": t_asm + t_sol,
        "K_red": K_red,
    }


def reproduce_table10():
    print("=" * 84)
    print("TABLE 10 REPRODUCTION: Full 10-Mode Spectrum on 81x81 Cloud (6,241 DoFs)")
    print("Parameters: beta=1.5, gamma=1, c_stab=4.0 -> s_div=120.00, alpha_stab=0.050")
    print("=" * 84)
    r0 = solve_mlpg_spectrum(81, 0.0)
    r25 = solve_mlpg_spectrum(81, 0.25)

    print(
        f"{'Mode':<6} {'kc_ref':>8} | {'0% kc':>8} {'Signed':>9} {'R_div':>7} | "
        f"{'25% kc':>8} {'Signed':>9} {'R_div':>7}"
    )
    print("-" * 84)
    abs0, abs25, rd0, rd25 = [], [], [], []
    for idx, (label, n, m) in enumerate(MODE_LIST):
        kc_ref = float(np.hypot(n, m))
        k0, d0 = r0["modes"][idx]
        k25, d25 = r25["modes"][idx]
        e0 = (k0 - kc_ref) / kc_ref * 100.0
        e25 = (k25 - kc_ref) / kc_ref * 100.0
        abs0.append(abs(e0))
        abs25.append(abs(e25))
        rd0.append(d0)
        rd25.append(d25)
        print(
            f"{label:<6} {kc_ref:8.4f} | {k0:8.4f} {e0:+8.2f}% {d0:6.2f}% | "
            f"{k25:8.4f} {e25:+8.2f}% {d25:6.2f}%"
        )
    print("-" * 84)
    print(
        f"Mean Abs Error : {np.mean(abs0):.2f}% (0% jitter)  |  {np.mean(abs25):.2f}% (25% jitter)"
    )
    print(
        f"Max Abs / R_div: {np.max(abs0):.2f}% / {np.max(rd0):.2f}% (0%)  |  "
        f"{np.max(abs25):.2f}% / {np.max(rd25):.2f}% (25%)"
    )
    print(
        f"Spurious in Band (<=10): {r0['spur_in_band']} (0%)  |  {r25['spur_in_band']} (25%)"
    )


def reproduce_49x49_end_to_end_runtime():
    print("\n" + "=" * 84)
    print("STAGE 6 BENCHMARK: 49x49 End-to-End Runtime & Conditioning (k=25, sigma=-0.05)")
    print("=" * 84)
    N = 49
    for jit in [0.0, 0.25]:
        coords, vectors, is_bnd = generate_cavity_mesh(
            Nx=N, Ny=N, Lx=LX, Ly=LX, interior_type="random", jitter_frac=jit, seed=42
        )

        t0 = time.perf_counter()
        Kg, Mg = assemble_vnmm_matrices_2d(
            coords=coords, vectors=vectors, basis="P1", s_div=6.0,
            Ncx=N - 1, Ncy=N - 1, points_per_dir=3, support_mode="gauss_point"
        )
        Kg_red, Mg_red, _ = apply_pec_boundary_conditions(Kg, Mg, is_bnd)
        t_asm_g = time.perf_counter() - t0
        t1 = time.perf_counter()
        _ = eigsh(Kg_red, k=25, M=Mg_red, sigma=-0.05, which="LM", tol=1e-10)
        t_sol_g = time.perf_counter() - t1
        t_tot_g = t_asm_g + t_sol_g
        cond_g = compute_sparse_cond2(Kg_red.tocsr())

        s_div_m = float(15.0 * ((N - 1) / 20.0) ** 1.5)
        alpha_m = float(4.0 / (N - 1))
        t0 = time.perf_counter()
        subs = build_mlpg_subdomains_2d(coords, vectors, is_bnd, alpha_s=0.80)
        Km, Mm, _ = assemble_mlpg_matrices_2d(
            subs, s_div=s_div_m, alpha_stab=alpha_m, assembly_mode="subdomain_galerkin"
        )
        Km_red, Mm_red, _ = apply_pec_boundary_conditions(Km, Mm, is_bnd)
        t_asm_m = time.perf_counter() - t0
        t1 = time.perf_counter()
        _ = eigsh(Km_red, k=25, M=Mm_red, sigma=-0.05, which="LM", tol=1e-10)
        t_sol_m = time.perf_counter() - t1
        t_tot_m = t_asm_m + t_sol_m
        cond_m = compute_sparse_cond2(Km_red.tocsr())

        print(f"\n[Jitter = {int(jit * 100)}%]")
        print(
            f"  VNMM-Gauss : t_asm={t_asm_g:.3f}s, t_sol={t_sol_g:.3f}s, t_tot={t_tot_g:.3f}s, "
            f"NNZ={Kg_red.nnz:,}, cond2={cond_g:.2e}"
        )
        print(
            f"  VNMM-MLPG  : t_asm={t_asm_m:.3f}s, t_sol={t_sol_m:.3f}s, t_tot={t_tot_m:.3f}s, "
            f"NNZ={Km_red.nnz:,}, cond2={cond_m:.2e}"
        )
        print(
            f"  Speedup    : assembly {t_asm_g / t_asm_m:.2f}x, total {t_tot_g / t_tot_m:.2f}x, "
            f"NNZ diff {(Km_red.nnz - Kg_red.nnz) / Kg_red.nnz * 100.0:+.1f}%"
        )


def main():
    reproduce_table10()
    reproduce_49x49_end_to_end_runtime()


if __name__ == "__main__":
    main()
