#!/usr/bin/env python3
"""
Master Reproducibility Runner for:
"A Truly Cell-Free Vector Nodal Meshless Method via Meshless Local Petrov--Galerkin (VNMM-MLPG)
for 2D Electromagnetic Eigenvalue Problems"
International Journal for Numerical Methods in Engineering (Wiley)

Usage:
    python3 run_all.py --help
    python3 run_all.py --list
    python3 run_all.py --all
    python3 run_all.py --tests
    python3 run_all.py --figures
    python3 run_all.py --tables
    python3 run_all.py --fig 4
    python3 run_all.py --stage 5s
"""

import os
import sys
import subprocess
import argparse
import time

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(ROOT_DIR, "scripts")
TESTS_DIR = os.path.join(ROOT_DIR, "tests")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

sys.path.insert(0, ROOT_DIR)
sys.path.insert(0, os.path.join(ROOT_DIR, "src"))
sys.path.insert(0, SCRIPTS_DIR)

FIGURES = {
    "1": {
        "name": "Figure 1",
        "desc": "Schematic of MLPG nodal cloud, autonomous integration subdomains, and 6-node stencils",
        "script": "generate_mlpg_figures.py",
        "output": "figures/fig_mlpg_subdomains_schematic.pdf"
    },
    "2": {
        "name": "Figure 2",
        "desc": "Method of Manufactured Solutions (MMS) refinement convergence (Dirichlet & Mixed Neumann)",
        "script": "generate_mlpg_figures.py",
        "output": "figures/fig_mlpg_mms_convergence.pdf"
    },
    "3": {
        "name": "Figure 3",
        "desc": "Parametric sensitivity of subdomain radius scaling factor alpha_s (coverage vs containment)",
        "script": "generate_mlpg_figures.py",
        "output": "figures/fig_mlpg_alpha_s_sensitivity.pdf"
    },
    "4": {
        "name": "Figure 4",
                "desc": "Visual comparison of analytical vs global VNMM-MLPG modal fields (TE10, TE11, TE21) on 21x21 & 81x81",
        "script": "generate_mlpg_figures.py",
        "output": "figures/fig_mlpg_modal_fields.pdf"
    },
    "5": {
        "name": "Figure 5",
        "desc": "Unstabilized cavity cutoff wavenumber error under 0% and 25% nodal jitter (15x15 to 49x49)",
        "script": "generate_mlpg_figures.py",
        "output": "figures/fig_mlpg_cavity_jitter.pdf"
    },
    "6": {
        "name": "Figure 6",
        "desc": "Orthogonal energy stabilization convergence across 8 grids (15x15 to 81x81): unstabilized vs gamma=2 vs gamma=1",
        "script": "generate_mlpg_figures.py",
        "output": "figures/fig_mlpg_cavity_stabilization.pdf"
    },
    "7": {
        "name": "Figure 7",
        "desc": "Computational benchmark comparing stabilized VNMM-MLPG (15..81) vs VNMM-Gauss (15..49) under 0% and 25% jitter",
        "script": "generate_mlpg_figures.py",
        "output": "figures/fig_mlpg_benchmark_comparison.pdf"
    }
}

STAGES = {
    "2": {
        "name": "Stage 2: Patch Test",
        "desc": "Galerkin patch test across 7 affine vector modes to machine precision (Table 2)",
        "scripts": ["mlpg_patch_test_2d.py"]
    },
    "3": {
        "name": "Stage 3: MMS & Multiplicity M(x)",
        "desc": "Auxiliary local-equilibrium MMS verification (Tables 3, 4) + Sec. 3.5.7 M(x) Rayleigh shift",
        "scripts": ["mlpg_source_problem_2d.py"]
    },
    "4": {
        "name": "Stage 4: Subdomain Sizing",
        "desc": "Theoretical bounds and parametric sweep of alpha_s in [0.707, 1.15] (Fig. 3)",
        "scripts": ["parametric_analysis_alpha_s_mlpg_2d.py"]
    },
    "5": {
        "name": "Stage 5: Baseline Cavity Spectrum",
        "desc": "21x21 s_div sweep, 10-mode spectrum at s_div=15.0, and superlinear s_div scaling (Tables 5, 6, 7)",
        "scripts": ["mlpg_cavity_eigenvalues_2d.py"]
    },
    "5f": {
        "name": "Stage 5 (Fine): Coupled (beta, gamma, c_stab) Sweeps (15..81)",
        "desc": "Full 8-grid stabilization & parameter-space sweeps regenerating Tables 8, 9 JSON datasets",
        "scripts": ["study_fine_mesh_beta_cstab_2d.py", "study_alpha_gamma_and_beta2_2d.py"]
    },
    "5s": {
        "name": "Stage 5/6: 81x81 Spectrum & 49x49 End-to-End Runtime",
        "desc": "Table 10 (81x81 10-mode spectrum, 0% & 25% jitter) and Sec. 3.6 49x49 end-to-end speedup & cond2",
        "scripts": ["reproduce_81x81_spectrum_and_runtime.py"]
    },
    "6": {
        "name": "Stage 6: Assembly & Sparsity Benchmark",
        "desc": "Assembly speedup (2.0x-3.7x) and stiffness NNZ reduction (-20% to -24%) vs VNMM-Gauss (Tables 11, 12)",
        "scripts": ["compute_eabe_benchmark.py"]
    }
}


def run_tests():
    print("=" * 80)
    print("RUNNING AUTOMATED UNIT TESTS (20 MLPG TESTS)")
    print("=" * 80)
    cmd = [sys.executable, "-m", "unittest", "discover", "-s", TESTS_DIR, "-p", "test_*.py", "-v"]
    res = subprocess.run(cmd)
    return res.returncode == 0


def generate_all_figures():
    print("=" * 80)
    print("GENERATING ALL 7 PUBLICATION FIGURES (.PDF & .PNG)")
    print("=" * 80)
    fig_script = os.path.join(SCRIPTS_DIR, "generate_mlpg_figures.py")
    res = subprocess.run([sys.executable, fig_script])
    return res.returncode == 0


def run_single_figure(fig_key):
    if fig_key not in FIGURES:
        print(f"Error: Figure '{fig_key}' not found. Available: {list(FIGURES.keys())}")
        return False
    item = FIGURES[fig_key]
    print(f"Generating {item['name']}: {item['desc']}...")
    fig_script = os.path.join(SCRIPTS_DIR, item['script'])
    import importlib.util
    spec = importlib.util.spec_from_file_location("gen_fig", fig_script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    func_map = {
        "1": mod.plot_fig1_schematic,
        "2": mod.plot_fig2_mms_convergence,
        "3": mod.plot_fig3_alpha_s_sensitivity,
        "4": mod.plot_fig4_modal_fields,
        "5": mod.plot_fig5_cavity_jitter,
        "6": mod.plot_fig6_cavity_stabilization,
        "7": mod.plot_fig7_benchmark_comparison,
    }
    func_map[fig_key]()
    print(f"Saved: {os.path.join(ROOT_DIR, item['output'])}")
    return True


def run_stage(stage_key):
    if stage_key not in STAGES:
        print(f"Error: Stage '{stage_key}' not found. Available: {list(STAGES.keys())}")
        return False
    item = STAGES[stage_key]
    print("=" * 80)
    print(f"RUNNING {item['name'].upper()}: {item['desc']}")
    print("=" * 80)
    for script_name in item['scripts']:
        script_path = os.path.join(SCRIPTS_DIR, script_name)
        res = subprocess.run([sys.executable, script_path])
        if res.returncode != 0:
            return False
    return True


def list_items():
    print("=" * 80)
    print("REPRODUCIBILITY CATALOG FOR VNMM-MLPG (IJNME)")
    print("=" * 80)
    print("\n[Figures]")
    for k, v in sorted(FIGURES.items(), key=lambda kv: int(kv[0])):
        print(f"  --fig {k:1s} : {v['name']:15s} -> {v['desc']}")
        print(f"           Artifact: {v['output']}")
    print("\n[Validation Stages & Tables]")
    for k in ["2", "3", "4", "5", "5f", "5s", "6"]:
        v = STAGES[k]
        print(f"  --stage {k:2s} : {v['name']} -> {v['desc']}")
        print(f"              Scripts: {', '.join('scripts/' + s for s in v['scripts'])}")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description="Master reproducibility runner for the VNMM-MLPG IJNME paper."
    )
    parser.add_argument("--all", action="store_true", help="Run tests, generate all 7 figures, and run core stages")
    parser.add_argument("--tests", action="store_true", help="Run automated unit test suite (20 tests)")
    parser.add_argument("--figures", action="store_true", help="Generate all 7 publication-quality figures (.pdf & .png)")
    parser.add_argument("--fig", type=str, help="Generate a specific figure (--fig 1 .. --fig 7)")
    parser.add_argument("--tables", action="store_true", help="Run core validation stages (2, 3, 4, 5, 5s, 6); use --stage 5f to regenerate fine-mesh JSON sweeps")
    parser.add_argument("--stage", type=str, help="Run a specific validation stage (--stage 2, 3, 4, 5, 5f, 5s, 6)")
    parser.add_argument("--list", action="store_true", help="List all reproducible figures, stages, and tables")

    args = parser.parse_args()

    if len(sys.argv) == 1 or args.list:
        list_items()
        return

    if args.tests:
        ok = run_tests()
        sys.exit(0 if ok else 1)
    elif args.figures:
        ok = generate_all_figures()
        sys.exit(0 if ok else 1)
    elif args.fig:
        ok = run_single_figure(args.fig)
        sys.exit(0 if ok else 1)
    elif args.tables:
        for s in ["2", "3", "4", "5", "5s", "6"]:
            if not run_stage(s):
                sys.exit(1)
    elif args.stage:
        ok = run_stage(args.stage)
        sys.exit(0 if ok else 1)
    elif args.all:
        t0 = time.time()
        print("\n>>> STEP 1: Running unit tests...")
        if not run_tests():
            sys.exit(1)
        print("\n>>> STEP 2: Generating all 7 publication figures...")
        if not generate_all_figures():
            sys.exit(1)
        for s in ["2", "5s", "6"]:
            print(f"\n>>> Executing Stage {s}...")
            if not run_stage(s):
                sys.exit(1)
        print(f"\nAll executions completed in {time.time() - t0:.2f} seconds.")


if __name__ == "__main__":
    main()
