#!/usr/bin/env python3
"""
Master Reproducibility Runner for:
"A Truly Cell-Free Vector Nodal Meshless Method via Meshless Local Petrov-Galerkin (VNMM-MLPG)
for 2D Electromagnetic Problems"
International Journal for Numerical Methods in Engineering (Wiley)

Usage:
    python3 run_all.py --help
    python3 run_all.py --list
    python3 run_all.py --all
    python3 run_all.py --tests
    python3 run_all.py --figures
    python3 run_all.py --tables
    python3 run_all.py --fig 1
    python3 run_all.py --stage 2
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
        "flag": "fig1",
        "output": "figures/fig_mlpg_subdomains_schematic.pdf"
    },
    "2": {
        "name": "Figure 2",
        "desc": "Method of Manufactured Solutions (MMS) asymptotic refinement convergence (Dirichlet & Neumann)",
        "script": "generate_mlpg_figures.py",
        "flag": "fig2",
        "output": "figures/fig_mlpg_mms_convergence.pdf"
    },
    "3": {
        "name": "Figure 3",
        "desc": "Parametric sensitivity of subdomain radius scaling factor alpha_s (coverage vs containment)",
        "script": "generate_mlpg_figures.py",
        "flag": "fig3",
        "output": "figures/fig_mlpg_alpha_s_sensitivity.pdf"
    },
    "4": {
        "name": "Figure 4",
        "desc": "Resonant cavity cutoff wavenumber convergence under 25% nodal jitter",
        "script": "generate_mlpg_figures.py",
        "flag": "fig4",
        "output": "figures/fig_mlpg_cavity_jitter.pdf"
    },
    "5": {
        "name": "Figure 5",
        "desc": "Orthogonal energy stabilization parametric sweep (c_stab) and quadratic rate recovery",
        "script": "generate_mlpg_figures.py",
        "flag": "fig5",
        "output": "figures/fig_mlpg_cavity_stabilization.pdf"
    },
    "6": {
        "name": "Figure 6",
        "desc": "Computational benchmark comparing VNMM-MLPG against cell-based VNMM-Gauss (0% and 25% jitter)",
        "script": "generate_mlpg_figures.py",
        "flag": "fig6",
        "output": "figures/fig_mlpg_benchmark_comparison.pdf"
    }
}

STAGES = {
    "2": {
        "name": "Stage 2: Patch Test",
        "desc": "Galerkin patch test across 7 affine vector modes to machine precision (~10^-14)",
        "script": "mlpg_patch_test_2d.py"
    },
    "3": {
        "name": "Stage 3: MMS Convergence",
        "desc": "MMS refinement study (Tables 2, 3, 4): L2 field and broken H(curl) semi-norm convergence",
        "script": "mlpg_source_problem_2d.py"
    },
    "4": {
        "name": "Stage 4: Subdomain Sizing",
        "desc": "Theoretical bounds and parametric sweep of alpha_s in [0.707, 1.15] (Table 5)",
        "script": "parametric_analysis_alpha_s_mlpg_2d.py"
    },
    "5": {
        "name": "Stage 5: Cavity Eigenvalues",
        "desc": "Cavity spectrum, superlinear s_div scaling, and orthogonal stabilization c_stab (Tables 6, 7, 8)",
        "script": "mlpg_cavity_eigenvalues_2d.py"
    },
    "6": {
        "name": "Stage 6: Computational Benchmark",
        "desc": "End-to-end benchmark: assembly speedup (2.0x-3.7x) and sparsity reduction (Tables 9, 10)",
        "script": "compute_eabe_benchmark.py"
    }
}


def run_tests():
    print("=" * 80)
    print("RUNNING AUTOMATED UNIT TESTS (20 MLPG TESTS)")
    print("=" * 80)
    cmd = [sys.executable, "-m", "unittest", "discover", "-s", TESTS_DIR, "-p", "test_*.py"]
    res = subprocess.run(cmd)
    return res.returncode == 0


def generate_all_figures():
    print("=" * 80)
    print("GENERATING ALL 6 PUBLICATION FIGURES (.PDF & .PNG)")
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
    # Import and call specific plot function
    import importlib.util
    spec = importlib.util.spec_from_file_location("gen_fig", fig_script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    func_map = {
        "1": mod.plot_fig1_schematic,
        "2": mod.plot_fig2_mms_convergence,
        "3": mod.plot_fig3_alpha_s_sensitivity,
        "4": mod.plot_fig4_cavity_jitter,
        "5": mod.plot_fig5_cavity_stabilization,
        "6": mod.plot_fig6_benchmark_comparison,
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
    script_path = os.path.join(SCRIPTS_DIR, item['script'])
    res = subprocess.run([sys.executable, script_path])
    return res.returncode == 0


def list_items():
    print("=" * 80)
    print("REPRODUCIBILITY CATALOG FOR VNMM-MLPG (IJNME)")
    print("=" * 80)
    print("\n[Figures]")
    for k, v in sorted(FIGURES.items()):
        print(f"  --fig {k:1s} : {v['name']:15s} -> {v['desc']}")
        print(f"           Artifact: {v['output']}")
    print("\n[Validation Stages & Tables]")
    for k, v in sorted(STAGES.items()):
        print(f"  --stage {k:1s} : {v['name']:30s} -> {v['desc']}")
        print(f"             Script: scripts/{v['script']}")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description="Master reproducibility runner for the VNMM-MLPG IJNME paper."
    )
    parser.add_argument("--all", action="store_true", help="Run tests, generate all figures, and run all stages")
    parser.add_argument("--tests", action="store_true", help="Run automated unit test suite (20 tests)")
    parser.add_argument("--figures", action="store_true", help="Generate all publication-quality figures (.pdf & .png)")
    parser.add_argument("--fig", type=str, help="Generate a specific figure (e.g. --fig 1, --fig 2, ..., --fig 6)")
    parser.add_argument("--tables", action="store_true", help="Run all validation stages to reproduce tables")
    parser.add_argument("--stage", type=str, help="Run a specific validation stage (e.g. --stage 2, --stage 3, ..., --stage 6)")
    parser.add_argument("--list", action="store_true", help="List all reproducible figures, stages, and tables")

    args = parser.parse_args()

    if len(sys.argv) == 1 or args.list:
        list_items()
        return

    if args.tests:
        run_tests()
    elif args.figures:
        generate_all_figures()
    elif args.fig:
        run_single_figure(args.fig)
    elif args.tables:
        for s in sorted(STAGES.keys()):
            run_stage(s)
    elif args.stage:
        run_stage(args.stage)
    elif args.all:
        t0 = time.time()
        print("\n>>> STEP 1: Running unit tests...")
        run_tests()
        print("\n>>> STEP 2: Generating all publication figures...")
        generate_all_figures()
        print("\n>>> STEP 3: Executing Stage 2 (Patch Test)...")
        run_stage("2")
        print("\n>>> STEP 4: Executing Stage 6 (Benchmark Comparison)...")
        run_stage("6")
        print(f"\nAll executions completed in {time.time() - t0:.2f} seconds.")


if __name__ == "__main__":
    main()

