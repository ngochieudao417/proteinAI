#!/usr/bin/env python3
"""
scripts/run_quandela_cloud.py
=============================
Quandela Cloud Execution Runner.
Executes live quantum simulation/QPU jobs on Quandela Cloud or local SLOS,
captures execution telemetry, and outputs folding scores to data/quantum_folding_scores.csv.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import perceval as pcvl

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from quantum_folding.circuit import build_ansatz_circuit
from quantum_folding.pipeline import QuantumFoldingPipeline
from quantum_folding.quandela_backend import QuandelaBackend


def setup_logger(log_file: Path) -> logging.Logger:
    """Setup dual file and console logger."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("QuandelaCloudRunner")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Console Handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File Handler
    fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger


def run_demo_circuit(backend: QuandelaBackend, logger: logging.Logger, shots: int = 100) -> None:
    """Execute live Hong-Ou-Mandel quantum bunching & 4-mode folding ansatz demo."""
    logger.info("=" * 70)
    logger.info("RUNNING LIVE QUANTUM DEMONSTRATION ON %s", backend.platform.upper())
    logger.info("=" * 70)

    # Demo 1: 2-mode Hong-Ou-Mandel Quantum Interference
    logger.info("[Demo 1] Submitting 2-mode Hong-Ou-Mandel circuit (|1, 1> -> 50:50 Beam Splitter)...")
    bs_circuit = pcvl.Circuit(2) // pcvl.BS()
    in_state = pcvl.BasicState([1, 1])

    t0 = time.time()
    res = backend.run(circuit=bs_circuit, input_state=in_state, shots=shots)
    dur = time.time() - t0

    counts = res.counts
    n_20 = counts.get("|2,0>", 0)
    n_02 = counts.get("|0,2>", 0)
    n_11 = counts.get("|1,1>", 0)
    total = n_20 + n_02 + n_11 or 1
    bunching_ratio = (n_20 + n_02) / total

    logger.info("Job ID: %s", res.job_id)
    logger.info("Platform: %s | Status: %s | Duration: %.2fs", res.platform, res.status, dur)
    logger.info("Raw Photon Counts: %s", counts)
    logger.info("HOM Quantum Bunching Ratio: %.4f (Theoretical: 1.0000)", bunching_ratio)
    logger.info("Coincidence Clicks |1,1>: %d (Theoretical: 0)", n_11)
    logger.info("Quantum Interference Fidelity: %.2f%%", bunching_ratio * 100.0)

    # Demo 2: 4-mode Parameterized Folding Circuit Ansatz
    logger.info("-" * 70)
    logger.info("[Demo 2] Submitting 4-mode parameterized folding ansatz circuit...")
    ansatz, params = build_ansatz_circuit(num_modes=4)
    # Set sample phase angles
    for idx, p in enumerate(params):
        p.set_value(0.15 * (idx + 1))

    ansatz_in = pcvl.BasicState([1, 0, 1, 0])  # Dual-rail 2-qubit state
    t0_ansatz = time.time()
    res_ansatz = backend.run(circuit=ansatz, input_state=ansatz_in, shots=shots)
    dur_ansatz = time.time() - t0_ansatz

    logger.info("Ansatz Job ID: %s", res_ansatz.job_id)
    logger.info("Ansatz Platform: %s | Status: %s | Duration: %.2fs", res_ansatz.platform, res_ansatz.status, dur_ansatz)
    logger.info("Ansatz Output Distribution: %s", res_ansatz.probabilities)
    logger.info("Demo execution completed successfully.")
    logger.info("=" * 70)


def main() -> int:
    parser = argparse.ArgumentParser(description="Quandela Cloud Quantum Execution Runner")
    parser.add_argument(
        "--platform",
        default="sim:slos",
        choices=["sim:slos", "sim:ascella", "qpu:ascella", "local:slos"],
        help="Target platform (default: sim:slos)",
    )
    parser.add_argument("--shots", type=int, default=100, help="Number of measurement shots (default: 100)")
    parser.add_argument("--token", default=None, help="Quandela API token (requires '_T_' prefix)")
    parser.add_argument("--demo", action="store_true", help="Run quick HOM & ansatz quantum demo")
    parser.add_argument(
        "--candidates-csv",
        default="notebooks/candidate_pool_for_screening.csv",
        help="Path to candidate pool CSV",
    )
    parser.add_argument(
        "--max-candidates",
        type=int,
        default=5,
        help="Maximum candidates to evaluate (default: 5)",
    )
    parser.add_argument(
        "--output-log",
        default="logs/quandela_cloud_execution.log",
        help="Path to execution log file",
    )
    parser.add_argument(
        "--output-csv",
        default="data/quantum_folding_scores.csv",
        help="Path to output scores CSV",
    )
    parser.add_argument("--local", action="store_true", help="Force local simulation mode")
    parser.add_argument("--vqe-iters", type=int, default=10, help="Number of VQE iterations (default: 10)")

    args = parser.parse_args()
    log_path = Path(args.output_log)
    logger = setup_logger(log_path)

    logger.info("=" * 70)
    logger.info("STARTING QUANDELA CLOUD QUANTUM ENGINE RUNNER")
    logger.info("Timestamp: %s", datetime.now(timezone.utc).isoformat())
    logger.info("Platform requested: %s | Mode: %s", args.platform, "Local" if args.local else "Cloud")
    logger.info("=" * 70)

    # Load .env if present
    env_file = Path(__file__).resolve().parent.parent / ".env"
    if env_file.exists():
        try:
            with open(env_file, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#") and "=" in line:
                        k, v = line.split("=", 1)
                        os.environ.setdefault(k.strip(), v.strip())
        except Exception:
            pass

    # Resolve token
    token = args.token or os.getenv("QUANDELA_TOKEN") or os.getenv("QUANDELA_API_KEY", "")

    platform = "local:slos" if args.local else args.platform
    backend = QuandelaBackend(token=token, platform=platform, fallback_to_sim=True, local_fallback=True)

    # Check live platform availability
    avail = backend.check_platform_availability()
    logger.info("Platform Availability Check: %s -> %s (Type: %s)", avail.get("name"), avail.get("status"), avail.get("type"))

    if args.demo:
        run_demo_circuit(backend, logger, shots=args.shots)

    # Candidate Screening
    candidates_path = Path(args.candidates_csv)
    if candidates_path.exists():
        logger.info("-" * 70)
        logger.info("Running Quantum Folding Screening on %s (max %d)...", candidates_path, args.max_candidates)
        pipeline = QuantumFoldingPipeline(backend=backend, window_size=6, vqe_max_iter=args.vqe_iters)
        df_scores = pipeline.evaluate_pool(
            candidates=candidates_path,
            max_candidates=args.max_candidates,
            output_csv=args.output_csv,
        )
        logger.info("Successfully evaluated %d candidates.", len(df_scores))
        logger.info("Top candidates by quantum folding stability score:")
        top_df = df_scores.sort_values(by="quantum_folding_stability_score", ascending=False).head(3)
        for _, row in top_df.iterrows():
            logger.info(
                "  Candidate %s (%s): Score=%.4f, VQE Energy=%.2f, Baseline=%.2f, Gap=%.2f, Platform=%s, Job ID=%s",
                row["name"],
                row["target"],
                row["quantum_folding_stability_score"],
                row["vqe_energy"],
                row["random_baseline_mean"],
                row["energy_gap"],
                row["execution_platform"],
                row["cloud_job_id"],
            )
    else:
        logger.warning("Candidate CSV not found at %s. Skipping screening.", candidates_path)

    logger.info("=" * 70)
    logger.info("Quandela Cloud Runner finished execution cleanly.")
    logger.info("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
