"""
quantum_folding.quandela_backend
================================
Quandela Cloud Backend and Execution Client.
Supports Quandela Cloud RemoteProcessor (sim:slos, sim:ascella, qpu:ascella)
and local Perceval SLOS simulation fallback.
"""

from __future__ import annotations

import logging
import os
import time
import uuid
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import requests
import perceval as pcvl
from perceval.algorithm import Sampler
from perceval.runtime.remote_job import RemoteJob

logger = logging.getLogger("QuandelaBackend")


@dataclass
class ExecutionResult:
    """Standardized result container for quantum circuit execution."""
    job_id: str
    platform: str
    status: str
    duration: float
    counts: Dict[str, int]
    probabilities: Dict[str, float]
    raw_results: Any
    global_perf: float
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "job_id": self.job_id,
            "platform": self.platform,
            "status": self.status,
            "duration": self.duration,
            "counts": self.counts,
            "probabilities": self.probabilities,
            "global_perf": self.global_perf,
            "metadata": self.metadata,
        }


class QuandelaBackend:
    """Client for Quandela Cloud and local Perceval photonic simulation."""

    DEFAULT_BASE_URL = "https://api.cloud.quandela.com"
    SUPPORTED_CLOUD_PLATFORMS = ("sim:slos", "sim:ascella", "qpu:ascella", "qpu:belenos")

    def __init__(
        self,
        token: Optional[str] = None,
        platform: str = "sim:slos",
        base_url: str = DEFAULT_BASE_URL,
        fallback_to_sim: bool = True,
        local_fallback: bool = True,
        timeout_seconds: int = 120,
        poll_interval: float = 0.5,
    ) -> None:
        """
        Initialize QuandelaBackend client.

        Args:
            token: Quandela API key with mandatory literal '_T_' prefix.
                   If None, defaults to QUANDELA_TOKEN environment variable.
            platform: Platform name ('sim:slos', 'sim:ascella', 'qpu:ascella', 'local:slos').
            base_url: Quandela Cloud API base URL.
            fallback_to_sim: If True, falls back to 'sim:slos' if target QPU is in maintenance.
            local_fallback: If True, falls back to local Perceval simulation if offline/unreachable.
            timeout_seconds: Timeout for cloud job execution.
            poll_interval: Polling frequency in seconds.
        """
        self.base_url = base_url.rstrip("/")
        self.fallback_to_sim = fallback_to_sim
        self.local_fallback = local_fallback
        self.timeout_seconds = timeout_seconds
        self.poll_interval = poll_interval
        self._jobs: Dict[str, Any] = {}

        # Resolve token
        resolved_token = token if token is not None else os.getenv("QUANDELA_TOKEN")
        self.is_local = (platform in ("local:slos", "local"))

        if self.is_local:
            self.token = None
            self.platform = "local:slos"
            logger.info("QuandelaBackend configured in local SLOS simulation mode.")
        else:
            if not resolved_token:
                if self.local_fallback:
                    logger.warning("No Quandela token provided. Falling back to local SLOS simulation.")
                    self.is_local = True
                    self.token = None
                    self.platform = "local:slos"
                else:
                    raise ValueError("API token must be a non-empty string starting with '_T_'.")
            else:
                self._validate_token(resolved_token)
                self.token = resolved_token
                self.platform = platform

        self._headers = {"Authorization": f"Bearer {self.token}"} if self.token else {}

    @staticmethod
    def _validate_token(token: str) -> None:
        """Enforces literal '_T_' prefix and non-empty string."""
        if not token or not isinstance(token, str):
            raise ValueError("API token must be a non-empty string.")
        if not token.startswith("_T_"):
            prefix = token[:5] if len(token) >= 5 else token
            raise ValueError(
                f"Invalid Quandela API token: must start with literal '_T_' prefix (got '{prefix}...')."
            )

    def check_platform_availability(self, platform_name: Optional[str] = None) -> Dict[str, Any]:
        """Query platform availability via GET /api/platform/{name}."""
        if self.is_local:
            return {"name": "local:slos", "status": "available", "type": "local_simulator"}

        target = platform_name or self.platform
        url = f"{self.base_url}/api/platform/{target}"
        try:
            r = requests.get(url, headers=self._headers, timeout=10)
            if r.status_code == 200:
                data = r.json()
                return {
                    "name": target,
                    "status": data.get("status", "unknown"),
                    "type": data.get("type", "unknown"),
                    "specs": data.get("specs", {}),
                    "raw": data,
                }
            return {"name": target, "status": "error", "code": r.status_code, "text": r.text}
        except Exception as e:
            logger.warning("Failed to reach platform availability endpoint: %s", e)
            return {"name": target, "status": "unreachable", "error": str(e)}

    def check_job_availability(self) -> Dict[str, Any]:
        """Query job queue availability via GET /api/jobs/availability."""
        if self.is_local:
            return {"max_concurrent_jobs": 1, "num_concurrent_jobs": 0}

        url = f"{self.base_url}/api/jobs/availability"
        try:
            r = requests.get(url, headers=self._headers, timeout=10)
            if r.status_code == 200:
                return r.json()
            return {"error": r.text, "code": r.status_code}
        except Exception as e:
            return {"status": "unreachable", "error": str(e)}

    def submit_job(
        self,
        circuit: pcvl.Circuit,
        input_state: pcvl.BasicState,
        shots: int = 100,
        min_detected_photons: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Asynchronously submit a circuit execution job."""
        if self.is_local:
            job_id = f"local-slos-{uuid.uuid4().hex[:12]}"
            self._jobs[job_id] = {
                "type": "local",
                "circuit": circuit,
                "input_state": input_state,
                "shots": shots,
                "status": "completed",
            }
            return {"job_id": job_id, "platform": "local:slos", "status": "completed"}

        # Validate platform status & fallback if necessary
        if self.fallback_to_sim and self.platform.startswith("qpu:"):
            info = self.check_platform_availability(self.platform)
            if info.get("status") != "available":
                logger.warning(
                    "Platform '%s' is %s. Falling back to 'sim:slos'.",
                    self.platform,
                    info.get("status"),
                )
                self.platform = "sim:slos"

        proc = pcvl.RemoteProcessor(self.platform, token=self.token)
        proc.set_circuit(circuit)
        proc.with_input(input_state)
        n_filter = min_detected_photons or input_state.n
        proc.min_detected_photons_filter(n_filter)

        sampler = Sampler(proc, max_shots_per_call=1000)
        remote_job = sampler.sample_count.execute_async(shots)
        job_id = remote_job.id
        self._jobs[job_id] = remote_job

        return {"job_id": job_id, "platform": self.platform, "status": str(remote_job.status())}

    def get_job_status(self, job_id: str) -> str:
        """Query status of a job ('waiting', 'running', 'completed', 'failed')."""
        if job_id not in self._jobs:
            if not self.is_local and self.token:
                try:
                    proc = pcvl.RemoteProcessor(self.platform, token=self.token)
                    job = RemoteJob.from_id(job_id, proc.get_rpc_handler())
                    self._jobs[job_id] = job
                    return str(job.status()).lower()
                except Exception:
                    pass
            return "unknown"

        job_obj = self._jobs[job_id]
        if isinstance(job_obj, dict):
            return job_obj.get("status", "unknown")
        return str(job_obj.status()).lower()

    def get_job_results(self, job_id: str) -> Dict[str, Any]:
        """Retrieve results of a completed job."""
        if job_id not in self._jobs:
            if not self.is_local and self.token:
                proc = pcvl.RemoteProcessor(self.platform, token=self.token)
                job = RemoteJob.from_id(job_id, proc.get_rpc_handler())
                return job.get_results()
            raise KeyError(f"Job ID '{job_id}' not found.")

        job_obj = self._jobs[job_id]
        if isinstance(job_obj, dict):
            return self._execute_local_sync(
                job_obj["circuit"], job_obj["input_state"], job_obj["shots"]
            )
        return job_obj.get_results()

    def run(
        self,
        circuit: pcvl.Circuit,
        input_state: pcvl.BasicState,
        shots: int = 100,
        command: str = "sample_count",
        min_detected_photons: Optional[int] = None,
        wait: bool = True,
        timeout: Optional[int] = None,
    ) -> ExecutionResult:
        """
        Execute a circuit with synchronous wait and telemetry capture.

        Args:
            circuit: Perceval Circuit instance.
            input_state: Input BasicState (Fock state).
            shots: Number of sample shots.
            command: 'sample_count' or 'probs'.
            min_detected_photons: Post-selection threshold (defaults to input_state.n).
            wait: If True, blocks until job completion.
            timeout: Maximum wait time in seconds.
        """
        t0 = time.time()
        timeout_limit = timeout or self.timeout_seconds

        if self.is_local:
            return self._run_local(circuit, input_state, shots, command)

        # Check platform availability & fallback
        if self.fallback_to_sim and self.platform.startswith("qpu:"):
            info = self.check_platform_availability(self.platform)
            if info.get("status") != "available":
                logger.warning(
                    "Target '%s' status is '%s'. Falling back to 'sim:slos'.",
                    self.platform,
                    info.get("status"),
                )
                self.platform = "sim:slos"

        try:
            proc = pcvl.RemoteProcessor(self.platform, token=self.token)
            proc.set_circuit(circuit)
            proc.with_input(input_state)
            n_filter = min_detected_photons or input_state.n
            proc.min_detected_photons_filter(n_filter)

            sampler = Sampler(proc, max_shots_per_call=1000)

            if command == "sample_count":
                job = sampler.sample_count.execute_async(shots)
            else:
                job = sampler.probs.execute_async()

            job_id = job.id
            self._jobs[job_id] = job

            if not wait:
                return ExecutionResult(
                    job_id=job_id,
                    platform=self.platform,
                    status=str(job.status()),
                    duration=time.time() - t0,
                    counts={},
                    probabilities={},
                    raw_results=None,
                    global_perf=1.0,
                )

            # Polling loop
            while not job.is_complete:
                if time.time() - t0 > timeout_limit:
                    raise TimeoutError(f"Cloud job {job_id} timed out after {timeout_limit}s.")
                time.sleep(self.poll_interval)

            if job.is_failed:
                raise RuntimeError(f"Cloud job {job_id} failed with status {job.status()}.")

            raw_results = job.get_results()
            counts, probs = self._decode_results(raw_results, command, shots)

            return ExecutionResult(
                job_id=job_id,
                platform=self.platform,
                status="SUCCESS",
                duration=time.time() - t0,
                counts=counts,
                probabilities=probs,
                raw_results=raw_results,
                global_perf=raw_results.get("global_perf", 1.0),
            )

        except Exception as e:
            if self.local_fallback:
                logger.warning("Cloud execution failed (%s). Falling back to local simulation.", e)
                return self._run_local(circuit, input_state, shots, command)
            raise

    def run_circuit(
        self,
        circuit: pcvl.Circuit,
        input_state: pcvl.BasicState,
        shots: int = 100,
    ) -> Dict[str, Any]:
        """
        Unified execution interface adhering to SCOPE.md interface contract.
        
        Returns a dictionary containing:
        - 'results': raw dictionary mapping BasicState or FockState to counts/probabilities
        - 'counts': string-keyed count map
        - 'probabilities': string-keyed normalized probability map
        - 'job_id': execution job identifier
        - 'platform': execution platform
        - 'status': execution status
        """
        res = self.run(circuit=circuit, input_state=input_state, shots=shots, command="sample_count")
        raw_res = res.raw_results
        results_map = raw_res.get("results", {}) if isinstance(raw_res, dict) else res.counts
        return {
            "results": results_map,
            "counts": res.counts,
            "probabilities": res.probabilities,
            "job_id": res.job_id,
            "platform": res.platform,
            "status": res.status,
            "duration": res.duration,
        }

    def _run_local(
        self,
        circuit: pcvl.Circuit,
        input_state: pcvl.BasicState,
        shots: int = 100,
        command: str = "sample_count",
    ) -> ExecutionResult:
        """Execute on local Perceval Processor('SLOS')."""
        t0 = time.time()
        raw = self._execute_local_sync(circuit, input_state, shots, command)
        counts, probs = self._decode_results(raw, command, shots)
        return ExecutionResult(
            job_id=f"local-slos-{uuid.uuid4().hex[:12]}",
            platform="local:slos",
            status="SUCCESS",
            duration=time.time() - t0,
            counts=counts,
            probabilities=probs,
            raw_results=raw,
            global_perf=raw.get("global_perf", 1.0),
        )

    @staticmethod
    def _execute_local_sync(
        circuit: pcvl.Circuit,
        input_state: pcvl.BasicState,
        shots: int = 100,
        command: str = "sample_count",
    ) -> Dict[str, Any]:
        """Synchronously execute circuit on local SLOS processor."""
        proc = pcvl.Processor("SLOS", circuit.m)
        proc.add(0, circuit)
        proc.with_input(input_state)
        sampler = Sampler(proc)
        if command == "sample_count":
            return sampler.sample_count.execute_sync(shots)
        return sampler.probs.execute_sync()

    @staticmethod
    def _decode_results(
        raw_results: Dict[str, Any], command: str, shots: int
    ) -> Tuple[Dict[str, int], Dict[str, float]]:
        """Normalize BSCount or BSDistribution into clean string-keyed dicts."""
        res_map = raw_results.get("results", {}) if isinstance(raw_results, dict) else {}
        counts: Dict[str, int] = {}
        probs: Dict[str, float] = {}

        if command == "sample_count":
            for k, v in res_map.items():
                counts[str(k)] = int(v)
            total = sum(counts.values()) or 1
            for k, v in counts.items():
                probs[k] = v / total
        else:
            for k, v in res_map.items():
                probs[str(k)] = float(v)
                counts[str(k)] = int(round(float(v) * shots))

        return counts, probs
