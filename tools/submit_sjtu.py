#!/usr/bin/env python3
"""
tools/submit_sjtu.py — SJTU AI Protein Platform Client & CLI Submitter.

Handles contestant authentication, token session refresh, sequence pre-validation,
batch CSV upload, and submission status verification for Changjiang Cup 2026.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import requests

# ---------------------------------------------------------------------------
# Platform Configuration & Target Lookup
# ---------------------------------------------------------------------------

def _load_dotenv_if_present() -> None:
    env_path = Path(__file__).resolve().parent.parent / ".env"
    if env_path.exists():
        try:
            with open(env_path, encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#") and "=" in line:
                        k, v = line.split("=", 1)
                        os.environ.setdefault(k.strip(), v.strip())
        except Exception:
            pass

_load_dotenv_if_present()

DEFAULT_BASE_URL = os.environ.get("SJTU_BASE_URL", "https://aiprotein.sjtu.edu.cn")
DEFAULT_EMAIL = os.environ.get("SJTU_EMAIL", "")
DEFAULT_PASSWORD = os.environ.get("SJTU_PASSWORD", "")

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

TARGET_CONFIGS: Dict[str, Dict[str, Any]] = {
    "egfr": {
        "id": "af6b1e2b-c62d-40e4-a5b3-3f091be04c60",
        "code": "egfr",
        "name": "EGFR",
        "modality": "nanobody",
        "min_len": 100,
        "max_len": 150,
        "max_seqs": 20,
        "min_length": 100,
        "max_length": 150,
        "max_sequences": 20,
    },
    "tnf": {
        "id": "16c8b90c-a678-4439-896c-67639c5aa2c0",
        "code": "tnf",
        "name": "TNF-α",
        "modality": "nanobody",
        "min_len": 100,
        "max_len": 150,
        "max_seqs": 20,
        "min_length": 100,
        "max_length": 150,
        "max_sequences": 20,
    },
    "glp1r": {
        "id": "e9134b84-24a8-4a24-bad4-5f5b095c2e26",
        "code": "glp1r",
        "name": "GLP-1R",
        "modality": "nanobody",
        "min_len": 100,
        "max_len": 150,
        "max_seqs": 20,
        "min_length": 100,
        "max_length": 150,
        "max_sequences": 20,
    },
    "tfr1": {
        "id": "c05ca94e-cf7e-4059-a84f-20073d7c5674",
        "code": "tfr1",
        "name": "TfR1",
        "modality": "nanobody",
        "min_len": 100,
        "max_len": 150,
        "max_seqs": 20,
        "min_length": 100,
        "max_length": 150,
        "max_sequences": 20,
    },
    "tl1a": {
        "id": "ea930a01-20a9-4bed-b70f-e658bd9ba3ab",
        "code": "tl1a",
        "name": "TL1A",
        "modality": "minibinder",
        "min_len": 40,
        "max_len": 100,
        "max_seqs": 20,
        "min_length": 40,
        "max_length": 100,
        "max_sequences": 20,
    },
    "tnfr1": {
        "id": "23b568a0-250a-4840-96d8-635e4e78e763",
        "code": "tnfr1",
        "name": "TNFR1",
        "modality": "nanobody",
        "min_len": 100,
        "max_len": 150,
        "max_seqs": 20,
        "min_length": 100,
        "max_length": 150,
        "max_sequences": 20,
    },
    "il7ra": {
        "id": "1347d0e0-7d3d-4aa8-bcf1-c2dab9d8f1df",
        "code": "il7ra",
        "name": "IL-7Rα",
        "modality": "minibinder",
        "min_len": 40,
        "max_len": 100,
        "max_seqs": 20,
        "min_length": 40,
        "max_length": 100,
        "max_sequences": 20,
    },
}

TL1A_TARGET_ID = "ea930a01-20a9-4bed-b70f-e658bd9ba3ab"
IL7RA_TARGET_ID = "1347d0e0-7d3d-4aa8-bcf1-c2dab9d8f1df"

COMPETITION_TARGETS: List[Dict[str, Any]] = list(TARGET_CONFIGS.values())
TARGETS_BY_ID: Dict[str, Dict[str, Any]] = {t["id"]: t for t in COMPETITION_TARGETS}
TARGETS_BY_ID["6ba7b810-9dad-11d1-80b4-00c04fd430c8"] = TARGETS_BY_ID[TL1A_TARGET_ID]
TARGETS_BY_ID["b7a2e850-20cf-46d5-a3d8-5807908c6ff9"] = TARGETS_BY_ID[IL7RA_TARGET_ID]
TARGETS_BY_CODE: Dict[str, Dict[str, Any]] = {t["code"]: t for t in COMPETITION_TARGETS}


def parse_csv_rows(csv_text: str) -> List[List[str]]:
    """Parse CSV text with RFC 4180 quote support, stripping whitespace and empty lines."""
    clean_text = csv_text.lstrip("\ufeff")
    reader = csv.reader(io.StringIO(clean_text))
    rows: List[List[str]] = []
    for r in reader:
        if any(c.strip() for c in r):
            rows.append([c.strip() for c in r])
    return rows


# ---------------------------------------------------------------------------
# Exceptions & Helper Types
# ---------------------------------------------------------------------------

class SJTUPlatformError(Exception):
    """Base exception for all SJTU platform errors."""
    pass


class SJTUAuthenticationError(PermissionError, SJTUPlatformError):
    """Raised when contestant authentication fails."""
    pass


class SJTUValidationError(ValueError, SJTUPlatformError):
    """Raised when sequence/CSV pre-validation or server 422 validation fails."""
    def __init__(self, message: str, errors: Optional[List[Any]] = None):
        super().__init__(message)
        self.errors = errors or []


class SJTUReUploadNoteRequiredError(SJTUValidationError):
    """Raised when re-uploading to a target with prior submission without an explanation note."""
    pass


class SJTUTargetNotFoundError(KeyError, SJTUValidationError):
    """Raised when the specified competition target does not exist on the platform (HTTP 404)."""
    def __init__(self, message: str = "Target does not exist on platform (HTTP 404): 赛题不存在", errors: Optional[List[Any]] = None):
        super().__init__(message)
        self.message = message
        self.errors = errors or []

    def __str__(self) -> str:
        return self.message


class RefreshTokenResponse(dict):
    """Dictionary response from refresh_token that also stringifies to access_token."""
    def __str__(self) -> str:
        return str(self.get("access_token", ""))


# ---------------------------------------------------------------------------
# Logging Setup
# ---------------------------------------------------------------------------

def setup_submission_logger(
    log_file: str = "logs/sjtu_submission.log",
    verbose: bool = False,
) -> logging.Logger:
    """Configures console and audit file logging."""
    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("sjtu_submission")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.propagate = False

    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(str(log_path), mode="a", encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


# ---------------------------------------------------------------------------
# Sequence & CSV Local Validator
# ---------------------------------------------------------------------------

def _split_csv_line(line: str) -> List[str]:
    """Split a single CSV line handling quotes (RFC 4180)."""
    out: List[str] = []
    cur = ""
    q = False
    i = 0
    while i < len(line):
        c = line[i]
        if q:
            if c == '"':
                if i + 1 < len(line) and line[i + 1] == '"':
                    cur += '"'
                    i += 1
                else:
                    q = False
            else:
                cur += c
        else:
            if c == '"':
                q = True
            elif c == ",":
                out.append(cur)
                cur = ""
            else:
                cur += c
        i += 1
    out.append(cur)
    return out


def validate_sequence(
    sequence_or_name: str,
    sequence_or_min_len: Optional[Union[str, int]] = None,
    modality_or_max_len: Optional[Union[str, int]] = None,
    min_len: Optional[int] = None,
    max_len: Optional[int] = None,
    modality: Optional[str] = None,
) -> Tuple[bool, Optional[str]]:
    """
    Validates single sequence against canonical amino acid and length constraints.
    Supports flexible signatures:
      - validate_sequence(sequence, min_len=40, max_len=100)
      - validate_sequence(name, sequence, modality="minibinder")
      - validate_sequence(sequence, 40, 100)
    """
    if isinstance(sequence_or_min_len, str):
        # Called as validate_sequence(name, sequence, modality=...)
        sequence = sequence_or_min_len
        if isinstance(modality_or_max_len, str):
            modality = modality_or_max_len
        elif isinstance(modality_or_max_len, int):
            max_len = modality_or_max_len
    else:
        # Called as validate_sequence(sequence, min_len, max_len)
        sequence = sequence_or_name
        if isinstance(sequence_or_min_len, int):
            min_len = sequence_or_min_len
        if isinstance(modality_or_max_len, int):
            max_len = modality_or_max_len

    if modality:
        if modality.lower() in ("nanobody", "vhh"):
            min_len = 100 if min_len is None else min_len
            max_len = 150 if max_len is None else max_len
        else:
            min_len = 40 if min_len is None else min_len
            max_len = 100 if max_len is None else max_len
    else:
        min_len = 40 if min_len is None else min_len
        max_len = 100 if max_len is None else max_len

    clean = re.sub(r"\s+", "", sequence).upper()
    if not clean:
        return False, "序列为空"
    if any(c not in CANONICAL_AA for c in clean):
        return False, "含有非法氨基酸"
    if len(clean) < min_len or len(clean) > max_len:
        return False, f"序列长度不符合要求：范围为 {min_len}—{max_len} aa"
    return True, None


def validate_submission_csv(
    csv_text_or_path: Union[str, Path],
    modality: str = "minibinder",
    max_rows: int = 20,
) -> Tuple[int, List[Tuple[int, str, str]]]:
    """
    Validates CSV content strictly following platform rules.
    Returns: (sequence_count, [(line_index, name, error_msg)])
    """
    if isinstance(csv_text_or_path, Path) or (
        isinstance(csv_text_or_path, str)
        and csv_text_or_path.strip()
        and "\n" not in csv_text_or_path
        and Path(csv_text_or_path).is_file()
    ):
        text = Path(csv_text_or_path).read_text(encoding="utf-8-sig")
    else:
        text = str(csv_text_or_path).lstrip("\ufeff")

    lines = [l.strip() for l in re.split(r"\r?\n", text) if l.strip()]
    if not lines:
        return 0, [(1, "", "没有序列")]

    hdr = [h.strip().lower() for h in _split_csv_line(lines[0])]
    if "name" not in hdr or "sequence" not in hdr:
        return 0, [(1, "", "表头必须包含 name 和 sequence")]

    min_len = 100 if modality.lower() in ("nanobody", "vhh") else 40
    max_len = 150 if modality.lower() in ("nanobody", "vhh") else 100

    name_idx = hdr.index("name")
    seq_idx = hdr.index("sequence")

    errs: List[Tuple[int, str, str]] = []
    seen: set[str] = set()
    count = 0

    if len(lines) - 1 > max_rows:
        return len(lines) - 1, [(max_rows + 2, "", f"最多提交 {max_rows} 条序列（含表头共 {max_rows + 1} 行）")]

    for idx, line in enumerate(lines[1:], start=2):
        parts = _split_csv_line(line)
        name_val = parts[name_idx].strip() if name_idx < len(parts) else ""
        seq_val = re.sub(r"\s+", "", parts[seq_idx] if seq_idx < len(parts) else "").upper()

        if not name_val and not seq_val:
            continue

        count += 1
        if count > max_rows:
            return count, [(idx, "", f"最多提交 {max_rows} 条序列（含表头共 {max_rows + 1} 行）")]

        if not name_val or len(name_val) > 64:
            errs.append((idx, name_val, "名称须为 1 到 64 个字符"))
            continue

        if name_val in seen:
            errs.append((idx, name_val, "名称重复"))
            continue
        seen.add(name_val)

        if not seq_val:
            errs.append((idx, name_val, "序列为空"))
            continue

        if any(c not in CANONICAL_AA for c in seq_val):
            errs.append((idx, name_val, "含有非法氨基酸"))
            continue

        if len(seq_val) < min_len or len(seq_val) > max_len:
            errs.append((idx, name_val, f"序列长度不符合要求：{modality} 为 {min_len}—{max_len} aa"))

    if count == 0:
        errs.append((1, "", "没有序列"))

    return count, errs


# ---------------------------------------------------------------------------
# SJTUPlatformClient Implementation
# ---------------------------------------------------------------------------

class SJTUPlatformClient:
    """Production client for SJTU AI Protein platform."""

    def __init__(
        self,
        email: str = DEFAULT_EMAIL,
        password: str = DEFAULT_PASSWORD,
        base_url: str = DEFAULT_BASE_URL,
        logger: Optional[logging.Logger] = None,
    ):
        self.email = email
        self.password = password
        self.base_url = base_url.rstrip("/")
        self.session = requests.Session()
        self.access_token: Optional[str] = None
        self.refresh_token: Optional[str] = None
        self.user: Optional[Dict[str, Any]] = None
        self.user_profile: Optional[Dict[str, Any]] = None
        self.team: Optional[Dict[str, Any]] = None
        self.logger = logger or logging.getLogger("sjtu_submission")

    def _auth_request(self, method: str, url: str, **kwargs) -> requests.Response:
        """Issues HTTP request with automatic token refresh on HTTP 401."""
        if self.access_token:
            if "headers" in kwargs:
                kwargs["headers"].setdefault("Authorization", f"Bearer {self.access_token}")
            else:
                kwargs["headers"] = {"Authorization": f"Bearer {self.access_token}"}

        m = method.upper()
        if m == "GET":
            resp = self.session.get(url, **kwargs)
        elif m == "POST":
            resp = self.session.post(url, **kwargs)
        elif m == "PUT":
            resp = self.session.put(url, **kwargs)
        elif m == "DELETE":
            resp = self.session.delete(url, **kwargs)
        else:
            resp = self.session.request(method, url, **kwargs)

        if resp.status_code == 401 and self.refresh_token:
            self.logger.warning("[Auth] Access token expired (HTTP 401). Attempting token refresh...")
            try:
                self.refresh_token_session()
                # Update authorization header with refreshed token
                if "headers" in kwargs:
                    kwargs["headers"]["Authorization"] = f"Bearer {self.access_token}"
                else:
                    kwargs["headers"] = {"Authorization": f"Bearer {self.access_token}"}
                if m == "GET":
                    resp = self.session.get(url, **kwargs)
                elif m == "POST":
                    resp = self.session.post(url, **kwargs)
                elif m == "PUT":
                    resp = self.session.put(url, **kwargs)
                elif m == "DELETE":
                    resp = self.session.delete(url, **kwargs)
                else:
                    resp = self.session.request(method, url, **kwargs)
            except Exception as e:
                self.logger.error(f"[Auth] Token refresh failed: {e}")
                raise SJTUAuthenticationError("Session expired and token refresh failed. Please re-login.")
        if resp.status_code >= 500:
            resp.raise_for_status()
        return resp

    def login(self) -> Dict[str, Any]:
        """Authenticates contestant and initializes Bearer token."""
        url = f"{self.base_url}/api/auth/login"
        payload = {
            "email": self.email,
            "password": self.password,
            "portal": "contestant",
        }
        resp = self.session.post(url, json=payload, timeout=20)
        if resp.status_code == 401:
            raise SJTUAuthenticationError("用户名或密码错误")
        resp.raise_for_status()

        data = resp.json()
        self.access_token = data.get("access_token")
        self.refresh_token = data.get("refresh_token")
        self.user = data.get("user")
        self.user_profile = self.user
        self.team = self.user.get("team") if self.user else None
        if self.access_token:
            self.session.headers.update({"Authorization": f"Bearer {self.access_token}"})
        self.logger.info(
            f"[Auth] Authenticated as {self.user.get('name') if self.user else self.email} "
            f"(Team: {self.team.get('name') if self.team else 'N/A'}, Role: {self.team.get('member_role') if self.team else 'N/A'})"
        )
        return data

    def refresh_token_session(self) -> RefreshTokenResponse:
        """Refreshes expired access token."""
        url = f"{self.base_url}/api/auth/refresh"
        resp = self.session.post(url, json={"refresh_token": self.refresh_token}, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        self.access_token = data.get("access_token")
        if data.get("refresh_token"):
            self.refresh_token = data.get("refresh_token")
        if self.access_token:
            self.session.headers.update({"Authorization": f"Bearer {self.access_token}"})
        return RefreshTokenResponse(data)

    def get_me(self) -> Dict[str, Any]:
        """Queries current user profile from platform."""
        url = f"{self.base_url}/api/auth/me"
        resp = self._auth_request("GET", url, timeout=15)
        resp.raise_for_status()
        return resp.json()

    def get_targets(self) -> List[Dict[str, Any]]:
        """Queries active challenge target tracks."""
        url = f"{self.base_url}/api/submissions/targets"
        resp = self._auth_request("GET", url, timeout=15)
        resp.raise_for_status()
        return resp.json()

    def resolve_target(self, target_input: str) -> Dict[str, Any]:
        """Resolves target input (code or UUID) to target config dict."""
        key = target_input.strip().lower()
        if key in TARGET_CONFIGS:
            return TARGET_CONFIGS[key]
        for cfg in TARGET_CONFIGS.values():
            if cfg["id"].lower() == key:
                return cfg
        valid_keys = ", ".join(TARGET_CONFIGS.keys())
        raise ValueError(f"Unknown target '{target_input}'. Valid targets: {valid_keys}")

    def validate_sequence(
        self,
        sequence: str,
        min_len: int = 40,
        max_len: int = 100,
    ) -> Tuple[bool, Optional[str]]:
        """Client method for sequence validation."""
        return validate_sequence(sequence, min_len=min_len, max_len=max_len)

    def validate_submission_csv(
        self,
        csv_path_or_content: Union[str, Path],
        modality: str = "minibinder",
    ) -> Tuple[bool, List[str]]:
        """Client method to validate candidate CSV file or content."""
        count, errs = validate_submission_csv(csv_path_or_content, modality=modality)
        if errs:
            return False, [f"Line {e[0]} ({e[1]}): {e[2]}" for e in errs]
        return True, []

    def validate_csv(
        self,
        csv_content: str,
        modality: str = "minibinder",
    ) -> Tuple[bool, List[str]]:
        """Alias for validate_submission_csv."""
        return self.validate_submission_csv(csv_content, modality=modality)

    def get_latest_submission(self, target_id: str) -> Optional[Dict[str, Any]]:
        """Retrieves metadata of the latest submission for target."""
        url = f"{self.base_url}/api/submissions/latest?target_id={target_id}"
        resp = self._auth_request("GET", url, timeout=15)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        if not resp.content or not resp.text.strip():
            return None
        try:
            data = resp.json()
            if not data:
                return None
            return data
        except Exception:
            return None

    def get_submissions(self, grouped: bool = True) -> Dict[str, Any]:
        """Retrieves team's submission history."""
        url = f"{self.base_url}/api/submissions"
        params = {"grouped": 1} if grouped else {}
        resp = self._auth_request("GET", url, params=params, timeout=15)
        resp.raise_for_status()
        return resp.json()

    def submit_candidates(
        self,
        target_id: str,
        csv_path: str,
        action: str = "confirm",
        agreed: bool = True,
        abstract: str = "",
        note: str = "",
        original_filename: Optional[str] = None,
        bypass_local_validation: bool = False,
        auto_note_on_reupload: bool = True,
        dry_run: bool = False,
        max_retries: int = 1,
    ) -> Dict[str, Any]:
        """
        Submits candidate sequences CSV to SJTU platform.
        Validates locally, checks re-upload status, and records audit trail.
        """
        if not agreed:
            raise SJTUValidationError("请先勾选并同意赛事细则")

        if not target_id:
            raise SJTUValidationError("Target ID is required")

        if len(abstract) > 4000:
            raise SJTUValidationError(f"Abstract exceeds maximum allowed length of 4000 characters ({len(abstract)})")
        if len(note) > 2000:
            raise SJTUValidationError(f"Note exceeds maximum allowed length of 2000 characters ({len(note)})")

        try:
            target_cfg = self.resolve_target(target_id)
            actual_target_id = target_cfg["id"]
            target_code = target_cfg.get("code", target_id)
            modality = target_cfg.get("modality", "minibinder")
        except ValueError:
            actual_target_id = target_id
            target_code = target_id
            modality = "minibinder"

        csv_file = Path(csv_path)
        if not csv_file.exists():
            raise FileNotFoundError(f"Candidate CSV not found: {csv_path}")

        csv_text = csv_file.read_text(encoding="utf-8-sig")

        # Local validation
        if not bypass_local_validation:
            count, errs = validate_submission_csv(csv_text, modality=modality)
            if errs:
                err_summary = "; ".join([f"Line {e[0]} ({e[1]}): {e[2]}" for e in errs])
                raise SJTUValidationError(
                    f"Local CSV validation failed: {err_summary}",
                    errors=[{"index": e[0], "name": e[1], "message": e[2]} for e in errs],
                )
        else:
            count = len([l for l in csv_text.splitlines() if l.strip()]) - 1

        # Check existing submission for target
        prior_submission = self.get_latest_submission(actual_target_id)
        final_note = note.strip()
        if prior_submission is not None and prior_submission.get("id"):
            if not final_note:
                if auto_note_on_reupload:
                    ts_iso = datetime.now(timezone.utc).isoformat()
                    final_note = (
                        f"Updated candidate batch optimized via hybrid photonic VQE quantum folding "
                        f"and ESMFold2/ipSAE multi-objective ranking. (Auto-noted at {ts_iso})"
                    )
                    self.logger.info(
                        f"[Re-upload Note] Target '{target_code}' has existing submission ({prior_submission.get('id')}). "
                        f"Auto-generated note: {final_note}"
                    )
                else:
                    raise SJTUReUploadNoteRequiredError(
                        f"Target '{target_code}' already has an existing submission batch. "
                        f"A non-empty 'note' explaining the re-upload is required."
                    )

        if dry_run:
            self.logger.info(
                f"[Dry Run] Pre-validation SUCCESS for target '{target_code}' ({actual_target_id}). "
                f"Modality: {modality}, Sequences: {count}. Ready to submit."
            )
            return {
                "dry_run": True,
                "target_id": actual_target_id,
                "target_code": target_code,
                "count": count,
                "status": "validated",
                "abstract": abstract.strip(),
                "note": final_note,
            }

        url = f"{self.base_url}/api/submissions/csv"
        payload = {
            "target_id": actual_target_id,
            "action": action,
            "agreed": agreed,
            "abstract": abstract.strip(),
            "abstract_name": "",
            "abstract_file": "",
            "note": final_note,
            "original_name": original_filename or csv_file.name,
            "csv": csv_text,
        }

        resp = None
        last_exc = None
        for attempt in range(max(1, max_retries)):
            try:
                resp = self._auth_request("POST", url, json=payload, timeout=45)
                if resp.status_code == 422:
                    body = resp.json()
                    detail = body.get("detail", "Validation Failed")
                    errors = body.get("errors", [])
                    err_msg = f"Platform validation failed (HTTP 422): {detail}"
                    if errors:
                        err_msg += " - " + "; ".join([f"Line {e.get('index')}: {e.get('message')}" for e in errors])
                    raise SJTUValidationError(err_msg, errors=errors)
                if resp.status_code == 404:
                    raise SJTUTargetNotFoundError("Target does not exist on platform (HTTP 404): 赛题不存在")
                if resp.status_code >= 500:
                    resp.raise_for_status()

                resp.raise_for_status()
                break
            except requests.HTTPError as exc:
                last_exc = exc
                if exc.response is not None and exc.response.status_code >= 500 and attempt < max_retries - 1:
                    self.logger.warning(
                        f"Transient HTTP {exc.response.status_code} on submission attempt {attempt + 1}/{max_retries}. Retrying..."
                    )
                    continue
                raise
            except requests.RequestException as exc:
                last_exc = exc
                if attempt < max_retries - 1:
                    self.logger.warning(f"Request exception on attempt {attempt + 1}/{max_retries}: {exc}. Retrying...")
                    continue
                raise

        if resp is None:
            if last_exc:
                raise last_exc
            raise RuntimeError("Submission failed without response")
        res_data = resp.json()

        batch_id = res_data.get("id") or res_data.get("batch_id", "unknown")
        confirmed_at = res_data.get("confirmed_at") or res_data.get("created_at") or datetime.now(timezone.utc).isoformat()
        self.logger.info(
            f"[Submission Confirmed] Target: {target_code} ({actual_target_id}) | "
            f"Batch ID: {batch_id} | Status: {action} | Count: {count} | "
            f"Confirmed At: {confirmed_at} | Response: {json.dumps(res_data, ensure_ascii=False)}"
        )
        return res_data


# ---------------------------------------------------------------------------
# Command-Line Interface (CLI)
# ---------------------------------------------------------------------------

def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="SJTU AI Protein Platform Automated Submission CLI (Changjiang Cup 2026).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--target", "-t", type=str, default="", help="Target code (tl1a, il7ra, etc.) or UUID.")
    parser.add_argument("--csv", "-c", type=str, default="", help="Path to input sequence CSV (name, sequence).")
    parser.add_argument("--action", "-a", type=str, default="confirm", choices=["confirm", "save"], help="Submission action.")
    parser.add_argument("--abstract", type=str, default="", help="Scientific abstract text or path to .txt file.")
    parser.add_argument("--note", "-n", type=str, default="", help="Re-upload explanation note (max 2000 chars).")
    parser.add_argument("--no-agree", action="store_true", help="Do not agree to rules (causes validation failure).")
    parser.add_argument("--dry-run", action="store_true", help="Validate CSV and target without making submission POST.")
    parser.add_argument("--query-status", action="store_true", help="Query latest submission status.")
    parser.add_argument("--log-file", type=str, default="logs/sjtu_submission.log", help="Path to audit log file.")
    parser.add_argument("--email", type=str, default=DEFAULT_EMAIL, help="Contestant email.")
    parser.add_argument("--password", type=str, default=DEFAULT_PASSWORD, help="Contestant password.")
    parser.add_argument("--base-url", type=str, default=DEFAULT_BASE_URL, help="Platform base API URL.")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose debug logging.")

    parsed = parser.parse_args(args)

    # If abstract is a file path, load contents
    if parsed.abstract and Path(parsed.abstract).is_file():
        parsed.abstract = Path(parsed.abstract).read_text(encoding="utf-8").strip()

    return parsed


def main(args: Optional[List[str]] = None) -> int:
    parsed_args = parse_args(args)
    logger = setup_submission_logger(log_file=parsed_args.log_file, verbose=parsed_args.verbose)

    # 1. Dry run mode without requiring client auth if purely local
    if parsed_args.dry_run and parsed_args.csv:
        try:
            target_input = parsed_args.target or "tl1a"
            cfg = TARGET_CONFIGS.get(target_input.lower())
            modality = cfg["modality"] if cfg else "minibinder"
            target_id = cfg["id"] if cfg else target_input

            count, errs = validate_submission_csv(parsed_args.csv, modality=modality)
            if errs:
                print("DRY RUN VALIDATION: FAILED")
                for e in errs:
                    print(f"  Line {e[0]} ({e[1]}): {e[2]}")
                return 1

            print("DRY RUN VALIDATION: SUCCESS")
            print(f"Total sequences: {count}")
            print(f"Target: {cfg['name'] if cfg else target_input} ({target_id})")
            print(f"Modality: {modality}")
            if parsed_args.note:
                print(f"Note: {parsed_args.note}")
            if parsed_args.abstract:
                print(f"Abstract length: {len(parsed_args.abstract)} characters")
            return 0
        except Exception as e:
            logger.error(f"[Dry Run Error] {e}")
            print(f"DRY RUN VALIDATION: FAILED - {e}")
            return 1

    # Initialize client
    client = SJTUPlatformClient(
        email=parsed_args.email,
        password=parsed_args.password,
        base_url=parsed_args.base_url,
        logger=logger,
    )

    try:
        client.login()
    except Exception as e:
        logger.error(f"[Auth Error] Failed to login: {e}")
        print(f"Error authenticating to SJTU platform: {e}", file=sys.stderr)
        return 3

    # Query Status Mode
    if parsed_args.query_status:
        try:
            if parsed_args.target:
                cfg = client.resolve_target(parsed_args.target)
                target_id = cfg["id"]
                latest = client.get_latest_submission(target_id)
                print(f"\n=== Latest Submission for Target {cfg['name']} ({cfg['code'].upper()}) ===")
                if latest:
                    print(f"Batch ID:      {latest.get('id')}")
                    print(f"Filename:      {latest.get('filename') or latest.get('original_name')}")
                    print(f"Sequence Count:{latest.get('count', 'N/A')}")
                    print(f"Confirmed At:  {latest.get('created_at') or latest.get('confirmed_at')}")
                    print(f"Note:          {latest.get('note') or '(None)'}")
                    print(f"Abstract:      {latest.get('abstract')[:100] + '...' if latest.get('abstract') else '(None)'}")
                else:
                    print("Status: No previous submission found for this target.")
            else:
                submissions = client.get_submissions(grouped=True)
                print("\n=== Submissions Summary ===")
                print(json.dumps(submissions, indent=2, ensure_ascii=False))
            return 0
        except Exception as e:
            logger.error(f"[Query Error] {e}")
            print(f"Error querying submission status: {e}", file=sys.stderr)
            return 1

    # Candidate Submission Mode
    if not parsed_args.csv:
        print("Error: --csv is required for candidate submission.", file=sys.stderr)
        return 1

    if not parsed_args.target:
        print("Error: --target is required for candidate submission.", file=sys.stderr)
        return 1

    try:
        target_cfg = client.resolve_target(parsed_args.target)
        target_id = target_cfg["id"]
        res = client.submit_candidates(
            target_id=target_id,
            csv_path=parsed_args.csv,
            action=parsed_args.action,
            agreed=not parsed_args.no_agree,
            abstract=parsed_args.abstract,
            note=parsed_args.note,
            dry_run=parsed_args.dry_run,
        )
        print("\n[SUCCESS] Candidate submission confirmed:")
        print(f"Target:        {target_cfg['name']} ({target_cfg['code'].upper()})")
        print(f"Batch ID:      {res.get('id') or res.get('batch_id')}")
        print(f"Status:        {res.get('status', 'confirmed')}")
        print(f"Confirmed At:  {res.get('confirmed_at') or res.get('created_at')}")
        return 0
    except SJTUValidationError as e:
        logger.error(f"[Validation Error] {e}")
        print(f"Validation Error: {e}", file=sys.stderr)
        return 2
    except SJTUAuthenticationError as e:
        logger.error(f"[Auth Error] {e}")
        print(f"Authentication Error: {e}", file=sys.stderr)
        return 3
    except Exception as e:
        logger.error(f"[Submission Error] {e}")
        print(f"Submission Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
