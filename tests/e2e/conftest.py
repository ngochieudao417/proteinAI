"""
tests/e2e/conftest.py - Shared E2E Test Fixtures, Mock Adapters, and Network Isolation.
Clean architecture without sys.modules virtual injection or dummy implementations.
"""

from __future__ import annotations

import csv
import io
import json
import logging
import os
import re
import sys
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import pytest
import requests
from requests.adapters import HTTPAdapter

# ---------------------------------------------------------------------------
# 1. System Path Configuration & Genuine Module Imports
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import perceval as pcvl
import quantum_folding
from quantum_folding.quandela_backend import QuandelaBackend, QuandelaCloudClient
from tools.submit_sjtu import (
    COMPETITION_TARGETS,
    TARGETS_BY_CODE,
    TARGETS_BY_ID,
    SJTUPlatformClient,
    parse_csv_rows,
    validate_sequence,
    validate_submission_csv,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# 2. Authoritative Constants
# ---------------------------------------------------------------------------
_env_path = Path(__file__).resolve().parent.parent.parent / ".env"
if _env_path.exists():
    try:
        with open(_env_path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    k, v = line.split("=", 1)
                    os.environ.setdefault(k.strip(), v.strip())
    except Exception:
        pass

VALID_QUANDELA_TOKEN = os.environ.get(
    "QUANDELA_TOKEN",
    "_T_eyJhbGciOiJIUzUxMiIsInR5cCI6IkpXVCJ9.dummy_test_payload.dummy_signature"
)

SJTU_EMAIL = os.environ.get("SJTU_EMAIL", "test_contestant@example.com")
SJTU_PASSWORD = os.environ.get("SJTU_PASSWORD", "test_contestant_password")
SJTU_BASE_URL = os.environ.get("SJTU_BASE_URL", "https://aiprotein.sjtu.edu.cn")

TL1A_TARGET_ID = "ea930a01-20a9-4bed-b70f-e658bd9ba3ab"
IL7RA_TARGET_ID = "1347d0e0-7d3d-4aa8-bcf1-c2dab9d8f1df"

CANONICAL_AMINO_ACIDS: Set[str] = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC_RESIDUES: Set[str] = set("ACFILMVWY")
POLAR_RESIDUES: Set[str] = set("DEGHIKNPQRST")
MJ_AA_ORDER: str = "ACDEFGHIKLMNPQRSTVWY"


# ---------------------------------------------------------------------------
# 3. In-Process Mock HTTPAdapter for SJTU Platform API
# ---------------------------------------------------------------------------
class MockSJTUPlatformAdapter(HTTPAdapter):
    """
    In-process HTTPAdapter simulating https://aiprotein.sjtu.edu.cn API.
    Enforces competition rules: sequence counts, modality length bounds,
    canonical amino acids, quote-safe CSV parsing, and metadata size limits.
    """

    def __init__(self):
        super().__init__()
        self.submissions_store: Dict[str, Dict[str, Any]] = {}
        self.auth_tokens: Set[str] = set()

    def send(self, request, **kwargs):
        url = request.url
        method = request.method
        resp = requests.Response()
        resp.request = request
        resp.url = url
        resp.encoding = "utf-8"

        def _json_resp(status_code: int, data: Any):
            resp.status_code = status_code
            resp.headers["Content-Type"] = "application/json; charset=utf-8"
            resp._content = json.dumps(data, ensure_ascii=False).encode("utf-8")
            return resp

        # 1. Login endpoint
        if "/api/auth/login" in url and method == "POST":
            body = json.loads(request.body.decode("utf-8") if request.body else "{}")
            if body.get("email") == SJTU_EMAIL and body.get("password") == SJTU_PASSWORD:
                tok = "mock-jwt-" + str(uuid.uuid4())
                ref = "mock-ref-" + str(uuid.uuid4())
                self.auth_tokens.add(tok)
                return _json_resp(200, {
                    "access_token": tok,
                    "refresh_token": ref,
                    "token_type": "bearer",
                    "user": {
                        "id": "07fe1517-b035-4805-a505-9a3fc68547fd",
                        "email": body["email"],
                        "name": "Mai Duc Khang",
                        "role": "contestant",
                        "team": {
                            "id": "66d2100d-6eb1-42ba-9dcd-0589b0654046",
                            "name": "QAAA Lab",
                            "member_role": "captain",
                        },
                    },
                })
            else:
                return _json_resp(401, {"detail": "用户名或密码错误"})

        # 2. Token refresh endpoint
        if "/api/auth/refresh" in url and method == "POST":
            new_tok = "refreshed-jwt-" + str(uuid.uuid4())
            self.auth_tokens.add(new_tok)
            return _json_resp(200, {
                "access_token": new_tok,
                "refresh_token": "new-refresh-tok",
                "token_type": "bearer",
            })

        # 3. Targets endpoint
        if "/api/submissions/targets" in url and method == "GET":
            return _json_resp(200, COMPETITION_TARGETS)

        # 4. CSV Submission endpoint
        if "/api/submissions/csv" in url and method == "POST":
            body = json.loads(request.body.decode("utf-8") if request.body else "{}")
            if not body.get("agreed", False):
                return _json_resp(422, {"detail": "请先勾选并同意赛事细则"})

            target_id = body.get("target_id")
            if target_id not in TARGETS_BY_ID and target_id not in TARGETS_BY_CODE:
                return _json_resp(404, {"detail": "赛题不存在"})

            target_info = TARGETS_BY_ID.get(target_id) or TARGETS_BY_CODE.get(target_id)
            min_len = target_info["min_length"]
            max_len = target_info["max_length"]
            modality = target_info["modality"]

            abstract = body.get("abstract", "")
            if abstract and len(abstract) > 4000:
                return _json_resp(422, {"detail": "摘要超过最大长度限制（最多4000字符）"})

            note = body.get("note", "")
            if note and len(note) > 2000:
                return _json_resp(422, {"detail": "备注超过最大长度限制（最多2000字符）"})

            csv_text = body.get("csv", "")
            rows = parse_csv_rows(csv_text)
            if not rows:
                return _json_resp(422, {
                    "detail": "批量校验未通过",
                    "errors": [{"index": 1, "name": "", "message": "没有序列"}],
                })

            hdr = [h.lower() for h in rows[0]]
            if "name" not in hdr or "sequence" not in hdr:
                return _json_resp(422, {
                    "detail": "批量校验未通过",
                    "errors": [{"index": 1, "name": "", "message": "表头必须包含 name 和 sequence"}],
                })

            seq_count = len(rows) - 1
            if seq_count == 0:
                return _json_resp(422, {
                    "detail": "批量校验未通过",
                    "errors": [{"index": 1, "name": "", "message": "没有序列"}],
                })

            if seq_count > 20:
                return _json_resp(422, {
                    "detail": "批量校验未通过",
                    "errors": [{"index": 22, "name": "", "message": "最多提交 20 条序列（含表头共 21 行）"}],
                })

            name_idx = hdr.index("name")
            seq_idx = hdr.index("sequence")
            seen_names = set()

            for row_idx, parts in enumerate(rows[1:], start=2):
                name_val = parts[name_idx] if name_idx < len(parts) else ""
                seq_val = re.sub(r"\s+", "", parts[seq_idx] if seq_idx < len(parts) else "").upper()

                if not name_val or len(name_val) > 64:
                    return _json_resp(422, {
                        "detail": "批量校验未通过",
                        "errors": [{"index": row_idx, "name": name_val, "message": "名称须为 1 到 64 个字符"}],
                    })

                if name_val in seen_names:
                    return _json_resp(422, {
                        "detail": "批量校验未通过",
                        "errors": [{"index": row_idx, "name": name_val, "message": "名称重复"}],
                    })
                seen_names.add(name_val)

                if any(c not in CANONICAL_AMINO_ACIDS for c in seq_val):
                    return _json_resp(422, {
                        "detail": "批量校验未通过",
                        "errors": [{"index": row_idx, "name": name_val, "message": "含有非法氨基酸"}],
                    })

                # Strict sequence length enforcement against target track
                if len(seq_val) < min_len or len(seq_val) > max_len:
                    return _json_resp(422, {
                        "detail": "批量校验未通过",
                        "errors": [{
                            "index": row_idx,
                            "name": name_val,
                            "message": f"序列长度不符合要求：{modality} 为 {min_len}—{max_len} aa",
                        }],
                    })

            batch_id = "batch-" + str(uuid.uuid4())
            record = {
                "id": batch_id,
                "target_id": target_id,
                "target_code": target_info["code"],
                "count": seq_count,
                "original_name": body.get("original_name", "submission.csv"),
                "confirmed_at": "2026-11-11T16:00:00+08:00",
                "status": "confirmed",
                "abstract": abstract,
                "note": note,
            }
            self.submissions_store[target_id] = record
            return _json_resp(200, record)

        # 5. Latest submission query
        if "/api/submissions/latest" in url and method == "GET":
            m = re.search(r"target_id=([^&]+)", url)
            target_id = m.group(1) if m else None
            record = self.submissions_store.get(target_id)
            return _json_resp(200, record)

        # 6. All submissions query
        if "/api/submissions" in url and method == "GET":
            items = list(self.submissions_store.values())
            return _json_resp(200, {
                "groups": [{"status": "submitted", "total": len(items), "items": items}]
            })

        return _json_resp(404, {"detail": "Not Found"})


# ---------------------------------------------------------------------------
# 4. Shared Pytest Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def mock_quandela_token() -> str:
    """Valid Quandela Cloud API token with mandatory '_T_' prefix."""
    return VALID_QUANDELA_TOKEN


@pytest.fixture
def mock_cloud_client() -> QuandelaBackend:
    """QuandelaBackend configured for deterministic simulation with fallback."""
    client = QuandelaBackend(token=VALID_QUANDELA_TOKEN, platform="sim:slos", local_fallback=True)
    client._simulate_local = True
    return client


@pytest.fixture(autouse=True)
def mock_sjtu_api(monkeypatch):
    """
    Mounts MockSJTUPlatformAdapter on requests.Session unless TEST_SJTU_LIVE=1 is set.
    """
    if os.environ.get("TEST_SJTU_LIVE") == "1":
        yield None
        return

    adapter = MockSJTUPlatformAdapter()
    orig_init = requests.Session.__init__

    def patched_init(self, *args, **kwargs):
        orig_init(self, *args, **kwargs)
        self.mount("https://aiprotein.sjtu.edu.cn", adapter)
        self.mount("http://aiprotein.sjtu.edu.cn", adapter)

    monkeypatch.setattr(requests.Session, "__init__", patched_init)
    yield adapter


@pytest.fixture
def sample_peptides() -> Dict[str, str]:
    """Standardized representative peptide sequences for testing."""
    return {
        "short_hp": "HPPH",
        "short_mj": "ACDEFY",
        "tl1a_sample": "AALEELLKALEELLKALEELLKALEPWHETVMSINLFMIEFIIAHPSPKEIQITWELKDFLKKNFEPVSKSFHDRLFPFGDYKWRSPLEIGK",
        "il7ra_sample": "AALEELLKALEELLKALEKLAPKRTQHADWTNSQIIQERVHEQDYRSDINIIFNIQYWPGLFGEVAQANVENDTTYHY",
        "homo_ala": "AAAAAAAA",
        "homo_gly": "GGGGGGGG",
        "hydrophobic": "LLLLLL",
        "polar": "KKKKKK",
        "boundary_40": "A" * 40,
        "boundary_85": "A" * 85,
        "boundary_100": "A" * 100,
        "boundary_150": "A" * 150,
        "underflow_39": "A" * 39,
        "overflow_101": "A" * 101,
        "invalid_aa": "ACDEFGHIKLMNPQRSTVWX",
    }


@pytest.fixture
def temp_candidate_csv(tmp_path: Path) -> str:
    """Generate a valid candidate pool CSV for testing."""
    csv_file = tmp_path / "candidate_pool.csv"
    csv_file.write_text(
        "name,sequence,target,length,epitope_group,generator,backbone_id\n"
        "CAND_01,AALEELLKALEELLKALEELLKALEPWHETVMSINLFMIEFIIAHPSPKEIQITWELKDFLKKNFEPVSKSFHDRLFPFGDYKWRSPLEIGK,TL1A,93,A,gen,bb_1\n"
        "CAND_02,AALEELLKALEELLKALEKLAPKRTQHADWTNSQIIQERVHEQDYRSDINIIFNIQYWPGLFGEVAQANVENDTTYHY,IL7RA,78,B,gen,bb_2\n"
        "CAND_03,AALEELLKALEELLKALEELLKALEPW,TL1A,26,A,gen,bb_3\n"
        "CAND_04,EKLRKLEKQARKQGDEVLVMLARMVLE,TL1A,26,A,gen,bb_4\n"
        "CAND_05,AALEELLKALEELLKALEKLAPKRTQH,IL7RA,26,B,gen,bb_5\n"
    )
    return str(csv_file)


@pytest.fixture
def temp_output_dir(tmp_path: Path) -> str:
    """Create a temporary output directory."""
    out_dir = tmp_path / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    return str(out_dir)
