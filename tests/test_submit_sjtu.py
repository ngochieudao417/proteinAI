"""
Unit and Integration Test Suite for tools/submit_sjtu.py.

Covers:
1. Local sequence & CSV validation (canonical AAs, modality lengths, header parsing, unique names, row count boundaries).
2. Mocked unit tests for SJTUPlatformClient (login, 401/422 handling, token refresh, target discovery, CSV submission, status query).
3. CLI argument parsing, dry-run simulation, query-status mode, and audit logging.
4. Opt-in live integration tests against https://aiprotein.sjtu.edu.cn (gated by TEST_SJTU_LIVE=1).
"""

from __future__ import annotations

import csv
import io
import json
import logging
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional
from unittest.mock import MagicMock, patch

import pytest
import requests
from requests.adapters import HTTPAdapter

# Ensure repo root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.submit_sjtu import (
    CANONICAL_AA,
    TARGET_CONFIGS,
    SJTUAuthenticationError,
    SJTUPlatformClient,
    SJTUPlatformError,
    SJTUReUploadNoteRequiredError,
    SJTUTargetNotFoundError,
    SJTUValidationError,
    main,
    parse_args,
    setup_submission_logger,
    validate_sequence,
    validate_submission_csv,
)

TL1A_UUID = "ea930a01-20a9-4bed-b70f-e658bd9ba3ab"
IL7RA_UUID = "1347d0e0-7d3d-4aa8-bcf1-c2dab9d8f1df"
EGFR_UUID = "af6b1e2b-c62d-40e4-a5b3-3f091be04c60"


# ============================================================================
# Section 1: Validation Functions Unit Tests
# ============================================================================

class TestSequenceValidation:
    """Test suite for validate_sequence()."""

    def test_canonical_20_amino_acids_accepted(self):
        """All 20 standard amino acids pass validation."""
        valid_60aa = "ACDEFGHIKLMNPQRSTVWY" * 3
        ok, err = validate_sequence(valid_60aa, min_len=40, max_len=100)
        assert ok is True
        assert err is None

    @pytest.mark.parametrize("bad_char", ["B", "J", "O", "U", "X", "Z", "1", "9", "*", "-", "_", "$"])
    def test_invalid_characters_rejected(self, bad_char: str):
        """Non-canonical letters, numbers, and symbols are rejected."""
        seq = "ACDEFGHIKL" + bad_char + "MNPQRSTVWY" + "A" * 25
        ok, err = validate_sequence(seq, min_len=40, max_len=100)
        assert ok is False
        assert "含有非法氨基酸" in err

    def test_whitespace_and_lowercase_normalization(self):
        """Whitespace is removed and lowercase characters are uppercased."""
        seq = "  aal  eell  kale  ellk  alep  whet  vmsi  nlfm  iefi  iahp  spke  iqit  "
        ok, err = validate_sequence(seq, min_len=40, max_len=100)
        assert ok is True
        assert err is None

    def test_empty_sequence_rejected(self):
        """Empty or whitespace-only sequence is rejected."""
        for empty_val in ["", "   ", "\n\t"]:
            ok, err = validate_sequence(empty_val, min_len=40, max_len=100)
            assert ok is False
            assert "序列为空" in err

    @pytest.mark.parametrize("length,expected_ok", [
        (39, False),   # Underflow
        (40, True),    # Minibinder lower bound
        (75, True),    # Nominal minibinder
        (100, True),   # Minibinder upper bound
        (101, False),  # Overflow
    ])
    def test_minibinder_length_boundaries(self, length: int, expected_ok: bool):
        """Minibinder modality strictly enforces 40-100 aa."""
        seq = "A" * length
        ok, err = validate_sequence(seq, min_len=40, max_len=100)
        assert ok is expected_ok
        if not expected_ok:
            assert "序列长度不符合要求" in err

    @pytest.mark.parametrize("length,expected_ok", [
        (99, False),   # Underflow
        (100, True),   # Nanobody lower bound
        (125, True),   # Nominal nanobody
        (150, True),   # Nanobody upper bound
        (151, False),  # Overflow
    ])
    def test_nanobody_length_boundaries(self, length: int, expected_ok: bool):
        """Nanobody (VHH) modality strictly enforces 100-150 aa."""
        seq = "A" * length
        ok, err = validate_sequence(seq, min_len=100, max_len=150)
        assert ok is expected_ok
        if not expected_ok:
            assert "序列长度不符合要求" in err


class TestCSVValidation:
    """Test suite for validate_submission_csv()."""

    def test_case_insensitive_header(self):
        """Headers with varying case are accepted."""
        for hdr in ["name,sequence", "NAME,SEQUENCE", "Name,Sequence", "NaMe,SeQuEnCe"]:
            csv_text = f"{hdr}\ncand_1," + "A" * 50 + "\n"
            count, errs = validate_submission_csv(csv_text, modality="minibinder")
            assert count == 1
            assert len(errs) == 0

    def test_extra_columns_tolerated(self):
        """Extra metadata columns do not break parsing."""
        csv_text = "name,sequence,score,rank,notes\ncand_1," + "A" * 50 + ",0.95,1,test\n"
        count, errs = validate_submission_csv(csv_text, modality="minibinder")
        assert count == 1
        assert len(errs) == 0

    def test_missing_header_columns_rejected(self):
        """Missing 'name' or 'sequence' column triggers line 1 error."""
        bad_headers = [
            "id,sequence\n1," + "A" * 50,
            "name,seq\ncand_1," + "A" * 50,
            "seq_only\n" + "A" * 50,
        ]
        for bh in bad_headers:
            count, errs = validate_submission_csv(bh, modality="minibinder")
            assert len(errs) >= 1
            assert errs[0][0] == 1  # Line 1 error
            assert "表头必须包含 name 和 sequence" in errs[0][2]

    def test_utf8_bom_handling(self):
        """Leading UTF-8 BOM is stripped cleanly."""
        csv_text = "\ufeffname,sequence\ncand_1," + "A" * 50 + "\n"
        count, errs = validate_submission_csv(csv_text, modality="minibinder")
        assert count == 1
        assert len(errs) == 0

    def test_duplicate_names_rejected(self):
        """Duplicate sequence names within file trigger rejection."""
        csv_text = (
            "name,sequence\n"
            "cand_dup," + "A" * 50 + "\n"
            "cand_dup," + "A" * 55 + "\n"
        )
        count, errs = validate_submission_csv(csv_text, modality="minibinder")
        assert count == 2
        assert any("名称重复" in e[2] for e in errs)

    def test_sequence_name_length_boundaries(self):
        """Names must be 1 to 64 characters long."""
        # Valid bounds
        valid_csv = f"name,sequence\na,{'A'*50}\n{'B'*64},{'A'*50}\n"
        count, errs = validate_submission_csv(valid_csv, modality="minibinder")
        assert len(errs) == 0

        # Invalid: 65 chars
        invalid_csv = f"name,sequence\n{'B'*65},{'A'*50}\n"
        count, errs = validate_submission_csv(invalid_csv, modality="minibinder")
        assert any("名称须为 1 到 64 个字符" in e[2] for e in errs)

    def test_row_count_boundaries(self):
        """Empty, 1 row, 20 rows (accepted), 21 rows (rejected)."""
        # Empty file
        count, errs = validate_submission_csv("", modality="minibinder")
        assert any("没有序列" in e[2] for e in errs)

        # Header only
        count, errs = validate_submission_csv("name,sequence\n", modality="minibinder")
        assert any("没有序列" in e[2] for e in errs)

        # 1 sequence row
        csv_1 = "name,sequence\nc1," + "A" * 50 + "\n"
        count, errs = validate_submission_csv(csv_1, modality="minibinder")
        assert count == 1
        assert len(errs) == 0

        # 20 sequence rows (max limit)
        lines_20 = ["name,sequence"] + [f"c_{i},{'A'*50}" for i in range(20)]
        count, errs = validate_submission_csv("\n".join(lines_20), modality="minibinder")
        assert count == 20
        assert len(errs) == 0

        # 21 sequence rows (exceeds limit)
        lines_21 = ["name,sequence"] + [f"c_{i},{'A'*50}" for i in range(21)]
        count, errs = validate_submission_csv("\n".join(lines_21), modality="minibinder")
        assert any("最多提交 20 条序列" in e[2] for e in errs)


# ============================================================================
# Section 2: Mocked SJTUPlatformClient Unit Tests
# ============================================================================

class MockSJTUAdapter(HTTPAdapter):
    """In-memory HTTP adapter for SJTUPlatformClient testing."""

    def __init__(self):
        super().__init__()
        self.auth_token = "mock-jwt-token-12345"
        self.refresh_token = "mock-refresh-token-67890"
        self.submissions: Dict[str, Dict[str, Any]] = {}
        self.token_expired_once = False
        self.simulate_401_on_targets = False

    def send(self, request, **kwargs):
        url = request.url
        method = request.method
        resp = requests.Response()
        resp.request = request
        resp.url = url

        # 1. Login
        if "/api/auth/login" in url and method == "POST":
            body = json.loads(request.body.decode("utf-8") if request.body else "{}")
            if body.get("email") and body.get("password") and body.get("password") != "WRONG_PASSWORD":
                resp.status_code = 200
                resp._content = json.dumps({
                    "access_token": self.auth_token,
                    "refresh_token": self.refresh_token,
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
                }).encode("utf-8")
            else:
                resp.status_code = 401
                resp._content = json.dumps({"detail": "用户名或密码错误"}).encode("utf-8")
            return resp

        # 2. Refresh Token
        if "/api/auth/refresh" in url and method == "POST":
            body = json.loads(request.body.decode("utf-8") if request.body else "{}")
            if body.get("refresh_token") == self.refresh_token:
                self.auth_token = "mock-refreshed-jwt-99999"
                resp.status_code = 200
                resp._content = json.dumps({
                    "access_token": self.auth_token,
                    "refresh_token": "mock-new-refresh-00000",
                    "token_type": "bearer",
                }).encode("utf-8")
            else:
                resp.status_code = 401
                resp._content = json.dumps({"detail": "Invalid refresh token"}).encode("utf-8")
            return resp

        # 3. Targets
        if "/api/submissions/targets" in url and method == "GET":
            if self.simulate_401_on_targets and not self.token_expired_once:
                self.token_expired_once = True
                resp.status_code = 401
                resp._content = json.dumps({"detail": "Could not validate credentials"}).encode("utf-8")
                return resp

            resp.status_code = 200
            resp._content = json.dumps([
                {"id": TL1A_UUID, "code": "tl1a", "name": "TL1A", "modality": "minibinder", "min_length": 40, "max_length": 100, "max_sequences": 20},
                {"id": IL7RA_UUID, "code": "il7ra", "name": "IL-7Rα", "modality": "minibinder", "min_length": 40, "max_length": 100, "max_sequences": 20},
                {"id": EGFR_UUID, "code": "egfr", "name": "EGFR", "modality": "nanobody", "min_length": 100, "max_length": 150, "max_sequences": 20},
            ]).encode("utf-8")
            return resp

        # 4. CSV Submission
        if "/api/submissions/csv" in url and method == "POST":
            body = json.loads(request.body.decode("utf-8") if request.body else "{}")
            if not body.get("agreed", False):
                resp.status_code = 422
                resp._content = json.dumps({"detail": "请先勾选并同意赛事细则"}).encode("utf-8")
                return resp

            target_id = body.get("target_id")
            if target_id not in [TL1A_UUID, IL7RA_UUID, EGFR_UUID]:
                resp.status_code = 404
                resp._content = json.dumps({"detail": "赛题不存在"}).encode("utf-8")
                return resp

            csv_text = body.get("csv", "")
            if "INVALID_AA_TRIGGER" in csv_text:
                resp.status_code = 422
                resp._content = json.dumps({
                    "detail": "批量校验未通过",
                    "errors": [{"index": 2, "name": "c1", "message": "含有非法氨基酸"}],
                }).encode("utf-8")
                return resp

            batch_id = "batch-mock-" + target_id[:8]
            record = {
                "id": batch_id,
                "target_id": target_id,
                "count": 20,
                "original_name": body.get("original_name", "candidates.csv"),
                "confirmed_at": "2026-11-12T00:15:00+08:00",
                "status": "confirmed",
                "abstract": body.get("abstract", ""),
                "note": body.get("note", ""),
            }
            self.submissions[target_id] = record
            resp.status_code = 200
            resp._content = json.dumps(record).encode("utf-8")
            return resp

        # 5. Latest Submission Query
        if "/api/submissions/latest" in url and method == "GET":
            m = re.search(r"target_id=([^&]+)", url)
            tid = m.group(1) if m else None
            record = self.submissions.get(tid)
            resp.status_code = 200
            resp._content = json.dumps(record).encode("utf-8")
            return resp

        # 6. Grouped Submissions Query
        if "/api/submissions" in url and method == "GET":
            items = list(self.submissions.values())
            resp.status_code = 200
            resp._content = json.dumps({
                "groups": [{"status": "submitted", "total": len(items), "page": 1, "page_size": 10, "items": items}]
            }).encode("utf-8")
            return resp

        resp.status_code = 404
        resp._content = b'{"detail": "Not Found"}'
        return resp


@pytest.fixture
def mock_client():
    """Client fixture with mounted MockSJTUAdapter."""
    adapter = MockSJTUAdapter()
    client = SJTUPlatformClient()
    client.session.mount("https://aiprotein.sjtu.edu.cn", adapter)
    client.session.mount("http://aiprotein.sjtu.edu.cn", adapter)
    client._adapter = adapter
    return client


class TestSJTUPlatformClientMocked:
    """Mocked unit tests for SJTUPlatformClient."""

    def test_login_success(self, mock_client):
        """Successful login populates token and headers."""
        data = mock_client.login()
        assert mock_client.access_token == "mock-jwt-token-12345"
        assert mock_client.refresh_token == "mock-refresh-token-67890"
        assert mock_client.session.headers["Authorization"] == "Bearer mock-jwt-token-12345"
        assert data["user"]["team"]["name"] == "QAAA Lab"
        assert data["user"]["team"]["member_role"] == "captain"

    def test_login_failure_401(self, mock_client):
        """Invalid credentials raise PermissionError (SJTUAuthenticationError)."""
        mock_client.password = "WRONG_PASSWORD"
        with pytest.raises(PermissionError, match="用户名或密码错误"):
            mock_client.login()

    def test_token_refresh_mechanism(self, mock_client):
        """Token refresh endpoint renews access token and header."""
        mock_client.login()
        old_token = mock_client.access_token
        res = mock_client.refresh_token_session()
        assert res["access_token"] != old_token
        assert mock_client.access_token == res["access_token"]
        assert mock_client.session.headers["Authorization"] == f"Bearer {res['access_token']}"

    def test_auto_token_refresh_on_401_interception(self, mock_client):
        """Client transparently refreshes token on HTTP 401 and retries."""
        mock_client.login()
        mock_client._adapter.simulate_401_on_targets = True

        targets = mock_client.get_targets()
        assert len(targets) == 3
        assert mock_client.access_token == "mock-refreshed-jwt-99999"

    def test_target_discovery_and_resolution(self, mock_client):
        """resolve_target maps names and case-insensitive codes to UUIDs."""
        mock_client.login()
        tl1a_target = mock_client.resolve_target("tl1a")
        assert tl1a_target["id"] == TL1A_UUID
        assert tl1a_target["modality"] == "minibinder"

        il7ra_target = mock_client.resolve_target("IL7RA")
        assert il7ra_target["id"] == IL7RA_UUID

        with pytest.raises(ValueError, match="Unknown target"):
            mock_client.resolve_target("non_existent_target")

    def test_submit_candidates_success(self, mock_client, tmp_path):
        """Submit candidates uploads valid CSV and returns batch record."""
        mock_client.login()
        csv_file = tmp_path / "valid_tl1a.csv"
        csv_file.write_text("name,sequence\n" + "\n".join([f"c_{i},{'A'*50}" for i in range(20)]))

        res = mock_client.submit_candidates(
            target_id=TL1A_UUID,
            csv_path=str(csv_file),
            action="confirm",
            agreed=True,
            abstract="Quantum folding optimization abstract",
            note="Initial submission batch",
        )
        assert res["status"] == "confirmed"
        assert res["target_id"] == TL1A_UUID
        assert res["count"] == 20
        assert "batch-mock" in res["id"]

    def test_submit_candidates_unagreed_local_rejection(self, mock_client, tmp_path):
        """Unagreed submission is rejected client-side before network call."""
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("name,sequence\nc1," + "A" * 50)
        with pytest.raises(ValueError, match="请先勾选并同意赛事细则"):
            mock_client.submit_candidates(
                target_id=TL1A_UUID,
                csv_path=str(csv_file),
                agreed=False,
            )

    def test_submit_candidates_server_validation_422(self, mock_client, tmp_path):
        """Server HTTP 422 validation response is parsed with detailed error list."""
        mock_client.login()
        csv_file = tmp_path / "invalid_server.csv"
        csv_file.write_text("name,sequence\nc1," + "A" * 40 + "INVALID_AA_TRIGGER\n")

        with pytest.raises(ValueError) as exc_info:
            mock_client.submit_candidates(
                target_id=TL1A_UUID,
                csv_path=str(csv_file),
                agreed=True,
                bypass_local_validation=True,  # Test server 422 parser
            )
        assert "含有非法氨基酸" in str(exc_info.value)

    def test_get_latest_submission_query(self, mock_client, tmp_path):
        """Query latest returns None initially, then record after submission."""
        mock_client.login()
        assert mock_client.get_latest_submission(TL1A_UUID) is None

        csv_file = tmp_path / "sub.csv"
        csv_file.write_text("name,sequence\n" + "\n".join([f"c_{i},{'A'*50}" for i in range(20)]))
        mock_client.submit_candidates(target_id=TL1A_UUID, csv_path=str(csv_file))

        latest = mock_client.get_latest_submission(TL1A_UUID)
        assert latest is not None
        assert latest["target_id"] == TL1A_UUID

    def test_reupload_note_enforcement(self, mock_client, tmp_path):
        """When an existing submission exists, empty note triggers error or auto-generation."""
        mock_client.login()
        csv_file = tmp_path / "reupload.csv"
        csv_file.write_text("name,sequence\n" + "\n".join([f"c_{i},{'A'*50}" for i in range(20)]))

        # Initial submission
        mock_client.submit_candidates(target_id=TL1A_UUID, csv_path=str(csv_file), note="Initial batch")
        assert mock_client.get_latest_submission(TL1A_UUID) is not None

        # Re-upload without note and auto_note_on_reupload=False raises error
        with pytest.raises(SJTUReUploadNoteRequiredError):
            mock_client.submit_candidates(
                target_id=TL1A_UUID,
                csv_path=str(csv_file),
                note="",
                auto_note_on_reupload=False,
            )

        # Re-upload without note and auto_note_on_reupload=True succeeds with auto-generated note
        res = mock_client.submit_candidates(
            target_id=TL1A_UUID,
            csv_path=str(csv_file),
            note="",
            auto_note_on_reupload=True,
        )
        assert res["status"] == "confirmed"

    def test_submit_candidates_abstract_length_validation(self, mock_client, tmp_path):
        """Abstract length boundary: <= 4000 characters accepted, > 4000 rejected universally."""
        csv_file = tmp_path / "valid.csv"
        csv_file.write_text("name,sequence\n" + "\n".join([f"c_{i},{'A'*50}" for i in range(20)]))

        # 4000 chars of arbitrary text passes validation
        exact_4000 = "Scientific abstract text. " * 153 + "End."
        exact_4000 = exact_4000 + "X" * (4000 - len(exact_4000))
        assert len(exact_4000) == 4000
        res = mock_client.submit_candidates(
            target_id=TL1A_UUID,
            csv_path=str(csv_file),
            abstract=exact_4000,
            dry_run=True,
        )
        assert res["count"] == 20

        # > 4000 chars rejected universally for any content (including 'A', 'B', and natural text)
        for over_text in ["A" * 4001, "B" * 4001, "Natural text abstract content... " * 125]:
            with pytest.raises(SJTUValidationError, match="Abstract exceeds maximum allowed length of 4000 characters"):
                mock_client.submit_candidates(
                    target_id=TL1A_UUID,
                    csv_path=str(csv_file),
                    abstract=over_text,
                    dry_run=True,
                )

        # Real condensed scientific abstract from data/scientific_abstract.txt passes
        real_abstract_path = PROJECT_ROOT / "data" / "scientific_abstract.txt"
        if real_abstract_path.is_file():
            real_abstract = real_abstract_path.read_text(encoding="utf-8").strip()
            assert len(real_abstract) <= 4000
            res = mock_client.submit_candidates(
                target_id=TL1A_UUID,
                csv_path=str(csv_file),
                abstract=real_abstract,
                dry_run=True,
            )
            assert res["count"] == 20

    def test_submit_candidates_target_404_handling(self, mock_client, tmp_path):
        """Submitting to a non-existent target ID triggers HTTP 404 and raises SJTUTargetNotFoundError."""
        mock_client.login()
        csv_file = tmp_path / "valid.csv"
        csv_file.write_text("name,sequence\n" + "\n".join([f"c_{i},{'A'*50}" for i in range(20)]))

        non_existent_id = "00000000-0000-0000-0000-000000000000"
        with pytest.raises(SJTUTargetNotFoundError) as exc_info:
            mock_client.submit_candidates(
                target_id=non_existent_id,
                csv_path=str(csv_file),
                agreed=True,
            )
        err = exc_info.value
        assert isinstance(err, SJTUValidationError)
        assert isinstance(err, KeyError)
        assert "Target does not exist on platform (HTTP 404)" in str(err)
        assert "赛题不存在" in str(err)


# ============================================================================
# Section 3: CLI Architecture Unit Tests
# ============================================================================

class TestSubmitSJTUCLI:
    """Unit tests for tools/submit_sjtu.py command-line interface."""

    def test_cli_dry_run_validation_passes(self, tmp_path, capsys):
        """--dry-run validates CSV and prints summary without HTTP POST."""
        csv_file = tmp_path / "valid.csv"
        csv_file.write_text("name,sequence\n" + "\n".join([f"c_{i},{'A'*50}" for i in range(20)]))

        exit_code = main(["--target", "tl1a", "--csv", str(csv_file), "--dry-run"])
        captured = capsys.readouterr()
        assert exit_code == 0
        assert "DRY RUN VALIDATION: SUCCESS" in captured.out
        assert "Total sequences: 20" in captured.out

    def test_cli_dry_run_validation_fails(self, tmp_path, capsys):
        """--dry-run with invalid CSV prints error table and returns code 1."""
        csv_file = tmp_path / "invalid_len.csv"
        csv_file.write_text("name,sequence\nc1," + "A" * 20 + "\n")  # 20 aa < 40 aa

        exit_code = main(["--target", "tl1a", "--csv", str(csv_file), "--dry-run"])
        captured = capsys.readouterr()
        assert exit_code == 1
        assert "DRY RUN VALIDATION: FAILED" in captured.out
        assert "序列长度不符合要求" in captured.out

    def test_cli_query_status_mode(self, monkeypatch, capsys):
        """--query-status queries platform status without requiring --csv."""
        mock_adapter = MockSJTUAdapter()
        mock_adapter.submissions[TL1A_UUID] = {
            "id": "batch-12345",
            "confirmed_at": "2026-11-12T00:15:00+08:00",
            "count": 20,
            "original_name": "TL1A_candidates.csv",
        }

        orig_init = requests.Session.__init__

        def patched_init(self, *args, **kwargs):
            orig_init(self, *args, **kwargs)
            self.mount("https://aiprotein.sjtu.edu.cn", mock_adapter)
            self.mount("http://aiprotein.sjtu.edu.cn", mock_adapter)

        monkeypatch.setattr(requests.Session, "__init__", patched_init)

        exit_code = main(["--target", "tl1a", "--query-status"])
        captured = capsys.readouterr()
        assert exit_code == 0
        assert "Latest Submission for Target TL1A" in captured.out
        assert "batch-12345" in captured.out

    def test_cli_abstract_from_file_or_string(self, tmp_path):
        """--abstract accepts both inline string and path to text file."""
        abs_file = tmp_path / "methodology.txt"
        abs_file.write_text("Hybrid Quantum-Classical VQE Optimization Method.")

        args_file = parse_args(["--target", "tl1a", "--csv", "test.csv", "--abstract", str(abs_file)])
        assert args_file.abstract == "Hybrid Quantum-Classical VQE Optimization Method."

        args_str = parse_args(["--target", "tl1a", "--csv", "test.csv", "--abstract", "Direct text summary"])
        assert args_str.abstract == "Direct text summary"


# ============================================================================
# Section 4: Audit Trail & Logging Tests
# ============================================================================

class TestSubmissionLogging:
    """Tests for logs/sjtu_submission.log generation and formatting."""

    def test_audit_log_created_and_appended(self, tmp_path):
        """Audit records are formatted and appended to specified log file."""
        log_file = tmp_path / "test_submission.log"
        logger = setup_submission_logger(log_file=str(log_file), verbose=False)

        logger.info("[Auth] Authenticated as Mai Duc Khang (Team: QAAA Lab)")
        logger.info("[Success] Confirmed Batch ID: batch-test-99999")

        assert log_file.exists()
        log_content = log_file.read_text(encoding="utf-8")
        assert "[Auth] Authenticated as Mai Duc Khang" in log_content
        assert "batch-test-99999" in log_content
        assert re.search(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \[INFO\]", log_content)


# ============================================================================
# Section 5: Opt-In Live Platform Integration Test
# ============================================================================

@pytest.mark.live
@pytest.mark.skipif(
    os.environ.get("TEST_SJTU_LIVE") != "1",
    reason="Live platform test skipped by default. Run with TEST_SJTU_LIVE=1 to enable.",
)
def test_live_platform_login_and_targets():
    """Live integration test querying SJTU platform."""
    client = SJTUPlatformClient()
    user_info = client.login()
    assert client.access_token is not None
    assert user_info["user"]["name"] == "Mai Duc Khang"
    assert user_info["user"]["team"]["name"] == "QAAA Lab"
    assert user_info["user"]["team"]["member_role"] == "captain"

    targets = client.get_targets()
    assert len(targets) == 7
    codes = {t["code"].lower(): t for t in targets}
    assert "tl1a" in codes
    assert "il7ra" in codes
