#!/usr/bin/env python3
"""
Unit tests for tools/ingest_matvenus_batch.py
"""

import os
import tempfile
import unittest

from tools.ingest_matvenus_batch import (
    compute_gravy,
    estimate_pi,
    check_developability_motifs,
    parse_raw_sequences,
    evaluate_candidate,
    ingest_batch
)


class TestIngestMatvenusBatch(unittest.TestCase):

    def test_compute_gravy(self):
        # Hydrophobic sequence (e.g. poly-leucine)
        poly_l = "LLLLLLLLLL"
        self.assertAlmostEqual(compute_gravy(poly_l), 3.8, places=2)
        # Hydrophilic sequence (e.g. poly-arginine)
        poly_r = "RRRRRRRRRR"
        self.assertAlmostEqual(compute_gravy(poly_r), -4.5, places=2)

    def test_estimate_pi(self):
        # Basic sequence
        basic_seq = "M" + "AK" * 40
        pi_basic = estimate_pi(basic_seq)
        self.assertGreater(pi_basic, 9.0)

        # Acidic sequence
        acidic_seq = "M" + "AD" * 40
        pi_acidic = estimate_pi(acidic_seq)
        self.assertLess(pi_acidic, 5.0)

    def test_check_developability_motifs(self):
        seq_clean = "M" + "AEKL" * 21
        self.assertEqual(len(check_developability_motifs(seq_clean)), 0)

        seq_with_motifs = "M" + "AEKL" * 10 + "NG" + "AEKL" * 5 + "DP"
        motifs = check_developability_motifs(seq_with_motifs)
        self.assertEqual(len(motifs), 2)
        self.assertTrue(any("NG" in m for m in motifs))
        self.assertTrue(any("DP" in m for m in motifs))

    def test_parse_raw_sequences(self):
        fasta_str = (
            ">IL7RA_A8_s1 len=92\n"
            "AALEELLKALEELLKALEELLKALEPWHETVMSINLFMIEFIIAHPSPKEIQITWELKDFLKKNFEPVSKSFHDRLFPFGDYKWRSPLEIGK\n"
            ">IL7RA_A8_s2 len=85\n"
            "AALEELLKALEELLKALEELLKALEKLAPKRTQHADWTNSQIIQERVHEQDYRSDINIIFNIQYWPGLFGEVAQANVENDTTYHY\n"
        )
        parsed = parse_raw_sequences(fasta_str)
        self.assertEqual(len(parsed), 2)
        self.assertEqual(parsed[0]["name"], "IL7RA_A8_s1")
        self.assertEqual(len(parsed[0]["sequence"]), 92)

    def test_evaluate_candidate(self):
        # Valid candidate
        cand_valid = {
            "name": "cand_good",
            "sequence": "M" + "AEKL" * 21  # 85 aa
        }
        res_valid = evaluate_candidate(cand_valid)
        self.assertTrue(res_valid["passed_qc"])
        self.assertEqual(res_valid["cys_count"], 0)

        # Candidate with Cys
        cand_cys = {
            "name": "cand_cys",
            "sequence": "M" + "AEKL" * 20 + "CAAA"
        }
        res_cys = evaluate_candidate(cand_cys)
        self.assertFalse(res_cys["passed_qc"])
        self.assertTrue(any("Contains 1 cysteine" in i for i in res_cys["hard_issues"]))

    def test_ingest_batch(self):
        fasta_input = (
            ">IL7RA_A8_s1\n"
            "AALEELLKALEELLKALEELLKALEPWHETVMSINLFMIEFIIAHPSPKEIQITWELKDFLKKNFEPVSKSFHDRLFPFGDYKWRSPLEIGK\n"
            ">IL7RA_A8_s2_bad_len\n"
            "AALEELLKALEELLKALE\n"  # Too short
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            out_csv = os.path.join(tmpdir, "ingested_manifest.csv")
            summary = ingest_batch(fasta_input, out_csv, target="IL7RA")
            self.assertEqual(summary["total_parsed"], 2)
            self.assertEqual(summary["total_passed_qc"], 1)
            self.assertEqual(summary["total_failed_qc"], 1)
            self.assertTrue(os.path.isfile(out_csv))


if __name__ == "__main__":
    unittest.main()
