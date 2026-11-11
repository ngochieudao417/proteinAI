#!/usr/bin/env python3
"""
Unit tests for tools/candidate_generator.py (pure unittest, standard library).
"""

import os
import json
import tempfile
import unittest

from tools.candidate_generator import (
    clean_pdb,
    generate_bindcraft_config,
    validate_candidate_sequence,
    create_screening_manifest,
    TARGET_SPECS,
    MINIBINDER_MIN_LEN,
    MINIBINDER_MAX_LEN
)


class TestCandidateGenerator(unittest.TestCase):

    def test_target_specs(self):
        self.assertIn("tl1a", TARGET_SPECS)
        self.assertIn("il7ra", TARGET_SPECS)
        
        tl1a = TARGET_SPECS["tl1a"]
        self.assertEqual(tl1a["length"], 180)
        self.assertIn("D104", tl1a["hotspots_group_a"])
        self.assertIn("W137", tl1a["hotspots_group_a"])
        self.assertIn("N62", tl1a["avoid_residues"])

        il7ra = TARGET_SPECS["il7ra"]
        self.assertEqual(il7ra["length"], 219)
        self.assertIn("K77", il7ra["hotspots_group_a"])
        self.assertIn("F79", il7ra["hotspots_group_a"])
        self.assertIn("N29", il7ra["avoid_residues"])

    def test_validate_candidate_sequence(self):
        # Valid 85 aa sequence
        valid_85 = "M" + "AEKL" * 21  # 1 + 84 = 85 aa
        ok, errs = validate_candidate_sequence(valid_85)
        self.assertTrue(ok)
        self.assertEqual(len(errs), 0)

        # Valid 100 aa sequence
        valid_100 = "M" + "AEKL" * 24 + "AAA"  # 1 + 96 + 3 = 100 aa
        ok, errs = validate_candidate_sequence(valid_100)
        self.assertTrue(ok)
        self.assertEqual(len(errs), 0)

        # Too short (< 85 aa)
        short_70 = "M" + "AEKL" * 17 + "A"  # 70 aa
        ok, errs = validate_candidate_sequence(short_70)
        self.assertFalse(ok)
        self.assertTrue(any("violates [85, 100]" in e for e in errs))

        # Too long (> 100 aa)
        long_105 = "M" + "AEKL" * 26  # 105 aa
        ok, errs = validate_candidate_sequence(long_105)
        self.assertFalse(ok)
        self.assertTrue(any("violates [85, 100]" in e for e in errs))

        # Invalid characters (e.g. B, Z, X)
        invalid_char = "M" + "AEKX" * 21
        ok, errs = validate_candidate_sequence(invalid_char)
        self.assertFalse(ok)
        self.assertTrue(any("invalid characters" in e for e in errs))

        # Free Cys (odd number of Cysteines)
        free_cys = "M" + "AEKL" * 20 + "CAAA"  # 1 Cys
        ok, errs = validate_candidate_sequence(free_cys)
        self.assertFalse(ok)
        self.assertTrue(any("free unpaired Cys" in e for e in errs))

    def test_generate_bindcraft_config(self):
        config = generate_bindcraft_config(
            target_key="tl1a",
            target_pdb="data/structures/2RE9_clean.pdb",
            target_chains=["A", "B", "C"],
            epitope_group="A",
            num_designs=50
        )
        self.assertEqual(config["binder_specs"]["min_length"], 85)
        self.assertEqual(config["binder_specs"]["max_length"], 100)
        self.assertIn("D104", config["target"]["target_hotspots"])
        self.assertEqual(config["backbone_generation"]["num_designs"], 50)

        config_il7ra = generate_bindcraft_config(
            target_key="il7ra",
            target_pdb="data/structures/7OPB_clean.pdb",
            target_chains=["A"],
            epitope_group="A"
        )
        self.assertEqual(config_il7ra["binder_specs"]["min_length"], 85)
        self.assertIn("K77", config_il7ra["target"]["target_hotspots"])

    def test_clean_pdb_with_sample(self):
        sample_pdb = (
            "HEADER    SAMPLE TEST PDB\n"
            "ATOM      1  N   LEU A  17      -9.174  53.580  -5.556  1.00103.51           N  \n"
            "ATOM      2  CA  LEU A  17     -10.560  54.097  -5.755  1.00103.74           C  \n"
            "ATOM      3  N   GLY A   0      -5.000  50.000  -5.000  1.00 50.00           N  \n"
            "ATOM      4  N   LEU B  17     -19.174  63.580 -15.556  1.00103.51           N  \n"
            "HETATM   10  O   HOH A 101       0.000   0.000   0.000  1.00 20.00           O  \n"
            "END\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "in.pdb")
            out_path = os.path.join(tmpdir, "out.pdb")
            with open(in_path, "w") as f:
                f.write(sample_pdb)

            res = clean_pdb(
                input_pdb_path=in_path,
                output_pdb_path=out_path,
                target_chains=["A"],
                remove_hetatm=True
            )
            self.assertEqual(res["total_atoms"], 2)  # LEU A 17 (N, CA); GLY A 0 stripped, HOH stripped, B stripped
            self.assertTrue(os.path.isfile(out_path))
            with open(out_path, "r") as f:
                content = f.read()
                self.assertNotIn("GLY A   0", content)
                self.assertNotIn("HOH", content)
                self.assertNotIn("LEU B", content)
                self.assertIn("LEU A  17", content)

    def test_create_screening_manifest(self):
        candidates = [
            {
                "name": "cand_tl1a_01",
                "sequence": "M" + "AEKL" * 21,  # 85 aa
                "target": "TL1A",
                "epitope_group": "A",
                "generator": "bindcraft",
                "backbone_id": "bb_001"
            },
            {
                "name": "cand_il7ra_01",
                "sequence": "M" + "AEKL" * 23,  # 93 aa
                "target": "IL7RA",
                "epitope_group": "A",
                "generator": "bindcraft",
                "backbone_id": "bb_002"
            },
            {
                "name": "cand_short_bad",
                "sequence": "M" + "AEKL" * 10,  # 41 aa, invalid
                "target": "TL1A"
            }
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            out_csv = os.path.join(tmpdir, "manifest.csv")
            result = create_screening_manifest(candidates, out_csv)
            self.assertEqual(result["total_input"], 3)
            self.assertEqual(result["total_valid"], 2)
            self.assertEqual(result["total_rejected"], 1)
            self.assertTrue(os.path.isfile(out_csv))

            with open(out_csv, "r") as f:
                lines = [line.strip() for line in f if line.strip()]
                self.assertEqual(len(lines), 3)  # Header + 2 valid
                self.assertEqual(lines[0], "name,sequence,target,length,epitope_group,generator,backbone_id")
                self.assertIn("cand_tl1a_01", lines[1])
                self.assertIn("cand_il7ra_01", lines[2])


if __name__ == "__main__":
    unittest.main()
