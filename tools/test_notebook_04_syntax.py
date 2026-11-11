#!/usr/bin/env python3
"""
tools/test_notebook_04_syntax.py — Verify notebook JSON validity and cell Python syntax.
"""

import json
import ast
import os
import unittest


class TestNotebook04Syntax(unittest.TestCase):

    def test_notebook_structure_and_syntax(self):
        nb_path = "notebooks/04_candidate_generation_bindcraft.ipynb"
        self.assertTrue(os.path.isfile(nb_path), f"File not found: {nb_path}")

        with open(nb_path, "r", encoding="utf-8") as f:
            nb = json.load(f)

        self.assertIn("cells", nb)
        self.assertGreaterEqual(len(nb["cells"]), 7)

        # Check all code cells parse as valid python
        for idx, cell in enumerate(nb["cells"]):
            cell_type = cell.get("cell_type")
            source = "".join(cell.get("source", []))
            if cell_type == "code":
                # Strip out IPython magic or title lines (e.g. #@title)
                clean_lines = []
                for line in source.splitlines():
                    if line.strip().startswith("%") or line.strip().startswith("!"):
                        continue
                    clean_lines.append(line)
                code_to_parse = "\n".join(clean_lines)

                try:
                    ast.parse(code_to_parse)
                except SyntaxError as e:
                    self.fail(f"SyntaxError in cell {idx}: {e}")

        # Check metadata
        self.assertEqual(nb.get("metadata", {}).get("accelerator"), "GPU")
        self.assertEqual(nb.get("metadata", {}).get("colab", {}).get("gpuType"), "L4")


if __name__ == "__main__":
    unittest.main()
