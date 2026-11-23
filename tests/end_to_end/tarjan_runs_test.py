#!/usr/bin/env python3
"""Check Tarjan run budgets against a built MtKaHyPar executable.

Run: python3 tests/end_to_end/tarjan_runs_test.py [path/to/MtKaHyPar]
"""

import argparse
import collections
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "release/mt-kahypar/application/MtKaHyPar"


class TarjanRunsTest(unittest.TestCase):
    def run_ip(self, *options, runs=3, algorithms=("random", "tarjan"), config=None):
        preset = ["--config", str(config)] if config else ["--preset-type=default"]
        command = [
            str(BINARY), "-h", str(ROOT / "tests/instances/karate_club.graph.hgr"),
            *preset, "-t", "1", "-k", "2", "-e", "0.03", "-o", "km1",
            "--i-enabled-ip-algos", *algorithms, "--ip-enabled-algos-tarjan", "bfs",
            "--i-runs", str(runs), "--i-use-adaptive-ip-runs=false",
            "--i-fm-refinement-rounds=0", "--i-perform-refinement-on-best-partitions=false",
            *options,
        ]
        return subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, timeout=30).stdout

    def assert_runs(self, output, random_runs, tarjan_runs):
        self.assertIn("RESULT algorithm=", output, output[-2000:])
        counts = collections.Counter(re.findall(r"Algorithm:\s+(\w+)\s+Quality:", output))
        self.assertEqual(counts["random"], random_runs)
        self.assertEqual(counts["tarjan"], tarjan_runs)
        # Each Tarjan attempt runs its separately selected inner BFS once.
        self.assertEqual(counts["bfs"], tarjan_runs)
        self.assertIn(f"initial_partitioning_tarjan_runs={tarjan_runs}", output)

    def test_inherits_common_budget(self):
        self.assert_runs(self.run_ip(), 3, 3)

    def test_independent_smaller_larger_and_zero_budgets(self):
        for runs in (1, 7, 0):
            with self.subTest(tarjan_runs=runs):
                self.assert_runs(self.run_ip(f"--i-tarjan-runs={runs}"), 3, runs)

    def test_tarjan_can_run_when_common_budget_is_zero(self):
        self.assert_runs(self.run_ip("--i-tarjan-runs=2", runs=0), 0, 2)

    def test_ini_override(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / "preset.ini"
            config.write_text((ROOT / "config/default_preset.ini").read_text()
                              + "\ni-tarjan-runs=2\n")
            self.assert_runs(self.run_ip(config=config), 3, 2)

    def test_rejects_negative_budget(self):
        output = self.run_ip("--i-tarjan-runs=-1")
        self.assertNotIn("RESULT algorithm=", output)
        self.assertIn("--i-tarjan-runs", output)

    def test_rejects_empty_task_pool(self):
        output = self.run_ip("--i-tarjan-runs=0", algorithms=("tarjan",))
        self.assertNotIn("RESULT algorithm=", output)
        self.assertIn("No initial partitioning runs are enabled", output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", nargs="?", type=Path, default=BINARY)
    args = parser.parse_args()
    BINARY = args.binary.resolve()
    unittest.main(argv=[__file__], verbosity=2)
