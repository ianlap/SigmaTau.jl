"""Checks for the web explorer's compiled analysis core (sigmatau-cli).

Runs the core as a subprocess and checks the JSON it prints against the
contract in app/web/README.md, the Stable32 reference values in
test/fixtures/validation, and (with SIGMATAU_PARITY=1) the same command run
through plain Julia. Standard library only.

    SIGMATAU_CLI=app/web/cli/build/bin/sigmatau-cli python3 -m unittest discover -s app/web/cli/test -v

Without SIGMATAU_CLI the tests look for app/web/cli/build/bin/sigmatau-cli.
"""
import csv
import json
import math
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
CLI_PROJECT = HERE.parent
REPO = CLI_PROJECT.parent.parent.parent
FIXTURES = REPO / "test" / "fixtures" / "validation"
ALL_DEVS = ["adev", "mdev", "tdev", "hdev", "mhdev", "htdev",
            "totdev", "mtotdev", "ttotdev", "htotdev", "mhtotdev",
            "theo1", "theoh", "pdev", "mtie", "tierms"]
SIGMA_X_DEVS = {"tdev", "htdev", "ttotdev", "mtie", "tierms"}


def core_path():
    env = os.environ.get("SIGMATAU_CLI")
    if env:
        return env
    exe = "sigmatau-cli.exe" if os.name == "nt" else "sigmatau-cli"
    return str(CLI_PROJECT / "build" / "bin" / exe)


def run(*args, cmd=None, timeout=300):
    """Run the core and return (exit code, parsed JSON or None, raw stdout)."""
    proc = subprocess.run((cmd or [core_path()]) + [str(a) for a in args],
                          capture_output=True, timeout=timeout)
    out = proc.stdout.decode()
    try:
        doc = json.loads(out)
    except json.JSONDecodeError:
        doc = None
    return proc.returncode, doc, out


def ok(*args):
    code, doc, out = run(*args)
    if code != 0 or doc is None:
        raise AssertionError(f"core failed (exit {code}) for {args}: {out[:500]}")
    return doc


def finite(v):
    return v is not None and math.isfinite(v)


@unittest.skipUnless(os.path.exists(core_path()), f"no core at {core_path()}")
class CoreTests(unittest.TestCase):

    def check_result(self, r):
        n = len(r["tau"])
        self.assertGreater(n, 0, r["dev"])
        self.assertEqual(len(r["sigma"]), n)
        self.assertEqual(r["axis"], "sigma_x" if r["dev"] in SIGMA_X_DEVS else "sigma_y")
        self.assertEqual(r["tau"], sorted(r["tau"]), r["dev"])
        self.assertTrue(all(t > 0 for t in r["tau"]), r["dev"])
        self.assertTrue(any(finite(s) and s > 0 for s in r["sigma"]), r["dev"])

    # ── Record and deviations ────────────────────────────────────────────

    def test_synthetic_record_schema(self):
        doc = ok("--synth", 4096, "--seed", 7)
        for key in ("sigmatau", "record", "options", "results", "elapsed_ms"):
            self.assertIn(key, doc)
        self.assertNotIn("spectrum", doc)
        self.assertNotIn("dynamic", doc)
        rec = doc["record"]
        self.assertEqual(rec["n"], 4096)
        self.assertEqual(rec["tau0"], 1.0)
        self.assertEqual(rec["kind"], "phase")
        self.assertEqual(rec["source"], "synthetic")
        self.assertEqual(len(rec["series"]["t"]), len(rec["series"]["x"]))
        chk = rec["check"]
        self.assertEqual(chk["samples"], 4096)
        self.assertEqual(chk["gaps"], 0)
        self.assertEqual(len(chk["outlier_t"]), min(chk["outliers"], 500))
        self.assertEqual(doc["options"]["grid"], "octave")
        self.assertTrue(doc["options"]["ci"])
        self.assertEqual([r["dev"] for r in doc["results"]], ["adev", "mdev", "hdev", "mhdev"])

    def test_same_seed_same_output(self):
        a = ok("--synth", 2048, "--seed", 3, "--devs", "adev")["results"]
        b = ok("--synth", 2048, "--seed", 3, "--devs", "adev")["results"]
        c = ok("--synth", 2048, "--seed", 4, "--devs", "adev")["results"]
        self.assertEqual(a, b)
        self.assertNotEqual(a[0]["sigma"], c[0]["sigma"])

    def test_every_estimator(self):
        doc = ok("--synth", 2048, "--devs", ",".join(ALL_DEVS))
        self.assertEqual([r["dev"] for r in doc["results"]], ALL_DEVS)
        for r in doc["results"]:
            with self.subTest(dev=r["dev"]):
                self.check_result(r)

    def test_confidence_intervals(self):
        doc = ok("--synth", 4096, "--devs", "adev,mdev,totdev")
        for r in doc["results"]:
            n = len(r["tau"])
            for key in ("ci_lo", "ci_hi", "edf", "noise"):
                self.assertEqual(len(r[key]), n, (r["dev"], key))
            for lo, s, hi in zip(r["ci_lo"], r["sigma"], r["ci_hi"]):
                if finite(lo) and finite(hi):
                    self.assertLessEqual(lo, s, r["dev"])
                    self.assertLessEqual(s, hi, r["dev"])
            self.assertTrue(all(isinstance(k, str) and k for k in r["noise"]), r["dev"])

    def test_no_ci(self):
        doc = ok("--synth", 2048, "--devs", "adev", "--no-ci")
        self.assertFalse(doc["options"]["ci"])
        r = doc["results"][0]
        for key in ("ci_lo", "ci_hi", "edf", "noise"):
            self.assertEqual(r[key], [], key)

    def test_grids(self):
        counts = {}
        for grid in ("decade", "octave", "half-octave", "all"):
            r = ok("--synth", 1024, "--devs", "adev", "--grid", grid)["results"][0]
            counts[grid] = len(r["tau"])
        self.assertLess(counts["decade"], counts["octave"])
        self.assertLess(counts["octave"], counts["half-octave"])
        self.assertLess(counts["half-octave"], counts["all"])

    def test_devs_none_reads_record_only(self):
        doc = ok("--synth", 64, "--devs", "none")
        self.assertEqual(doc["results"], [])
        self.assertEqual(doc["record"]["n"], 64)

    def test_series_cap(self):
        doc = ok("--synth", 8192, "--devs", "none", "--series", 500)
        self.assertLessEqual(len(doc["record"]["series"]["x"]), 500)

    # ── Files ───────────────────────────────────────────────────────────

    def test_stable32_file_matches_stable32(self):
        """The Stable32 .DAT fixture, against Stable32's own numbers."""
        dat = FIXTURES / "stable32gen.DAT"
        ref = {}
        with open(FIXTURES / "stable32out" / "stable32_data_full.csv", newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                ref[(row["Type"], float(row["Tau"]))] = float(row["Sigma"])
        names = {"adev": "Overlapping Allan", "mdev": "Modified Allan",
                 "hdev": "Overlapping Hadamard", "tdev": "Time"}
        doc = ok("--input", dat, "--time-col", 0, "--value-col", 1,
                 "--devs", ",".join(names), "--grid", "octave", "--no-ci")
        self.assertEqual(doc["record"]["n"], 8192)
        self.assertEqual(doc["record"]["tau0"], 1.0)
        for r in doc["results"]:
            checked = 0
            for tau, sigma in zip(r["tau"], r["sigma"]):
                want = ref.get((names[r["dev"]], tau))
                if want is None or not finite(sigma):
                    continue
                with self.subTest(dev=r["dev"], tau=tau):
                    # Stable32 writes five significant figures.
                    self.assertAlmostEqual(sigma / want, 1.0, delta=2e-4)
                checked += 1
            self.assertGreaterEqual(checked, 8, r["dev"])

    def write(self, name, text):
        path = Path(self.tmp.name) / name
        path.write_text(text, encoding="utf-8")
        return str(path)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="sigmatau-core-test-")

    def tearDown(self):
        self.tmp.cleanup()

    def phase_text(self, n, gap_at=None, nan_at=None, comment=True):
        """A phase record with a time column. gap_at: ten time stamps go
        missing from there on; nan_at: five values from there are NaN."""
        lines = ["# synthetic phase record, ns"] if comment else []
        for i in range(n):
            t = float(i)
            if gap_at is not None and i >= gap_at:
                t += 10.0
            x = 1e-3 * math.sin(0.01 * i) + 1e-4 * ((i * 7919) % 101 - 50) / 50
            v = "NaN" if nan_at is not None and nan_at <= i < nan_at + 5 else f"{x:.9e}"
            lines.append(f"{t:.1f} {v}")
        return "\n".join(lines) + "\n"

    def test_text_file_with_scaling_and_detrend(self):
        path = self.write("rec.txt", self.phase_text(1024))
        doc = ok("--input", path, "--time-col", 1, "--value-col", 2, "--scaling", 1e-9,
                 "--detrend", "linear", "--devs", "adev")
        self.assertEqual(doc["record"]["n"], 1024)
        self.check_result(doc["results"][0])
        self.assertLess(max(abs(x) for x in doc["record"]["series"]["x"] if x is not None), 1e-11)

    def test_frequency_record(self):
        lines = "\n".join(f"{1e-11 * ((i * 104729) % 97 - 48) / 48:.9e}" for i in range(2048))
        path = self.write("freq.txt", lines + "\n")
        doc = ok("--input", path, "--kind", "frequency", "--time-col", 0, "--value-col", 1,
                 "--tau0", 1, "--devs", "adev")
        self.assertEqual(doc["record"]["kind"], "frequency")
        self.check_result(doc["results"][0])

    def test_reference_column(self):
        rows = [f"{i} {1e-9 * i:.6e} {2e-9 * i:.6e}" for i in range(512)]
        path = self.write("two.txt", "\n".join(rows) + "\n")
        doc = ok("--input", path, "--time-col", 1, "--value-col", 3, "--reference-col", 2,
                 "--devs", "none", "--series", 0)
        self.assertEqual(doc["record"]["n"], 512)

    def test_fillgaps(self):
        """Gap filling goes through FFTW, the trim-sensitive path."""
        path = self.write("gappy.txt", self.phase_text(1024, gap_at=600, nan_at=200))
        raw = ok("--input", path, "--devs", "none")
        self.assertEqual(raw["record"]["n"], 1024)
        self.assertEqual(raw["record"]["check"]["gaps"], 5)
        filled = ok("--input", path, "--fillgaps", "--devs", "adev")
        self.assertTrue(filled["options"]["fillgaps"])
        self.assertEqual(filled["record"]["n"], 1034)
        self.assertEqual(filled["record"]["check"]["gaps"], 0)
        self.check_result(filled["results"][0])

    def test_remove_outliers(self):
        text = self.phase_text(2048, comment=False).splitlines()
        t, x = text[1000].split()
        text[1000] = f"{t} {float(x) + 5.0:.9e}"
        path = self.write("spike.txt", "\n".join(text) + "\n")
        before = ok("--input", path, "--devs", "none")
        self.assertGreaterEqual(before["record"]["check"]["outliers"], 1)
        after = ok("--input", path, "--remove-outliers", "--devs", "adev")
        self.assertTrue(after["options"]["remove_outliers"])
        self.check_result(after["results"][0])

    # ── Spectrum and dynamic ────────────────────────────────────────────

    def test_spectrum(self):
        doc = ok("--synth", 8192, "--devs", "none", "--spectrum")
        sp = doc["spectrum"]
        self.assertEqual(sp["nperseg"], 2048)
        n = len(sp["f"])
        self.assertGreater(n, 100)
        self.assertEqual(len(sp["sy"]), n)
        self.assertEqual(len(sp["sx"]), n)
        self.assertGreater(sp["f"][0], 0)
        self.assertEqual(sp["f"], sorted(sp["f"]))
        self.assertTrue(all(finite(v) and v > 0 for v in sp["sy"] + sp["sx"]))

    def test_spectrum_nperseg(self):
        sp = ok("--synth", 8192, "--devs", "none", "--spectrum", "--nperseg", 512)["spectrum"]
        self.assertEqual(sp["nperseg"], 512)
        self.assertEqual(len(sp["f"]), 256)

    def test_dynamic(self):
        for dev in ("dadev", "dhdev"):
            with self.subTest(dev=dev):
                dyn = ok("--synth", 4096, "--devs", "none", "--dynamic", dev,
                         "--window", 512, "--step", 256)["dynamic"]
                self.assertEqual(dyn["dev"], dev)
                self.assertEqual(dyn["window"], 512)
                self.assertEqual(dyn["step"], 256)
                self.assertGreater(len(dyn["t"]), 5)
                self.assertEqual(len(dyn["map"]), len(dyn["t"]))
                for row in dyn["map"]:
                    self.assertEqual(len(row), len(dyn["tau"]))
                self.assertTrue(any(finite(v) and v > 0 for row in dyn["map"] for v in row))

    def test_dynamic_aliases(self):
        for alias, dev in (("adev", "dadev"), ("hdev", "dhdev")):
            dyn = ok("--synth", 2048, "--devs", "none", "--dynamic", alias)["dynamic"]
            self.assertEqual(dyn["dev"], dev)

    def test_dynamic_defaults(self):
        dyn = ok("--synth", 4096, "--devs", "none", "--dynamic", "dadev")["dynamic"]
        self.assertEqual(dyn["window"], 512)
        self.assertEqual(dyn["step"], 256)

    # ── Errors ──────────────────────────────────────────────────────────

    def assert_error(self, *args, contains=None):
        code, doc, out = run(*args)
        self.assertEqual(code, 2, out)
        self.assertIsNotNone(doc, out)
        self.assertEqual(list(doc), ["error"])
        self.assertTrue(doc["error"])
        if contains:
            self.assertIn(contains, doc["error"])

    def test_errors(self):
        self.assert_error("--synth", 100, "--bogus", contains="unknown option --bogus")
        self.assert_error("--synth", "many", contains="expects an integer")
        self.assert_error("--synth", 100, "--tau0", "x", contains="expects a number")
        self.assert_error("--synth", contains="needs a value")
        self.assert_error("--synth", 100, "--devs", "nodev", contains="unknown deviation nodev")
        self.assert_error("--synth", 4, contains="only 4 samples")
        self.assert_error("--synth", 100, "--dynamic", "nodev", contains="dynamic")
        self.assert_error("--input", Path(self.tmp.name) / "missing.txt")

    def test_spectrum_rejects_gaps(self):
        path = self.write("nan.txt", self.phase_text(1024, nan_at=600))
        self.assert_error("--input", path, "--devs", "none", "--spectrum", contains="without gaps")

    def test_help(self):
        code, _, out = run("--help")
        self.assertEqual(code, 0)
        self.assertIn("usage", out)
        code, _, _ = run()
        self.assertEqual(code, 1)


@unittest.skipUnless(os.environ.get("SIGMATAU_PARITY") == "1" and os.path.exists(core_path())
                     and shutil.which("julia"), "set SIGMATAU_PARITY=1 with julia on PATH")
class ParityTests(unittest.TestCase):
    """The trimmed binary prints the same numbers as plain Julia."""

    def test_same_numbers_as_julia(self):
        args = ["--synth", 4096, "--seed", 11, "--devs", "adev,mdev,hdev,totdev,mtie",
                "--spectrum", "--dynamic", "dhdev"]
        julia = [shutil.which("julia"), f"--project={CLI_PROJECT}", "--startup-file=no",
                 "-m", "SigmaTauCLI"]
        _, compiled, _ = run(*args)
        code, reference, out = run(*args, cmd=julia, timeout=1800)
        self.assertEqual(code, 0, out[:500])
        for doc in (compiled, reference):
            doc.pop("elapsed_ms")
        self.assertEqual(compiled, reference)


if __name__ == "__main__":
    unittest.main()
