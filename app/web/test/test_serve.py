"""Request tests for app/web/serve.py.

Starts the server in a thread on a free port and sends it the requests the
page sends. Standard library only. Set SIGMATAU_CLI to the compiled binary
(cli/build/bin/sigmatau); without it serve.py falls back to the built binary
or to `julia --project=cli`.

    SIGMATAU_CLI=cli/build/bin/sigmatau python3 -m unittest discover -s app/web/test -v
"""
import importlib.util
import json
import os
import threading
import time
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

SERVE = Path(__file__).resolve().parents[1] / "serve.py"


def load_serve():
    spec = importlib.util.spec_from_file_location("serve", SERVE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


serve = load_serve()


class RequestToArgsTests(unittest.TestCase):
    def test_file_goes_after_the_command(self):
        with __import__("tempfile").TemporaryDirectory() as d:
            args = serve.request_to_args(
                {"args": ["analyze", "--devs", "adev"], "file": {"name": "../../x.txt", "content": "1\n"}}, d)
            self.assertEqual(args[0], "analyze")
            self.assertEqual(Path(args[1]).parent, Path(d))     # no path escape
            self.assertEqual(Path(args[1]).name, "x.txt")
            self.assertEqual(args[2:], ["--devs", "adev"])

    def test_rejects(self):
        for req in ({"args": []}, {"args": "analyze"}, {"args": ["rm"]},
                    {"args": ["analyze", "--out", "/tmp/x"]}, {"args": ["analyze", "--out=/tmp/x"]},
                    {"args": ["analyze", 3]}):
            with self.assertRaises(ValueError, msg=req):
                serve.request_to_args(req, "/nonexistent")

    def test_live_rejects(self):
        for args in (None, "live", [], ["analyze"], ["live", "sim", "--out", "x"], ["live", 3]):
            with self.assertRaises(ValueError, msg=args):
                serve.Live(["true"], args)


@unittest.skipUnless(serve.cli_command(), "no sigmatau binary and no julia on PATH")
class ServeTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), serve.Handler)
        cls.base = f"http://127.0.0.1:{cls.server.server_address[1]}"
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()

    def request(self, path, body=None, raw=None):
        data = raw if raw is not None else (json.dumps(body).encode() if body is not None else None)
        req = urllib.request.Request(self.base + path, data=data,
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=600) as resp:
                return resp.status, resp.read()
        except urllib.error.HTTPError as e:
            return e.code, e.read()

    def post(self, body=None, raw=None, path="/api/analyze"):
        status, out = self.request(path, body=body, raw=raw)
        return status, json.loads(out)

    def test_health(self):
        status, out = self.request("/api/health")
        self.assertEqual(status, 200)
        doc = json.loads(out)
        self.assertTrue(doc["ok"])
        self.assertEqual(doc["core"], serve.cli_command()[0])

    def test_page_served(self):
        status, out = self.request("/")
        self.assertEqual(status, 200)
        self.assertIn(b"<html", out.lower())

    def test_synthetic_analyze(self):
        status, doc = self.post({"args": [
            "analyze", "--synth", "4096", "--seed", "5", "--devs", "adev,hdev",
            "--grid", "octave", "--confidence", "0.95", "--spectrum", "--series", "1000",
            "--format", "json"]})
        self.assertEqual(status, 200, doc)
        self.assertEqual(doc["schema_version"], 1)
        self.assertEqual([r["deviation"] for r in doc["results"]], ["adev", "hdev"])
        self.assertEqual(doc["analysis"]["confidence"], 0.95)
        self.assertIn("spectrum", doc)

    def test_file_check(self):
        content = "# t x\n" + "".join(f"{i} {1e-9 * ((i * 31) % 17):.3e}\n" for i in range(1024))
        status, doc = self.post({
            "args": ["check", "--time-col", "1", "--value-col", "2", "--detrend", "mean",
                     "--format", "json"],
            "file": {"name": "../../escape.txt", "content": content}})
        self.assertEqual(status, 200, doc)
        self.assertEqual(doc["record"]["n"], 1024)

    def test_dynamic(self):
        status, doc = self.post({"args": ["dynamic", "--synth", "4096", "--dev", "dhdev",
                                          "--format", "json"]})
        self.assertEqual(status, 200, doc)
        self.assertIn("dynamic", doc)

    def test_live(self):
        status, doc = self.post({"args": ["live", "sim", "--fast", "--count", "20", "--every", "10", "--m", "1"]},
                                path="/api/live/start")
        self.assertEqual(status, 200, doc)
        lines, since = [], 0
        for _ in range(200):
            status, out = self.request(f"/api/live/poll?id={doc['id']}&since={since}")
            poll = json.loads(out)
            lines += poll["lines"]
            since = poll["next"]
            if not poll["running"] and not poll["lines"]:
                break
            time.sleep(0.05)
        self.assertIsNone(poll["error"])
        self.assertEqual([json.loads(l)["n"] for l in lines], [10, 20])
        self.assertTrue(json.loads(lines[-1])["done"])
        status, _ = self.post({}, path="/api/live/stop")
        self.assertEqual(status, 200)
        status, _ = self.request(f"/api/live/poll?id={doc['id']}")
        self.assertEqual(status, 404)

    def test_live_error(self):
        status, doc = self.post({"args": ["live", "tcp://127.0.0.1:1"]}, path="/api/live/start")
        self.assertEqual(status, 200, doc)
        for _ in range(100):
            poll = json.loads(self.request(f"/api/live/poll?id={doc['id']}")[1])
            if not poll["running"]:
                break
            time.sleep(0.05)
        self.assertIn("could not connect", poll["error"])

    def test_core_error_is_422(self):
        status, doc = self.post({"args": ["analyze", "--synth", "4", "--format", "json"]})
        self.assertEqual(status, 422)
        self.assertIn("samples", doc["error"])

    def test_bad_requests(self):
        status, doc = self.post(raw=b"{not json")
        self.assertEqual(status, 400)
        self.assertIn("bad JSON", doc["error"])
        status, doc = self.post({"args": ["rm", "-rf"]})
        self.assertEqual(status, 400)
        status, doc = self.post({"args": ["analyze", "--out", "x"]})
        self.assertEqual(status, 400)
        status, doc = self.post({"args": ["analyze"]}, path="/api/nothing")
        self.assertEqual(status, 404)


if __name__ == "__main__":
    unittest.main()
