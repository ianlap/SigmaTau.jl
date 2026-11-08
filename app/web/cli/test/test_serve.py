"""Request tests for app/web/serve.py against the compiled core.

Starts the server in a thread on a free port and sends it the requests the
page sends. Standard library only; see test_core.py for how the core is found.
"""
import importlib.util
import json
import os
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

from test_core import CLI_PROJECT, core_path

SERVE = CLI_PROJECT.parent / "serve.py"


def load_serve():
    spec = importlib.util.spec_from_file_location("serve", SERVE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@unittest.skipUnless(os.path.exists(core_path()), f"no core at {core_path()}")
class ServeTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._env = os.environ.get("SIGMATAU_CLI")
        os.environ["SIGMATAU_CLI"] = core_path()
        serve = load_serve()
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), serve.Handler)
        cls.base = f"http://127.0.0.1:{cls.server.server_address[1]}"
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        if cls._env is None:
            os.environ.pop("SIGMATAU_CLI", None)
        else:
            os.environ["SIGMATAU_CLI"] = cls._env

    def request(self, path, body=None, raw=None):
        data = raw if raw is not None else (json.dumps(body).encode() if body is not None else None)
        req = urllib.request.Request(self.base + path, data=data,
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
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
        self.assertEqual(doc["core"], core_path())

    def test_page_served(self):
        status, out = self.request("/")
        self.assertEqual(status, 200)
        self.assertIn(b"<html", out.lower())

    def test_synthetic_request(self):
        status, doc = self.post({
            "source": {"type": "synth", "n": 4096, "seed": 5},
            "devs": ["adev", "hdev"], "grid": "octave", "confidence": 0.95,
            "bias_correction": True, "series": 1000,
            "spectrum": True, "dynamic": "dadev",
        })
        self.assertEqual(status, 200, doc)
        self.assertEqual([r["dev"] for r in doc["results"]], ["adev", "hdev"])
        self.assertEqual(doc["options"]["confidence"], 0.95)
        self.assertIn("spectrum", doc)
        self.assertEqual(doc["dynamic"]["dev"], "dadev")

    def test_file_request(self):
        content = "# t x\n" + "".join(f"{i} {1e-9 * ((i * 31) % 17):.3e}\n" for i in range(1024))
        status, doc = self.post({
            "source": {"type": "file", "name": "../../escape.txt", "content": content},
            "kind": "phase", "time_col": 1, "value_col": 2, "scaling": 1.0,
            "detrend": "mean", "devs": ["mdev"], "confidence": None,
        })
        self.assertEqual(status, 200, doc)
        self.assertEqual(doc["record"]["n"], 1024)
        self.assertFalse(doc["options"]["ci"])
        self.assertEqual(doc["results"][0]["ci_lo"], [])

    def test_core_error_is_422(self):
        status, doc = self.post({"source": {"type": "synth", "n": 4}})
        self.assertEqual(status, 422)
        self.assertIn("samples", doc["error"])

    def test_bad_requests(self):
        status, doc = self.post(raw=b"{not json")
        self.assertEqual(status, 400)
        self.assertIn("bad JSON", doc["error"])
        status, doc = self.post({"source": {"type": "url"}})
        self.assertEqual(status, 400)
        status, doc = self.post({"source": {"type": "synth"}}, path="/api/nothing")
        self.assertEqual(status, 404)


if __name__ == "__main__":
    unittest.main()
