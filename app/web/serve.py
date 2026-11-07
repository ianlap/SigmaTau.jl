#!/usr/bin/env python3
"""Local server for the SigmaTau web explorer (development shell).

Serves app/web/ui/ and answers POST /api/analyze by running sigmatau-cli.
Standard library only.

    python3 app/web/serve.py               # http://127.0.0.1:9385
    python3 app/web/serve.py --port 8000 --open
    python3 app/web/serve.py --fetch-plotly   # save Plotly into ui/vendor for offline use

Which analysis core it runs, in order:
  1. $SIGMATAU_CLI, a path to a compiled sigmatau-cli
  2. app/web/cli/build/bin/sigmatau-cli(.exe), the juliac bundle
  3. julia --project=app/web/cli -m SigmaTauCLI (slow first call, no build)

The request and response formats are in app/web/README.md.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HERE = Path(__file__).resolve().parent
UI = HERE / "ui"
CLI_PROJECT = HERE / "cli"
MAX_BODY = 512 * 1024 * 1024
TIMEOUT_S = 600
PLOTLY_URL = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"


def cli_command():
    env = os.environ.get("SIGMATAU_CLI")
    if env:
        return [env]
    exe = "sigmatau-cli.exe" if os.name == "nt" else "sigmatau-cli"
    built = CLI_PROJECT / "build" / "bin" / exe
    if built.exists():
        return [str(built)]
    julia = shutil.which("julia")
    if julia:
        return [julia, f"--project={CLI_PROJECT}", "--startup-file=no", "-m", "SigmaTauCLI"]
    return None


def request_to_args(req, workdir):
    """Translate a JSON request into sigmatau-cli flags."""
    args = []
    src = req.get("source") or {}
    if src.get("type") == "synth":
        args += ["--synth", str(int(src.get("n", 8192))), "--seed", str(int(src.get("seed", 42)))]
    elif src.get("type") == "file":
        name = Path(str(src.get("name") or "record.txt")).name
        path = Path(workdir) / name
        path.write_text(str(src.get("content", "")), encoding="utf-8")
        args += ["--input", str(path)]
    else:
        raise ValueError("source.type must be 'file' or 'synth'")

    def opt(key, flag, conv=str):
        if req.get(key) is not None:
            args.extend([flag, str(conv(req[key]))])

    opt("kind", "--kind")
    opt("tau0", "--tau0", float)
    opt("time_col", "--time-col", int)
    opt("value_col", "--value-col", int)
    opt("reference_col", "--reference-col", int)
    opt("header", "--header", int)
    opt("scaling", "--scaling", float)
    opt("detrend", "--detrend")
    opt("grid", "--grid")
    opt("series", "--series", int)
    devs = req.get("devs")
    if devs is not None:
        args += ["--devs", ",".join(str(d) for d in devs) or "none"]
    if req.get("fillgaps"):
        args.append("--fillgaps")
    if req.get("remove_outliers"):
        args.append("--remove-outliers")
    if req.get("confidence") is None and "confidence" in req:
        args.append("--no-ci")
    else:
        opt("confidence", "--confidence", float)
    if req.get("bias_correction") is False:
        args.append("--no-bias-correction")
    return args


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=str(UI), **kw)

    def log_message(self, fmt, *a):
        if os.environ.get("SIGMATAU_WEB_LOG"):
            super().log_message(fmt, *a)

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def send_json(self, status, payload):
        body = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == "/api/health":
            cmd = cli_command()
            return self.send_json(200, {"ok": cmd is not None, "core": cmd[0] if cmd else None})
        return super().do_GET()

    def do_POST(self):
        if self.path != "/api/analyze":
            return self.send_json(404, {"error": "unknown endpoint"})
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0 or length > MAX_BODY:
            return self.send_json(413, {"error": "request body missing or too large"})
        try:
            req = json.loads(self.rfile.read(length))
        except json.JSONDecodeError as e:
            return self.send_json(400, {"error": f"bad JSON: {e}"})
        cmd = cli_command()
        if cmd is None:
            return self.send_json(503, {"error": "no analysis core: build sigmatau-cli or put julia on PATH"})
        with tempfile.TemporaryDirectory(prefix="sigmatau-web-") as work:
            try:
                args = request_to_args(req, work)
            except (ValueError, TypeError) as e:
                return self.send_json(400, {"error": str(e)})
            try:
                proc = subprocess.run(cmd + args, capture_output=True, timeout=TIMEOUT_S)
            except subprocess.TimeoutExpired:
                return self.send_json(504, {"error": "analysis timed out"})
        out = proc.stdout.strip()
        if out.startswith(b"{"):
            return self.send_json(200 if proc.returncode == 0 else 422, out)
        err = proc.stderr.decode(errors="replace").strip().splitlines()
        return self.send_json(500, {"error": err[-1] if err else f"core exited with {proc.returncode}"})


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--port", type=int, default=9385)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--open", action="store_true", help="open the page in the default browser")
    ap.add_argument("--fetch-plotly", action="store_true", help="download Plotly into ui/vendor and exit")
    ns = ap.parse_args()
    if ns.fetch_plotly:
        import urllib.request
        dest = UI / "vendor" / "plotly.min.js"
        dest.parent.mkdir(exist_ok=True)
        urllib.request.urlretrieve(PLOTLY_URL, dest)
        print(f"saved {dest} ({dest.stat().st_size // 1024} KiB)")
        return 0
    server = ThreadingHTTPServer((ns.host, ns.port), Handler)
    url = f"http://{ns.host}:{ns.port}/"
    core = cli_command()
    print(f"SigmaTau web explorer on {url}  (core: {' '.join(core) if core else 'none found'})", flush=True)
    if ns.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    sys.exit(main())
