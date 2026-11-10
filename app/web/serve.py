#!/usr/bin/env python3
"""Local server for the SigmaTau web explorer (development shell).

Serves app/web/ui/ and answers POST /api/analyze by running the sigmatau binary,
and runs `sigmatau live` for the Live page (/api/live/*). Standard library only.

    python3 app/web/serve.py               # http://127.0.0.1:9385
    python3 app/web/serve.py --port 8000 --open
    python3 app/web/serve.py --fetch-plotly   # save Plotly into ui/vendor for offline use

Which analysis core it runs, in order:
  1. $SIGMATAU_CLI, a path to a compiled sigmatau
  2. cli/build/bin/sigmatau(.exe), the juliac bundle
  3. julia --project=cli -m SigmaTauCLI (slow first call, no build)

The request and response formats are in app/web/README.md.
"""
import argparse
import json
import os
import select
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

HERE = Path(__file__).resolve().parent
UI = HERE / "ui"
CLI_PROJECT = HERE.parent.parent / "cli"
MAX_BODY = 512 * 1024 * 1024
TIMEOUT_S = 600
PLOTLY_URL = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"


def cli_command():
    env = os.environ.get("SIGMATAU_CLI")
    if env:
        return [env]
    exe = "sigmatau.exe" if os.name == "nt" else "sigmatau"
    built = CLI_PROJECT / "build" / "bin" / exe
    if built.exists():
        return [str(built)]
    julia = shutil.which("julia")
    if julia:
        return [julia, f"--project={CLI_PROJECT}", "--startup-file=no", "-m", "SigmaTauCLI"]
    return None


COMMANDS = {"check", "analyze", "spectrum", "dynamic", "drift", "timescale", "model"}


def request_to_args(req, workdir):
    """Check a request {args: [command, flags...], file: {name, content} | null}
    and return the sigmatau arguments, with the file saved under workdir."""
    args = req.get("args")
    if not isinstance(args, list) or not args or not all(isinstance(a, str) for a in args):
        raise ValueError("args must be a non-empty list of strings")
    if args[0] not in COMMANDS:
        raise ValueError(f"unknown command {args[0]!r}")
    if any(a == "--out" or a.startswith("--out=") for a in args):
        raise ValueError("--out is not allowed here")
    f = req.get("file")
    if f is None:
        # The record is the uploaded file or a synthetic one, never a path on this machine.
        if not any(a == "--synth" or a.startswith("--synth=") for a in args):
            raise ValueError("a request without a file must use --synth")
        return list(args)
    name = Path(str(f.get("name") or "record.txt")).name or "record.txt"
    path = Path(workdir) / name
    path.write_text(str(f.get("content", "")), encoding="utf-8")
    return [args[0], str(path), *args[1:]]


def hung_up(sock):
    """True once the page has closed its end of the connection."""
    try:
        if not select.select([sock], [], [], 0)[0]:
            return False
        return sock.recv(1, socket.MSG_PEEK) == b""
    except BlockingIOError:
        return False
    except OSError:
        return True


def run_core(argv, env, sock=None, timeout=TIMEOUT_S):
    """Run the core and return (returncode, stdout, stderr). If the page hangs
    up first (it cancelled the run), stop the core and return None, so a
    cancelled analysis does not keep the processor busy. Raises
    subprocess.TimeoutExpired after `timeout` seconds."""
    proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    done = {}
    worker = threading.Thread(target=lambda: done.update(out=proc.communicate()), daemon=True)
    worker.start()
    deadline = time.monotonic() + timeout
    while worker.is_alive():
        worker.join(0.1)
        if sock is not None and worker.is_alive() and hung_up(sock):
            proc.kill()
            worker.join()
            return None
        if worker.is_alive() and time.monotonic() > deadline:
            proc.kill()
            worker.join()
            raise subprocess.TimeoutExpired(argv, timeout)
    out, err = done["out"]
    return proc.returncode, out, err


def core_error(stderr, code):
    """The message of a failed core, from its stderr. A Julia failure to load
    ends in a trailer ("in expression starting at stdin:5"); the cause is the
    `ERROR:` line, up to the stack trace. Otherwise the last line."""
    lines = stderr.decode(errors="replace").strip().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("ERROR:"):
            msg = []
            for l in lines[i:]:
                if l.startswith("Stacktrace:") or (msg and l.startswith("ERROR:")):
                    break
                msg.append(l)
            return " ".join(msg).strip()[:600]
    return lines[-1] if lines else f"core exited with {code}"


class Live:
    """One `sigmatau live` process: its JSON lines are kept as they arrive
    and handed out by index, so a page can poll for the ones it has not seen."""

    def __init__(self, cmd, args):
        if not isinstance(args, list) or not args or args[0] != "live":
            raise ValueError("a live request starts with the live command")
        if not all(isinstance(a, str) for a in args) or any(a == "--out" or a.startswith("--out=") for a in args):
            raise ValueError("args must be strings, without --out")
        self.lines = []
        self.error = None
        self.lock = threading.Lock()
        self.proc = subprocess.Popen(cmd + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.reader = threading.Thread(target=self._read, daemon=True)
        self.reader.start()

    def _read(self):
        with self.proc.stdout:
            for line in self.proc.stdout:
                with self.lock:
                    self.lines.append(line.decode(errors="replace").strip())
        with self.proc.stderr:
            stderr = self.proc.stderr.read()
        if self.proc.wait():
            self.error = core_error(stderr, self.proc.returncode)

    def poll(self, since):
        running = self.reader.is_alive()
        with self.lock:
            lines = self.lines[since:]
        return {"lines": lines, "next": since + len(lines), "running": running, "error": self.error}

    def stop(self):
        if self.proc.poll() is None:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.proc.kill()


LOOPBACK = ("127.0.0.1", "::1", "localhost")
LIVE = {}            # id -> Live; one measurement at a time
LIVE_IDS = iter(range(1, 1 << 62))


def stop_live():
    for live in LIVE.values():
        live.stop()
    LIVE.clear()


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=str(UI), **kw)

    def log_message(self, fmt, *a):
        if os.environ.get("SIGMATAU_WEB_LOG"):
            super().log_message(fmt, *a)

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def foreign(self):
        """A request the page on this machine did not send: on loopback, one
        addressed to another host name (DNS rebinding); a POST that is not
        JSON (a cross-site form, which skips the browser's preflight)."""
        host = urlparse("//" + (self.headers.get("Host") or "")).hostname
        if self.server.server_address[0] in LOOPBACK and host not in LOOPBACK:
            return True
        return self.command == "POST" and self.headers.get("Content-Type", "").split(";")[0].strip() != "application/json"

    def send_json(self, status, payload):
        body = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        url = urlparse(self.path)
        if url.path.startswith("/api/") and self.foreign():
            return self.send_json(403, {"error": "only the explorer page on this computer may call the API"})
        if url.path == "/api/health":
            cmd = cli_command()
            return self.send_json(200, {"ok": cmd is not None, "core": cmd[0] if cmd else None})
        if url.path == "/api/live/poll":
            q = parse_qs(url.query)
            try:
                live = LIVE[int(q["id"][0])]
                since = int(q.get("since", ["0"])[0])
            except (KeyError, ValueError):
                return self.send_json(404, {"error": "no such live measurement"})
            return self.send_json(200, live.poll(since))
        return super().do_GET()

    def do_POST(self):
        if self.path not in ("/api/analyze", "/api/live/start", "/api/live/stop"):
            return self.send_json(404, {"error": "unknown endpoint"})
        if self.foreign():
            return self.send_json(403, {"error": "only the explorer page on this computer may call the API"})
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0 or length > MAX_BODY:
            return self.send_json(413, {"error": "request body missing or too large"})
        try:
            req = json.loads(self.rfile.read(length))
        except json.JSONDecodeError as e:
            return self.send_json(400, {"error": f"bad JSON: {e}"})
        if self.path == "/api/live/stop":
            stop_live()
            return self.send_json(200, {"ok": True})
        cmd = cli_command()
        if cmd is None:
            return self.send_json(503, {"error": "no analysis core: build cli/ (see cli/README.md) or put julia on PATH"})
        if self.path == "/api/live/start":
            args = req.get("args")
            # Device and file paths read this machine: only for a page on this machine.
            local = self.server.server_address[0] in LOOPBACK
            if not local and isinstance(args, list) and len(args) > 1 and args[1] != "sim" and not str(args[1]).startswith("tcp://"):
                return self.send_json(403, {"error": "device and file sources are only served on 127.0.0.1"})
            stop_live()
            try:
                live = Live(cmd, args)
            except (ValueError, TypeError) as e:
                return self.send_json(400, {"error": str(e)})
            live_id = next(LIVE_IDS)
            LIVE[live_id] = live
            return self.send_json(200, {"id": live_id})
        with tempfile.TemporaryDirectory(prefix="sigmatau-web-") as work:
            try:
                args = request_to_args(req, work)
            except (ValueError, TypeError) as e:
                return self.send_json(400, {"error": str(e)})
            try:
                # The compiled core runs on one thread unless told otherwise.
                env = dict(os.environ)
                env.setdefault("JULIA_NUM_THREADS", "auto")
                result = run_core(cmd + args, env, self.connection)
            except subprocess.TimeoutExpired:
                return self.send_json(504, {"error": "analysis timed out"})
        if result is None:
            return          # the page cancelled the run and hung up; nobody is listening
        code, stdout, stderr = result
        out = stdout.strip()
        if out.startswith(b"{"):
            return self.send_json(200 if code == 0 else 422, out)
        return self.send_json(500, {"error": core_error(stderr, code)})


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
    finally:
        stop_live()


if __name__ == "__main__":
    sys.exit(main())
