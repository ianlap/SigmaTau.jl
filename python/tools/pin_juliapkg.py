"""Pin python/src/sigmatau/juliapkg.json to the registered SigmaTau.jl release.

In the repo the file dev-installs SigmaTau.jl from this checkout ("path":
"../../..", "dev": true), which only works from a checkout. A wheel instead
names the General-registry release whose version equals sigmatau.__version__:

    python python/tools/pin_juliapkg.py      # "version": "=<version>"

It exits non-zero if __version__ differs from Project.toml's version, or if
that version is not yet in the General registry.
"""

from __future__ import annotations

import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
JULIAPKG = ROOT / "python" / "src" / "sigmatau" / "juliapkg.json"
INIT = ROOT / "python" / "src" / "sigmatau" / "__init__.py"
PROJECT = ROOT / "Project.toml"
REGISTRY = "https://raw.githubusercontent.com/JuliaRegistries/General/master/S/SigmaTau/Versions.toml"


def read_version(path: Path, pattern: str) -> str:
    m = re.search(pattern, path.read_text(), re.MULTILINE)
    if m is None:
        sys.exit(f"no version found in {path}")
    return m.group(1)


def registered(version: str) -> bool:
    try:
        with urllib.request.urlopen(REGISTRY, timeout=30) as r:
            text = r.read().decode()
    except urllib.error.HTTPError as e:
        if e.code == 404:  # the package itself is not registered yet
            return False
        sys.exit(f"could not read the General registry ({e})")
    except OSError as e:
        sys.exit(f"could not read the General registry ({e})")
    return f'["{version}"]' in text


def main() -> None:
    py = read_version(INIT, r'^__version__ = "([^"]+)"')
    jl = read_version(PROJECT, r'^version = "([^"]+)"')
    if py != jl:
        sys.exit(f"sigmatau.__version__ {py} differs from Project.toml version {jl}")
    if not registered(jl):
        sys.exit(f"SigmaTau.jl {jl} is not in the General registry yet; register it first")
    spec = json.loads(JULIAPKG.read_text())
    entry = spec["packages"]["SigmaTau"]
    spec["packages"]["SigmaTau"] = {"uuid": entry["uuid"], "version": f"={jl}"}
    JULIAPKG.write_text(json.dumps(spec, indent=4) + "\n")
    print(json.dumps(spec, indent=4))


if __name__ == "__main__":
    main()
