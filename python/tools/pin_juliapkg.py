"""Rewrite python/src/sigmatau/juliapkg.json for a release build.

In the repo the file dev-installs SigmaTau.jl from this checkout ("path":
"../../..", "dev": true), which only works from a checkout. A wheel must name a
SigmaTau.jl that juliapkg can fetch on its own:

    python python/tools/pin_juliapkg.py git        # url + rev = v<version> (tag on GitHub)
    python python/tools/pin_juliapkg.py registry   # version = <version> (General registry)

<version> is sigmatau.__version__, which must equal Project.toml's version.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
JULIAPKG = ROOT / "python" / "src" / "sigmatau" / "juliapkg.json"
INIT = ROOT / "python" / "src" / "sigmatau" / "__init__.py"
PROJECT = ROOT / "Project.toml"
URL = "https://github.com/ianlap/SigmaTau.jl"


def read_version(path: Path, pattern: str) -> str:
    m = re.search(pattern, path.read_text(), re.MULTILINE)
    if m is None:
        sys.exit(f"no version found in {path}")
    return m.group(1)


def main() -> None:
    if len(sys.argv) != 2 or sys.argv[1] not in ("git", "registry"):
        sys.exit("usage: pin_juliapkg.py git|registry")
    py = read_version(INIT, r'^__version__ = "([^"]+)"')
    jl = read_version(PROJECT, r'^version = "([^"]+)"')
    if py != jl:
        sys.exit(f"sigmatau.__version__ {py} differs from Project.toml version {jl}")
    spec = json.loads(JULIAPKG.read_text())
    entry = spec["packages"]["SigmaTau"]
    pinned = {"uuid": entry["uuid"]}
    if sys.argv[1] == "git":
        pinned |= {"url": URL, "rev": f"v{jl}"}
    else:
        pinned |= {"version": f"={jl}"}
    spec["packages"]["SigmaTau"] = pinned
    JULIAPKG.write_text(json.dumps(spec, indent=4) + "\n")
    print(json.dumps(spec, indent=4))


if __name__ == "__main__":
    main()
