#!/usr/bin/env bash
# Build the desktop app and put the compiled core inside it.
#
# The core goes in after `cargo tauri build` rather than through Tauri's
# "resources", because Tauri copies symlinks as files and the juliac bundle
# has about 60 of them (libLLVM, libjulia, ...), which would more than double
# the app's size. Build the core first (cli/README.md).
set -euo pipefail
cd "$(dirname "$0")"
CORE=../../cli/build
[ -x "$CORE/bin/sigmatau" ] || { echo "no core at cli/build: build it first (cli/README.md)" >&2; exit 1; }
cargo tauri build "$@"
case "$(uname -s)" in
  Darwin)
    APP="target/release/bundle/macos/SigmaTau Explorer.app"
    rsync -a --delete "$CORE/" "$APP/Contents/Resources/core/"
    chmod -R u+w "$APP/Contents/Resources/core"
    codesign --force --deep --sign - "$APP"
    echo "$APP: $(du -sh "$APP" | cut -f1)"
    ;;
  *)
    echo "Copy $CORE/ next to the executable as core/ (see README)." >&2
    ;;
esac
