#!/usr/bin/env bash
# Build the desktop app with the compiled core inside it.
#
#   app/desktop/bundle.sh            # macOS: .app, .dmg and .zip; Linux: .deb and .AppImage
#
# Build the core first (cli/README.md); the outputs land in app/desktop/dist/.
#
# macOS: the core goes in after `cargo tauri build` rather than through
# Tauri's "resources", because Tauri copies symlinks as files and the juliac
# bundle has about 60 of them (libLLVM, libjulia, ...), which would more than
# double the app's size. The app is then signed:
#   - with $APPLE_SIGNING_IDENTITY (a "Developer ID Application: ..." identity
#     in the keychain) and the hardened runtime, and notarized when
#     $APPLE_ID, $APPLE_PASSWORD (an app-specific password) and $APPLE_TEAM_ID
#     are set;
#   - otherwise ad hoc, which runs on this Mac and, on other Macs, after
#     right-click > Open (see README).
# Linux: Tauri ships the core as a resource (tauri.linux.conf.json).
set -euo pipefail
cd "$(dirname "$0")"
CORE=../../cli/build
[ -x "$CORE/bin/sigmatau" ] || { echo "no core at cli/build: build it first (cli/README.md)" >&2; exit 1; }
chmod -R u+w "$CORE"   # juliac leaves some dylibs read-only, which breaks re-copies
VERSION=$(python3 -c 'import json; print(json.load(open("tauri.conf.json"))["version"])')
ARCH=$(uname -m)
mkdir -p dist

# The app is signed below, after the core goes in. Keep Tauri from signing
# it first: in CI the signing variables are set but empty when the secrets
# are absent, and Tauri then fails looking for an identity named "".
env -u APPLE_SIGNING_IDENTITY -u APPLE_ID -u APPLE_PASSWORD -u APPLE_TEAM_ID \
  -u APPLE_CERTIFICATE -u APPLE_CERTIFICATE_PASSWORD cargo tauri build "$@"

case "$(uname -s)" in
  Darwin)
    APP="target/release/bundle/macos/SigmaTau Explorer.app"
    rsync -a --delete "$CORE/" "$APP/Contents/Resources/core/"
    chmod -R u+w "$APP/Contents/Resources/core"
    if [ -n "${APPLE_SIGNING_IDENTITY:-}" ]; then
      sign() { codesign --force --timestamp --options runtime --entitlements entitlements.plist --sign "$APPLE_SIGNING_IDENTITY" "$@"; }
      # Inside out: every dylib and the core, then the app.
      find "$APP/Contents/Resources/core" -type f \( -name '*.dylib' -o -perm -u+x \) -print0 |
        while IFS= read -r -d '' f; do file -b "$f" | grep -q Mach-O && sign "$f"; done
      sign "$APP"
    else
      codesign --force --deep --sign - "$APP"
    fi
    codesign --verify --deep --strict "$APP"

    NAME="SigmaTau-Explorer-$VERSION-macos-$ARCH"
    rm -f "dist/$NAME.zip" "dist/$NAME.dmg"
    ditto -c -k --keepParent "$APP" "dist/$NAME.zip"
    if [ -n "${APPLE_SIGNING_IDENTITY:-}" ] && [ -n "${APPLE_ID:-}" ]; then
      xcrun notarytool submit "dist/$NAME.zip" --apple-id "$APPLE_ID" \
        --password "$APPLE_PASSWORD" --team-id "$APPLE_TEAM_ID" --wait
      xcrun stapler staple "$APP"
      rm "dist/$NAME.zip"
      ditto -c -k --keepParent "$APP" "dist/$NAME.zip"
    fi
    STAGE=$(mktemp -d)
    cp -R "$APP" "$STAGE/"
    ln -s /Applications "$STAGE/Applications"
    hdiutil create -quiet -volname "SigmaTau Explorer" -srcfolder "$STAGE" -format UDZO "dist/$NAME.dmg"
    rm -rf "$STAGE"
    [ -n "${APPLE_SIGNING_IDENTITY:-}" ] && codesign --force --timestamp --sign "$APPLE_SIGNING_IDENTITY" "dist/$NAME.dmg"
    echo "$APP: $(du -sh "$APP" | cut -f1)"
    ;;
  Linux)
    cp target/release/bundle/deb/*.deb dist/
    cp target/release/bundle/appimage/*.AppImage dist/
    ;;
esac
ls -lh dist
