#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SHELL_DIR="$ROOT/tools/ketcher-shell"
OUT="$ROOT/web/vendor/ketcher"
VER="3.7.0"

cd "$SHELL_DIR"
if [[ ! -d node_modules ]]; then
  npm install
fi
npm run build

rm -rf "$OUT"
mkdir -p "$OUT"
cp -R "$SHELL_DIR/dist/." "$OUT/"
echo "$VER" > "$OUT/VERSION"

# ketcher-core / indigo deps expect Node globals in the browser.
# Source shell index.html ships the shim; fail loudly if Vite drops it.
if ! grep -q 'window.global' "$OUT/index.html"; then
  echo "error: $OUT/index.html missing window.global shim (ketcher will fail with 'global is not defined')" >&2
  exit 1
fi

echo "Wrote $OUT (ketcher $VER)"
