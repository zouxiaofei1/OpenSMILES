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
echo "Wrote $OUT (ketcher $VER)"
