#!/usr/bin/env bash
# Helper notes for distributing PIC++GUI. This does not produce signed installers;
# it prints the local binary path and dependency hints for packaging tools.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/build/bin/PIC++GUI"

if [[ ! -x "$BIN" ]]; then
  echo "Build the GUI first: ./scripts/build.sh" >&2
  exit 1
fi

echo "PIC++GUI binary: $BIN"
echo ""
echo "Packaging hints:"
echo "  macOS:   use macdeployqt after brew install qt"
echo "           macdeployqt \"$BIN\" -dmg"
echo "  Linux:   consider linuxdeploy / AppImage with Qt plugins bundled"
echo "  Windows: windeployqt next to PIC++GUI.exe"
echo ""
echo "Runtime also needs the repo's inputFiles/demo/ tree (or pass the repo root)."
