#!/usr/bin/env bash
#
# Verifies that the version string is identical across the Cargo workspace,
# the Python package, and the npm package. These three are bumped together by
# release.sh; this guard catches drift on every PR rather than only at tag push.
#
# Uses only grep/sed so it needs no extra tooling in CI.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Version from [workspace.package] in the root Cargo.toml.
cargo_version=$(
  sed -n '/^\[workspace\.package\]/,/^\[/p' "$root/Cargo.toml" \
    | grep -m1 '^[[:space:]]*version' | sed -E 's/.*"([^"]+)".*/\1/'
)
# Version from [project] in pyproject.toml.
py_version=$(
  sed -n '/^\[project\]/,/^\[/p' "$root/python/pyproject.toml" \
    | grep -m1 '^[[:space:]]*version' | sed -E 's/.*"([^"]+)".*/\1/'
)
# "version": "x.y.z" in package.json.
npm_version=$(
  grep -m1 '"version"' "$root/wasm/package.json" | sed -E 's/.*"version"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/'
)

echo "cargo (workspace): $cargo_version"
echo "python (pyproject): $py_version"
echo "npm (package.json): $npm_version"

if [[ -z "$cargo_version" ]]; then
  echo "ERROR: could not read workspace version from Cargo.toml" >&2
  exit 1
fi

if [[ "$cargo_version" != "$py_version" || "$cargo_version" != "$npm_version" ]]; then
  echo "ERROR: version mismatch — run 'make bump-patch' (or minor/major) to resync." >&2
  exit 1
fi

echo "OK: all package versions match ($cargo_version)"
