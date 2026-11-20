#!/usr/bin/env bash
set -euxo pipefail

kind="$1"

case "$kind" in
  patch|minor|major) ;;
  *)
    echo "Usage: release.sh {patch|minor|major}"
    exit 1
    ;;
esac

# Refuse to release without changelog entries (before touching any version).
bash scripts/changelog.sh check

cargo workspaces version "$kind" --yes --no-git-commit

if git diff --quiet Cargo.toml Cargo.lock */Cargo.toml; then
  exit 0
fi

VERSION=$(toml get Cargo.toml workspace.package.version | tr -d '"')

toml set python/pyproject.toml project.version "$VERSION" | sponge python/pyproject.toml
jq --arg v "$VERSION" '.version = $v' wasm/package.json | sponge wasm/package.json

# Sanity-check that all three package manifests now agree before tagging.
bash scripts/check-versions.sh

# Turn [Unreleased] into this version's section so the release commit includes it.
bash scripts/changelog.sh rotate "$VERSION"

# Create a fresh release commit (never --amend: that would rewrite the previous,
# possibly already-pushed, commit).
git commit -am "Release v$VERSION"
git tag "v$VERSION"