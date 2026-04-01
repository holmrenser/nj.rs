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

cargo workspaces version "$kind" --yes --no-git-commit

if git diff --quiet Cargo.toml Cargo.lock */Cargo.toml; then
  exit 0
fi

VERSION=$(toml get Cargo.toml workspace.package.version | tr -d '"')

toml set python/pyproject.toml project.version "$VERSION" | sponge python/pyproject.toml
jq --arg v "$VERSION" '.version = $v' wasm/package.json | sponge wasm/package.json

git commit --amend -am "Release v$VERSION"
git tag "v$VERSION"