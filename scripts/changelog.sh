#!/usr/bin/env bash
#
# Helpers for keeping CHANGELOG.md in sync with releases.
#
#   changelog.sh check            fail unless [Unreleased] has content
#   changelog.sh rotate X.Y.Z     turn [Unreleased] into [X.Y.Z] - <today>, and
#                                 start a fresh empty [Unreleased] above it
#   changelog.sh section X.Y.Z    print the body of the [X.Y.Z] section
#
# Uses only awk/grep so it needs no extra tooling in CI.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
file="${CHANGELOG_FILE:-$root/CHANGELOG.md}"

# Prints the body of the section whose heading starts with "## [$1]".
section() {
  awk -v head="## [$1]" '
    index($0, head) == 1 { found = 1; next }
    found && /^## \[/    { exit }
    found                { print }
  ' "$file"
}

# True when the text on stdin contains a non-blank, non-heading line.
has_content() { grep -qE '^[[:space:]]*[^[:space:]#]' ; }

cmd="${1:-}"
case "$cmd" in
  check)
    if ! section Unreleased | has_content; then
      echo "CHANGELOG.md: the [Unreleased] section is empty; document the release first." >&2
      exit 1
    fi
    ;;
  rotate)
    version="${2:?usage: changelog.sh rotate X.Y.Z}"
    "$0" check
    if grep -q "^## \[$version\]" "$file"; then
      echo "CHANGELOG.md already has a section for $version" >&2
      exit 1
    fi
    tmp="$(mktemp)"
    awk -v v="$version" -v d="$(date +%F)" '
      !done && /^## \[Unreleased\]/ {
        print "## [Unreleased]\n\n## [" v "] - " d
        done = 1; next
      }
      { print }
    ' "$file" > "$tmp"
    mv "$tmp" "$file"
    ;;
  section)
    version="${2:?usage: changelog.sh section X.Y.Z}"
    body="$(section "$version")"
    if ! printf '%s\n' "$body" | has_content; then
      echo "CHANGELOG.md has no entry for $version" >&2
      exit 1
    fi
    # Strip leading and trailing blank lines.
    printf '%s\n' "$body" | awk 'NF { started = 1 } started { lines[++n] = $0 } END { while (n > 0 && lines[n] == "") n--; for (i = 1; i <= n; i++) print lines[i] }'
    ;;
  *)
    echo "Usage: changelog.sh {check|rotate X.Y.Z|section X.Y.Z}" >&2
    exit 2
    ;;
esac
