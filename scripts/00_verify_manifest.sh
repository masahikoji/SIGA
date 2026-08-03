#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
if command -v sha256sum >/dev/null 2>&1; then
  HASH='sha256sum'
elif command -v shasum >/dev/null 2>&1; then
  HASH='shasum -a 256'
else
  echo 'Neither sha256sum nor shasum is available.' >&2
  exit 1
fi
while IFS='  ' read -r expected file; do
  [ -z "$expected" ] && continue
  actual=$($HASH "$file" | awk '{print $1}')
  if [ "$actual" != "$expected" ]; then
    echo "Hash mismatch: $file" >&2
    exit 1
  fi
done < MANIFEST_SHA256.txt
echo 'PASS: every file listed in MANIFEST_SHA256.txt matches.'
