#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# DEGORA quickstart (repository-root forwarder).
#
# The real script lives at scripts/degora_quickstart.sh. This forwarder exists
# so the command works from the root of a checkout too, which is where a reader
# who has just run `git clone && cd DEGORA` actually is:
#
#   bash degora_quickstart.sh --ref v0.4.17
#
# Every option is passed through unchanged; see scripts/degora_quickstart.sh
# or run `bash degora_quickstart.sh --help`.
#
# Downloaded on its own, outside a checkout, this file cannot do the work: it
# says which file to download instead rather than failing with "No such file".
# ---------------------------------------------------------------------------
set -euo pipefail

here=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
real="$here/scripts/degora_quickstart.sh"

if [ ! -f "$real" ]; then
  printf 'ERROR: %s\n' "this is only the repository-root forwarder; scripts/degora_quickstart.sh is missing." >&2
  printf '%s\n' "To run DEGORA standalone, download the real script instead:" >&2
  printf '%s\n' "  curl -fsSLO https://raw.githubusercontent.com/kangk1204/DEGORA/main/scripts/degora_quickstart.sh" >&2
  printf '%s\n' "  bash degora_quickstart.sh" >&2
  exit 1
fi

exec bash "$real" "$@"
