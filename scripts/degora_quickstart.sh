#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# DEGORA quickstart: prepare an environment, build the demo, open the browser.
#
#   Inside a checkout:   bash scripts/degora_quickstart.sh
#   Standalone:          bash degora_quickstart.sh          (clones into ./DEGORA)
#
# Supported: Ubuntu (and other Linux), macOS (Intel and Apple silicon),
# Windows 11 via WSL2 Ubuntu. Re-running is safe.
#
# Options:
#   --port N        Preferred browser/API port (default 8765; the next free
#                   port is used when it is taken).
#   --dir PATH      Where to place or find the checkout when run standalone.
#   --ref NAME      Branch or tag to check out (default: the repository default).
#   --config PATH   Serve this DEGORA config instead of the bundled demo.
#   --update        git pull an existing checkout before installing.
#   --no-browser    Do not try to open a browser (headless or remote shells).
#   --no-demo       Skip demo creation; serve an existing database.
#   -h, --help      Show this help.
# ---------------------------------------------------------------------------
set -euo pipefail

PORT=8765
TARGET_DIR=""
GIT_REF=""
CONFIG_PATH=""
DO_UPDATE=0
OPEN_BROWSER=1
BUILD_DEMO=1
REPO_URL="https://github.com/kangk1204/DEGORA.git"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
note() { printf '==> %s\n' "$*"; }

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; }

while [ "$#" -gt 0 ]; do
  case "$1" in
    --port) [ "$#" -ge 2 ] || die "--port needs a value"; PORT="$2"; shift 2 ;;
    --dir) [ "$#" -ge 2 ] || die "--dir needs a value"; TARGET_DIR="$2"; shift 2 ;;
    --ref) [ "$#" -ge 2 ] || die "--ref needs a value"; GIT_REF="$2"; shift 2 ;;
    --config) [ "$#" -ge 2 ] || die "--config needs a value"; CONFIG_PATH="$2"; shift 2 ;;
    --update) DO_UPDATE=1; shift ;;
    --no-browser) OPEN_BROWSER=0; shift ;;
    --no-demo) BUILD_DEMO=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option: $1 (try --help)" ;;
  esac
done

case "$PORT" in
  ''|*[!0-9]*) die "--port must be a whole number" ;;
esac
[ "$PORT" -ge 1 ] && [ "$PORT" -le 65535 ] || die "--port must be between 1 and 65535"

# --- platform -------------------------------------------------------------
case "$(uname -s)" in
  Darwin) PLATFORM=macos ;;
  Linux)
    if grep -qiE '(microsoft|wsl)' /proc/version 2>/dev/null; then PLATFORM=wsl; else PLATFORM=linux; fi
    ;;
  *) PLATFORM=other ;;
esac
note "Platform: $PLATFORM ($(uname -s) $(uname -m))"

command -v git >/dev/null 2>&1 || die "git is required. Install it (macOS: xcode-select --install; Ubuntu: sudo apt install git) and re-run."

# --- interpreter ----------------------------------------------------------
# macOS /usr/bin/python3 is often 3.9, which the package metadata rejects.
supported_python() {
  "$1" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' >/dev/null 2>&1
}

PY=""
for candidate in python3.12 python3.13 python3.11 python3.10 python3 python; do
  if command -v "$candidate" >/dev/null 2>&1 && supported_python "$candidate"; then
    PY="$candidate"
    break
  fi
done
if [ -z "$PY" ]; then
  printf '\n'
  printf 'ERROR: no Python 3.10 or newer was found.\n' >&2
  case "$PLATFORM" in
    macos) printf 'Install one with:  brew install python@3.12\n' >&2 ;;
    linux|wsl) printf 'Install one with:  sudo apt update && sudo apt install python3.12 python3.12-venv\n' >&2 ;;
    *) printf 'Install Python 3.10 or newer, then re-run this script.\n' >&2 ;;
  esac
  exit 1
fi
note "Python: $PY ($("$PY" --version 2>&1))"

# --- locate or create the checkout ---------------------------------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT=""
for candidate in "$SCRIPT_DIR/.." "$SCRIPT_DIR"; do
  if [ -f "$candidate/pyproject.toml" ] && grep -q '^name = "degora"' "$candidate/pyproject.toml" 2>/dev/null; then
    REPO_ROOT="$(cd -- "$candidate" && pwd -P)"
    break
  fi
done

if [ -n "$REPO_ROOT" ]; then
  note "Using this checkout: $REPO_ROOT"
else
  BASE_DIR="${TARGET_DIR:-$PWD}"
  mkdir -p "$BASE_DIR"
  BASE_DIR="$(cd -- "$BASE_DIR" && pwd -P)"
  if [ -d "$BASE_DIR/DEGORA/.git" ]; then
    REPO_ROOT="$BASE_DIR/DEGORA"
    note "Found an existing checkout: $REPO_ROOT"
  else
    note "Cloning DEGORA into $BASE_DIR/DEGORA"
    git clone --quiet "$REPO_URL" "$BASE_DIR/DEGORA"
    REPO_ROOT="$BASE_DIR/DEGORA"
  fi
fi
cd "$REPO_ROOT"

if [ -n "$GIT_REF" ]; then
  note "Checking out $GIT_REF"
  git fetch --quiet origin "$GIT_REF" \
    || die "could not fetch '$GIT_REF' from origin; check the branch or tag name."
  REF_COMMIT="$(git rev-parse --verify --quiet "FETCH_HEAD^{commit}")" \
    || die "'$GIT_REF' does not point at a commit."
  if git checkout --quiet "$GIT_REF" 2>/dev/null; then
    # A branch or tag of that name already existed locally; make sure it is
    # not an old copy, otherwise the run would silently serve stale code.
    if [ "$(git rev-parse HEAD)" != "$REF_COMMIT" ]; then
      git merge --ff-only --quiet "$REF_COMMIT" 2>/dev/null \
        || die "local '$GIT_REF' has diverged from origin; rename or delete it and re-run."
    fi
  else
    git checkout --quiet -b "$GIT_REF" "$REF_COMMIT" || die "could not check out '$GIT_REF'."
  fi
elif [ "$DO_UPDATE" -eq 1 ]; then
  note "Updating the checkout"
  git pull --ff-only --quiet || die "git pull failed; resolve local changes and re-run without --update."
fi
if git rev-parse --git-dir >/dev/null 2>&1; then
  note "Commit: $(git log --oneline -1 2>/dev/null || echo unknown)"
fi

# --- virtual environment --------------------------------------------------
VENV_DIR="$REPO_ROOT/.venv"
if [ -d "$VENV_DIR" ] && ! supported_python "$VENV_DIR/bin/python"; then
  note "Existing .venv uses an unsupported Python; recreating it"
  rm -rf "$VENV_DIR"
fi
if [ ! -d "$VENV_DIR" ]; then
  note "Creating the virtual environment"
  if ! "$PY" -m venv "$VENV_DIR" 2>/tmp/degora_venv_error.$$; then
    printf '\n' >&2
    cat /tmp/degora_venv_error.$$ >&2 || true
    rm -f /tmp/degora_venv_error.$$
    printf '\n' >&2
    case "$PLATFORM" in
      linux|wsl) printf 'Debian/Ubuntu ships venv separately. Install it with:\n  sudo apt install %s-venv\n' "$PY" >&2 ;;
      *) printf 'Could not create a virtual environment with %s.\n' "$PY" >&2 ;;
    esac
    exit 1
  fi
  rm -f /tmp/degora_venv_error.$$
fi
# shellcheck source=/dev/null
. "$VENV_DIR/bin/activate"

"$VENV_DIR/bin/python" -m pip install --upgrade pip --quiet
note "Installing DEGORA and dependencies (the first run takes about a minute)"
"$VENV_DIR/bin/python" -m pip install -e . --quiet
note "Installed: $(degora --version)"

# --- workspace ------------------------------------------------------------
if [ -n "$CONFIG_PATH" ]; then
  [ -f "$CONFIG_PATH" ] || die "config not found: $CONFIG_PATH"
  note "Running your config: $CONFIG_PATH"
  RUN_LOG="$(mktemp)"
  degora run "$CONFIG_PATH" 2>&1 | tee "$RUN_LOG"
  # `degora run` prints "- Database: <path>" on success.
  DB_PATH="$(sed -n 's/^- Database: //p' "$RUN_LOG" | tail -1)"
  rm -f "$RUN_LOG"
  if [ -z "$DB_PATH" ] || [ ! -f "$DB_PATH" ]; then
    die "Could not locate the score database for $CONFIG_PATH. Run 'degora serve <output_dir>/degora_scores.db' directly."
  fi
elif [ "$BUILD_DEMO" -eq 1 ]; then
  note "Building the demo workspace"
  rm -rf degora-demo
  degora demo degora-demo >/dev/null
  degora run degora-demo/degora_demo_config.xlsx
  DB_PATH="$REPO_ROOT/degora-demo/results/degora_scores.db"
else
  DB_PATH="$REPO_ROOT/degora-demo/results/degora_scores.db"
  [ -f "$DB_PATH" ] || die "No database at $DB_PATH. Re-run without --no-demo."
fi

# --- choose a free port so the printed URL is always correct --------------
PORT="$(
  "$VENV_DIR/bin/python" - "$PORT" <<'PY'
import socket
import sys

start = int(sys.argv[1])
for port in range(start, min(start + 50, 65536)):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        try:
            probe.bind(("127.0.0.1", port))
        except OSError:
            continue
    print(port)
    break
else:
    print(start)
PY
)"

URL="http://127.0.0.1:${PORT}"
cat <<EOF

==================================================
 Dashboard:  ${URL}
 Database:   ${DB_PATH}
 Checkout:   ${REPO_ROOT}
 Stop:       press Ctrl+C in this window
==================================================

EOF

# --- open a browser once the server is listening --------------------------
if [ "$OPEN_BROWSER" -eq 1 ]; then
  (
    sleep 3
    case "$PLATFORM" in
      macos) command -v open >/dev/null 2>&1 && open "$URL" ;;
      wsl)
        if command -v wslview >/dev/null 2>&1; then wslview "$URL"
        elif command -v explorer.exe >/dev/null 2>&1; then explorer.exe "$URL" || true
        fi
        ;;
      linux) command -v xdg-open >/dev/null 2>&1 && xdg-open "$URL" ;;
      *) : ;;
    esac
  ) >/dev/null 2>&1 &
fi

exec degora serve "$DB_PATH" --port "$PORT"
