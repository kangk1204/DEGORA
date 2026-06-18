#!/usr/bin/env bash
# Download and unpack the archived DEGORA reproduction data (DEG tables + locked
# gold panels) so the configs in reproducibility/datasets/ can run.
#
# The configs reference data/deg/raw/... paths; this script restores those files
# at the repository root. Run it ONCE from the repo root:
#
#     bash reproducibility/fetch_reproduction_data.sh
#
# The data archive lives on Zenodo (see the Data Availability statement in the
# paper). Set DEGORA_REPRO_DATA_URL to override the download URL, or point
# DEGORA_REPRO_DATA_ZIP at an already-downloaded zip to skip the download.
set -euo pipefail

VERSION="v1"
BUNDLE="degora_reproduction_data_${VERSION}"
# Zenodo record for the v1 reproduction-data archive (DOI 10.5281/zenodo.20739978).
DEFAULT_URL="https://zenodo.org/records/20739978/files/${BUNDLE}.zip"
URL="${DEGORA_REPRO_DATA_URL:-$DEFAULT_URL}"
EXPECTED_SHA256="8e757e96319e215d18f9071d0bd84ac21ff393fc3881562de8174af946d52b82"

# Locate the repo root (the folder that has both outputs/ and reproducibility/).
here="$(cd "$(dirname "$0")/.." && pwd)"
if [ ! -d "$here/outputs" ] || [ ! -d "$here/reproducibility" ]; then
  echo "error: run this from the DEGORA repo root (could not find outputs/ + reproducibility/)." >&2
  exit 1
fi
cd "$here"

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
zip_path="${DEGORA_REPRO_DATA_ZIP:-$work/${BUNDLE}.zip}"

if [ -n "${DEGORA_REPRO_DATA_ZIP:-}" ]; then
  echo "Using local bundle: $zip_path"
else
  if [[ "$URL" == *PLACEHOLDER* ]]; then
    cat >&2 <<'MSG'
error: the Zenodo download URL is still a PLACEHOLDER.
  Either:
    - edit DEFAULT_URL in this script with the published Zenodo record id, or
    - download the bundle manually and re-run with:
        DEGORA_REPRO_DATA_ZIP=/path/to/degora_reproduction_data_v1.zip \
          bash reproducibility/fetch_reproduction_data.sh
MSG
    exit 1
  fi
  echo "Downloading reproduction data:"
  echo "  $URL"
  curl -L --fail -o "$zip_path" "$URL"
fi

# Verify the archive checksum (MANDATORY): a corrupted, truncated, or wrong-version
# bundle must not silently reproduce. python3 is already required (used below to
# unpack), so this check always runs regardless of whether sha256sum is installed.
echo "Verifying bundle checksum..."
got="$(python3 - "$zip_path" <<'PY'
import hashlib, sys
h = hashlib.sha256()
with open(sys.argv[1], "rb") as fh:
    for chunk in iter(lambda: fh.read(1 << 20), b""):
        h.update(chunk)
print(h.hexdigest())
PY
)"
if [ "$got" = "$EXPECTED_SHA256" ]; then
  echo "checksum OK ($got)"
elif [ "${DEGORA_REPRO_DATA_SKIP_CHECKSUM:-0}" = "1" ]; then
  echo "warning: SHA-256 mismatch, but DEGORA_REPRO_DATA_SKIP_CHECKSUM=1 is set; continuing." >&2
  echo "  expected $EXPECTED_SHA256 ; got $got" >&2
else
  echo "error: bundle SHA-256 mismatch -- refusing to unpack a non-canonical archive." >&2
  echo "  expected $EXPECTED_SHA256" >&2
  echo "  got      $got" >&2
  echo "  The download may be corrupted, truncated, or a different bundle version." >&2
  echo "  (To use a different bundle on purpose, re-run with DEGORA_REPRO_DATA_SKIP_CHECKSUM=1.)" >&2
  exit 1
fi

echo "Unpacking..."
python3 -m zipfile -e "$zip_path" "$work/unpacked"

# The archive holds <BUNDLE>/data/... ; merge that data/ into the repo root.
src="$work/unpacked/${BUNDLE}/data"
if [ ! -d "$src" ]; then
  echo "error: expected data/ inside the archive at $src" >&2
  exit 1
fi
mkdir -p data
cp -r "$src/." data/

n=$(find data/deg/raw -type f 2>/dev/null | wc -l | tr -d ' ')
echo ""
echo "Done. Restored $n DEG-table files under data/deg/raw/ (+ gold panels under data/studies/gold/)."
echo "You can now run, e.g.:"
echo "  degora run reproducibility/datasets/01_ifn_rnaseq/config.xlsx"
