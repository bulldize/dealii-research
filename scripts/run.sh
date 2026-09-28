#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
# Override for another machine; no changes to the existing deal.II installation.
dealii_dir="${DEAL_II_DIR:-$repo_dir/../Dealii/install-local}"
export DYLD_LIBRARY_PATH="$dealii_dir/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
exec "$repo_dir/build/giesekus" "${1:-$repo_dir/configs/manufactured.prm}" "${2:-$repo_dir/results/manufactured}"
