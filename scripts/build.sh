#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
dealii_dir="${DEAL_II_DIR:-$repo_dir/../Dealii/install-local}"
cmake -S "$repo_dir" -B "$repo_dir/build" -DCMAKE_BUILD_TYPE=Release -DDEAL_II_DIR="$dealii_dir"
cmake --build "$repo_dir/build" -j "${BUILD_JOBS:-4}"
ctest --test-dir "$repo_dir/build" --output-on-failure
