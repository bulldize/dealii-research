#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_dir"
run_root="${RUN_ROOT:-$repo_dir/results/full-paper}"
mkdir -p "$run_root"
exec 9>"$run_root/pipeline.lock"
flock -n 9 || { echo 'A pipeline is already running'; exit 2; }
echo "$$" > "$run_root/pipeline.pid"
phase(){ printf '%s %s\n' "$(date -Is)" "$1" > "$run_root/phase.txt"; echo "$(cat "$run_root/phase.txt")"; }
trap 'rc=$?; if [ "$rc" -ne 0 ]; then phase "FAILED exit=$rc (see pipeline.log)"; fi' EXIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DEAL_II_NUM_THREADS=1
export DEAL_II_DIR="${DEAL_II_DIR:-$HOME/opt/dealii-9.7.1}"
numlib="$HOME/opt/numerics/usr/lib/x86_64-linux-gnu"
export LD_LIBRARY_PATH="$numlib:$numlib/openblas-pthread:$numlib/blas:$numlib/lapack:$DEAL_II_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
phase 'Preparing user-local numerical libraries'
if [ ! -f "$numlib/openblas-pthread/libopenblas.so.0" ]; then
  mkdir -p "$HOME/opt/numerics/debs"
  (cd "$HOME/opt/numerics/debs"
   apt-get download libopenblas0-pthread libgfortran5 libquadmath0
   for package in ./*.deb; do dpkg-deb -x "$package" "$HOME/opt/numerics"; done)
fi
phase 'Building deal.II 9.7.1 in user directory'
if [ ! -f "$DEAL_II_DIR/lib/cmake/deal.II/deal.IIConfig.cmake" ]; then
  mkdir -p "$HOME/opt/dealii-9.7.1-source"
  if [ -f "$HOME/opt/dealii-9.7.1-minimal.tar.xz" ] && xz -t "$HOME/opt/dealii-9.7.1-minimal.tar.xz"; then
    tar -xJf "$HOME/opt/dealii-9.7.1-minimal.tar.xz" -C "$HOME/opt/dealii-9.7.1-source"
  else
    curl --fail --location --retry 5 --connect-timeout 30 https://codeload.github.com/dealii/dealii/tar.gz/refs/tags/v9.7.1 -o "$HOME/opt/dealii-9.7.1-github.tar.gz"
    tar -xzf "$HOME/opt/dealii-9.7.1-github.tar.gz" --strip-components=1 -C "$HOME/opt/dealii-9.7.1-source"
  fi
  cmake -G Ninja -S "$HOME/opt/dealii-9.7.1-source" -B "$HOME/opt/dealii-9.7.1-build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$DEAL_II_DIR" \
    -DDEAL_II_ALLOW_AUTODETECTION=OFF -DDEAL_II_ALLOW_BUNDLED=ON \
    -DDEAL_II_WITH_LAPACK=ON -DLAPACK_LIBRARIES="$numlib/openblas-pthread/libopenblas.so.0" \
    -DBLAS_LIBRARIES="$numlib/openblas-pthread/libopenblas.so.0" -DDEAL_II_WITH_UMFPACK=ON \
    -DDEAL_II_FORCE_BUNDLED_UMFPACK=ON -DDEAL_II_COMPONENT_EXAMPLES=OFF \
    -DDEAL_II_COMPONENT_DOCUMENTATION=OFF -DDEAL_II_COMPONENT_PYTHON_BINDINGS=OFF \
    > "$run_root/dealii-configure.log" 2>&1
  cmake --build "$HOME/opt/dealii-9.7.1-build" -j "${BUILD_JOBS:-6}" > "$run_root/dealii-build.log" 2>&1
  cmake --install "$HOME/opt/dealii-9.7.1-build" > "$run_root/dealii-install.log" 2>&1
fi
phase 'Building and testing research library'
BUILD_JOBS="${BUILD_JOBS:-6}" scripts/build.sh > "$run_root/research-build.log" 2>&1
phase 'Running manufactured and contraction smoke checks'
[ -f "$run_root/smoke-mms/summary.json" ] || scripts/run.sh configs/manufactured.prm "$run_root/smoke-mms" > "$run_root/smoke-mms.log" 2>&1
python3 - "$run_root/smoke-mms" <<'CHECK'
import csv,json,sys
from pathlib import Path
p=Path(sys.argv[1]);r=list(csv.DictReader((p/'history.csv').open()))[1:]
assert json.loads((p/'summary.json').read_text())['completed']
assert max(abs(float(x['energy_balance_residual'])) for x in r)<1e-9
assert max(abs(float(x['weak_divergence_inf'])) for x in r)<1e-10
print('Manufactured smoke validation passed')
CHECK
[ -f "$run_root/smoke-contraction/summary.json" ] || scripts/run.sh configs/smoke/contraction.prm "$run_root/smoke-contraction" > "$run_root/smoke-contraction.log" 2>&1
phase 'Measuring full-size memory and CPU throughput'
python3 scripts/calibrate_server.py "$run_root/calibration" > "$run_root/calibration.log" 2>&1
source "$run_root/calibration/settings.sh"
export OPENBLAS_NUM_THREADS=$CONVERGENCE_THREADS OMP_NUM_THREADS=$CONVERGENCE_THREADS
phase 'Running 42 convergence cases'
python3 scripts/run_suite.py --suite convergence --jobs "${CONVERGENCE_JOBS:-2}" --output "$run_root/convergence" > "$run_root/convergence.log" 2>&1 || echo 'Some convergence cases failed; retained their logs.'
export OPENBLAS_NUM_THREADS=$CONTRACTION_THREADS OMP_NUM_THREADS=$CONTRACTION_THREADS
phase 'Running 24 contraction cases'
python3 scripts/run_suite.py --suite contraction --jobs "${CONTRACTION_JOBS:-1}" --output "$run_root/contraction" > "$run_root/contraction.log" 2>&1 || echo 'Some contraction cases failed; retained their logs.'
phase 'Preparing Python analysis dependencies'
python3 -m pip --version >/dev/null 2>&1 || { curl --fail --location --retry 3 https://bootstrap.pypa.io/get-pip.py -o "$HOME/opt/get-pip.py"; python3 "$HOME/opt/get-pip.py" --user; }
python3 -m pip install --user 'matplotlib>=3.8,<4' 'scipy>=1.10,<2' > "$run_root/python-install.log" 2>&1
phase 'Generating figures and completion audit'
python3 scripts/plot_suite.py "$run_root" > "$run_root/plotting.log" 2>&1
if python3 -c 'import json,sys;sys.exit(not all(x["status"]=="completed" for x in json.load(open(sys.argv[1]))))' "$run_root/audit.json"; then
 phase 'FINISHED: all 66 cases completed; figures ready for scientific review'
else
 phase 'FINISHED WITH INCOMPLETE CASES: inspect REPORT.md and audit.json'
fi
