#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
: "${NVCC:=nvcc}"
: "${CUDA_ARCH:=80}"
case "${1:-}" in
  RE) source=re_native_parallel.cu; binary=re ;;
  DE) source=de_native_simulation_parallel.cu; binary=de ;;
  *) echo "Usage: bash src/cuda/build.sh RE|DE" >&2; exit 2 ;;
esac
mkdir -p ../../build/cuda
"$NVCC" -O3 -std=c++17 --default-stream per-thread --fmad=false \
  --prec-div=true --prec-sqrt=true -Xcompiler=-pthread \
  -Xcompiler=-ffp-contract=off -Xcompiler=-fno-fast-math \
  -Xcompiler=-march=native -arch="sm_${CUDA_ARCH}" "$source" \
  -o "../../build/cuda/$binary"
