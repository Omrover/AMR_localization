#!/usr/bin/env bash
# Pruebas en host de amr_odom (cinemática omni + EKF).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
LIB="$(dirname "$HERE")"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT
${CC:-gcc} -std=c99 -D_DEFAULT_SOURCE -Wall -Wextra -Werror -Wshadow \
  -fsanitize=address,undefined -g -I"$LIB" "$LIB/omni.c" "$LIB/ekf.c" "$HERE/host_test.c" \
  -lm -o "$OUT/host_test"
"$OUT/host_test"
