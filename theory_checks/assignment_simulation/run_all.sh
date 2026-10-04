#!/bin/sh
set -eu
cd "$(dirname "$0")"
CXX=${CXX:-c++}
if "$CXX" -O3 -std=c++17 -fopenmp simulate_sign_products.cpp -o simulate_sign_products 2>compile.log; then
  printf '%s\n' 'Compiled with OpenMP.'
else
  printf '%s\n' 'OpenMP unavailable; compiling the deterministic serial version.'
  "$CXX" -O3 -std=c++17 simulate_sign_products.cpp -o simulate_sign_products 2>>compile.log
fi
./simulate_sign_products --self-test > self_test.log
sh run_new.sh > run.log 2>&1
python3 analyze_split.py > analysis.log 2>&1
printf '%s\n' 'Wrote split_results.csv and directions_and_covariances.json.'
