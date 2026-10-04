#!/bin/sh
set -eu
cd "$(dirname "$0")"
./simulate_sign_products 5 400 200000 0.80 926210031 4 train_n400.bin
./simulate_sign_products 5 400 100000 0.80 926210057 4 eval_n400.bin
./simulate_sign_products 5 2000 200000 0.80 926210083 4 train_n2000.bin
./simulate_sign_products 5 2000 100000 0.80 926210109 4 eval_n2000.bin
./simulate_sign_products 5 400 200000 0.50 926210135 4 train_control.bin
./simulate_sign_products 5 400 100000 0.50 926210161 4 eval_control.bin
