#!/bin/sh
set -eu
cd "$(dirname "$0")"
python audit_certificate_independent.py range_certificate_results.json.gz --output independent_audit.json
python verify_correlated_witnesses.py
python verify_product_boxes.py
python audit_product_boxes_independent.py
python verify_production_bounds.py
python verify_permuted_block_covariance.py
python verify_sampling_variance_relation.py

python verify_new_results.py
