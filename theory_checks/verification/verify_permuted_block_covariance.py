#!/usr/bin/env python3
"""Exact finite-block checks for the added stratified permuted-block theorem.

Only Python's standard library is needed. The checks enumerate all balanced
sign vectors, then all independent pairs/triples for block sizes 2, 4 and 6.
They verify the covariance identities used in the proof, not an asymptotic
claim by simulation. Run this file from any working directory.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
import json


def balanced_vectors(length: int) -> list[tuple[int, ...]]:
    if length < 2 or length % 2:
        raise ValueError('Block length must be an even integer >= 2.')
    out = []
    for plus in combinations(range(length), length // 2):
        indices = set(plus)
        out.append(tuple(1 if i in indices else -1 for i in range(length)))
    return out


def check(length: int) -> dict:
    vectors = balanced_vectors(length)
    count = len(vectors)
    for i in range(length):
        assert sum(v[i] for v in vectors) == 0
        for j in range(i + 1, length):
            assert Q(sum(v[i] * v[j] for v in vectors), count) == -Q(1, length-1)
    totals = [0] * 7
    for a, b, c in product(vectors, repeat=3):
        ab = sum(x*y for x, y in zip(a, b))
        ac = sum(x*y for x, y in zip(a, c))
        bc = sum(x*y for x, y in zip(b, c))
        abc = sum(x*y*z for x, y, z in zip(a, b, c))
        increments = (ab, ab*ab, ab*ac, ab*bc, abc, abc*abc, ab*abc)
        totals = [x+y for x, y in zip(totals, increments)]
    means = [Q(x, count**3) for x in totals]
    expected_pair = Q(length**2, length-1)
    expected_triple = Q(length) - Q(length, (length-1)**2)
    assert means == [Q(0), expected_pair, Q(0), Q(0), Q(0), expected_triple, Q(0)]
    return {
        'block_length': length,
        'balanced_vectors': count,
        'independent_triples_enumerated': count**3,
        'pair_reward_variance': str(expected_pair),
        'pair_variance_per_observation': str(expected_pair / length),
        'triple_reward_variance': str(expected_triple),
        'distinct_pair_cross_moments': '0',
        'pair_triple_cross_moment': '0',
        'all_exact_checks_passed': True,
    }


def main() -> None:
    result = {'method': 'exhaustive enumeration with exact rational arithmetic',
              'cases': [check(length) for length in (2, 4, 6)]}
    output = Path(__file__).with_name('permuted_block_exact_checks.json')
    output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
