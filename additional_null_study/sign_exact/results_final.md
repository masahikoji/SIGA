# Values used in the manuscript (Table 3, Supplementary Tables 11-12, Proposition 4, Corollary 4)

Restricted ratio r_n = v'Var(G_n)v / (n v'Pi v), G_n = sum_i h(S_i) Z_i^(1) Z_i^(2), on pi'v = 0.
Deterministic forward recursion (sparsedp.py); "discarded" = probability mass of pruned states;
error bound for every entry: discarded x n / min_s pi_s.  Reproduce with reproduce_sign_tables.py.

## One factor (tolerance 1e-16; discarded <= 1.4e-10; error < 1e-7)
| (pi+, p) | n=10 | n=20 | n=40 | n=60 | n=80 | 1/n extrap. |
|---|---|---|---|---|---|---|
| (.5,.6)  | 1.0208 | 1.0208 | 1.0188 | 1.0168 | 1.0158 | 1.013 |
| (.5,.7)  | 1.0634 | 1.0585 | 1.0501 | 1.0455 | 1.0430 | 1.036 |
| (.5,.8)  | 1.0977 | 1.0899 | 1.0815 | 1.0781 | 1.0764 | 1.071 |
| (.3,.8)  | 1.0996 | 1.0956 | 1.0855 | 1.0803 | 1.0778 | 1.070 |
| (.1,.8)  | 1.0762 | 1.1018 | 1.1058 | 1.0994 | 1.0936 | 1.081 |
| (.5,.95) | 1.0597 | 1.0595 | 1.0593 | 1.0592 | 1.0592 | 1.059 |
| (.3,.95) | 1.0604 | 1.0583 | 1.0558 | 1.0549 | 1.0544 | 1.053 |
| (.1,.95) | 1.0505 | 1.0598 | 1.0541 | 1.0485 | 1.0450 | 1.036 |
| (.5,.99) | 1.0146 | 1.0148 | 1.0150 | 1.0150 | 1.0150 | 1.015 |
| (.5,.999)| 1.0015 | 1.0016 | 1.0016 | 1.0016 | 1.0016 | 1.0016 |
| any, p=1 | 1 exactly for every n (Corollary 4; exact_p1_F1.py: Psi_D = Pi + (1+q-q^2)/(q(1-q)) pi pi') |

## Two factors
(.5,.4), p=.8 (tol 1e-13; discarded 4.5e-7 at n=40): n=10: 1.0586 1.0836 1.0837; n=20: 1.0604 1.0790 1.0796;
  n=30: 1.0573 1.0743 1.0746; n=40: 1.0544 1.0710 1.0710
(.5,.4), p=.95 (tol 1e-14; discarded 3.9e-9): n=10: 1.0594 1.0622 1.1894; n=20: 1.0640 1.0666 1.2050; n=30: 1.0653 1.0674 1.2054
uniform, p=.8, n=30 (tol 1e-13): Walsh |A|=1: 1.0734 (twice), |A|=2: 1.0594

## Three factors
uniform, p=.8: |A|=3 for n=3..10: .9904 .9960 .9949 .9964 .9962 .9969 .9969 .9974 (exact for n<=9; n=10 discarded 7.3e-6)
               |A|=2 for n=3,4,6,8,10: 1.0192 1.0160 1.0285 1.0356 1.0398 ; |A|=1: 1.0480 1.0514 1.0643 1.0683 1.0694
(.5,.4,.3), p=.8 (tol 1e-12): smallest eigenvalue n=3,4,6,8: .9921 .9964 .9966 .9970

## n = 3 closed forms (uniform law, full interaction; exact_bruteforce.py)
F=1: 1 + 4/3 pq eta^2 ; F=2: 1 + 1/6 p^2 eta^2 ; F=3: 1 - 1/6 pq eta^2 (= 619/625 at p = 4/5) ; F=4: 1 - 1/32 p^2 eta^2
