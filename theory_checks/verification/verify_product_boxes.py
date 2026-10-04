#!/usr/bin/env python3
"""Verify every corner of two five-factor prevalence boxes by exact arithmetic.
Regression coefficients for a product binary law are multi-affine in the
marginal probabilities. Every interior coefficient vector is consequently a
convex combination of the corner vectors; the base functional is concave.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json
from verify_range_certificate import product_probabilities, certificate

BOXES={
 'main_model_box': [('0.499','0.501'),('0.399','0.401'),('0.299','0.301'),('0.199','0.201'),('0.099','0.101')],
 'production_trial_box': [('0.465','0.467'),('0.591','0.593'),('0.286','0.288'),('0.153','0.155'),('0.299','0.301')]
}
results={}
for name,box in BOXES.items():
    all_results=[]
    for i,corner in enumerate(product(*box)):
        r=certificate(product_probabilities(corner))
        assert r['certificate_pass'],(name,corner,r)
        assert r['all_zero_rays_structural_constant_sign']
        r['marginals']=list(corner);all_results.append(r)
    c0=min(Q(r['minimum_base']['value']) for r in all_results)
    ray=min(Q(r['minimum_ray']['value']) for r in all_results)
    c1=min(Q(r['minimum_zero_overall_slope']['value']) for r in all_results)
    results[name]={'box':box,'corners_verified':len(all_results),'minimum_c0':str(c0),
                   'minimum_ray':str(ray),'minimum_c1':str(c1),
                   'minimum_c1_decimal':float(c1),'all_corners_pass':True,'corners':all_results}
    print(name,'corners=',len(all_results),'c0=',c0,'ray=',ray,'c1=',c1,'decimal=',float(c1),flush=True)
Path(__file__).with_name('product_box_certificates.json').write_text(json.dumps(results,indent=2)+'\n')
