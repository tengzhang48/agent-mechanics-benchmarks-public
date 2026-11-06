"""Controlled alternative to assess whether a temporal lag improves robustness."""
import json
from common import ROOT
from global_checks import run

for filename in ['damage_problems/3.5e-08.json','saturation_problem.json','extreme_problem.json']:
    problem=json.loads((ROOT/'verification'/filename).read_text())
    name='lag_'+problem['id']
    run(name,n=8,spc=40,cycles=5,distorted=True,problem=problem,lag=True)
