"""Explore a progressively more damaging loading without prescribed damage."""
import json
from common import ROOT,write_json
from global_checks import run

for gc in [3.5e-7,3.5e-8,3.5e-9]:
    prob=json.loads((ROOT/'host'/'dev_problem'/'d0.json').read_text())
    prob['id']=f'D0_Gc_{gc:g}'
    prob['parameters']['G_c']=gc
    write_json(ROOT/'verification'/'damage_problems'/f'{gc:g}.json',prob)
    run(prob['id'],n=8,spc=40,cycles=5,distorted=True,problem=prob)
