"""Frozen-source runtime eta sweep, preserving all failed and accepted states."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from common import ROOT,compiled
from run_case import run

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--eta',type=float,nargs='+',default=[.1,.5,1.,2.]);ap.add_argument('--restart',action='store_true');a=ap.parse_args()
 module,contract=compiled();summaries=[]
 source=ROOT.parent/'lumped_psilag_candidate/spatial_floor_distorted16_s3_c2_spc80'
 state=np.load(source/'last_state.npz');last=list(csv.DictReader((source/'increments.csv').open()))[-1]
 for eta in a.eta:
  cc={**contract,'props':contract['props'].copy()};cc['props'][12]=eta
  tag=str(eta).replace('.','p')
  kwargs={}
  if a.restart:
   kwargs['restart']=dict(U=state['U'],svars=state['svars'],t_start=float(last['t']),duration=float(last['dt']),steps=1)
  out=ROOT/(('restart_' if a.restart else 'spatial_floor_')+'eta'+tag)
  result=run(ROOT/'problems/notch_spatial_floor.json',mesh='distorted',n=16,seed=3,cycles=2,spc=80,out=out,module_contract=(module,cc),**kwargs)
  summaries.append(dict(eta=eta,**result))
 (ROOT/'verification'/('restart_eta_sweep.json' if a.restart else 'spatial_eta_sweep.json')).write_text(json.dumps(summaries,indent=1)+'\n')
if __name__=='__main__':main()
