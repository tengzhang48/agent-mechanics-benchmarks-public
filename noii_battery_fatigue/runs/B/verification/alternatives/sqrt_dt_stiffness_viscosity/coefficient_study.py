"""Final frozen-source runtime C sweep for eta_effective=C/sqrt(DTIME)."""
import argparse,json
from common import ROOT,compiled
from run_case import run

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--coefficient',type=float,nargs='+',default=[.1,.3,1.]);a=ap.parse_args()
 module,contract=compiled();summaries=[]
 for coefficient in a.coefficient:
  cc={**contract,'props':contract['props'].copy()};cc['props'][12]=coefficient
  tag=str(coefficient).replace('.','p')
  out=ROOT/('spatial_floor_C'+tag)
  result=run(ROOT/'problems/notch_spatial_floor.json',mesh='distorted',n=16,seed=3,cycles=2,spc=80,out=out,module_contract=(module,cc))
  summaries.append(dict(C=coefficient,**result))
 (ROOT/'verification'/('spatial_C'+tag+'.json')).write_text(json.dumps(summaries,indent=1)+'\n')
if __name__=='__main__':main()
