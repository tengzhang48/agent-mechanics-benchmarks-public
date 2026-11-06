"""Compare declared PNEWDT accuracy controls for accepted-history coupling.

Sources and compiler work are isolated under tmp/.  Every case uses the
unchanged host solver and retains the physical equations/history update.
"""
from pathlib import Path
import argparse
import hashlib
import importlib
import json
import shutil
import sys

import numpy as np

from common import ROOT, fehost
from run_case import run


CONTROL = '''    ! Declared accuracy requests for the accepted-history phase split.
    pred_gap_old = 1.0d0/(1.0d0+(1.0d0-PROPS(5))*SVARS(i))
    pred_gap_new = 1.0d0/(1.0d0+(1.0d0-PROPS(5))*trial_state(i))
    if (PROPS(12) > 0.0d0) then
      if (pred_gap_old-pred_gap_new > PROPS(12)) goto 900
    end if
    if (PROPS(11) > 0.0d0 .and. (PROPS(13) > 0.0d0 .or. PROPS(14) > 0.0d0)) then
      if (PROPS(11)*(log(max(SVARS(i+1),PROPS(10)))-log(PROPS(10))) < log(10.0d0)) then
        alpha_limit = max(PROPS(13)*PROPS(10),PROPS(14)*SVARS(i+1))
        if (trial_state(i+1)-SVARS(i+1) > alpha_limit) goto 900
      end if
    end if
'''


def prepare():
    source = ROOT/'verification/alternatives/accepted_h_lag/element'
    scratch = ROOT/'tmp/time_control_candidate'
    (scratch/'element').mkdir(parents=True, exist_ok=True)
    for name in ['generated_model.for', 'uel.f90', 'model.py', 'lump_reaction.py', 'generate.py', 'uel_contract.json']:
        shutil.copyfile(source/name, scratch/'element'/name)
    wrapper = (source/'uel.f90').read_text()
    wrapper = wrapper.replace('geom_scale, p_trial', 'geom_scale, p_trial, pred_gap_old, pred_gap_new, alpha_limit')
    wrapper = wrapper.replace('NPROPS /= 11', 'NPROPS /= 14')
    anchor = '    if (trial_state(i) < SVARS(i) .or. trial_state(i+1) < SVARS(i+1)) goto 900\n'
    assert wrapper.count(anchor) == 1
    wrapper = wrapper.replace(anchor, anchor+CONTROL)
    (scratch/'element/uel.f90').write_text(wrapper)
    contract = fehost.load_contract(str(scratch))
    contract['prop_names'] += ['history_damage_limit', 'alpha_absolute_limit', 'alpha_relative_limit']
    contract['props'] += [.05, 1., .1]
    (scratch/'element/uel_contract.json').write_text(json.dumps(contract, indent=2)+'\n')
    manifest = {f'element/{p.name}':hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (scratch/'element').iterdir() if p.is_file()}
    digest = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()[:12]
    name = 'uel_time_control_'+digest
    build = scratch/'tmp'/digest
    if list(build.glob(name+'*.so')):
        sys.path.insert(0, str(build))
        module = importlib.import_module(name)
    else:
        module = fehost.build_module(str(scratch), contract, str(build), name)
    return module, contract, manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variant', nargs='+', default=['H_only', 'absolute1', 'relative01', 'loose', 'gap02'])
    parser.add_argument('--spc', nargs='+', type=int, default=[20, 80])
    parser.add_argument('--cycles', type=int, default=2)
    args = parser.parse_args()
    module, contract, manifest = prepare()
    variants = dict(H_only=[.05, 0., 0.], absolute1=[.05, 1., 0.],
                    relative01=[.05, 1., .1], loose=[.05, 5., .2], gap02=[.02, 1., .1],
                    off=[0., 0., 0.])
    out = ROOT/'verification/time_control_study'
    out.mkdir(exist_ok=True)
    rows = []
    for label in args.variant:
        controls = variants[label]
        for spc in args.spc:
            cc = {**contract, 'props':contract['props'][:11]+controls}
            path = out/f'{label}_spc{spc}'
            summary = run(ROOT/'verification/problems/uniform_fatigue_floor.json', n=8,
                          cycles=args.cycles, spc=spc, out=path, uniform_strain=.03,
                          module_contract=(module, cc), source_manifest=manifest)
            if (path/'increments.csv').exists():
                increments = np.atleast_1d(np.genfromtxt(path/'increments.csv', delimiter=',', names=True))
                floor = increments[increments['H_max']>1e14]
                summary['floor_onset_H_gt_1e14'] = float(floor['t'][0]) if len(floor) else None
                summary['final_alpha'] = float(increments['alpha_max'][-1])
                summary['accepted_subincrements'] = len(increments)
            row = dict(variant=label, controls=controls, spc=spc, **summary)
            rows.append(row)
            (out/'results.json').write_text(json.dumps(rows, indent=2)+'\n')
            print(json.dumps({k:v for k,v in row.items() if k not in ['problem','element_source_sha256']}), flush=True)


if __name__ == '__main__':
    main()
