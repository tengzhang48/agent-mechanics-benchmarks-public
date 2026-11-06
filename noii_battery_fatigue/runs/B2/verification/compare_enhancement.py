"""Matched accepted-history D0 comparison: enhanced versus ordinary Q1 strain."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    core=(ROOT/'element/generated_kernel.for').read_text()
    for hook in (
        '      DOUBLE PRECISION :: fcorr(3,3,4), kcorr(24,24)\n',
        '      CALL enhanced_modes(coords_2d,U,PROPS,fcorr,kcorr,PNEWDT)\n      IF (PNEWDT .LT. 1.0D0) RETURN\n\n',
        '        F = F + fcorr(:,:,kk)\n','      AMATRX = AMATRX - kcorr\n\n'):
        assert core.count(hook)==1;core=core.replace(hook,'')
    archive=ROOT/'verification/exploration/basic_split_corrected'
    assert core==(archive/'generated_kernel.for').read_text()
    for name in ('accepted_damage.f90','uel.f90'):
        assert (ROOT/'element'/name).read_bytes()==(archive/name).read_bytes()
    enhanced=ROOT/'verification/runs/D0_enhanced_accepted_H_structured_16_s0_spc40_cy20'
    basic=ROOT/'verification/exploration/basic16_matched40'
    if not all((p/'run.json').exists() for p in (basic,enhanced)):
        print('Matched runs still pending');return
    bs=json.loads((basic/'run.json').read_text());es=json.loads((enhanced/'run.json').read_text())
    assert bs['status']==es['status']=='ok'
    for k in ('parameters','n','mesh','seed','cycles','steps_per_cycle','dt'):
        assert bs[k]==es[k],k
    bc=json.loads((basic/'contract.json').read_text());ec=json.loads((enhanced/'contract.json').read_text())
    for k in ('dof_names','props','prop_names','history_slots','residual_scale','lflags','time_scheme','damage_reaction'):
        assert bc[k]==ec[k],k
    keys=('d_max','D_int','C_int','E_el')
    def load(path):
        return np.atleast_1d(np.genfromtxt(path/'observables.csv',delimiter=',',names=True))
    b,e=load(basic),load(enhanced)
    assert np.array_equal(b['t'],e['t'])
    result=dict(status='passed',scope='n=16, 20 cycles x 40 steps, identical history split/reaction/properties/tolerances',
        core_identical_after_removing_enhancement_hooks=True,reaction_and_guard_byte_identical=True,
        basic=dict(path=str(basic.relative_to(ROOT)),final={k:float(b[k][-1]) for k in keys}),
        enhanced=dict(path=str(enhanced.relative_to(ROOT)),final={k:float(e[k][-1]) for k in keys}))
    ref=ROOT/'verification/reference_runs/D0_Q4_structured_8_spc40_cy20_acceptedH'
    if (ref/'run.json').exists():
        r=json.loads((ref/'run.json').read_text());assert r['status']=='ok'
        data=load(ref);assert np.array_equal(b['t'],data['t'])
        result['independent_Q4_acceptedH_reference']=dict(path=str(ref.relative_to(ROOT)),
            final={k:float(data[k][-1]) for k in keys},note='Same accepted-history time split; independent spatial/constitutive implementation and consistent high-order reaction')
        for label,row in (('basic',b),('enhanced',e)):
            result[label]['final_relative_difference_from_reference']={k:float((row[k][-1]-data[k][-1])/data[k][-1]) for k in keys}
    (ROOT/'verification/enhancement_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
