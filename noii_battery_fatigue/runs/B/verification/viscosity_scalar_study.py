"""Independent accuracy study of bounded algorithmic damage relaxation."""
import argparse
import json
from pathlib import Path

import numpy as np

from uniform_time_reference import solve


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--form',default='history',choices=['history','accepted_stiffness','phase_stiffness','sqrt_dt_phase_stiffness'])
    args=ap.parse_args()
    reference = solve(163840)['final_alpha']
    records = []
    coefficients=(0., .1, .3, 1.) if args.form=='sqrt_dt_phase_stiffness' else (0., .1, .5, 1., 2.)
    for eta in coefficients:
        for spc in (20, 80, 320, 1280, 10240):
            result = solve(spc, implicit='psilag', predictor_floor=.01, relaxation=eta,
                           relaxation_form=args.form)
            data = result.pop('rows')
            result['alpha_relative_difference'] = result['final_alpha']/reference-1
            work = data[data[:, 0] < .8, 3]
            direction = np.sign(np.diff(work))
            result['W_direction_reversals_before_t08'] = int(np.sum(direction[:-1]*direction[1:] < 0))
            records.append(result)
            print(json.dumps(result), flush=True)
    name='viscosity_study' if args.form=='history' else 'viscosity_'+args.form
    out = Path(__file__).resolve().parent/f'uniform_time_reference/{name}.json'
    out.write_text(json.dumps(dict(reference_alpha=reference, results=records), indent=2)+'\n')


if __name__ == '__main__':
    main()
