"""Refine the independent uniform scalar problem for energy-driver splitting."""
import json
from pathlib import Path

import numpy as np

from uniform_time_reference import solve


def main():
    reference = solve(163840)['final_alpha']
    rows = []
    for theta in (0., .001, .01, .1):
        for spc in (20, 40, 80, 160, 320, 1280, 10240, 40960):
            result = solve(spc, implicit='psilag', predictor_floor=theta)
            data = result.pop('rows')
            result['alpha_rel_error'] = result['final_alpha']/reference-1
            work = data[data[:, 0] < .8, 3]
            direction = np.sign(np.diff(work))
            result['W_direction_reversals_before_t08'] = int(np.sum(direction[:-1]*direction[1:] < 0))
            rows.append(result)
            print(json.dumps(result), flush=True)
    out = Path(__file__).resolve().parent/'uniform_time_reference/psilag_study.json'
    out.write_text(json.dumps(dict(reference_alpha=reference, results=rows), indent=2)+'\n')


if __name__ == '__main__':
    main()
