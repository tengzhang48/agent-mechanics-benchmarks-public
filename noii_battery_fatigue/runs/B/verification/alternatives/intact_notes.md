2026-10-20 19:49:47 UTC: Created isolated intact-fraction candidate by copy of main source; no main or tool edits. Sixth Q1 field vfrac satisfies L2[(d-1)+vfrac]=0. Degradation, W and H use vfrac accurately below ulp(d); damage reaction is d-(1-k)*vfrac*H.

2026-10-20 19:51:51 UTC: Candidate tests complete. All element source files frozen; no main element files changed by this agent.
Generic material verification needed cancellation-aware states for the intact-fraction candidate: d derivative1 cannot be resolved by default FD when H*vfrac is1e5 or1e18; q subtraction also hides tiny derivative. Used frozen vfrac1e-5,q0 and fatiguefloor vfrac1e-18,q0; generation precheck uses verified nonzero active state. No tolerances changed or package edits.
D0 n8 spc40 onecycle: ok, 2.621606s, 0cutbacks. Floor n8 spc80 twocycles: ok, 5.224881s, 0cutbacks.
Comparison and floor regime details:
```json
{
  "D0_final_primary_q_max_absolute_difference": {
    "ux": 1.8188575645616822e-15,
    "uy": 1.27675647831893e-15,
    "c": 5.029310301551959e-14,
    "d": 2.799982468104645e-13,
    "q": 3.51745000314041e-15
  },
  "D0_all_output_observables_max_absolute_difference": {
    "d_max": 6.376577014349273e-12,
    "d_min": 2.928144854949348e-13,
    "D_int": 9.36855424509865e-13,
    "C_int": 9.98590099499097e-12,
    "E_el": 2.087758205714029e-17,
    "A_damage": 0.0
  },
  "D0_final_history_max_absolute_difference": 3.0238034298690764e-12,
  "floor_cycles": [
    {
      "file": "fields_cycle_001.npz",
      "d_min": 1.0,
      "d_max": 1.0000000000000002,
      "H_min": 9.350649350649228e+17,
      "H_max": 9.350649350649437e+17,
      "alpha_min": 1.4420968814036077e-07,
      "alpha_max": 1.4420968814047374e-07,
      "fatigue_min": 1e-15,
      "aux_min": 1.0694551389958249e-18,
      "aux_max": 1.0694551389958465e-18,
      "vfrac_constraint_max": 2.2311406006402716e-16
    },
    {
      "file": "fields_cycle_002.npz",
      "d_min": 1.0,
      "d_max": 1.0000000000000002,
      "H_min": 9.350649350649228e+17,
      "H_max": 9.350649350649437e+17,
      "alpha_min": 2.805733245039964e-07,
      "alpha_max": 2.8057332450411146e-07,
      "fatigue_min": 1e-15,
      "aux_min": 1.0694551389958247e-18,
      "aux_max": 1.0694551389958465e-18,
      "vfrac_constraint_max": 2.2311406006402716e-16
    }
  ],
  "largest_direct_FD_row_relative_error": 1.4786587213524347e-09
}
```
Source SHA256 at freeze:
```json
{
  "element/generate.py": "a5eddf55ba3917f60f46d8a6e2c78ddb772a3e400f34561f69e7813dfcbcb44c",
  "element/generated_model.for": "744b5303512e26f261cd5d3158f85e2d8492fb20712bc2a901cb66062812e1ca",
  "element/model.py": "aa5f851698b7e14d0eaadfbcb20c840306544724ded7c73dcbf961bec11401f5",
  "element/uel.f90": "5608c8db5306f92725c8192fd4f23dfc9d0ab8f459577d4de37d48b48a05dcd9",
  "element/uel_contract.json": "f11dbe86781ca6096d80491d9f265472f11d51b1772b3fb78d3370623c979736"
}
```
Unverified: broader mesh/time convergence and hidden organizer problems. These use unmodified supplied host, not Abaqus.
