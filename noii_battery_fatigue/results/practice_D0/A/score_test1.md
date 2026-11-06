# Test 1 score: A

Problem D0, 10 cycles. Tolerances: `{"obs": {"d_max": 0.04, "D_int": 0.04, "C_int": 0.005}, "field_rms": 0.01, "drop": 1.5, "converged": 0.005, "history": -1e-12, "d_range": [-0.01, 1.01], "d_reached": 0.9}`

| check | pass | detail |
|---|---|---|
| C1_key_fields | **no** | ux 2.33e-03, uy 4.93e-03, c 2.32e-03, d 4.97e-02 |
| C1_key_observables | **no** | d_max 1.15e-01, D_int 1.14e-01, C_int 2.47e-03, E_el 1.63e-01 |
| C2_converges_to_key | yes | d_max 3.6e-01->1.2e-01, D_int 4.3e-01->1.1e-01, C_int 9.9e-03->2.5e-03, field_ux 8.8e-03->2.3e-03, field_uy 2.0e-02->4.9e-03, field_c 9.1e-03->2.3e-03, field_d 1.9e-01->5.0e-02 |
| C3_distortion | yes | dist0: d_max 2.44e-02, D_int 2.22e-02, C_int 5.15e-04, E_el 2.65e-02; dist1: d_max 1.60e-02, D_int 1.89e-02, C_int 4.68e-04, E_el 2.41e-02 |
| R0_runs_complete | yes |  |
| R1_history | yes | declared True, worst change 0.0 |
| R2_hidden_state | yes | test=hidden_state, max_abs_dU=0.0, max_abs_dSVARS=0.0 |
| R3_dtime_zero | yes | test=dtime_zero, element=63, finite=True, pnewdt=0.5 |
| R4_large_step | yes | complete=True, cutbacks=0 |
| R5_damage_near_1 | not run | not run: no stress variant given |

Info (not scored):

```
{
 "E_el_vs_key_n32": 0.1633660716657926,
 "same_dt_T40_vs_key_n32": {
  "d_max": 0.13263993611454508,
  "D_int": 0.11751637892285405,
  "C_int": 0.002502213141862583,
  "E_el": 0.15978482149287007
 }
}
```
