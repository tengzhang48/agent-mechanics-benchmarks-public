# Test 1 score: B2

Problem D0, 10 cycles. Tolerances: `{"obs": {"d_max": 0.04, "D_int": 0.04, "C_int": 0.005}, "field_rms": 0.01, "drop": 1.5, "converged": 0.005, "history": -1e-12, "d_range": [-0.01, 1.01], "d_reached": 0.9}`

| check | pass | detail |
|---|---|---|
| C1_key_fields | yes | ux 2.06e-04, uy 1.49e-04, c 1.34e-04, d 2.80e-03 |
| C1_key_observables | yes | d_max 9.02e-03, D_int 6.98e-03, C_int 8.84e-05, E_el 1.33e-01 |
| C2_converges_to_key | yes | d_max 2.2e-02->9.0e-03, D_int 2.6e-02->7.0e-03, C_int 1.9e-04->8.8e-05, field_ux 7.7e-04->2.1e-04, field_uy 3.5e-04->1.5e-04, field_c 4.5e-04->1.3e-04, field_d 1.0e-02->2.8e-03 |
| C3_distortion | yes | dist0: d_max 1.19e-02, D_int 1.34e-02, C_int 3.05e-04, E_el 2.82e-02; dist1: d_max 1.50e-02, D_int 1.30e-02, C_int 3.02e-04, E_el 2.67e-02 |
| R0_runs_complete | yes |  |
| R1_history | yes | declared True, worst change 0.0 |
| R2_hidden_state | yes | test=hidden_state, max_abs_dU=0.0, max_abs_dSVARS=0.0 |
| R3_dtime_zero | yes | test=dtime_zero, element=60, finite=True, pnewdt=0.5 |
| R4_large_step | yes | complete=True, cutbacks=0 |
| R5_damage_near_1 | not run | not run: no stress variant given |

Info (not scored):

```
{
 "E_el_vs_key_n32": 0.13343230839571865,
 "same_dt_T40_vs_key_n32": {
  "d_max": 0.006579228833001982,
  "D_int": 0.006308610706395438,
  "C_int": 5.7422852610749824e-05,
  "E_el": 0.12882215694600044
 }
}
```
