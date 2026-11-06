# Test 1 score: B

Problem D0, 10 cycles. Tolerances: `{"obs": {"d_max": 0.04, "D_int": 0.04, "C_int": 0.005}, "field_rms": 0.01, "drop": 1.5, "converged": 0.005, "history": -1e-12, "d_range": [-0.01, 1.01], "d_reached": 0.9}`

| check | pass | detail |
|---|---|---|
| C1_key_fields | **no** | ux 2.96e-03, uy 6.48e-03, c 3.05e-03, d 6.21e-02 |
| C1_key_observables | **no** | d_max 1.82e-01, D_int 1.37e-01, C_int 3.29e-03, E_el 1.75e-01 |
| C2_converges_to_key | yes | d_max 5.5e-01->1.8e-01, D_int 5.1e-01->1.4e-01, C_int 1.3e-02->3.3e-03, field_ux 1.1e-02->3.0e-03, field_uy 2.5e-02->6.5e-03, field_c 1.2e-02->3.1e-03, field_d 2.3e-01->6.2e-02 |
| C3_distortion | yes | dist0: d_max 2.46e-02, D_int 2.30e-02, C_int 5.48e-04, E_el 2.65e-02; dist1: d_max 1.85e-02, D_int 2.00e-02, C_int 5.19e-04, E_el 2.44e-02 |
| R0_runs_complete | yes |  |
| R1_history | yes | declared True, worst change 0.0 |
| R2_hidden_state | yes | test=hidden_state, max_abs_dU=0.0, max_abs_dSVARS=0.0 |
| R3_dtime_zero | yes | test=dtime_zero, element=63, finite=True, pnewdt=0.5 |
| R4_large_step | yes | complete=True, cutbacks=0 |
| R5_damage_near_1 | not run | not run: no stress variant given |

Info (not scored):

```
{
 "E_el_vs_key_n32": 0.17466755921572227,
 "same_dt_T40_vs_key_n32": {
  "d_max": 0.1800011878189059,
  "D_int": 0.13522673771814026,
  "C_int": 0.003087051960022719,
  "E_el": 0.16747927710645238
 }
}
```
