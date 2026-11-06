2026-10-20 19:48:35 UTC: Created isolated candidate with copied model/generator/wrapper/contract. Added Q1 L2 log-history projection and exact positive nodal damage-row scaling including its z derivative. No main element, host or tool modifications. No git metadata present.

2026-10-20 19:51:51 UTC: Candidate tests complete. All element source files frozen; no main element files changed by this agent.
D0 n8 spc40 onecycle: ok, 2.959253s, 0cutbacks. Floor n8 spc80 twocycles: ok, 22.295650s, 8cutbacks.
Comparison and floor regime details:
```json
{
  "D0_final_primary_q_max_absolute_difference": {
    "ux": 2.6354749719376125e-10,
    "uy": 1.5440037356290182e-10,
    "c": 1.610790723027833e-08,
    "d": 3.436497501008162e-07,
    "q": 1.9630907164593092e-09
  },
  "D0_all_output_observables_max_absolute_difference": {
    "d_max": 3.4623240074971307e-07,
    "d_min": 8.258518605086529e-09,
    "D_int": 4.4687377936947215e-08,
    "C_int": 4.167255540643566e-09,
    "E_el": 1.471854025987475e-13,
    "A_damage": 0.0
  },
  "D0_final_history_max_absolute_difference": 3.3060614332836025e-06,
  "floor_cycles": [
    {
      "file": "fields_cycle_001.npz",
      "d_min": 1.0,
      "d_max": 1.0000000000000002,
      "H_min": 9.350649350649228e+17,
      "H_max": 9.350649350649437e+17,
      "alpha_min": 1.4420968813471826e-07,
      "alpha_max": 1.4420968813483123e-07,
      "fatigue_min": 1e-15,
      "aux_min": 41.37938237100518,
      "aux_max": 41.3793823710052
    },
    {
      "file": "fields_cycle_002.npz",
      "d_min": 1.0,
      "d_max": 1.0000000000000002,
      "H_min": 9.350649350649228e+17,
      "H_max": 9.350649350649437e+17,
      "alpha_min": 2.8057332449835386e-07,
      "alpha_max": 2.8057332449846895e-07,
      "fatigue_min": 1e-15,
      "aux_min": 41.379382371005185,
      "aux_max": 41.3793823710052
    }
  ],
  "largest_direct_FD_row_relative_error": 1.4474444980082133e-08
}
```
Source SHA256 at freeze:
```json
{
  "element/generate.py": "d886bc5fe3eb8ca546fe9e4655b23e4f2ad936a5c575886a14561de609beb01f",
  "element/generated_model.for": "6420fc8c8097ef1202b7d738a52caa22cf0712cce702a25e6f15350661a11db4",
  "element/model.py": "14553f592cd2d6e6f4c3405054e3f72410a2b3aa36d18b300deae36635dfe756",
  "element/uel.f90": "e1c0c9eaf8d870d5addc2d6ef1a5c16d2497f1ac5f6c0ab7a00a77dfb8466e08",
  "element/uel_contract.json": "03914279814bb409476095f51b2b3eaa4f4d60ac6842842595aef5c012ab5f3a"
}
```
Unverified: broader mesh/time convergence and hidden organizer problems. These use unmodified supplied host, not Abaqus.
