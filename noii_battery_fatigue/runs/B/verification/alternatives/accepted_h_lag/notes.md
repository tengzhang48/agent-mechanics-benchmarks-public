2026-10-20 20:27:04 UTC: Createdisolatedlumped+accepted-H-lagcandidateafterimplicitlumpednotch32failedstep9. OnlyphaseHcarrierchangesfromcurrentHtohhist_old; allnodalreactionquadrature, gradient, physicsandend-statehistoryformulasunchanged. No production/host/toolchanges.

2026-10-20 20:28:41 UTC: Basiclagspatialfloor n16distseed3alphaT1e-10c2 passed160steps0cutbacks28.16s; fullspatialfatigue/Hextents verification/spatial_floor_extents.json. Parenttestscontrolleralternativesbeforechoosingcontrolledspatialfloorvariant.

2026-10-20 20:30:45 UTC: Basiclump+lag notch32distseed3firstcycle80spc passed80steps0cutbacks100.26s,max3observedlinearsolves/attempt. Independent24x24JacobianFD maxrowrelative1.2385e-9; independentstateerror2.6e-14; bitwisehiddenmemorydiscardchecksandstate-preservingDTIME0passed. RawtoolFortran3969lines,3984afterpostgenextension,89linecheckedwrapper.
```json
{
  "d_min": 0.16443779082141516,
  "d_max": 1.0,
  "v_min": 3.372440674686497e-17,
  "v_max": 0.835562209178585,
  "H_min": 0.014231502157650617,
  "H_max": 220488.01797785456,
  "alpha_min": 9.01494058669591e-08,
  "alpha_max": 6.46768543163558e-05,
  "c_min": 0.35,
  "c_max": 0.8179107684837159,
  "all_increment_d_max": 1.0,
  "intact_constraint_max": 5.551115123125783e-16
}
```
Formal accepted-history lag is first-order on smooth regular paths. Fatigue-jump onset and irreversible alpha overshoot remain accuracy risks; controller variants are under separate parent investigation. No claims of Abaqus execution or hidden-case success.
