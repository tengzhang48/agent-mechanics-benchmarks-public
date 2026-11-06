2026-10-20 20:18:34 UTC: Isolated active-tie production snapshot and lumped-reaction candidate created. Material phase_storage returns H carrier only; deterministic output extension replaces phase reaction and its complete chain rule. All other physics and history remain unchanged. No production, host or tool changes.

2026-10-20 20:21:26 UTC: Deterministic regeneration reproduced Fortran SHA256 c01af3958986091358dcb17fc061dc2f9db929f2949671e7898e63447fef4d0b. Generation bookkeeping separates unmodified-tool output, 5 assembly lines rewritten as10lines, added nodal-diagonal loop, and handwritten checked wrapper. Actual generated WeakForm Python phase_storage is an H carrier; only compiled postprocessed element is the lumped physical phase residual. problem.verify alone does not test lumping; independent24x24FD+assembly does.
```json
{
  "generator": "abaqus_ufl",
  "pinned_commit": "a02ed7c",
  "package_modifications": [],
  "output_transformations": [
    "Rename UEL entry to CHEM_UEL_CORE for checked wrapper",
    "Lump phase reaction to nodal d/vfrac with GP H; multiply generated H chain-rule tangent by -a*vfrac_node(row) and add nodal diagonals"
  ],
  "python_tangent_checks": {
    "active_undiminished": true,
    "active_fatigue": true,
    "unloading_frozen": true,
    "fatigue_floor": true,
    "compression": true
  },
  "sha256": "c01af3958986091358dcb17fc061dc2f9db929f2949671e7898e63447fef4d0b",
  "generated_lines_before_output_extension": 3985,
  "postprocessed_total_lines": 4000,
  "assembly_extension_original_lines_replaced": 5,
  "assembly_extension_replacement_lines": 10,
  "assembly_extension_added_lines": 10,
  "handwritten_wrapper_lines": 89
}
```

2026-10-20 20:26:12 UTC: Prepared candidate-local spatial fatigue-floor case fromnotch_sides withalphaT=1e-10, unchangedkfg=.5/Gc1.5e-6/floor1e-15, two80-stepcycles, free coupledmechanics/transport on n16distortedseed3. Execution willfollowpendingnotch32result. No productionproblemfile changes.
