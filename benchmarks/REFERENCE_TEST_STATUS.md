# Reference Test Status

This table records whether each benchmark is currently a smoke test or a
numerical comparison against a standalone engine reference solution.

- **Engine reference** means the driver output is compared with a committed,
  normalized result produced by the corresponding standalone chemistry engine.
- **Smoke** means CTest verifies only that the driver runs successfully. The
  final column records what remains before a reference comparison can be
  enabled.

The names below are the CTest names reported by `ctest -N`. Tests can also be
selected by classification with `ctest -L smoke` or
`ctest -L engine-reference`.

## Batch chemistry

| CTest name | Engine | Status | Why no engine-reference comparison yet |
|---|---|---|---|
| `batch_chem_calcite_short_ac_pflotran__smoke` | PFLOTRAN | Smoke | The calcite-family time integration and output mapping must be reconciled with the standalone result. |
| `batch_chem_calcite_short_pc_pflotran__smoke` | PFLOTRAN | Smoke | The calcite-family time integration and output mapping must be reconciled with the standalone result. |
| `batch_chem_calcite_vf_pc_pflotran__smoke` | PFLOTRAN | Smoke | The current standalone and Alquimia results have material differences, including concentration differences of about 14%; these must be understood before accepting a reference. |
| `batch_chem_general_reaction_ac_pflotran__engine_reference` | PFLOTRAN | Engine reference | — |
| `batch_chem_general_reaction_pc_pflotran__smoke` | PFLOTRAN | Smoke | The AC variant is covered; this alternate driver configuration still needs its mapping and tolerances validated independently. |
| `batch_chem_ion_exchange_valocchi_ac_pflotran__smoke` | PFLOTRAN | Smoke | The PC/native-condition variant is covered; this alternate driver configuration still needs a separately validated mapping. |
| `batch_chem_ion_exchange_valocchi_pc_pflotran__engine_reference` | PFLOTRAN | Engine reference | — |
| `batch_chem_isotherms_ac_pflotran__engine_reference` | PFLOTRAN | Engine reference | — |
| `batch_chem_isotherms_pc_pflotran__smoke` | PFLOTRAN | Smoke | The AC variant is covered; this alternate driver configuration has not yet been calibrated against the same native result. |
| `batch_chem_surface_complexation_2_ac_pflotran__engine_reference` | PFLOTRAN | Engine reference | — |
| `batch_chem_surface_complexation_2_pc_pflotran__smoke` | PFLOTRAN | Smoke | The AC variant is covered; this alternate driver configuration has not yet been calibrated against the same native result. |
| `batch_chem_microbial_pflotran__smoke` | PFLOTRAN | Smoke | A standalone reference, comparable output fields, and native-to-driver column mappings have not yet been added. |
| `batch_chem_calcite_short_ac_crunch__smoke` | CrunchFlow | Smoke | The calcite-family time integration and output mapping must be reconciled with the standalone result. |
| `batch_chem_calcite_short_cc_crunch__smoke` | CrunchFlow | Smoke | The calcite-family time integration and output mapping must be reconciled with the standalone result. |
| `batch_chem_calcite_vf_ac_crunch__smoke` | CrunchFlow | Smoke | The calcite-family time integration and output mapping must be reconciled with the standalone result. |
| `batch_chem_calcite_vf_cc_crunch__smoke` | CrunchFlow | Smoke | The current standalone and Alquimia results have material pH and primary-concentration differences that must be understood first. |
| `batch_chem_general_reaction_ac_crunch__engine_reference` | CrunchFlow | Engine reference | — |
| `batch_chem_general_reaction_cc_crunch__smoke` | CrunchFlow | Smoke | The AC variant is covered; this alternate driver configuration still needs its mapping and tolerances validated independently. |
| `batch_chem_ion_exchange_valocchi_ac_crunch__smoke` | CrunchFlow | Smoke | The CC/native-condition variant is covered; this alternate driver configuration still needs a separately validated mapping. |
| `batch_chem_ion_exchange_valocchi_cc_crunch__engine_reference` | CrunchFlow | Engine reference | — |
| `batch_chem_isotherms_ac_crunch__engine_reference` | CrunchFlow | Engine reference | — |
| `batch_chem_isotherms_cc_crunch__smoke` | CrunchFlow | Smoke | The AC variant is covered; this alternate driver configuration has not yet been calibrated against the same native result. |
| `batch_chem_surface_complexation_2_ac_crunch__smoke` | CrunchFlow | Smoke | The standalone/driver path currently produces non-finite values; comparison is deferred until that behavior is resolved. |
| `batch_chem_surface_complexation_2_cc_crunch__smoke` | CrunchFlow | Smoke | Reference work is deferred until the non-finite surface-complexation behavior is resolved and the CC mapping is validated. |
| `batch_chem_calcite_co2_crunch_cc__smoke` | CrunchFlow | Smoke | A normalized reference and mappings for the gas and mineral outputs have not yet been added. |
| `batch_chem_ch4_o2_crunch_cc__smoke` | CrunchFlow | Smoke | A normalized reference and mappings for the gas and aqueous-kinetics outputs have not yet been added. |

## Transport

Transport reference testing is planned for future work. All transport cases
currently remain smoke tests because the reference generator and comparator
handle zero-dimensional batch time series only. Transport comparisons will
also need spatial field selection, mesh/cell alignment, time alignment, and
normalization of the native PFLOTRAN and CrunchFlow output formats.

| CTest name | Engine | Status | Why no engine-reference comparison yet |
|---|---|---|---|
| `transport_tracer_1d_pflotran__smoke` | PFLOTRAN | Smoke | Spatial concentration fields and time levels are not yet normalized or aligned. |
| `transport_tritium_1d_pflotran__smoke` | PFLOTRAN | Smoke | Spatial concentration fields and time levels are not yet normalized or aligned. |
| `transport_calcite_1d_pflotran__smoke` | PFLOTRAN | Smoke | Spatial aqueous and mineral fields need mappings plus mesh and time alignment. |
| `transport_isotherms_1d_pflotran__smoke` | PFLOTRAN | Smoke | Spatial aqueous and sorbed fields need mappings plus mesh and time alignment. |
| `transport_ion_exchange_1d_pflotran__smoke` | PFLOTRAN | Smoke | Spatial aqueous and ion-exchange fields need mappings plus mesh and time alignment. |
| `transport_farea_1d_pflotran__smoke` | PFLOTRAN | Smoke | Spatial aqueous and mineral surface-area fields need mappings plus mesh and time alignment. |
| `transport_tracer_1d_crunch__smoke` | CrunchFlow | Smoke | Spatial concentration fields and time levels are not yet normalized or aligned. |
| `transport_tritium_1d_crunch__smoke` | CrunchFlow | Smoke | Spatial concentration fields and time levels are not yet normalized or aligned. |
| `transport_calcite_1d_crunch__smoke` | CrunchFlow | Smoke | Spatial aqueous and mineral fields need mappings plus mesh and time alignment. |
| `transport_isotherms_1d_crunch__smoke` | CrunchFlow | Smoke | Spatial aqueous and sorbed fields need mappings plus mesh and time alignment. |
| `transport_ion_exchange_1d_crunch__smoke` | CrunchFlow | Smoke | Spatial aqueous and ion-exchange fields need mappings plus mesh and time alignment. |

## Coverage summary

| Group | Engine reference | Smoke | Total |
|---|---:|---:|---:|
| Batch chemistry | 7 | 19 | 26 |
| Transport | 0 | 11 | 11 |
| **All benchmarks** | **7** | **30** | **37** |

Unit tests are not included in this table; they are classified separately with
the CTest label `unit`.
