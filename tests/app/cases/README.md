# fepx executable error cases

Configurations run by the `test_fepx_error_<case>_np<N>` tests (tests/CMakeLists.txt) on the
mesh `tests/regression/cases/meshes/sx2.msh` (one phase), each in its own working directory on 1
and 2 ranks. Each must end every rank with the fatal error listed there:

| Case | Error |
|---|---|
| `phase_count` | number_of_phases 2 on a one-phase mesh |
| `rate_away` | strain rate away from the strain target (step schedule) |
| `no_target` | a preset step without a strain target or a time target (configuration) |
| `scys_vertex_rank` | `initialization isotropic_viscoplastic` with hcp basal + prismatic only: the active systems span 4 deviatoric dimensions, no yield-surface vertex (step 6.7) |
| `opt_both_forms` | `simulation.opt` (`<case>.opt`) with both `$ElsetCrss` and `$ElementCrss` (step 7.1, D59) |
| `opt_count` | `$ElsetCrss` with 2 values for an fcc elset (one, or 12) |
| `opt_tau_s_low` | `$ElsetCrssSat` below the initial strength of the elset's elements |
| `km_opt_both_forms` | Kocks–Mecking phase (step 7.2, D60) with both `$ElsetRho0` and `$ElementRho0` |
| `km_opt_count` | `$ElsetRho0` with 2 values for an elset (one) |
| `km_opt_zero` | `$ElsetRho0` value 0 (the initial density must be positive) |
| `load_rate_away` | a load target (step 8.1) that the strain rate moves away from, from zero load: detected before the run |
| `stress_rate_away_run` | after a tensile strain step, a stress target of 0 approached with a positive strain rate: detected at the start of step 2 during the run |
| `stress_reached_start` | a stress target within its `stress_tol` of the value reached by the previous step: already reached |
| `subincr_depth` | `subincrementation on` (step 7.3, Q35, D63) with `point_max_iters_state 1` and `subincrement_max_depth 2`: every point solve fails at any length; increment 1 is rejected at dt/4 (limit L1) |
| `subincr_successive` | `subincrementation on`, `point_max_iters_state 3`, `subincrement_max_successive 1`: increment 2 is subdivided, increment 3's first attempt is rejected (global iteration failed) and the run stops (limit L2) |

Run cases (must complete; the log lines named are checked):
- `opt_warnings_run` (`test_fepx_run_opt_warnings_np<N>`, step 7.1): two phases on the two
  elsets of sx2 (`opt_warnings_run.phase`), phase 1 crystal plasticity with
  `precipitate_strengthening on`, phase 2 elastic; `opt_warnings_run.opt` gives `$ElsetCrss` for
  both: the warnings "precipitate_strengthening is ignored" (phase 1) and "plasticity parameters
  for phase 2, which is elastic: ignored" are logged.
- `km_opt_warnings_run` (`test_fepx_run_km_opt_warnings_np<N>`, step 7.2): two phases on the
  two elsets of sx2, phase 1 `hardening kocks_mecking`, phase 2 `saturation_isotropic`;
  `km_opt_warnings_run.opt` gives `$ElsetCrss` and `$ElsetRho0` for both: Hill's μ₀ of phase 1
  and the warnings "Crss, CrssSat or HardCoeff values for phase 1, whose hardening is
  kocks_mecking" and "Rho0 values for phase 2, whose hardening is not kocks_mecking" are logged.
- `general_strain_target_run` (`test_fepx_run_general_strain_target_np<N>`): `general`
  conditions with a compressive strain target: the set_bc velocities are reversed (step plan).
- `crystal_plasticity_run` (`test_fepx_run_crystal_plasticity_np<N>`, step 6.2): an hcp
  crystal-plasticity phase (two active families, given out of canonical order): the
  initialization elastic step of 0.03 tau_0 / (largest row sum of [C]) / rate (ROADMAP D49) and
  the first increment shortened by it are logged, and the run completes.
- `subincr_run` (`test_fepx_run_subincr_np<N>`, step 7.3, Q35, D63): `subincrementation on`
  with `point_max_iters_state 3`: increments 2-4 are rejected at full length and solved in parts
  (the rejection and the per-increment summary lines are checked, identical on 1 and 2 ranks);
  `print convergence` with the attempt columns; the run completes.
- `multi_step_run` (`test_fepx_run_multi_step_np<N>`, step 5.2): two steps, the second with a
  shortened last increment (step plan, last increment).
- `load_unload_run` (`test_fepx_run_load_unload_np<N>`, step 8.1): elastic hcp, a stress target
  of 40 MPa, unloading to a load target of 0 N, then a strain target planned during the run from
  the strain reached; the plan and end lines of each step are checked.
- `necking_run` (`test_fepx_run_necking_np<N>`, step 8.1, v2-fixed C54): the v2-fixed necking case
  started as v2; the run stops at increment 23 (as v2-fixed) with the necking error, and writes
  step 1 (`stress.step1`).
- `landing_run` (`test_fepx_run_landing_np<N>`, step 8.2a, addendum_triaxial_load_control.md
  §3.6): fcc crystal to 400 N near yield and back to 0 N with `dtime_factor 2.0`: increment 3
  passes 400 N and lands in 5 false-position attempts, the unloading in 1 (log lines checked);
  `print convergence` with the attempt columns.
- `landing_max_bc_iter` (`test_fepx_run_landing_max_bc_iter_np<N>`, step 8.2a): the same with
  `max_bc_iter 3`: the run stops with the error "no landing on the target" after writing the
  last accepted increment as step 1 (`stress.step1`).

Periodic conditions (step 8.3b, addendum_periodic_mpc.md §5 V7; tests
`test_fepx_error_<case>_np<N>` of the periodic block, each on its own mesh):

| Case | Mesh | Error |
|---|---|---|
| `periodic_no_relations` | `sx2.msh` | `boundary_conditions periodic` on a mesh without `$Periodicity` (v2-fixed C69) |
| `periodic_axis_x` | `periodic_cube_x.msh` | a mesh periodic along x only (P8) |
| `periodic_bc_secondary` | `periodic_box_2x1x1.msh` | `set_bc vel` on node 29, the x drive pair's secondary, whose x velocity is constrained (11 imposed; Q-P6) |
| `periodic_free_rotation` | `periodic_cube.msh` | `deformation_gradient_rate` with only 11 imposed: the three free-rotation warnings (Q-P2), then the run diverges (as v2, finding P2) |

