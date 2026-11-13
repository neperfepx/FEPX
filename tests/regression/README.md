# Regression tests

Regression cases run FEPX on a case and compare every result value with FEPX's own stored
reference data, `reference/fepx/<case>/` (CLAUDE.md §6). Until step 9.13 (2026-10-25) most cases
were compared with the results of the FEPX v2-fixed baseline (`legacy/fepx2-fixed`, with its own
regression suite of v2 cases); the comparison with v2 was then complete (ROADMAP D90), and those
cases were converted to FEPX references generated from `ac211915` (tag `v2-final-comparison`, the
last commit compared with v2). What the comparisons showed is kept below as a record ("Record:
the comparison with v2-fixed") and in `dev/validation/fepx_v2_equivalence.md`; v2-fixed, its
cases and its reference data are in the tag.

## Running

```
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L regression -LE tutorial -j 8 --output-on-failure
ctest --test-dir build -L tutorial -j 8 --output-on-failure     # the tutorial cases, about 25 min
```

Needs Python 3 with NumPy.

Each case is one CTest test `regression_fepx_<case>` (labels `regression`, `fepx`; the tutorial
cases also `tutorial`). It copies the case inputs to `build/tests/regression/work/<case>/` (the
source tree is never written to), runs FEPX with MPI, and compares `simulation.sim/results` with
`reference/fepx/<case>/results` using `tests/tools/compare_sim.py`. Implementation: `run_case.cmake`,
case list: `CMakeLists.txt`. A case with `REFERENCE <other>` compares with another case's
reference (the same results from other inputs or another number of ranks); a case with
`SAME_AS <other>` compares with that case's results of the same build.

## Regenerating reference data

```
FEPX_REGRESSION_UPDATE=1 ctest --test-dir build -L regression -j 8
```

replaces `reference/fepx/<case>/` with the new results for every case that **fails** the
comparison (or has no reference) and leaves the others untouched, so that run-to-run noise below
the tolerance does not churn the reference data. `FEPX_REGRESSION_UPDATE=all` replaces every
reference unconditionally. Cases with `REFERENCE` or `SAME_AS` never write a reference. This is
allowed only under CLAUDE.md §2.5: a documented reason in the commit message and the ROADMAP,
and the user's approval when the change is not caused by a documented bug fix. Review `git diff
--stat` of `reference/` before committing: it must change exactly the cases expected to change.
References are generated with a Release build, Open MPI, on the case's number of ranks.

History of the references:

- `fepx_km_*`, `fepx_irr_*`, `fepx_mixed_unload_*`, `fepx_triaxial_mixed_*`,
  `iso_uniaxial_clr`, `poly_fcc_cyclic_subincr` (features v2 lacks): FEPX's own
  references since steps 7.2 to 8.2c (first generated 2026-10-12, gfortran 15.2.0, Open MPI
  5.0.10); `manual_tutorial_*` since step 8.5 (2026-10-24).
- The 61 former FEPX-vs-v2 cases (every other case with a reference of its own, the
  `fepx_tutorial_*` cases included): generated on 2026-10-25 from `ac211915` (step 9.13, D90,
  approved by the user), where every one of them passed its comparison with v2-fixed (full CI of
  2026-10-25), with a Release build, gfortran 15.2.0, Open MPI 5.0.10: the results of that run
  copied as they were. Reproduced on HEAD (`src/` unchanged since `ac211915`) and in the other CI
  builds as measured below.
- Step 9.2 (2026-10-27, ROADMAP D97, approved by the user): the default preconditioner became
  `block_jacobi` without constraints, so the 75 unconstrained cases with a reference of their own
  (every case but the `fepx_mpc_*`, `fepx_pbc_*` and periodic tutorial cases) were regenerated
  with `FEPX_REGRESSION_UPDATE=all` (Release, gfortran 15.2.0, Open MPI 5.0.10); the 9
  constrained cases (Jacobi, unchanged) were not run. Against the former (Jacobi) references only
  `triaxial_clr` and `triaxial_clr_dual` exceeded `fepx_own_tol` (CG-tolerance noise,
  `dev/benchmarks/preconditioner_study/README.md` M4); all were replaced so that each reference is
  the result of the default preconditioner.

## Comparison and tolerances

`tests/tools/compare_sim.py TEST_SIM REF_SIM` checks that:

1. both `results/` trees contain the same files (a missing or extra file fails);
2. each pair of files has the same number of lines and of values per line; comment lines
   (starting with `%`, e.g. the headers of `forces` files) are identical text;
3. **every value** on every line agrees:

   |test − ref| ≤ max(rtol·|ref|, ftol·scale, atol),  scale = max |ref| over the reference file.

A NaN or infinity always fails. The `ftol·scale` floor handles components that are tiny compared
with the rest of the field (e.g. shear stresses near zero in uniaxial loading), where a pure
relative test is meaningless.


Every case uses `--rtol 0 --ftol 1e-9` (`fepx_own_tol`); the constant-load-rate cases
`triaxial_clr*` add `--atol 1e-18`: they write fields that vanish up to rounding (the
velocities of a dwell, the plastic fields of the elastic step; field scale about 1e-20), whose
noise is reproduced by the optimized builds but not by the `-O0` build of the coverage job
(`triaxial_clr_dwell`, step-2 velocities: 2e-20 against an allowed 1.7e-29, found in step
9.13's full CI); the same floor as in their comparison with v2-fixed. Three further exceptions:
`iso_uniaxial_clr_triaxial` (`--ftol 1e-12`) and `km_sx_bcc_subincr_on` (bit for
bit) against results of the same build, and the restart's cycle 2 (its steps 3 and 4 only).
A reference does not depend on the build beyond the solver's iteration error and rounding: the
cases with tightened solver tolerances (`fepx_tight_*`, `fepx_km_*`, `*_tol`) by construction,
the others because the same iterations are taken in every build.

Measured 2026-10-25 (step 9.13; worst |test − ref| / allowed over every file of a case):

| Build | Cases (regression, not tutorial) | Bit for bit | Others (worst / allowed) |
|---|---|---|---|
| Release, Open MPI (also the 18 tutorial cases: all bit for bit) | 67 | 61 | `km_sx_bcc` 0.125, `poly_fcc_cyclic_subincr` 0.018, `dual_mixedhard_np1` 0.0056, `km_poly_fcc_rho_opt` 0.0034, `dual_bcc_fcc_np1` 0.0014, `irr_sx_bcc` 0.0006 |
| Checked, Open MPI | 67 | 61 | the same values |
| Checked and Release, MPICH | 67 | 62 | the same values (`irr_sx_bcc` bit for bit) |

The 61 former FEPX-vs-v2 references, generated at `ac211915`, are reproduced bit for bit on 2
ranks in every build (and on 4 ranks for the tutorial cases, Release); the `_np1` cases (1 rank
against the 2-rank reference) and the older FEPX references (generated in steps 7.2 to 7.3 with
earlier builds) differ by at most 0.125 of the tolerance.

`tight_poly_fcc_part2_block_jacobi` (step 9.2) runs the tight case with
`cg_preconditioner block_jacobi` against the (Jacobi) reference of `tight_poly_fcc_part2`,
with `fepx_own_tol`: the two preconditioners take different CG iterations (296 against 353 in
the first SA solve), so the results agree to the solver tolerances only. Measured 2026-10-27
(Release, Open MPI): at most 8e-14 of field scale (8e-5 of the tolerance).

Per-variable overrides (`--var-tol VAR=RTOL,FTOL[,ATOL]`) and per-case options (`COMPARE_ARGS` in
`CMakeLists.txt`) are allowed only with a documented reason here.

### Manual tutorials (label `tutorial`, step 8.5)

The inputs of the manual's Tutorials page, `tutorials/<name>` (user, 2026-10-24: the
manual shows FEPX's own inputs, the same examples as the v2 tutorials, with FEPX's defaults: no
`initialization isotropic_viscoplastic`, no `midnode_linearization on`, the default `load_tol`
and solver tolerances), on the same meshes, 4 ranks, against FEPX's own references
`reference/fepx/manual_tutorial_*` (`regression_fepx_manual_tutorial_*`, regenerated by
`FEPX_REGRESSION_UPDATE=1`), so that every input the manual shows is run. Tolerance
`fepx_own_tol` (`--rtol 0 --ftol 1e-9`), as the other own-reference cases at default solver
tolerances. Cycle 2 of the restart compares its steps 3 and 4 only. First generated 2026-10-24
(Release, Open MPI 5.0.10, 4 ranks; 77 MB).

Step 9.11 (ROADMAP D95): the cases are the lines of tier `test` of
`tutorials/tutorials.txt` (read by `CMakeLists.txt` and by
`tutorials/tools/make_tutorial_pages.sh`, which makes the tutorial pages' figures from the same inputs).
Added 2026-10-26: `mixed_targets` (4 ranks) and `hcp_basal`, `hcp_prismatic`, `hcp_pyramidal`
(1 rank, mesh `sx2.msh`, orientation files in `tutorials/hcp_single_crystal/`), new
references. The inputs of `multiaxial_csr`, `multiaxial_clr`, `restart_cycle1`, `restart_cycle2`
and `overwriting` print `stress` and `strain` (and the restart cycles `forces`) for the pages'
stress-strain curves (`neper -S -resmesh`): their references gained these files only, the files
already there unchanged (byte-identical; `FEPX_REGRESSION_UPDATE=1`). The crystal types of the
inputs are written in upper case (`neper -S` aborts on lower case, `dev/neper_findings.md` A3;
only `results/` is compared, so no reference changed). The cruciform specimen (tier `manual`,
mesh `tutorial_cruciform.msh`, about 1.5 h on 4 ranks) is run by the page script only.

Test meshes of Neper 5.0.1-713-devel (step 9.18, ROADMAP D110, 2026-11-02, approved by the user):
the meshes with Neper commands were regenerated, their geometry changed, and the 57 cases on them
have new references (`FEPX_REGRESSION_UPDATE=1`, Release, Open MPI). The differences from the former
references are in `dev/validation/reference_change_neper_713.txt`. The cases whose geometry did not
change (the tutorial meshes of v2) kept their references, on Neper's partition tables.

Test meshes of Neper 5.0.1-731-devel (step 9.20, ROADMAP D113, 2026-11-05, approved by the user):
every mesh with Neper commands regenerated with the same commands, so that all come from one Neper
version. The tessellations are those of 713; the meshes are coarser at the same `-rcl` (sx2 263 → 145
tets, poly 615 → 454, poly_periodic 363 → 239, tutorial_overwriting 5904 → 3835) and stiffer: forces
up to 8 % of their peak higher (`km_poly_fcc_rho_opt`; single crystals 2 to 6 %). The 56 cases on
them have new references (`FEPX_REGRESSION_UPDATE=1`, Release, Open MPI), listed with their
differences in `dev/validation/reference_change_neper_731.txt`; `pbc_box_uniaxial` changed before
(the regularized box, 9.20 item 5). The cases on v2's meshes, whose re-partitioned copies changed
only their tables, kept their references.


### Verification: single-point time integration against RK4 (label `verification`)

`hardening_integration_run` (driver `tests/material/drive_hardening_point.f90`) and the checks
`hardening_integration_{check,spin_check,modifier_check,km_limit_check}` (cases 1–6, 7–10,
11–14, 15–18 of `tests/verification/check_hardening_integration.py`; addendum_hardening.md §5
test 3 and addendum_point_verification.md). Relative errors of the finest run (dt 0.05 s)
against the RK4 reference; first order checked by the observed orders.

| Tolerance | Default | Rationale (measured 2026-10-17, step 7.4a) |
|---|---|---|
| `--tol` | 2e-4 | Cases 1–5 (step 6.4/7.2): worst 1.2e-4. |
| `--tol-irr` | 1e-3 | Case 6 (step 7.3): worst 4.0e-4 (softening with η = 100). |
| `--tol-spin` | 1e-3 | Cases 7–18, stress, strengths and orientation (the last only where the lattice has rotated more than 0.1°): worst 4.6e-4 (hcp stress at 4 s, in the elastic–plastic transition); at 40 s at most 1.1e-4 (cyclic: 4.2e-4). |
| `--tol-extrap` | 5e-7 | Check 5, Richardson extrapolation 2 y(0.05) − y(0.1), stress and strengths, all cases but the cyclic one: worst 1.2e-7 (hcp). |
| `--tol-extrap-ori` | 5e-6 | Check 5, orientation: worst 1.5e-6. Mutation: without the term −½[E, Ė] of the plastic spin (P1) in the reference, cases 7–9 give 2.8e-5–3.4e-5 (and stress 5.6e-7–1.5e-6) and fail; without the strain terms of the plastic spin, or with the slip spin's sign reversed, case 8 fails the first-order check. |
| reference | 0.01 × tol | RK4 against twice its step: at most 1e-8 (cyclic case 11 with its 5× finer step: 1.2e-6). |


## Cases

Inputs in `cases/<case>/` (FEPX syntax), meshes in `cases/meshes/` (`generate_meshes.sh`
records the Neper commands; `triaxial_cube.msh` and `dual_bcc_fcc.msh` are the meshes of v2's
cases `LEGACY_triaxial_csr` / `LEGACY_triaxial_clr` and `05_dual_bcc_fcc`, copied in step 9.13),
the manual's tutorials in `tutorials/` (inputs; pages in `doc/tutorials/`). What each case covers is listed with it in
`CMakeLists.txt`; the record below describes the former FEPX-vs-v2 cases in detail (their names
there are those of the v2-fixed cases, e.g. `sx_fcc` for `sx_fcc`).

## Record: the comparison with v2-fixed (steps 0.3 to 9.4)

Kept as written until step 9.13. The v2-fixed cases (`regression_v2_*`), their inputs
(`cases/<case>` without the `fepx_` prefix, `legacy/fepx2-fixed/tests`), their references
(`reference/v2-fixed/`), `dev/tools/sim_from_h5.py` and the `FEPX_BUILD_V2` build were removed in
step 9.13; they are in the tag `v2-final-comparison`. Tolerances of v2-fixed against its own
references: rtol 1e-5 (v2 writes 7 significant digits; 1, 2 and 4 ranks at most 2.5 % of it),
ftol 1e-6, atol 0, `h5diff -p 1e-5` for the HDF5 case 70.

### FEPX against v2-fixed


FEPX results are not expected to equal v2's to the last digit (ROADMAP D43): v2 always runs
its elasto-viscoplastic model (the elastic cases suppress slip with g_0 = 1e10), and the
formulations differ (v2's deviatoric 5-vector plus volumetric split against FEPX's recoupled
Mandel 6-vectors). Tests `regression_fepx_<case>` (label `fepx`) run FEPX on a case and
compare with a v2-fixed reference.

**What agreement to expect, and why** (measured in step 5.4, 2026-10-10, at the user's
request):

1. *Elastic, with v2's HCP constraint (c33 = c11 + c12 − c13, so that v2's
   deviatoric/volumetric split is exact): the two codes compute the same solution to first and
   second order in the strain.* Measured on `elastic_sx_hcp` at three strain levels, both codes
   converged tightly (nl_tol_strict 1e-9, cg_tol 1e-12) and v2's values read at full precision
   from its HDF5 output: max |FEPX − v2| / field scale for the stress 2.2e-13, 2.2e-10 and
   2.2e-7 at final strains 1e-4, 1e-3 and 1e-2 (elastic strain alike; displacement ≤ 2e-8;
   orientation ≤ 5e-10). The difference grows as ε³ (1000 times per 10 times the strain, and
   about 2³ from step 1 to step 2 of each run): a third-order term whose source is not
   identified. Ruled out: β = det(I + e^e) (v2 also uses the exact determinant,
   `solveit_evp_mod4.f90::determinant`, stresses divided in `finalize_res_mod1.f90`), the
   volumetric update (v2's e_kk += dt d_kk equals the trace of FEPX's equation, the spin
   term being traceless), and v2's single pass of the point state loop (FEPX made to stop
   after one pass, toler_state 1e10, differs from v2 by 4.0e-7 instead of 2.2e-7 at 1e-2).
   Candidates not checked: v2's element assembly (separate deviatoric and volumetric
   stiffness, its {h} terms) and its 5-vector rotation of the stiffness.
2. *What the regression comparisons see is v2's output precision.* v2 writes `e13.7` (7
   significant digits): a rounding of up to half a unit in the 7th digit, 5e-8 to 5e-7 of the
   value (e.g. 0.1412194E+01 ± 5e-7, relative 3.5e-7); forces `e15.5` (5 digits); coordinates
   absolute 5e-8. The measured 1–5e-7 of the multi-step cases is this rounding: it does not
   change with the strain (1e-4 to 1e-2) and appears in the orientation too.
3. *Iteration error at FEPX's default tolerances* (nl_tol_strict 5e-4): up to about 1e-5 of
   field scale after a single first increment (FEPX starts from the prescribed velocities, v2
   from its viscoplastic estimate); negligible after further increments (the previous
   velocities are a good start). Visible in the one-increment cases only.
4. *Crystal plasticity (step 6.6, measured 2026-10-11).* With `stress_init scys_vertex_v2` (named `scys_vertex` in 6.6; v2's
   increment-1 estimate from its tabulated yield-surface vertices) and every tolerance of both
   codes tightened (global 1e-10, CG 1e-13, point Newton 1e-11, state loop and strength update
   1e-12; v2-fixed's `point_tol_state`, `point_tol_hardening`, CHANGES.md A1), FEPX and
   v2-fixed's full-precision HDF5 output agree to **1e-11 of each field's scale** for every
   output both write (fcc, bcc, hcp with v2's c33 constraint, single crystals and an 8-grain
   polycrystal, through 1.2 % strain with hardening; orientations to 1e-14): the two codes
   solve the same discretized equations. At the default tolerances they agreed to at most
   3.4e-4 of field scale in 6.6 (crss 2.9e-4 after 10 increments), attributed then to the solver
   tolerances; most of it was a deviation of FEPX's strength update from v2 (it tested the
   start estimate before any Newton step), corrected in step 7.1: now at most 3.2e-5 (crss
   1.6e-6; `dev/theory/addendum_hardening_modifiers.md` §7). (The elastic ε³ term of item
   1 does not appear: 1.2e-11 for hcp at 0.6 % strain, where ε³ scaling predicts 5e-8; it is
   specific to v2's elastic-only cases, whose state loop does not iterate, see the note on
   τ₀ = 1e10 below.) Without v2's HCP constraint (c33 independent, D3) the two codes solve
   different material problems and are not compared.

Tolerances therefore sit above v2's output precision with a margin (multi-step cases), or above
FEPX's first-increment iteration error (one-increment cases), not at the formulation
difference, which the `.sim` files cannot resolve.

| Cases | Options | Measured (2026-10-09, step 4.4) |
|---|---|---|
| `fepx_elastic_{sx_hcp,poly_fcc,poly_hcp}_incr1` against `elastic_*_incr1` (one elastic increment, strain ε = 2.5e-4, `uniaxial_grip`; FEPX with `midnode_linearization on` and v2's HCP constraint c33 = c11 + c12 − c13) | `--rtol 0 --ftol 5e-5 --var-tol coo=0,2e-7 --var-tol forces=0,1e-4` (step 5.4; 3e-4 in step 4.4) | Re-measured 2026-10-10 (step 5.4, after 5.2's v2-style orientation update): max \|FEPX − v2\| / field scale: stress 1.1e-5, 7.5e-6, 7.2e-6; strain_el 9.1e-6, 6.8e-6, 5.3e-6; disp ≤ 7.4e-7; forces ≤ 2.3e-5; coo 5e-8 (v2's 7-digit output). These are **FEPX's iteration error** in a first increment at the default tolerances (nl_tol_strict 5e-4; it starts from the prescribed velocities, v2 from its viscoplastic estimate): with FEPX tight (nl_tol_strict 1e-9, cg_tol 1e-12) the stress difference to v2 is 1.2e-7 (to v2 also run tight: 9.9e-8; v2 tight vs default: 2e-7). The step-4.4 note ("the differences stay the same with FEPX tight: formulation differences") held before step 5.2, when FEPX kept the start-of-increment stiffness (D40). |
| `fepx_elastic_{sx_hcp,poly_hcp,poly_fcc}` against `elastic_*` (2 steps of 2 increments to 5e-4 and 1e-3), `elastic_sx_hcp_unload` (1e-3, then back to 5e-4), `fepx_elastic_general_{dstrain,dtime}` (`general` conditions, 10 + 5 increments to 2e-3 and 3e-3, D44), `elastic_sx_hcp_triaxial` (`triaxial` with `control_component 33`, 11 and 22 by rates, D45; v2: three `minimal` strain-rate conditions); v2's total strain (`elts/strain`) ignored (not written by FEPX yet) | `--rtol 0 --ftol 2e-6 --var-tol coo=0,2e-7 --var-tol forces=0,1e-4 --ignore elts/strain/*` | Measured 2026-10-10 (step 5.4), worst over steps and variables: stress 0.8–4.4e-7, strain_el 1.4–4.6e-7, disp 1.0–2.5e-7, ori 1.1–4.9e-7, vel 0.5–2.8e-7, forces ≤ 3.9e-5, coo 5e-8. FEPX tight gives the same values (the default iteration is converged after the first increment), and they do not change with the strain (1e-4 to 1e-2): they are v2's 7-digit output rounding (item 2 above; the formulations agree to 2e-10 at this strain, item 1); the face forces differ more and are written by v2 with 5 digits. Tolerances: about 4× the worst field difference, 2.6× the worst force difference. A mutation check (no geometry update in the app) fails all ten comparisons. |
| `fepx_{sx_fcc,sx_bcc,sx_hcp,poly_fcc_part2}` against the v2-fixed cases of the same names (crystal plasticity, saturation_isotropic, 2 steps of 5 increments to 0.2 % and 1 %; FEPX with `stress_init scys_vertex_v2`, `midnode_linearization on`, hcp with c33 = c11 + c12 − c13; default tolerances) | `--rtol 0 --ftol 1e-3 --var-tol coo=0,1e-6` | Step 6.6 (2026-10-11): at most 3.4e-4 of field scale (slip, plastic strain at step 1; crss 2.9e-4); coordinates 3.6e-7; face forces at most 0.13 of the tolerance. Step 7.1 (strength update corrected to v2's rule): at most 3.2e-5 (crss 1.6e-6). |
| `sx_hcp_scys` against `sx_hcp`; `fepx_tight_{sx_hcp,poly_fcc_part2}_scys` against `tight_*` (step 6.7, D57: the same cases started from the yield-surface vertices computed from the phase's slip systems and `tau_0`, `stress_init scys_vertex`) | as the cases they repeat (`--ftol 1e-3`, coordinates 1e-6; tight: `--ftol 1e-9`) | Step 6.7 (2026-10-11): against v2-fixed the same as with v2's tables (3.4e-4 at most at default tolerances, 1.3e-11 tight); against the `scys_vertex_v2` runs at most 1.0e-11 (default) and 3.5e-13 (tight) of field scale, hcp included: the vertices only start the Newton solve of v2's estimate (dev/theory/addendum_yield_vertices.md §1, §9). |
| `fepx_{sx_fcc_sat_evo,sx_fcc_precip,sx_hcp_precip_bowing,poly_fcc_cyclic,poly_fcc_elset_opt,poly_fcc_elt_opt}` against the v2-fixed cases of the same names (step 7.1, D59: saturation evolution, precipitates (cutting with both sets; bowing only, hcp), cyclic hardening with a load reversal at the same rate, per-elset and per-element values from `simulation.opt`; started and converged as the step 6.6 cases) | `--rtol 0 --ftol 1e-3 --var-tol coo=0,1e-6` (as step 6.6) | Step 7.1 (2026-10-11): at most 2.4e-5 of field scale (slip, plastic strain; poly_fcc_elt_opt), 1.1e-5 for poly_fcc_elset_opt, below 1e-5 for the others; coordinates at most 0.5 of their tolerance (cyclic). Mutation: without the reset of the accumulated slip at the reversal, `poly_fcc_cyclic` fails (crss of step 2, 1.5 times the tolerance). |
| `fepx_{sx_fcc_aniso,sx_bcc_aniso,dual_bcc_fcc}` against the v2-fixed cases of the same names, `dual_bcc_fcc_np1` (the same inputs on 1 rank) against `dual_bcc_fcc` (step 7.4a: v2's fully anisotropic hardening, v2 tests 32/33, 1 on the diagonal and 1.4 off it = FEPX's default matrix; two phases, fcc + bcc, v2 test 05's parameters and mesh, one phase per partition; started and converged as the step 6.6 cases) | `--rtol 0 --ftol 1e-3 --var-tol coo=0,1e-6` (as step 6.6) | Step 7.4a (2026-10-17): worst 0.050 (fcc), 0.43 (bcc), 0.060 (two phases, 2 ranks and 1 rank) of the tolerance (recorded in 7.4a as "5.0e-5, 4.3e-4, 6.0e-5 of field scale", which is this ratio times 1e-3; the bcc value is set by the coordinates, 4.3e-7 mm); outputs other than the coordinates at most 3.0e-5, 5.4e-5, 3.9e-5 of field scale (re-measured in step 7.4b). Mutation: FEPX with the isotropic law against the anisotropic references fails by 28 (fcc) and 56 (bcc) times the tolerance. The two-phase comparison found v2-fixed error C51 (references of `05_dual_bcc_fcc`, `46_mixedhard` regenerated; before the fix FEPX was 1.15e3 times the tolerance away). |
| `fepx_{sx_hcp_family,sx_bcc_112,sx_bct,sx_bcc_cyclic,sx_bcc_precip,dual_mixedhard}` against the v2-fixed cases of the same names, `dual_mixedhard_np1` (the same inputs on 1 rank) against `dual_mixedhard` (step 7.4b, one case per feature: hcp with per-family m and `tau_0`, three families (v2 test 10); bcc {110} + {112} (v2 test 44); bct, ten families (v2 test 04's parameters without saturation evolution); bcc cyclic hardening with a load reversal (v2 test 16, threshold a = 0.002); bcc precipitates, both sets, cutting governs (v2 test 40, a_p 50 MPa); two phases with v2 test 46's hardening matrices, fcc fully anisotropic and bcc v2's 2 × 2 block form written as `interaction_matrix` rows; started and converged as the step 6.6 cases) | `--rtol 0 --ftol 1e-3 --var-tol coo=0,1e-6` (as step 6.6) | Step 7.4b (2026-10-17): worst 0.058 (hcp), 0.42 ({112}), 0.39 (bct), 0.48 (cyclic), 0.32 (precipitates), 0.095 (two phases, 2 ranks and 1 rank) of the tolerance, set by the coordinates except for hcp; outputs other than the coordinates at most 3.9e-5, 4.0e-5, 1.6e-4 (slip), 5.4e-5, 6.5e-5, 7.0e-5 of field scale. Mutation (each feature changed in FEPX only): equal m for the three hcp families 15.8 times the tolerance, {112} strength equal to {110} 838, bct strengths all 22 MPa 3.2e3, cyclic hardening off 18.5, precipitates off 765, phase 2 with the default matrix 67. |
| `fepx_tight_{sx_fcc_aniso,sx_bcc_aniso}` against `tight_sx_{fcc,bcc}_aniso` (step 7.4a: the anisotropic cases with every tolerance tightened, as the step 6.6 `tight_*` cases) | `--rtol 0 --ftol 1e-9` (equivalence, D55) | Step 7.4a (2026-10-17): at most 7.2e-12 (fcc) and 2.0e-11 (bcc) of field scale. |
| `fepx_tight_{sx_fcc,sx_bcc,sx_hcp,poly_fcc_part2}` against v2-fixed cases `tight_*` (the same materials, 2 steps of 2 increments to 0.2 % and 0.6 %, every tolerance tightened in both codes, v2's HDF5 output converted at full precision by `dev/tools/sim_from_h5.py`; every output both codes write) | `--rtol 0 --ftol 1e-9` (equivalence, D55) | Step 6.6 (2026-10-11): at most 1.3e-11 of field scale (worst 0.013 of the tolerance). The v2 cases `tight_*` compare their HDF5 output with h5diff like `70_hdf5_output`. |
| `load_target_unload`, `load_target_general` against the v2-fixed cases of the same names with `_tol` (since step 8.2a: both codes with `load_tol 0.5`, above v2-fixed's largest overshoot, 0.35 N, so that FEPX's landing rule redoes no increment and both end each step at the first increment within 0.5 N of the target; step 8.1: load targets; fcc on `sx2.msh`, started and converged as the step 6.6 cases; presets: grip along z to 450 N and unloaded to 0 N; general conditions: x1 pulled to 400 N and back to 100 N, the velocities reversed; FEPX takes the signed strain rates, v2 reverses itself) | `--rtol 0 --ftol 1e-3 --var-tol coo=0,1e-6` (as step 6.6) | Step 8.1 (2026-10-17): worst 0.47 and 0.29 of the tolerance (Release, Open MPI, 2 ranks); the increments chosen by the load rule are the same as v2-fixed's (18 and 15). The comparison is insensitive to `dtime_factor` 1.001 → 1 (0.56, 0.41): that rule is checked by `verification_fepx_hcp_sx_load_targets` and `test_step_schedule`. Step 8.2a (2026-10-17, against `load_target_{unload,general}_tol`): 0.446 and 0.373 of the tolerance, the same increments as v2-fixed (9 + 6 and 8 + 4); no landing attempt. |
| `triaxial_csr` against `LEGACY_triaxial_csr`, `triaxial_csr_pert_floor` against `triaxial_csr_pert_floor` (step 8.2b, addendum_triaxial_load_control.md §5 V3: triaxial load control at a constant strain rate, v2 CSR; v2's hcp unit cube, compression along x to (−300, 0, 300) and (−600, 0, 600) N, 22 and 33 on each step's load path; `load_tol 10` (v2's tolerance below 1000 N), `min_pert_frac 0.01` in the second; started as the step 6.6 cases, `midnode_linearization on`: v2's solves linearize the mid-nodes, and since v2-fixed C65 (D70) its triaxial drivers keep them in the stored coordinates, as its uniaxial driver; finding F6) | `--rtol 0 --ftol 1e-4 --var-tol coo=0,1e-6`, the second with v2's force files (an extra length column) not compared | Measured 2026-10-19 after v2-fixed C65: 5.5e-5 and 5.2e-5 of field scale (vel and velgrad of step 2), the same increments and iterations. Measured 2026-10-18, before C65 (FEPX with `midnode_linearization off`) (Release, Open MPI, 2 ranks): the same increments (11) and boundary-velocity iterations (14, 3 with a correction) as v2-fixed; `LEGACY_triaxial_csr`: rotrate and velgrad of step 2 1.2e-3 of field scale, every other field within 1e-3, coordinates within 1e-6; `triaxial_csr_pert_floor`: 0.54 of the tolerance, face forces within 0.05 N of v2's 5-digit values. Mutations: the control's path start taken as 0 instead of the measured value fails both (and the verification cases); face forces of the trial solves on the start geometry fail both: v2's trial loads are on the end-of-increment geometry too (addendum §8, survey correction). |
| `triaxial_csr_tol` against v2-fixed case `triaxial_csr_tol` (step 8.2b, V3 with tightened tolerances: `LEGACY_triaxial_csr` with the load tolerance 0.05 N in both codes (above v2's largest overshoot, 0.03 N, so that both end each step at the same increment), `dtime_factor 1.0`, `dtime_min` 1e-6 (each step approached from below), global and CG tolerances 1e-7 and 1e-13 (1e-10 stalls near 2e-6 in v2's plastic increments of this case); v2's HDF5 output at full precision) | `--rtol 0 --ftol 2e-7` | Measured 2026-10-19 after v2-fixed C65 (D70; `midnode_linearization on` in both): the same 14 increments and iterations; at most 5.6e-8 of field scale (rotrate of step 2), stresses 1.9e-9, coordinates 2.3e-12. Measured 2026-10-18, before C65 (Q41(b): v2's accepted state update on unlinearized coordinates): state fields at most 2.8e-4 of field scale, the rates of step 2 (last increment 0.15 s) at most 2.75e-3: the order of the load tolerance (each code stops anywhere within 0.05 N; FEPX's path from the measured start values, Q-T4, v2's from the previous targets). With 1e-3 N instead, v2 passes step 1's target by 0.03 N, FEPX redoes and lands that increment, and the two runs' increments drift apart (rates 1.9e-2). |
| `triaxial_mixed_poly_fcc` against FEPX's own reference (`reference/fepx/`, step 8.2b, V5 partly: the 8-grain fcc polycrystal, triaxial with control 33: strain to 0.8 % with σ11 in time to −50 MPa, unloading to σ33 = 100 MPa with σ11 on the path to 0, reloading past yield to 375 MPa with σ11 to −100 MPa; `stress_tol 0.01`; default solver tolerances) | `--rtol 0 --ftol 1e-9` | Measured 2026-10-18 against the 2-rank Release Open MPI run: Debug 4.9e-12 of field scale; 1, 2 and 4 ranks within 6.8e-13 (stresses, `verification_fepx_rank_independence_triaxial`). |
| `triaxial_clr`, `triaxial_clr_dwell`, `triaxial_clr_dual`, `triaxial_clr_dwell_comp` against the v2-fixed cases `LEGACY_triaxial_clr`, `triaxial_clr_dwell`, `triaxial_clr_dual`, `triaxial_clr_dwell_comp` (step 8.2c, addendum_triaxial_load_control.md §5 V3: triaxial load control at a constant load rate, v2 CLR; v2's fcc unit cube (the two-phase mesh of 05_dual_bcc_fcc for `_dual`), 33 at 2 N/s, 11 and 22 in time; dwells as `load_rate 0` steps with `target_duration`; v2's equal load-step increments as `number_of_incrs`; `load_tol 5` (v2's `load_tol_abs`); `midnode_linearization on`, v2-fixed C65, D70) | `triaxial_clr`: `--rtol 0 --ftol 5e-3 --atol 1e-18 --var-tol coo=0,1e-6`; `_dwell`: `--ftol 5e-4`, vel 1e-2, slip and strain_pl 5e-3; `_dual`: `--ftol 1e-3`, vel 2e-3, absolute 2 MPa on stress, 1e-5 on strain_el, strain, disp and coo (the unloaded step 3); `_dwell_comp` (since step 8.5: `load_tol 0.1` in both inputs): `--ftol 1e-3`, coo 5e-6; zero-scale fields atol 1e-18; v2's force files (an extra length column) not compared | Measured 2026-10-19 after v2-fixed C65: `LEGACY_triaxial_clr` 2.6e-3 of field scale (disp of the unloaded step 3); `_dwell` 2.6e-4 on the state fields, 7.1e-3 on vel and 2.9e-3 on slip and strain_pl after the reload (Q-T4); `_dual` 3.7e-4 at steps 1 and 2, vel of step 3 1.6e-3, step 3 within the absolute tolerances (1.27 MPa); `_dwell_comp` 5.0e-2 (F5, Q-T9). Step 8.5 (2026-10-24, after v2-fixed C98, each step starting from the loads reached): `_dwell_comp` at `load_tol 5` 0.14 of the residual stress scale after the unload (1.39 × the then tolerance 0.1): each code holds the dwell loads at other points of the 5 N band (v2 F33 −376.1, −375.3, −378.8 N; FEPX −374.1, −372.4, −372.0 N) and the creep goes as the load to the power 50; with the same `load_tol` in both inputs the worst difference falls with it (5, 1, 0.5, 0.1 N: 0.14, 1.8e-2, 9.0e-3, 7.9e-4 of field scale), the same 36 increments; at 0.1 N (user, 2026-10-24; v2-fixed reference regenerated for this input change) worst 0.86 of the standard tolerance (strain_el, step 7), coo 3.6e-6 mm. Measured 2026-10-18 before C65 (FEPX with `midnode_linearization off`; Release, Open MPI, 2 ranks): the same increments as v2-fixed (8, 14, 8, 36), in `LEGACY_triaxial_clr` the same iterations (10 + 4) and, in step 1, the same loads at every iteration; `LEGACY_triaxial_clr` 7.9e-3 on the stresses of step 3 (unloaded, field scale 2 MPa), every other field within 5e-3; `_dwell` 2e-4 on the state fields, the reload after the second dwell (from the measured load, Q-T4) 7.6e-3 on vel and 2.5e-3 on slip and strain_pl; `_dual` 1.7e-3 at steps 1 and 2, the unloaded step 3 within the load tolerance (1.28 MPa); `_dwell_comp` 4.8e-2 (the −375 N step ends 1.4 N apart in the two codes, at 0.8 % plastic strain; the episode starts differ, Q-T9, finding F5), the creep velocities of its dwells by up to a factor 2 (slip rate as load^(1/m), 1/m = 50). Mutation: the values to reach taken at an attempt's start instead of its end fails `triaxial_clr` (and V1, V5). |
| `triaxial_clr_tol` against v2-fixed case `triaxial_clr_tol` (step 8.2c, V3 tightened: `LEGACY_triaxial_clr` with the load tolerance 1e-3 N in both codes, global and CG tolerances 1e-7 and 1e-13 and the material-point tolerances of the step 6.6 cases; v2's HDF5 output at full precision) | `--rtol 0 --ftol 1e-6 --atol 1e-18` | Measured 2026-10-19 after v2-fixed C65 (D70; `midnode_linearization on` in both): the same 8 increments and iterations; step 1 (elastic) 1e-12 of field scale; steps 2 and 3 at most 4.7e-7 (the rates and disp of the unloaded step 3), stresses 2.4e-9. Measured 2026-10-18 before C65: step 2 1.3e-4 on the stresses, 1.5e-3 on the plastic rates, step 3 7.6e-3; the cause (Q41(b), validation document §8c) was finding F6: v2's triaxial drivers solved on linearized mid-nodes but updated the accepted increment's state on the unlinearized coordinates (an earlier check with FEPX's `midnode_linearization on` matched neither half). |
| `iso_uniaxial_clr` against FEPX's own reference, `iso_uniaxial_clr_triaxial` against its results (step 8.2c, V2: the isotropic elastic cube under `uniaxial_symmetry` at a load rate, and the same loading with `triaxial` and only 33 load-controlled) | own: `--rtol 0 --ftol 1e-9`; the triaxial run: `--ftol 1e-12` | Measured 2026-10-18: the two runs identical (Release and Debug). |
| `triaxial_mixed_clr_sx_fcc` against FEPX's own reference (step 8.2c, V5: fcc single crystal, triaxial with control 33: strain to 0.4 % with σ11 in time to −20 MPa, unloading at −200 MPa/s to 0 with σ11 in time, a 0.5 s dwell with σ11 held by its rate 0, reloading at a strain rate to 300 MPa with σ11 on the path to −30 MPa; `stress_tol 0.05`; default solver tolerances) | `--rtol 0 --ftol 1e-9` | Measured 2026-10-18: 1, 2 and 4 ranks within 1.7e-13 (stresses, `verification_fepx_rank_independence_triaxial_clr`). |
| `mpc_plane_faces` against v2-fixed case `mpc_plane_faces`, `mpc_linear` against `mpc_linear` (step 8.3a, addendum_periodic_mpc.md §5 V6: multi-point constraints on the 8-grain bcc polycrystal `poly_part2.msh`, general conditions, uniaxial along z; `mpc_plane_faces`: v2 test 47's `set_mpc1` on x0 and x1 (x), crystal plasticity to 0.6 % strain, started as the step 6.6 cases; `mpc_linear`: elastic, `set_mpc1` on x1 (x) with a condition on its secondary node x1y1z1 (C66), `set_mpc` on y1 (y, c = 0.5, d = 2e-4, C67); `midnode_linearization on`) | `mpc_plane_faces`: the step 6.6 tolerances (`--rtol 0 --ftol 1e-3 --var-tol coo=0,1e-6`); `mpc_linear`: those of step 5.4 (`--ftol 2e-6`, coo 2e-7, forces 1e-4), the step-0 velocities not compared (FEPX writes the initial guess projected onto the constraints, offsets included; v2 the prescribed values only) | Measured 2026-10-21 (Release, Open MPI, 2 ranks): `mpc_plane_faces` at most 8.8e-5 of field scale (slip; strain_pl 5.7e-5, vel 2.4e-5, stress 6.4e-6, coo 5.7e-8); `mpc_linear` 3.6e-7 (vel; stress 1.3e-7, coo 5.0e-8). The case of `mpc_plane_faces` with even partitioning on 1, 2 and 4 ranks (`verification_fepx_rank_independence_mpc`, MPC groups across the partitions): stresses within the 1e-7 tolerance. |
| `fepx_pbc_{box_uniaxial,poly_full_state,poly_uniaxial,poly_velgrad_free,poly_plane_strain}` against the v2-fixed cases `pbc_*` (step 8.3b, addendum_periodic_mpc.md §5 V6: periodic conditions, crystal plasticity with v2 tests 49/50's bcc material, `initialization isotropic_viscoplastic_v2`, `midnode_linearization on`; the 2 × 1 × 1 box and the periodic 8-grain polycrystal; Cases below) | the step 6.6 tolerances with `--remove-mean coo,disp,vel` (FEPX fixes the rigid translation by the zero-mean fluctuation, v2 does not) and `--ignore macro` (FEPX only); `pbc_poly_full_state` also `--var-tol slip=0,1e-3,1e-8 --var-tol strain_pl=0,1e-3,1e-8` | Measured 2026-10-21 (Release, Open MPI, 2 ranks): at most 0.39 of the tolerance (3.9e-4 of field scale) in every case; `pbc_poly_full_state`'s step 1 (0.1 %) is nearly elastic: one element starts to slip (slip 1.6e-7, the field scale), where slip ∝ (τ/τ_crit)^50 turns a stress difference within its tolerance into 1.4 % of that slip (2.3e-9 slip, 9.6e-10 strain_pl): the absolute floor of 1e-8 is about 4 × that. Before D78 (`stress_init scys_vertex_v2`, no isotropic velocity start) these runs diverged (stresses up to 1.8e7 MPa): the vertex estimate started from the offset field of the periodicity relations, d_eq up to 8.7 1/s in the small face elements of Neper's periodic meshes. |
| Effect of the isotropic viscoplastic start (D78, measured 2026-10-21): every FEPX-vs-v2 case started as v2, before (`stress_init scys_vertex_v2`: increment 1 from the boundary-condition guess) and after (`initialization isotropic_viscoplastic_v2`: from v2's isotropic viscous solve) | unchanged | Worst diff / allowed unchanged to three digits in every default-tolerance case (`sx_*`, `poly_*`, `dual_*`, `load_target_*`, `mpc_plane_faces`, `triaxial_*`, `*_tol`): their worst field is set by v2's 7-digit output and the face forces, not by the start. With every tolerance tightened (`tight_*`, `--ftol 1e-9`) the start shows: `tight_sx_bcc` 0.0133 → 5.9e-4, `tight_sx_hcp` and `tight_sx_hcp_scys` 0.0123 → 0.0065, the others unchanged (0.0019–0.020). The FEPX references of `fepx_{irr_poly_bcc_opt,irr_sx_bcc,km_poly_fcc_rho_opt,km_sx_bcc,poly_fcc_cyclic_subincr}` were regenerated (user, 2026-10-21): the new start changes their results by about 1e-6 of the forces and 1e-11 of the fields (the iteration tolerance), up to 675 × their `--ftol 1e-9`. The norms of the three isotropic iterations equal v2-fixed's log (`sx_fcc`: 0.11768, 0.022962, 0.014255 against 0.1177, 0.02296, 0.01425). |
| Step 9.4 (2026-10-25, D56, D89): `isotropic_viscoplastic_v2` (v2's vertex tables) removed; every case that used it runs with `initialization isotropic_viscoplastic` (the computed vertices); `sx_hcp_scys` and `fepx_tight_{sx_hcp,poly_fcc_part2}_scys` removed (identical to `sx_hcp`, `fepx_tight_{sx_hcp,poly_fcc_part2}` after the switch) | unchanged | Measured (Release, Open MPI, 2 ranks): worst diff / allowed `tight_sx_fcc` 0.0019, `tight_sx_bcc` 5.9e-4, `tight_sx_hcp` 0.0065, `tight_sx_fcc_aniso` 0.0072, `tight_sx_bcc_aniso` 0.020, `tight_poly_fcc_part2` 0.0115; `sx_fcc` 0.05, `sx_bcc` 0.325, `sx_hcp` 0.050, `pbc_poly_full_state` 0.388; FEPX-own `poly_fcc_cyclic_subincr` 0.0177 (no reference regenerated). |
| Final FEPX-vs-v2 comparison (2026-10-25, ROADMAP D90): every FEPX-vs-v2 regression and tutorial case passed in the full CI on `ac211915`; from then on FEPX is compared only with itself (the cases converted to FEPX's own references from `ac211915`, step 9.13) | pass | `dev/tools/ci_local.sh` (every job) on `ac211915`. |
| `mixed_unload_sx_fcc`, `mixed_unload_poly_fcc` against FEPX's own references (step 8.1, addendum_load_targets.md §5 V4: mixed per-step targets, which v2 lacks; fcc single crystal and 8-grain polycrystal, strain to 1 %, unloaded to zero stress, reloaded to 1.2 % in a step planned during the run; every tolerance tightened in the case) | `--rtol 0 --ftol 1e-9` | Measured 2026-10-17 against the references (Release, Open MPI, 2 ranks): 1 and 4 ranks at most 4.4e-5 of the tolerance. Step 1 of the single crystal is bit-identical to a strain-only run. Step 8.2a (landing rule, default `stress_tol` 0.1 MPa): references unchanged, the unloading steps end within the tolerance (−0.0047 and −0.026 MPa) with no increment redone; the tests pass against the 8.1 references. |
| `km_sx_bcc`, `km_poly_fcc_rho_opt` against FEPX's own references (`reference/fepx/`, step 7.2, D60: Kocks–Mecking hardening, which v2 lacks; bcc iron with Barton's parameters on `sx2.msh`, 2 steps to 0.2 % and 2 %; the 8-grain fcc polycrystal with `$ElsetRho0` from 1e4 to 8e7 mm⁻², 2 steps to 0.2 % and 1 %; every tolerance tightened in the case) | `--rtol 0 --ftol 1e-9` | Measured 2026-10-12 against the reference (Release, Open MPI, 2 ranks): 1 and 4 ranks, Debug: at most 6.2e-14 of field scale; MPICH Release: identical. |
| `irr_sx_bcc` against FEPX's own reference (`reference/fepx/`, step 7.3, D61/D62: irradiation softening, which v2 lacks; bcc iron with Barton's parameters and loops 3.61e12 mm⁻³ on `sx2.msh`, elastic to 0.3 % in 3 increments, then 45 increments of 1e-4 through the peak stress into the stress drop, damped state loop; every tolerance tightened in the case) | `--rtol 0 --ftol 1e-9` | Measured 2026-10-12 against the reference (Release, Open MPI, 2 ranks): 1 and 4 ranks at most 4.6e-5 of the allowed difference (about 5e-14 of field scale); identical with the damping floor 0.01 instead of 0.05 (the floor is not reached in this case). |
| `irr_poly_bcc_opt` against FEPX's own reference (`reference/fepx/`, step 7.3, D61/D62: the 8-grain bcc iron polycrystal on `poly_part2.msh` with Barton's parameters, Kocks–Mecking from 2e7 mm⁻² with irradiation softening, `$ElsetIrrLoopDensity` 1.5e12 to 5.5e12 mm⁻³ (the cfg value replaced), the cfg loop diameter 100 b, ⟨111⟩ loops; elastic to 0.3 % in 3 increments, then 84 increments of 5e-4 to 4.5 %: peak 679 MPa at 0.55 %, minimum 590 MPa at 2.65 %, re-hardening to 609 MPa; **default solver tolerances** (user, 2026-10-12: tightened ones cost about 30 s per increment). Replaces `fepx_irr_poly_fcc_opt` (96bda01; user: the case must reach the re-hardening).) | `--rtol 0 --ftol 1e-9` | Measured 2026-10-12 against the 2-rank Release Open MPI run: 1 and 4 ranks at most 1.8e-2 of the allowed difference, Debug 1.3e-2, MPICH identical (about 2e-11 of field scale at most), so the default tolerances do not need a looser comparison here. Reference regenerated 2026-10-12 for D62 revision 3 (user-approved; the damped iterates change at default tolerances, about 9e-7 of field scale; σ₃₃ unchanged). |
| `poly_fcc_cyclic_subincr` against FEPX's own reference (`reference/fepx/`, step 7.3, Q35, D63: `poly_fcc_cyclic` with `subincrementation on`; increments 6, 7, 8 solved as 3 parts down to dt/4, increment 16 as 2, no point left failed; **default solver tolerances**: with every tolerance tightened one point of increment 6 still fails at dt/16 and the run stops at the depth limit) | `--rtol 0 --ftol 1e-9` | Measured 2026-10-12 against the 2-rank Release Open MPI run: 1 and 4 ranks 1.6e-2 of the allowed difference, Debug 1.2e-2, MPICH identical; the same subdivisions in every run. |
| `sx_fcc_explicit_options` against the v2-fixed case `sx_fcc` (step 9.6, D88: `sx_fcc` with `initialization isotropic_viscoplastic` and the keywords that no other run uses, at values that leave the results unchanged: `read_ori_from_file` with the mesh's orientation in `meshes/sx2.ori`, `cg_max_iters`, `global_nr_switch`, `global_nr_slow_ratio`, `point_max_iters_newton`, `point_max_iters_hardening` at their defaults, `max_strain`, `max_eqstrain`, `max_total_time` not reached; a quick test, so that the quick tests use every keyword) | `fepx_v2_cp_tol` (as `sx_fcc`) | Measured 2026-10-24 (Debug, Open MPI, 2 ranks): worst 0.05 of the tolerance, as `sx_fcc` (D78). |
| `km_sx_bcc_subincr_on` against the results of `km_sx_bcc` **in the same build** (its work directory; CTest fixture), step 7.3, D63: `subincrementation on` without a rejected attempt must not change any result | `--rtol 0 --ftol 0 --atol 0` | Bit for bit by construction (the same operations); a stored reference is not used, since the builds differ in the last digits. |
| Note on v2's `g_0` = 1e10 (τ₀; user question, step 5.4) | | Slip is exactly zero already (γ̇ = γ̇₀ \|τ/τ₀\|^(1/m) underflows: (1e-8)^100), so a larger τ₀ changes nothing in the slip; but v2's per-point state loop stops when the stress norm changes by less than toler_state × τ₀ = 1e5 MPa (`solve_evp_mod2.f90`), so v2 never iterates the coupling of the lattice spin and the stress (its Newton stress solve is converged), whatever τ₀ ≥ 1e10. FEPX iterates it (D40); the difference is part of the small values above. |


#### v2 tutorials (label `tutorial`, step 8.5, D83, D86)

FEPX v2's tutorials (all but the cruciform specimen, about 2.5 h per case on v2) in FEPX
syntax, `cases/tutorial_*` (not the manual's inputs: see "Manual tutorials" below), on v2's tutorial meshes at
msh 2.4 (`cases/meshes/tutorial_*`), 4 ranks, against the v2-fixed references `tutorial_*` (v2's
tutorial inputs, repaired in v2-fixed C89–C92, C96). Outside the default runs:
`ctest -L tutorial -j 8` (with `-DFEPX_BUILD_V2=ON`, Release; about 30 min), job `tutorial` of
`dev/tools/ci_local.sh`.

Tolerance `tutorial_tol`: `fepx_v2_cp_tol` (`--rtol 0 --ftol 1e-3 --var-tol coo=0,1e-6`) with
`--atol 1e-12` and, on `slip` and `sliprate`, 1e-15: v2 rounds its text output to an absolute
1e-12, slip-system fields to 1e-15 (`OUTPUT_PRECISION`, CHANGES.md C7 note; FEPX writes
unrounded values, D25), so that, for example, the plastic strain of an elastic step (4.7e-14 in
FEPX at a slip rate of 6.3e-13 /s, the same in both codes) is 0 in v2. v2's triaxial force files
(an extra length column) are not compared; the periodic case removes the mean of `coo` and `disp`
and does not compare `macro`.

| Case | Notes | Measured (2026-10-23), worst difference / allowed |
|---|---|---|
| `simple_simulation`, `simple_rve`, `simple_rve_grip`, `periodic_rve` | uniaxial and periodic | pass (step 8.5, 2026-10-22) |
| `multiaxial_csr` | triaxial CSR, strain rate doubled in step 2: FEPX scales the previous velocities by the ratio as v2 (D86); both start each step from the loads reached (v2-fixed C98) | 0.49 |
| `multiaxial_clr` | triaxial CLR with a dwell: every load held by `load_rate 0` (as v2-fixed C98) | 0.42 |
| `overwriting` | uniaxial load targets: `load_tol 2` in both inputs (v2's copy `cases/tutorial_overwriting`): at the tutorial's `dtime_min` v2 passes 425 N by 1.27 N after being 1.20 N short, so no FEPX landing tolerance alone reproduces its increments | 0.49 |
| `restart_cycle1`, `restart_cycle2` | triaxial CSR with cyclic hardening and a reversal; cycle 2 restarts from cycle 1 (v2 writes cycle 2's steps 3 and 4 only). **Both inputs with tightened solver tolerances** (global 5e-6, global_abs 1e-12, cg 1e-10, point 1e-6, point_state 1e-7, point_hardening 1e-6; v2's copies `cases/tutorial_restart_cycle*`; user decision, 2026-10-24). The 0.75 % difference before v2-fixed C100 (accumulated slip added at every trial solve) did not shrink with tighter tolerances. After C100, at the default tolerances, one value of one element (4045 of 5641, slip system 5, at the onset of slip) was 1.0e-2 of field scale, every other within 2.2e-4: against the converged answer (the tutorial's load_tol 10 kept, the solver tolerances tightened by 100 and 10 000) v2 was off by 1.15e-2 there (4.6 % of the value) and FEPX by 1.0e-3; v2 with only its point tolerances tightened is converged there, with only the global ones mostly, with only cg_tol unchanged. At the tightened tolerances the codes agree to 2.1e-7 of field scale everywhere (4e-9 at element 4045), and each code changes by less than 9e-7 at the tighter level | 0.27 (coo), 0.50 (cycle 2: defrate_pl_eq of step 3, 5e-13 against v2's rounding floor 1e-12); every field within 6e-4 of its tolerance otherwise |


### v2-fixed's cases

#### v2's test cases (74, from `legacy/fepx2-fixed/tests`)

Every v2 test case is included (ROADMAP D23), on 2 ranks as in v2's own suite:

| Group | Cases | Covers |
|---|---|---|
| Crystal types | 01–04, 05 | bcc, fcc, hcp, bct; dual phase bcc + fcc |
| Hardening | 06–19, 32–39, 46 | saturation evolution, anisotropic m, coplanar / fully anisotropic interaction, cyclic, anisotropic τ₀, mixed hardening |
| Orientations | 22–30 | rodrigues, euler-bunge, euler-kocks, axis-angle, quaternion; active and passive |
| Precipitates, bcc {112} | 40–45 | precipitation hardening; bcc {112} slip |
| Boundary conditions | 31, 47, 49, 50 | general BCs, multi-point constraints, periodic (triaxial, uniaxial) |
| Per-element/elset input | 51–60 | `.opt` files with per-element / per-elset τ₀, τ_s0, h₀ |
| Restart | 61, 62 (×3) | plain, CSR, CLR; part 2 continues from part 1's restart files (CTest fixture) |
| Loading targets | 63–69, LEGACY_* | target load, strain (dtime, dstrain, incr), time; legacy triaxial CSR/CLR input |
| HDF5 output | 70 | compared with `h5diff` |

Cases 05, 46, 67, 68 and 69 are disabled in v2's own suite (`test.cmake-disable`); they run
normally and are included. Not included: `tests/Testing` and `tests/test_mesh_gen` (not cases).

All v2 cases run on 2 ranks, as in v2's own suite. Two extra cases check rank-count invariance
against the 2-rank references (`REFERENCE` option): `02_fcc_np4` (4 ranks) and `03_hcp_np3` (3
ranks). Pristine v2 hung on these small meshes with 3 or more ranks (CHANGES.md C6, fixed).
A third, `dual_bcc_fcc_np1`, runs the two-phase case `dual_bcc_fcc` (below) on 1 rank, both
phases on one rank, against its 2-rank reference (CHANGES.md C39: v2 failed from the first
plastic increment; after the fix 7.3e-10 of field scale, 7e-4 of the tolerance).

#### New cases made with Neper (`cases/`)

Meshes: `cases/meshes/` (script `generate_meshes.sh` records the Neper and Gmsh commands and
versions). Material parameters are those of v2 cases 01 (bcc), 02 (fcc) and 03 (hcp). Uniaxial
grip loading along z at 0.01 s⁻¹.

| Case | Mesh | Ranks | Purpose |
|---|---|---|---|
| sx_fcc, sx_bcc, sx_hcp | single-crystal cube `sx2.msh` (2 cells, one orientation), 145 tets (263 with Neper 5.0.1-713-devel, D113) | 2 | single crystals of each slip family, elasto-plastic to 1 % strain |
| poly_fcc_part1/2/4 | 8-grain polycrystal, 454 tets (615 with 713), Neper partition tables for 1, 2, 4 partitions | 1, 2, 4 | partitioned meshes on the matching number of ranks |
| elastic_sx_hcp | single-crystal cube `sx2.msh` | 2 | elastic only (τ₀ = 1e10), HCP: affected by the kelas fix (CHANGES.md C1) |
| elastic_poly_fcc, elastic_poly_hcp | polycrystal, 2 partitions | 2 | elastic only (τ₀ = 1e10) |
| elastic_sx_hcp_incr1, elastic_poly_fcc_incr1, elastic_poly_hcp_incr1 | as the three elastic cases above | 2 | one increment (ε = 2.5e-4) of those cases, printing only what FEPX writes (coo, disp, forces, stress, strain_el): the v2-fixed references of the FEPX first-increment comparison (step 4.4, "FEPX against v2-fixed") |
| elastic_sx_hcp_unload | elastic_sx_hcp, targets 1e-3 then 5e-4 | 2 | unloading by strain reversal (v2's target sign); step 5.4 reference of the FEPX comparison |
| elastic_general_dstrain, elastic_general_dtime | v2 tests 65 / 68, BCC, `g_0 1e10`, on poly_part2.msh (their own meshes are msh 2.2.3, which FEPX does not read), strains scaled down 10× (targets 0.002, 0.003 with `dstrain 0.0002`; z1 velocity 0.1 with `dtime 0.002`; the same 10 + 5 increments) | 2 | `general` conditions with strain targets and dstrain, and with time targets and dtime (D44); step 5.4 references. Both end at z1 displacement 3e-3, mean σzz 586 MPa |
| sx_fcc_sat_evo, sx_fcc_precip, sx_hcp_precip_bowing | `sx2.msh`; sx_fcc / sx_hcp with saturation evolution (`m_prime 0.05`, `gammadot_s0 0.01`), precipitates with both sets (`a_p 50`, case 41 otherwise: cutting governs, τ_p = 56.6 MPa), precipitate bowing only (hcp, τ_p = 85.9 MPa) | 2 | step 7.1 references of the FEPX comparison (D59); v2's own cases 06–09 never harden (τ_s below τ₀ at their slip rates, CHANGES.md C50 note) and 40–43 have τ_p = 0.03 MPa |
| sx_fcc_aniso, sx_bcc_aniso, tight_sx_fcc_aniso, tight_sx_bcc_aniso | `sx2.msh`; sx_fcc / sx_bcc (and tight_sx_fcc / tight_sx_bcc) with `hardening anisotropic`, `interaction_matrix_parameters 1.0 1.4` (v2 tests 32/33) | 2 | step 7.4a references of the FEPX comparison; the tight ones with HDF5 output |
| dual_bcc_fcc | v2 test 05's mesh (2 grains, one phase and one partition each); fcc (phase 1) and bcc (phase 2) with case 05's parameters, loaded as sx_fcc | 2 | step 7.4a reference of the two-phase FEPX comparison (also on 1 rank); found C51 |
| sx_hcp_family, sx_bcc_112, sx_bct, sx_bcc_cyclic, sx_bcc_precip | `sx2.msh`; hcp with m 0.010 / 0.011 / 0.013 and g_0 390 / 468 / 663 MPa per family (v2 test 10); bcc with `g_0_bcc_112 240` (v2 test 44); bct with v2 test 04's parameters (no saturation evolution); bcc cyclic hardening (v2 test 16, `cyclic_a 0.002`), tension to 0.4 % and back to 0; bcc precipitates (v2 test 40, `a_p 50`) | 2 | step 7.4b references of the FEPX comparison, one per feature |
| dual_mixedhard, dual_mixedhard_np1 | v2 test 05's mesh; v2 test 46's two phases (fcc `interaction_matrix_parameters 1.0 1.4`, bcc the seven-parameter 2 × 2 block form), loaded as sx_fcc | 2, 1 | step 7.4b reference of the FEPX comparison; the 1-rank run against the 2-rank reference (C39) |
| load_target_unload | `sx2.msh`; sx_fcc's material, uniaxial grip along z at 0.01 /s, `target_load3 450 0` (load past yield, unload to zero load), `dtime 0.05`, `dtime_min 0.001 0.002` | 2 | v2-fixed C52 (first increment of a following load-target step) and C53 (`target_load2/3` read as 0); step 8.1 reference of the FEPX load-target comparison |
| load_target_general | `sx2.msh`; sx_fcc's material, `boundary_conditions general` (x0 held, x1 at 0.01 mm/s, both gripped), `target_load1 400 100`, `dtime 0.05`, `dtime_min 0.001 0.002` | 2 | step 8.1 reference of the FEPX load-target comparison with general conditions (velocities reversed by v2 for the second step) |
| load_target_tol_rel | load_target_unload with `load_tol_rel 20`, against load_target_unload's reference | 2 | v2-fixed C64 (`load_tol` and `load_tol_rel` read into the uniaxial load-target tolerance); step 8.2 |
| load_target_unload_tol | load_target_unload with `load_tol 0.5` (step 1 ends at 449.90 N, within 0.5 N short of the target, step 2 at −0.03 N) | 2 | step 8.2a: reference of `load_target_unload` (FEPX lands within `load_tol`; with 0.5 N it redoes no increment and ends where v2-fixed with the same `load_tol` ends) |
| load_target_general_tol | load_target_general with `load_tol 0.5` (399.87 N, 99.92 N) | 2 | step 8.2a: reference of `load_target_general` (as load_target_unload_tol) |
| dual_bcc_fcc_g0, dual_bcc_fcc_g0_opt | dual_bcc_fcc with g_0 190 / 230 MPa and g_s0 330 / 360 MPa (fcc / bcc): given in the phases (2 ranks), and from `simulation.opt` (one elset per phase; 1 rank, compared with `dual_bcc_fcc_g0` at the default tolerance) | 2, 1 | CHANGES.md C37 (strengths from file went to the wrong elements with several phases on a rank; without the fix 2.5e4 times the tolerance, after it identical) |
| poly_fcc_cyclic, poly_fcc_elset_opt, poly_fcc_elt_opt | `poly_part2.msh`; cyclic hardening (`cyclic_a 0.002`, `cyclic_c 3.5`) to 0.4 % and back to 0 at the same rate (v2 resets the accumulated slip); `g_0`, `g_s0`, `h_0 from_file` with per-elset values (8 elsets) and with per-element values (12 strengths per element) | 2 | step 7.1 references (D59) |
| triaxial_clr_dwell | LEGACY_triaxial_clr's unit cube and material, triaxial constant load rate along z (2 N/s), load steps to (−37.5, 0, 100), (−93.75, 0, 250), (−112.5, 0, 300) N with a dwell after steps 1 (20 s, printed) and 2 (10 s, not printed), then a reload | 2 | CHANGES.md C56 (dwell print flags swapped) and C57 (dwell placement and the load steps' data and print flags read with the merged step index); step 8.2 |
| triaxial_clr_dwell_comp | LEGACY_triaxial_clr's unit cube and material, triaxial constant load rate along z in compression (2 N/s): to (37.5, 0, −100) N, two 20 s dwells (below yield), to (140.625, 0, −375) N, two 20 s dwells (creep), unload to zero | 2 | CHANGES.md C60 (CLR perturbations without a lower bound, a zero control velocity in the scaling loop) and C61 (negative first dtime of a dwell after a dwell); step 8.2 |
| triaxial_csr_pert_floor | LEGACY_triaxial_csr (hcp unit cube, triaxial constant strain rate 1e-4 /s, compression along x to (−300, 0, 300) and (−600, 0, 600) N) with `min_pert_frac 0.01` | 2 | CHANGES.md C58 (the perturbations' lower bound was inactive for a negative control velocity); step 8.2 |
| triaxial_csr_tol | LEGACY_triaxial_csr with `load_tol_abs 0.05`, `load_tol_rel 0`, `dtime_factor 1.0`, `dt_min` 1e-6, global and CG tolerances 1e-7, 1e-13, HDF5 output | 2 | step 8.2b: reference of `triaxial_csr_tol` (tightened comparison, addendum §5 V3) |
| triaxial_clr_tol | LEGACY_triaxial_clr with `load_tol_abs 0.001`, global and CG tolerances 1e-7, 1e-13, the material-point tolerances of the step 6.6 cases (`sx_tol 1e-11`, `point_tol_state`, `point_tol_hardening` 1e-12, `sx_max_iters_state 400`), HDF5 output | 2 | step 8.2c: reference of `triaxial_clr_tol` (tightened comparison, addendum §5 V3) |
| triaxial_clr_dual | v2 test 05's two-grain unit cube and its fcc and bcc phases, LEGACY_triaxial_clr's loading (triaxial constant load rate along z, three steps) | 2 | CHANGES.md C59 (the initial face velocities' elastic estimate gave every phase phase 1's moduli); step 8.2 |
| mpc_plane_faces | v2 test 47_mpc_plane_face's bcc material and conditions (`set_bc vel` z0, z1 at 1e-2 mm/s, corners x0y0z0, x0y1z0; `set_mpc1 vel x0 x`, `set_mpc1 vel x1 x`) on `poly_part2.msh`, `target_strain33` 0.002, 0.004, 0.006, `dstrain 0.001` | 2 | step 8.3a: reference of `mpc_plane_faces` (v2's own mesh of test 47 is msh 2.2.3, which FEPX does not read, D33) |
| mpc_linear | the same conditions without the x0 MPC, elastic (`g_0 1e10`), with `set_bc vel x1y1z1 x 5e-4` on a secondary node of the x1 group (C66) and `set_mpc y1 y 0.5 2e-4` (C67), `target_strain33` 0.001, 0.002, `dstrain 0.0005` | 2 | step 8.3a: reference of `mpc_linear` |
| elastic_sx_hcp_triaxial | HCP, titanium constants, `g_0 1e10`, on `tests/verification/cases/meshes/sx2_rotx45.msh`; `set_bc strainrate 11 -0.25e-2`, `22 -0.3e-2`, `33 1e-2`, each `minimal`; `target_strain33 0.0005 0.0010` | 2 | three-component strain control (D45; after CHANGES.md C27 v2's "minimal" conditions equal the triaxial preset). v2 times the steps with the largest \|rate\| (here component 33). The rotated crystal couples normal and shear components: non-homogeneous field (σ spread ≈ 5 MPa about σzz ≈ 124 MPa), no analytic solution; step 5.4 reference |
| pbc_poly_full_state, pbc_poly_uniaxial, pbc_poly_velgrad_free, pbc_poly_plane_strain | the periodic 8-grain bcc polycrystal `meshes/poly_periodic_part2.msh` (Neper, `-periodicity all -reg 1`, 239 tets with Neper 5.0.1-731-devel (363 with 713; step 9.20, D113), 2 partitions; regularized in step 9.19, the references regenerated then: the earlier mesh, without `-reg`, 2149 tets, of the v2 comparison above, had sliver tets, median/min volume 2.8e5, that amplified rounding ~100x, so that the results differed between CPUs by up to 37 x the tolerance; v2's own mesh of tests 49/50 is msh 2.2.3); v2 tests 49/50's material; v2 test 49's six imposed rates (two steps), test 50's 11 imposed and the rest free, `velgrad` 11 = 1e-2 with 12 = 13 = 23 = 0 (C70's configuration), and the plane-strain compression of v2's periodic tutorial | 2 | step 8.3b references of the FEPX periodic comparison (V6) |
| pbc_box_uniaxial | the periodic 2 × 1 × 1 box `tests/verification/cases/meshes/periodic_box_2x1x1_part2.msh` (non-cubic, non-planar boundaries; C72; Neper `-periodicity all -reg 1`, 168 tets; regularized in step 9.20, the reference regenerated then (Release, Open MPI, 2 ranks; user's approval 2026-11-05): the earlier mesh, without `-reg`, 880 tets, of the v2 comparison above, had sliver tets, median/min volume 34363, on which `fepx_periodic_v1_full`'s nine components did not converge, and on which turning off glibc's AVX2 and FMA routines moved this case by 0.878 of its tolerance (0.0111 on the new mesh); the volume-averaged σ11 moved by 0.3 to 0.8 %, 508.0 to 506.4 MPa at the last increment); v2 test 50's loading | 2 | step 9.20 reference of `pbc_box_uniaxial` |

Two further cases use v2's 02_fcc inputs with a reduced output selection: `output_elts_only`
(only stress, strain, ori) and `output_nodes_only` (only coo, disp, vel); pristine v2 crashed when
no nodal or no elemental result was printed (CHANGES.md C9).

Notes: v2 ignores the Neper partition table and splits elements into equal contiguous blocks
itself, so the partitioned meshes do not change how v2 distributes work. Neper numbers elements by
partition, so results of the 1-, 2- and 4-partition meshes are not comparable element by element.
Sanity check of the elastic cases: elastic_poly_fcc carries 162 MPa at 0.1 % strain (apparent
modulus 162 GPa; the Voigt average for these constants is about 151 GPa, and the grips restrain the
ends laterally).


### Record: pristine v2 against v2's own committed reference data (2026-10-07)

v2's committed `ref.sim` data were compared with the pristine-v2 runs above using
`compare_sim.py` with the default tolerances (all values; v2's own `fepx --diff` compared only the
first value of each line, with looser tolerances, CHANGES.md C5). 60 of 73 cases agree; case 70
(HDF5) agrees under `h5diff -p 1e-5`. The 13 that do not:

| Case | Files differing (of total) | Nature |
|---|---|---|
| 02_fcc | 3 of 78 | `crss` files are written but not in the reference (the reference predates `print crss`); all values agree |
| 05_dual_bcc_fcc | 7 of 78 | values; `crss` has 12 values per element where the reference has 1 (disabled in v2's suite) |
| 46_mixedhard | 2 of 78 | values (disabled in v2's suite) |
| 49_pbc_triaxial_loading | 31 of 48 | values (periodic BCs) |
| 50_pbc_uniaxial_loading | 27 of 48 | values (periodic BCs) |
| 61_restart_clr_test_1 | 16 of 52 | values (CLR control) |
| 62_restart_test_2, 62_restart_csr_test_2, 62_restart_clr_test_2 | 50/50, 26/26, 52/52 | step 3 written but not in the reference (plain and CLR); values differ |
| 67, 68, 69 target_time | 52 of 78 each | values (disabled in v2's suite) |
| LEGACY_triaxial_clr | 63 of 104 | values (CLR control) |

Investigated in step 0.4 (CHANGES.md C8): the committed references are outdated. The restart chains
equal continuous runs bit for bit, the CLR runs reach their load targets, the target-time cases
equal the equivalent target-strain case, and 02/05/46/50 agree in substance. Case 49 applied shear
strain rates twice because its input gave both ij and ji components; v2 now rejects such input and
case 49's input was corrected to 6 components (CHANGES.md C15).
They do not affect the harness, which compares with reference data generated from v2-fixed itself.

