# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Verification cases: simulations checked against analytic solutions (CLAUDE.md §6), run by
# tests/regression/run_case.cmake (no reference data; CHECK_SCRIPT does the check). Description
# and tolerances: tests/verification/README.md.
#   verification_<case>     FEPX executable (label "verification"). The v2-fixed cases
#                           (verification_v2_*) were removed with v2 in step 9.13 (D90).

# Tier labels of this directory's tests, at its end (step 9.6, src/cmake/FEPXTesting.cmake).
fepx_defer_test_tiers()

# fepx_add_verification(<name> RANKS <n> MESH <file> CHECK_SCRIPT <script> CHECK_ARGS <arg>...
#                        [SOURCE <case>] [CASE_DIR <dir>])   (SOURCE: the inputs of cases/<case>,
#                        default <name>; CASE_DIR: the inputs of another directory, e.g. a
#                        regression case)
function(fepx_add_verification name)
  cmake_parse_arguments(ARG "" "RANKS;MESH;CHECK_SCRIPT;SOURCE;CASE_DIR" "CHECK_ARGS" ${ARGN})
  fepx_mesh_for_ranks(ARG_MESH ${ARG_MESH} ${ARG_RANKS})
  if(NOT ARG_SOURCE)
    set(ARG_SOURCE ${name})
  endif()
  if(NOT ARG_CASE_DIR)
    set(ARG_CASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cases/${ARG_SOURCE})
  endif()
  string(REPLACE ";" "|" check_args "${ARG_CHECK_ARGS}")
  add_test(NAME verification_${name}
           COMMAND ${CMAKE_COMMAND}
                   -DEXECUTABLE=$<TARGET_FILE:fepx>
                   -DCASE_SOURCE_DIR=${ARG_CASE_DIR}
                   -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/${name}
                   -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                   -DMPIEXEC_NUMPROC=${MPIEXEC_NUMPROC_FLAG}
                   -DNUM_RANKS=${ARG_RANKS}
                   "-DMPIEXEC_EXTRA=${FEPX_MPIEXEC_EXTRA_FLAGS_TEXT}"
                   -DMESH_FILE=${ARG_MESH}
                   -DCOMPARE=OFF
                   -DPYTHON=${Python3_EXECUTABLE}
                   -DCHECK_SCRIPT=${ARG_CHECK_SCRIPT}
                   "-DCHECK_ARGS=${check_args}"
                   -P ${FEPX_ROOT_DIR}/tests/regression/run_case.cmake)
  set_tests_properties(verification_${name} PROPERTIES
                       PROCESSORS ${ARG_RANKS}
                       TIMEOUT ${fepx_timeout_300}
                       LABELS "verification")
endfunction()

# --- FEPX -------------------------------------------------------------------------------------

# Kocks-Mecking hardening (step 7.2, D60, addendum_kocks_mecking.md §9 item 6): unirradiated bcc
# iron with Barton et al. (2013)'s parameters, slightly off [100], uniaxial at 100 /s to 4 % strain;
# exact checks (initial strengths with Hill's mu_0, crss from km_rho at every step, the network
# term not decreasing) and a report of it against Barton's Fig. 3 (approximate, not a pass/fail check).
fepx_add_verification(fepx_bcc_km_barton RANKS 2
  MESH ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/sx2_off100.msh
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_km_barton.py)

# Single-crystal uniaxial-stress elasticity against the analytic compliance (one increment to
# 1e-4 strain, uniaxial_minimal; elastic strain = total strain): HCP with c33 independent of the
# other constants (D3; titanium values, v2's constraint would give 185400), along z and along x,
# and BCT, on the 2-cell single crystal rotated 45 deg about x.
set(check_uniaxial ${CMAKE_CURRENT_SOURCE_DIR}/check_uniaxial_elastic.py)
set(mesh_rotx45 ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/sx2_rotx45.msh)
fepx_add_verification(fepx_hcp_sx_uniaxial_elastic RANKS 2 MESH ${mesh_rotx45}
  CHECK_SCRIPT ${check_uniaxial}
  CHECK_ARGS --step 1 --strain-var strain_el --c11 162400 --c12 92000 --c13 69000 --c33 180700
             --c44 46700 --axis 1 0 0 --angle 45)
# Step 9.2 (addendum_block_jacobi.md §6 T6): the same case with each preconditioner given
# explicitly (block_jacobi the default since the study, Q-B1).
foreach(precond IN ITEMS jacobi block_jacobi none)
  fepx_add_verification(fepx_hcp_sx_uniaxial_elastic_${precond} RANKS 2 MESH ${mesh_rotx45}
    CHECK_SCRIPT ${check_uniaxial}
    CHECK_ARGS --step 1 --strain-var strain_el --c11 162400 --c12 92000 --c13 69000 --c33 180700
               --c44 46700 --axis 1 0 0 --angle 45)
endforeach()
fepx_add_verification(fepx_hcp_sx_uniaxial_elastic_x RANKS 2 MESH ${mesh_rotx45}
  CHECK_SCRIPT ${check_uniaxial}
  CHECK_ARGS --step 1 --strain-var strain_el --c11 162400 --c12 92000 --c13 69000 --c33 180700
             --c44 46700 --axis 1 0 0 --angle 45 --load-dir 1)
# The copies of the input files (simulation.sim/inputs): with a tessellation simulation.tess in
# the run directory, copied with simulation.cfg and simulation.msh although FEPX does not read it.
fepx_add_verification(fepx_inputs_tess RANKS 2 MESH ${mesh_rotx45}
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_inputs.py
  CHECK_ARGS --expect simulation.cfg simulation.msh simulation.tess)
# Step 5.2: the HCP case loaded and unloaded over several increments (geometry and lattice
# updated after each); checked at the end of the unloading step.
fepx_add_verification(fepx_hcp_sx_uniaxial_elastic_multistep RANKS 2 MESH ${mesh_rotx45}
  CHECK_SCRIPT ${check_uniaxial}
  CHECK_ARGS --step 2 --strain-var strain_el --c11 162400 --c12 92000 --c13 69000 --c33 180700
             --c44 46700 --axis 1 0 0 --angle 45)
# Step 5.3: triaxial strain control (D45) on the HCP crystal with its c-axis along y (no
# normal-shear coupling in the sample frame, so the solution is homogeneous): all three
# components prescribed (11 by rates, 22 by targets, control 33, load and reverse), and a
# biaxial case with 22 free and 33 held.
# Step 8.1: load and stress targets on the elastic HCP crystal (addendum_load_targets.md §5
# V1-V2): stress 50 MPa, load -20 N, stress 0; completion, increment bounds, and the overshoot
# (f - 1)(v* - v_prev) of the aimed last increment (exact extrapolation of a linear response).
fepx_add_verification(fepx_hcp_sx_load_targets RANKS 2 MESH ${mesh_rotx45}
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_load_target.py
  CHECK_ARGS --targets stress:50 load:-20 stress:0 --factor 1.5 --dt-max 0.1 0.1 0.07 --dt-min 1e-4
             --stress-tol 2.5 --load-tol 1.0)
# Step 8.2a: landing on the targets (addendum_triaxial_load_control.md §3.6, §5 V7): the same
# case with the default tolerances 0.1 MPa and 0.1 N; each step's aimed last increment passes its
# target by more and is redone once by false position, landing on the target (elastic: linear in
# the increment length); Hooke's law at the landed step 1.
fepx_add_verification(fepx_hcp_sx_landing RANKS 2 MESH ${mesh_rotx45}
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_landing.py
  CHECK_ARGS --targets stress:50 load:-20 stress:0 --attempts 1 1 1 --dt-max 0.1 0.1 0.07)
fepx_add_verification(fepx_hcp_sx_landing_hooke RANKS 2 MESH ${mesh_rotx45}
  SOURCE fepx_hcp_sx_landing CHECK_SCRIPT ${check_uniaxial}
  CHECK_ARGS --step 1 --strain-var strain_el --c11 162400 --c12 92000 --c13 69000 --c33 180700
             --c44 46700 --axis 1 0 0 --angle 45)
# Step 8.2b: triaxial load control (addendum_triaxial_load_control.md §5 V1) on an isotropic
# elastic unit cube (lambda 1e5, mu 5e4 MPa): (1) constant strain rate on 11 to load targets, 22
# and 33 load-controlled on each step's load path (v2 CSR; LEGACY_triaxial_csr's targets), landing
# on each target; (2) 33 strain-controlled, 11 stress-controlled in time (Q-T8 (ii)), 22 free.
# Every increment on its values within the tolerance; Hooke's law and the face stresses
# (homogeneous solution) at both steps.
set(check_load_control ${CMAKE_CURRENT_SOURCE_DIR}/check_triaxial_load_control.py)
set(mesh_id ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/sx2_id.msh)
fepx_add_verification(fepx_iso_triaxial_csr RANKS 2 MESH ${mesh_id}
  CHECK_SCRIPT ${check_load_control}
  CHECK_ARGS --control 11 --modes path path --rows 11:load:-300,-600 22:load:0,0 33:load:300,600
             --load-tol 0.01 --lame 100000 50000 --hooke-steps 1 2)
fepx_add_verification(fepx_iso_triaxial_time RANKS 2 MESH ${mesh_id}
  CHECK_SCRIPT ${check_load_control}
  CHECK_ARGS --control 33 --modes time time --rows 33:load:0,0 11:stress:-50,0 --stress-tol 0.01
             --lame 100000 50000 --hooke-steps 1 2)
# Step 8.2c: triaxial load control at a load or stress rate (iteration (B), v2 CLR; addendum §5 V1,
# V2) on the same cube: (1) 33 at 100 N/s to 200 N and back to 0, 11 and 22 in time; (2) 33 at
# 50 MPa/s to 100 MPa, a 1 s dwell, back to 0 MPa, 11 by its stress rate, 22 in time (an elastic
# dwell: vanishing face velocities, v2-fixed C60); (3) the uniaxial preset at a load rate
# (regression_fepx_iso_uniaxial_clr, the same as triaxial with 11 and 22 free). Every
# increment's values (the control's too) in time within the tolerance; Hooke's law and the face
# stresses at every step.
fepx_add_verification(fepx_iso_triaxial_clr RANKS 2 MESH ${mesh_id}
  CHECK_SCRIPT ${check_load_control}
  CHECK_ARGS --control 33 --modes clr clr --rows 33:load:200,0 11:load:-100,0 22:load:50,0
             --load-tol 0.01 --lame 100000 50000 --hooke-steps 1 2 --atol 1e-9)
fepx_add_verification(fepx_iso_triaxial_clr_stress RANKS 2 MESH ${mesh_id}
  CHECK_SCRIPT ${check_load_control}
  CHECK_ARGS --control 33 --modes clr clr clr --rows 33:stress:100,h,0 11:stress:r-10,r0,r10
             22:stress:20,20,0 --stress-tol 0.01 --lame 100000 50000 --hooke-steps 1 2 3
             --atol 1e-9)
fepx_add_verification(fepx_iso_uniaxial_clr RANKS 2 MESH ${mesh_id}
  CASE_DIR ${FEPX_ROOT_DIR}/tests/regression/cases/iso_uniaxial_clr
  CHECK_SCRIPT ${check_load_control}
  CHECK_ARGS --control 33 --modes clr clr --rows 33:load:100,20 --load-tol 0.01
             --lame 100000 50000 --hooke-steps 1 2)
# Step 8.5 (D86, addendum_triaxial_load_control.md §3.4 revised): at a step change of triaxial
# control with the control component at a strain rate in both steps, the previous face velocities,
# the load-controlled ones included, times the ratio of its strain rates (v2 vel_factor): a strain
# rate jump (x 2) and a reversal (x -1) on triaxial_csr's hcp cube; the first
# boundary-velocity iteration of step 2 against the last of step 1.
foreach(variant_ratio jump:2 reversal:-1)
  string(REPLACE ":" ";" variant_ratio "${variant_ratio}")
  list(GET variant_ratio 0 variant)
  list(GET variant_ratio 1 ratio)
  fepx_add_verification(fepx_triaxial_csr_${variant} RANKS 2
    MESH ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/triaxial_cube.msh
    CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_step_change_velocities.py
    CHECK_ARGS --control 11 --load-comps 22 33 --ratio ${ratio})
endforeach()
# Step 8.2c: creep in a dwell above yield (addendum §5 V4): fcc single crystal, uniaxial_grip at
# 100 N/s to 440 N (past yield), dwells of 2 s and 1 s (max_strain_incr 2e-4), unloading: the load
# held within load_tol, the face creeping on, v2's dwell increments (13 + 4, strain-limited).
fepx_add_verification(fepx_fcc_sx_creep_dwell RANKS 2
  MESH ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/sx2.msh
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_dwell.py
  CHECK_ARGS --component 33 --dwells 2:0.5:6.4 3:0.5:7.4 --episode-starts 2
             --max-strain-incr 2e-4 --length0 1.0 --load-tol 0.1)
# ... and every control mode in one run (V5, regression_fepx_triaxial_mixed_clr_sx_fcc):
# stress 11 in time during a strain step, then 33 unloaded at a stress rate with 11 in time, a
# dwell with 11 held by its rate 0, and the reload at a strain rate on the load path.
fepx_add_verification(fepx_triaxial_mixed_clr_sx_fcc RANKS 2
  MESH ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/sx2.msh
  CASE_DIR ${FEPX_ROOT_DIR}/tests/regression/cases/triaxial_mixed_clr_sx_fcc
  CHECK_SCRIPT ${check_load_control}
  CHECK_ARGS --control 33 --modes time clr clr path --rows 33:stress:0,0,h,300
             11:stress:-20,0,r0,-30 --stress-tol 0.05 --lame 1 1)
# ... and the crystal-plasticity polycrystal of regression_fepx_triaxial_mixed_poly_fcc
# (V5, the parts available in 8.2b): stress 11 in time during a strain step, then on the load
# path of two stress-target steps (unloading, reloading past yield), within stress_tol.
fepx_add_verification(fepx_triaxial_mixed_poly_fcc RANKS 2
  MESH ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_part2.msh
  CASE_DIR ${FEPX_ROOT_DIR}/tests/regression/cases/triaxial_mixed_poly_fcc
  CHECK_SCRIPT ${check_load_control}
  CHECK_ARGS --control 33 --modes time path path --rows 33:stress:0,100,375
             11:stress:-50,0,-100 --stress-tol 0.01 --lame 1 1)
set(check_triaxial ${CMAKE_CURRENT_SOURCE_DIR}/check_triaxial_elastic.py)
set(mesh_rotx90 ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/sx2_rotx90.msh)
fepx_add_verification(fepx_hcp_sx_triaxial_elastic RANKS 2 MESH ${mesh_rotx90}
  CHECK_SCRIPT ${check_triaxial}
  CHECK_ARGS --step 2 --c11 162400 --c12 92000 --c13 69000 --c33 180700 --c44 46700
             --axis 1 0 0 --angle 90 --prescribed 11 22 33 --expect=-2e-4,-3e-4,1e-3)
fepx_add_verification(fepx_hcp_sx_triaxial_elastic_free RANKS 2 MESH ${mesh_rotx90}
  CHECK_SCRIPT ${check_triaxial}
  CHECK_ARGS --step 1 --c11 162400 --c12 92000 --c13 69000 --c33 180700 --c44 46700
             --axis 1 0 0 --angle 90 --prescribed 11 33 --expect=1.5e-3,0,0)
fepx_add_verification(fepx_bct_sx_uniaxial_elastic RANKS 2 MESH ${mesh_rotx45}
  CHECK_SCRIPT ${check_uniaxial}
  CHECK_ARGS --step 1 --strain-var strain_el --c11 72.3e3 --c12 59.4e3 --c13 35.8e3 --c33 88.4e3
             --c44 22.0e3 --c66 24.0e3 --axis 1 0 0 --angle 45)

# Step 6.2: crystal plasticity with constant strengths. fcc along [001] (sx2_id.msh): the 8
# systems with Schmid factor 1/sqrt(6) give the steady power-law flow stress at 1 % strain. hcp
# with one m per family and prismatic slip active (the FEPX form of the C11 case): the stress
# must stay uniaxial.
fepx_add_verification(fepx_fcc_sx001_flow_stress RANKS 2
  MESH ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/sx2_id.msh
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_flow_stress.py
  CHECK_ARGS --step 1 --crss 100 --m 0.02 --gammadot0 1 --rate 1e-3 --strain 0.01
             --schmid 0.40824829046386301 --num-systems 8 --tol 1e-5)
# Step 9.2 (addendum_block_jacobi.md §6 T6): the same case with the Jacobi and the block-Jacobi
# preconditioner given explicitly (block_jacobi the default since the study, Q-B1).
foreach(precond IN ITEMS jacobi block_jacobi)
  fepx_add_verification(fepx_fcc_sx001_flow_stress_${precond} RANKS 2
    MESH ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/sx2_id.msh
    CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_flow_stress.py
    CHECK_ARGS --step 1 --crss 100 --m 0.02 --gammadot0 1 --rate 1e-3 --strain 0.01
               --schmid 0.40824829046386301 --num-systems 8 --tol 1e-5)
endforeach()
# hcp with c33 independent (step 6.6; FEPX v2 cannot represent it, D3): c-axis loading, the
# elastic modulus 1/S33 and the steady flow stress of the 12 pyramidal <c+a> systems.
fepx_add_verification(fepx_hcp_sx_c_axis_c33 RANKS 2
  MESH ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/sx2_id.msh
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_hcp_c_axis.py
  CHECK_ARGS --c11 162.4e3 --c12 92.0e3 --c13 69.0e3 --c33 180.7e3 --c44 46.7e3 --c-over-a 1.587
             --elastic-step 1 --elastic-strain 1e-5 --flow-step 2 --flow-strain 0.01 --crss 300
             --m 0.02 --gammadot0 1 --rate 1e-3)
fepx_add_verification(fepx_hcp_sx_uniaxial_aniso_m RANKS 2 MESH ${mesh_rotx90}
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_uniaxial_stress.py
  CHECK_ARGS --step 1 --tol 1e-6)

# Step 6.5: the crystal-plasticity outputs of a polycrystal with isotropic hardening, consistent
# with each other (slip, plastic strain, slip rates, strengths, rotation rate), and the work of
# the homogeneous fcc [001] crystal against the trapezoidal rule of sigma_zz d_zz.
fepx_add_verification(fepx_fcc_poly_outputs RANKS 2
  MESH ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_part1.msh
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_cp_outputs.py
  CHECK_ARGS --times 0.2 0.3 0.5 0.8 1.2 --first-step=2 --m=0.05 --gammadot0=1 --h0=200
             --tau-s=330 --tau-0=210 --rotrate)
# Step 9.6 (D88): the Kocks-Mecking and irradiation outputs (km_k2, km_rhorate, km_rho, dmg,
# dmgrate) against the laws and the backward-Euler updates (measured 2026-10-24: the laws to
# 4e-16, the updates to 2.2e-5 of the increment at default tolerances, limit 1e-4; with the rates
# of the increment's start, i.e. forward Euler, the updates would differ by 6 to 62 %).
fepx_add_verification(fepx_km_irr_outputs RANKS 2
  MESH ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/sx2.msh
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_km_irr_outputs.py
  CHECK_ARGS --times 0.3 0.35 0.4 0.45 0.5 --first-step=2 --k1=4.5e5 --k20=14
             --gammadot-k0=1e10 --n-k=100)
fepx_add_verification(fepx_fcc_sx001_work RANKS 2
  MESH ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/sx2_id.msh
  CHECK_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/check_work.py
  CHECK_ARGS --times 2 4 7 10 --rate=1e-3 --first-step=2)

# --- Rank independence ------------------------------------------------------------------------

# FEPX rank independence with per-element values from simulation.opt (step 7.1, D59): the case of
# regression_fepx_poly_fcc_elt_opt (12 strengths per element, per-element tau_s and h_0) on
# 1, 2 and 4 ranks (even partitioning); the element stresses must agree to 1e-7 relative.
add_test(NAME verification_fepx_rank_independence_elt_opt
         COMMAND ${CMAKE_COMMAND}
                 -DEXECUTABLE=$<TARGET_FILE:fepx>
                 -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cases/fepx_fcc_elt_opt_ranks
                 -DMESH=${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_part2.msh
                 -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/fepx_rank_independence_elt_opt
                 -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                 -DMPIEXEC_NUMPROC=${MPIEXEC_NUMPROC_FLAG}
                 "-DMPIEXEC_EXTRA=${FEPX_MPIEXEC_EXTRA_FLAGS_TEXT}"
                 "-DRANKS=1|2|4"
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCHECK_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/check_rank_independence.py
                 -DRTOL=1e-7
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/rank_independence.cmake)
set_tests_properties(verification_fepx_rank_independence_elt_opt PROPERTIES
                     PROCESSORS 4
                     TIMEOUT ${fepx_timeout_600}
                     LABELS "verification;fepx")

# Periodic conditions (step 8.3b, addendum_periodic_mpc.md §5): isotropic elastic bodies on Neper
# periodic meshes, the 2 x 1 x 1 box (non-planar "raspberry" boundaries, non-cubic: V3) and the
# one-cell cube (planar faces), checked by check_periodic.py (relations on the velocities and
# coordinates, homogeneity and the volume-averaged stress, imposed and free components).
set(check_periodic ${CMAKE_CURRENT_SOURCE_DIR}/check_periodic.py)
set(mesh_periodic_box ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/periodic_box_2x1x1.msh)
set(mesh_periodic_cube ${CMAKE_CURRENT_SOURCE_DIR}/cases/meshes/periodic_cube.msh)
# V1: small strain, all six components imposed, two steps with different rates; C : e.
fepx_add_verification(fepx_periodic_v1_small RANKS 2 MESH ${mesh_periodic_box}
  CHECK_SCRIPT ${check_periodic}
  CHECK_ARGS --steps 1 2 --homogeneous 1e-9 --symmetric --elastic 200000 100000 --erel 2e-4)
# Step 9.6 (D88): shear control components (no run used one before). Isotropic elastic cube,
# two steps with a reversal; small_strain controlled by 23, 13, 12 (pure shear, the rest free),
# deformation_gradient_rate by 32, 31, 21 (simple shear: the transposed component and one of
# each other shear pair imposed at zero, the rest free). Checks: relations and homogeneity to
# round-off; the free components' nominal stress P = J sigma F^-T (conjugate to Fbar-dot; the
# Cauchy components of free normal directions are the shear stress times the strain, not 0) at
# most 1e-9 (measured 3e-16); stress = C : e within 2e-4 as V1 (measured 1.5e-5 to 5.3e-5: the
# finite-strain terms, O(strain) at 2e-5).
foreach(case IN ITEMS "23;small" "13;small" "12;small"
                      "32;full;23=0 21=0 31=0" "31;full;13=0 21=0 32=0" "21;full;12=0 31=0 32=0")
  list(GET case 0 comp)
  list(GET case 1 kinematics)
  # The free components: the three normal ones and the two other shear components.
  set(free_components 11 22 33 23 13 12)
  foreach(shear IN ITEMS 23 13 12)
    string(SUBSTRING ${comp} 1 1 second)
    string(SUBSTRING ${comp} 0 1 first)
    if(shear STREQUAL comp OR shear STREQUAL "${second}${first}")
      list(REMOVE_ITEM free_components ${shear})
    endif()
  endforeach()
  if(kinematics STREQUAL "small")
    set(kinematics_args --symmetric)
  else()
    list(GET case 2 zeros)
    separate_arguments(zeros)
    set(kinematics_args --imposed ${zeros})
  endif()
  fepx_add_verification(fepx_periodic_shear_${comp} RANKS 2 MESH ${mesh_periodic_cube}
    CHECK_SCRIPT ${check_periodic}
    CHECK_ARGS --steps 1 2 --homogeneous 1e-9 --free-nominal ${free_components} --stol 1e-9
               --elastic 200000 100000 --erel 2e-4 --rtol 1e-9 ${kinematics_args})
endforeach()
# V1: the nine components of Fbar-dot imposed; element velocity gradients = Fbar-dot Fbar^-1.
fepx_add_verification(fepx_periodic_v1_full RANKS 2 MESH ${mesh_periodic_cube}
  CHECK_SCRIPT ${check_periodic}
  CHECK_ARGS --steps 1 --homogeneous 1e-9 --velgrad --rtol 1e-9
             --imposed 11=1.0e-3 12=4.0e-4 13=-2.0e-4 21=-1.0e-4 22=-3.0e-4 23=2.5e-4
             31=3.0e-4 32=-1.5e-4 33=5.0e-4)
# The same on the 2 x 1 x 1 box (step 9.20): the multi-cell case, non-planar boundaries.
fepx_add_verification(fepx_periodic_v1_full_box RANKS 2 MESH ${mesh_periodic_box}
  SOURCE fepx_periodic_v1_full CHECK_SCRIPT ${check_periodic}
  CHECK_ARGS --steps 1 --homogeneous 1e-9 --velgrad --rtol 1e-9
             --imposed 11=1.0e-3 12=4.0e-4 13=-2.0e-4 21=-1.0e-4 22=-3.0e-4 23=2.5e-4
             31=3.0e-4 32=-1.5e-4 33=5.0e-4)
# V2: only Fbar-dot 11 imposed (small strain; and with 12 = 13 = 23 = 0 imposed with the nine
# components, the exact diagonal): uniaxial stress, free components zero, lateral rates -nu L11.
fepx_add_verification(fepx_periodic_v2_small RANKS 2 MESH ${mesh_periodic_box}
  CHECK_SCRIPT ${check_periodic}
  CHECK_ARGS --steps 1 --homogeneous 1e-7 --symmetric --free 22 33 23 13 12 --stol 1e-8
             --poisson 0.3333333333333333 --rtol 1e-9)
fepx_add_verification(fepx_periodic_v2_full RANKS 2 MESH ${mesh_periodic_cube}
  CHECK_SCRIPT ${check_periodic}
  CHECK_ARGS --steps 1 --homogeneous 1e-7 --free 22 33 23 13 12 --stol 1e-8
             --poisson 0.3333333333333333 --imposed 11=1.0e-3 --rtol 1e-9)
# Q-P6: a velocity condition on an independent dof (all three components of node 2).
fepx_add_verification(fepx_periodic_vel_condition RANKS 2 MESH ${mesh_periodic_box}
  CHECK_SCRIPT ${check_periodic}
  CHECK_ARGS --steps 1 --homogeneous 1e-9 --symmetric --rtol 1e-9)

# FEPX rank independence with periodic conditions (step 8.3b, V5): a periodic 8-grain
# polycrystal (Neper periodic tessellation) with crystal plasticity on 1, 2 and 4 ranks (even
# partitioning: the relations cross the partitions); the element stresses must agree to 1e-7.
# Each run's limit scales with the build (a Debug 1-rank run took about 400 s on the mesh before
# step 9.19; on the regularized mesh, 363 tets, the three runs took about 55 s in Debug; 239 tets since 9.20).
add_test(NAME verification_fepx_rank_independence_periodic
         COMMAND ${CMAKE_COMMAND}
                 -DEXECUTABLE=$<TARGET_FILE:fepx>
                 -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cases/fepx_periodic_ranks
                 -DMESH=${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_periodic_part2.msh
                 -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/fepx_rank_independence_periodic
                 -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                 -DMPIEXEC_NUMPROC=${MPIEXEC_NUMPROC_FLAG}
                 "-DMPIEXEC_EXTRA=${FEPX_MPIEXEC_EXTRA_FLAGS_TEXT}"
                 "-DRANKS=1|2|4"
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCHECK_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/check_rank_independence.py
                 -DRTOL=1e-7
                 -DRUN_TIMEOUT=${fepx_timeout_300}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/rank_independence.cmake)
set_tests_properties(verification_fepx_rank_independence_periodic PROPERTIES
                     PROCESSORS 4
                     TIMEOUT ${fepx_timeout_600}
                     LABELS "verification;fepx")

# FEPX rank independence with multi-point constraints (step 8.3a, V5): the case of
# regression_fepx_mpc_plane_faces (set_mpc1 on x0 and x1, crystal plasticity) on 1, 2 and 4
# ranks (even partitioning: the MPC groups cross the partitions); the element stresses must
# agree to 1e-7 relative.
add_test(NAME verification_fepx_rank_independence_mpc
         COMMAND ${CMAKE_COMMAND}
                 -DEXECUTABLE=$<TARGET_FILE:fepx>
                 -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cases/fepx_mpc_ranks
                 -DMESH=${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_part2.msh
                 -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/fepx_rank_independence_mpc
                 -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                 -DMPIEXEC_NUMPROC=${MPIEXEC_NUMPROC_FLAG}
                 "-DMPIEXEC_EXTRA=${FEPX_MPIEXEC_EXTRA_FLAGS_TEXT}"
                 "-DRANKS=1|2|4"
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCHECK_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/check_rank_independence.py
                 -DRTOL=1e-7
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/rank_independence.cmake)
set_tests_properties(verification_fepx_rank_independence_mpc PROPERTIES
                     PROCESSORS 4
                     TIMEOUT ${fepx_timeout_600}
                     LABELS "verification;fepx")

# Restart (step 8.4, dev/theory/addendum_restart.md §8): each case's whole history (full.cfg)
# and the same history interrupted and continued from restart files (part1.cfg, part2.cfg, ...,
# one after the other, each later part with restart on and only its own steps) must leave
# byte-identical simulation.sim directories (restart_chain.cmake; inputs/ and restart/ excluded):
# monotonic tension (restart_monotonic); a reversal at the restart with cyclic hardening and
# dtime/dstrain steps (restart_reversal); an unload to a stress target after the restart
# (restart_load_unload); three runs with step 1 not printed, Kocks-Mecking and
# sub-incrementation (restart_three_parts); a load rate with dwells, restarted within the dwell
# episode (restart_dwell); triaxial load control with stress rates and a dwell in three runs
# (restart_triaxial); periodic conditions reversed at the restart (restart_periodic); general
# conditions with MPCs and the isotropic viscoplastic start (restart_mpc).
function(fepx_add_restart_chain name)
  cmake_parse_arguments(ARG "" "RANKS;MESH" "" ${ARGN})
  fepx_mesh_for_ranks(ARG_MESH ${ARG_MESH} ${ARG_RANKS})
  add_test(NAME verification_fepx_${name}
           COMMAND ${CMAKE_COMMAND}
                   -DEXECUTABLE=$<TARGET_FILE:fepx>
                   -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cases/${name}
                   -DMESH=${ARG_MESH}
                   -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/fepx_${name}
                   -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                   -DMPIEXEC_NUMPROC=${MPIEXEC_NUMPROC_FLAG}
                   "-DMPIEXEC_EXTRA=${FEPX_MPIEXEC_EXTRA_FLAGS_TEXT}"
                   -DNUM_RANKS=${ARG_RANKS}
                   -DRUN_TIMEOUT=${fepx_timeout_300}
                   -DH5DIFF=${FEPX_H5DIFF}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/restart_chain.cmake)
  set_tests_properties(verification_fepx_${name} PROPERTIES
                       PROCESSORS ${ARG_RANKS}
                       TIMEOUT ${fepx_timeout_600}
                       LABELS "verification;fepx")
endfunction()
set(mesh_poly2 ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_part2.msh)
set(mesh_sx2 ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/sx2.msh)
fepx_add_restart_chain(restart_monotonic RANKS 2 MESH ${mesh_poly2})
fepx_add_restart_chain(restart_reversal RANKS 2 MESH ${mesh_poly2})
fepx_add_restart_chain(restart_load_unload RANKS 2 MESH ${mesh_poly2})
fepx_add_restart_chain(restart_three_parts RANKS 2 MESH ${mesh_poly2})
fepx_add_restart_chain(restart_dwell RANKS 2 MESH ${mesh_sx2})
fepx_add_restart_chain(restart_triaxial RANKS 2 MESH ${mesh_sx2})
fepx_add_restart_chain(restart_periodic RANKS 2
  MESH ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_periodic_part2.msh)
fepx_add_restart_chain(restart_mpc RANKS 2 MESH ${mesh_poly2})

# HDF5 output (step 9.1, ROADMAP D93; builds with FEPX_ENABLE_HDF5): every output variable in
# both formats on the fcc + hcp cube (12 and 18 slip systems, a step not printed, 4 ranks), the
# HDF5 file against the text files bit for bit, the mesh and the XDMF file
# (check_hdf5_output.py); the same with hdf5 alone, its HDF5 and XDMF files identical to those
# (h5diff); periodic conditions (the macroscopic table); restart chains with both formats
# (restart_monotonic's: forces, convergence) and with hdf5 alone (restart_periodic's: macro).
if(FEPX_HAVE_HDF5)
  set(mesh_dual ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/dual_bcc_fcc.msh)
  set(check_hdf5 ${CMAKE_CURRENT_SOURCE_DIR}/check_hdf5_output.py)
  fepx_add_verification(fepx_output_hdf5 RANKS 4 MESH ${mesh_dual}
    CHECK_SCRIPT ${check_hdf5} CHECK_ARGS --mesh ${mesh_dual})
  fepx_add_verification(fepx_output_hdf5_only RANKS 4 MESH ${mesh_dual}
    CHECK_SCRIPT ${check_hdf5}
    CHECK_ARGS --hdf5-only --reference ${CMAKE_CURRENT_BINARY_DIR}/work/fepx_output_hdf5/simulation.sim)
  set_tests_properties(verification_fepx_output_hdf5 PROPERTIES FIXTURES_SETUP output_hdf5)
  set_tests_properties(verification_fepx_output_hdf5_only PROPERTIES
                       FIXTURES_REQUIRED output_hdf5)
  fepx_add_verification(fepx_output_hdf5_periodic RANKS 2
    MESH ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_periodic_part2.msh
    CHECK_SCRIPT ${check_hdf5}
    CHECK_ARGS --mesh ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_periodic_part2.msh)
  fepx_add_restart_chain(restart_hdf5 RANKS 2 MESH ${mesh_poly2})
  fepx_add_restart_chain(restart_hdf5_periodic RANKS 2
    MESH ${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_periodic_part2.msh)
endif()

# FEPX rank independence with Kocks-Mecking hardening and per-elset initial densities from
# simulation.opt (step 7.2, D60): the case of regression_fepx_km_poly_fcc_rho_opt at
# default tolerances on 1, 2 and 4 ranks (even partitioning); the element stresses must agree to
# 1e-7 relative.
add_test(NAME verification_fepx_rank_independence_km
         COMMAND ${CMAKE_COMMAND}
                 -DEXECUTABLE=$<TARGET_FILE:fepx>
                 -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cases/fepx_fcc_km_ranks
                 -DMESH=${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_part2.msh
                 -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/fepx_rank_independence_km
                 -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                 -DMPIEXEC_NUMPROC=${MPIEXEC_NUMPROC_FLAG}
                 "-DMPIEXEC_EXTRA=${FEPX_MPIEXEC_EXTRA_FLAGS_TEXT}"
                 "-DRANKS=1|2|4"
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCHECK_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/check_rank_independence.py
                 -DRTOL=1e-7
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/rank_independence.cmake)
set_tests_properties(verification_fepx_rank_independence_km PROPERTIES
                     PROCESSORS 4
                     TIMEOUT ${fepx_timeout_600}
                     LABELS "verification;fepx")

# FEPX rank independence with triaxial load control (step 8.2b, V5): the case of
# regression_fepx_triaxial_mixed_poly_fcc on 1, 2 and 4 ranks (even partitioning); the
# element stresses must agree to 1e-6 relative (the boundary-velocity iteration stops anywhere
# within stress_tol 0.01 MPa of ~400 MPa: 2.5e-5 at most; measured 6.8e-13). Each run may take
# 300 s (Debug: 1200 s; the 1-rank Debug run exceeds the script's default 240 s).
add_test(NAME verification_fepx_rank_independence_triaxial
         COMMAND ${CMAKE_COMMAND}
                 -DEXECUTABLE=$<TARGET_FILE:fepx>
                 -DCASE_DIR=${FEPX_ROOT_DIR}/tests/regression/cases/triaxial_mixed_poly_fcc
                 -DMESH=${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_part2.msh
                 -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/fepx_rank_independence_triaxial
                 -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                 -DMPIEXEC_NUMPROC=${MPIEXEC_NUMPROC_FLAG}
                 "-DMPIEXEC_EXTRA=${FEPX_MPIEXEC_EXTRA_FLAGS_TEXT}"
                 "-DRANKS=1|2|4"
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCHECK_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/check_rank_independence.py
                 -DRTOL=1e-6
                 -DRUN_TIMEOUT=${fepx_timeout_300}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/rank_independence.cmake)
set_tests_properties(verification_fepx_rank_independence_triaxial PROPERTIES
                     PROCESSORS 4
                     TIMEOUT ${fepx_timeout_600}
                     LABELS "verification;fepx")

# ... and with iteration (B) (step 8.2c, V5): the case of regression_fepx_triaxial_mixed_clr_
# sx_fcc on 1, 2 and 4 ranks (even partitioning); the element stresses must agree to 1e-6 relative
# (each iteration stops anywhere within stress_tol 0.05 MPa of ~300 MPa: 1.7e-4 at most).
add_test(NAME verification_fepx_rank_independence_triaxial_clr
         COMMAND ${CMAKE_COMMAND}
                 -DEXECUTABLE=$<TARGET_FILE:fepx>
                 -DCASE_DIR=${FEPX_ROOT_DIR}/tests/regression/cases/triaxial_mixed_clr_sx_fcc
                 -DMESH=${FEPX_ROOT_DIR}/tests/regression/cases/meshes/sx2.msh
                 -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/fepx_rank_independence_triaxial_clr
                 -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                 -DMPIEXEC_NUMPROC=${MPIEXEC_NUMPROC_FLAG}
                 "-DMPIEXEC_EXTRA=${FEPX_MPIEXEC_EXTRA_FLAGS_TEXT}"
                 "-DRANKS=1|2|4"
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCHECK_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/check_rank_independence.py
                 -DRTOL=1e-6
                 -DRUN_TIMEOUT=${fepx_timeout_300}
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/rank_independence.cmake)
set_tests_properties(verification_fepx_rank_independence_triaxial_clr PROPERTIES
                     PROCESSORS 4
                     TIMEOUT ${fepx_timeout_600}
                     LABELS "verification;fepx")

# The terminal output (step 9.7, ROADMAP D94, dev/plan/terminal_output.md §3.8): the whole
# standard output of a short run on 2 ranks (fcc crystal plasticity, two steps of two
# increments, a mesh of 2 partitions) compared with its reference after masking
# the variable fields (version, build, wall times, the global-iteration lines and counts;
# check_terminal_output.py). Reference regenerated with FEPX_REGRESSION_UPDATE=1, only under
# CLAUDE.md §2.5.
fepx_mpiexec_command(terminal_output_command 2 $<TARGET_FILE:fepx>)
add_test(NAME verification_fepx_terminal_output
         COMMAND ${CMAKE_COMMAND}
                 "-DCOMMAND=${terminal_output_command}"
                 -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cases/terminal_output
                 -DMESH=${FEPX_ROOT_DIR}/tests/regression/cases/meshes/poly_part2.msh
                 -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/work/fepx_terminal_output
                 -DPYTHON=${Python3_EXECUTABLE}
                 -DCHECK_SCRIPT=${CMAKE_CURRENT_SOURCE_DIR}/check_terminal_output.py
                 -DREFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/cases/terminal_output/reference_stdout.txt
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/terminal_output.cmake)
set_tests_properties(verification_fepx_terminal_output PROPERTIES
                     PROCESSORS 2
                     TIMEOUT ${fepx_timeout_300}
                     LABELS "verification;fepx")
