! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> FEPX executable, step 5.2 scope: an elastic run of one or more loading steps, each of one
!> or more increments, from a configuration file and a mesh, with the results written to
!> simulation.sim (step 0: initial state; then the end of each printed step). Each phase has an
!> elastic model or a crystal plasticity model (step 6.2) with its hardening law (step 6.4, D51:
!> saturation_isotropic, saturation_anisotropic, or none for constant strengths) and its
!> modifiers (step 7.1, D59: saturation_evolution, cyclic_hardening, precipitate_strengthening).
!> With read_from_opt_file, per-elset or per-element initial strengths, saturation strengths
!> and hardening coefficients from simulation.opt replace the phase's values (step 7.1, D59;
!> strengths from the file are used without the precipitate strength, with a warning; values of
!> elastic phases are ignored, with a warning). Cyclic hardening's accumulated slip is reset at
!> the start of every step that reverses the control strain rate (D59).
!>
!> Initialization (ROADMAP D49, D54, D78: initialization elastic_step, the default): a run with a
!> crystal
!> plasticity phase starts with a small elastic step, part of step 1 and not written: the
!> strain that gives about 1/100 of the expected yield stress 3 tau_0 (the smallest tau_0 of
!> the crystal-plasticity phases, a Taylor factor of about 3) with the largest stiffness, bounded
!> by the largest absolute row sum of [C] (Gershgorin: an upper bound of its eigenvalues), at the
!> fastest strain rate of step 1 (general conditions: the largest |set_bc velocity| over the
!> smallest bounding-box length); if that is not shorter than the first increment, half the
!> first increment. The first increment is shortened by it; its velocities are the initial
!> guess of the first increment. initialization isotropic_viscoplastic (D78, v2's start,
!> dev/theory/addendum_isotropic_viscoplastic.md): no elastic step; after step 0 is written, v2's
!> isotropic viscoplastic solve (isotropic_start_mod: a linear incompressible viscous problem,
!> three augmented-Lagrangian iterations, with the step-1 conditions and constraints) gives the
!> velocity guess of the first increment, and every point solve of increment 1 starts from
!> 0.5 x v2's viscoplastic estimate (D54 (3)), started from the yield-surface vertices computed
!> at setup from each phase's active systems and tau_0 (yield_vertices_compute; a fatal error
!> when they span fewer than 5 deviatoric dimensions, D57).
!>
!> Files (working directory, as v2): simulation.cfg, simulation.msh, and, when the configuration
!> asks for them, simulation.ori (read_ori_from_file), simulation.phase (read_phase_from_file),
!> simulation.opt (read_from_opt_file). Output: simulation.sim (sim_writer_mod).
!>
!> Source : FEPX (ROADMAP steps 4.4, 5.1 and 5.2, D43, D44): the sequence of v2 fepx.f90 and
!>          driver_uniaxial_control_mod1.f90 (read input, distribute the mesh, set up the
!>          material, then the time-stepping loop: boundary conditions of the step, increment
!>          solve, geometry update x = x + dt v (v2 solveit_evp_mod2.f90::update_state_evp),
!>          output), with v2's loading velocity (strain rate times the original gage length along
!>          the loading axis; triaxial: along each prescribed axis, D45; general: the set_bc
!>          velocities times the step's velocity factor),
!>          the increments of the step schedule (step_schedule_mod), v2's initial velocity
!>          guess (the previous increment's velocities, scaled at a step change by the ratio of
!>          the step strain rates, v2 vel_factor), and the run limits of D44 (v2-fixed C32).
!>          Contains no numerical method of its own.
!>
!> Steps (step 8.1, dev/theory/addendum_load_targets.md §3): each step ends at a strain, load
!> (force on the face x1/y1/z1 of the control component), stress (that force over the face's
!> current area) or time target. Steps are planned from their actual start (step_schedule_mod):
!> before the run up to the first load or stress target, the later ones when they start. A load
!> or stress target step starts with its direction checked from the current value
!> (step_schedule_load_start), takes increments of at most dtime (dstrain), the next one chosen
!> by v2's linear extrapolation within [dtime_min, dtime], and ends after the first increment
!> that reaches its target (v2's uniaxial load target, v2-fixed C52), landing within load_tol
!> (stress_tol) of it (step 8.2a, dev/theory/addendum_triaxial_load_control.md §3.6, Q-T13,
!> Q-T14): an increment that passes the target by more is discarded (coordinates, velocities
!> and state restored, as a rejected sub-incrementation attempt) and redone with the length
!> given by false position (Press et al. 2007, Numerical Recipes 3rd ed., §9.2;
!> step_load_false_position) on the bracket of the lengths and values reached, at most
!> max_bc_iter times, then a fatal error after writing the last accepted increment as the
!> current step; dtime_min and dtime_factor do not apply to these attempts, and only the
!> accepted one counts in the time, the increments and the output. With check_necking on,
!> after every increment, a control-face force that changes against the step's direction
!> stops the run after writing that increment as the current step (v2-fixed C54).
!>
!> Multi-point constraints (step 8.3a, D73, dev/theory/addendum_periodic_mpc.md §3.2): the
!> set_mpc1 and set_mpc lines give, for every node of a node set but the first, v = c v_first + d
!> in a direction (constraint_mod, built once after the mesh). At every step start the velocity
!> conditions are carried through the groups (constraint_propagate_bcs, v2-fixed C66) and the
!> initial velocity guess is projected onto the constraints (v2's uniaxial driver); every
!> increment solves the reduced system (increment_mod, pcg_mod with the constraints).
!>
!> Periodic conditions (step 8.3b, D71-D73, dev/theory/addendum_periodic_mpc.md §3.3): with
!> boundary_conditions periodic the periodicity relations of the mesh give v = T v_hat + g
!> (periodic_mod, built once after the mesh: drive pairs, lengths, the replicated table); at every
!> step start the set_bc vel conditions are applied as given (an error on a constrained dof), the
!> offsets of the step's imposed rates of Fbar-dot are set, and the initial guess is projected
!> onto the relations and shifted to the zero-mean fluctuation; every increment solves the
!> reduced system with the zero-mean shift of each correction (increment_mod). The steps are
!> those of the control component (its strain is a component of Fbar - I); after every
!> increment Fbar, Fbar-dot and the volume-averaged stress are written to results/macro.
!>
!> Triaxial load control (step 8.2b, dev/theory/addendum_triaxial_load_control.md §3.1, §3.4,
!> iteration (A); decisions Q-T4, Q-T8 (ii), Q-T10, Q-T12, Q-T18; ROADMAP D66): a component other
!> than the control component with a target_load or target_stress row in a step is
!> load-controlled: its face x1/y1/z1 moves at a uniform normal velocity found in every increment
!> attempt (and sub-increment) by the boundary-velocity iteration of v2's triaxial CSR driver
!> (bc_iterate; the pure parts in bc_iteration_mod): a solve at the current face velocities,
!> its load (stress) compared with the value to reach (on the step's load path at the control's
!> progress when the control has a load or stress target, v2 CSR; interpolated in time when it
!> has a strain or time target, Q-T8 (ii)), within load_tol (stress_tol) the solve is the
!> attempt; otherwise one solve with each load-controlled face velocity perturbed, v2's linear
!> system for the correction, and again, at most max_bc_iter iterations (then a fatal error).
!> Every solve starts from the start of the attempt; its face forces are taken on the
!> end-of-increment geometry (Q-T10); a failed solve rejects the attempt as a failed increment
!> solve does (sub-incrementation). Face-velocity changes are spread over the body as v2's
!> homogeneous field Delta V (x - x_min)/L along the component. The face velocities are carried
!> from increment to increment and from step to step; a run whose step 1 has load-controlled
!> components starts from v2's elastic estimate (Voigt E and nu of the phases by volume
!> fraction, Q-T12). Landing (step 8.2a) redoes whole attempts, the iteration included.
!> Rows given as rates (load_rate, stress_rate on a component other than the control component,
!> step 8.2c) are load-controlled too: their value to reach is the start value plus the rate
!> times the time since the step's start (bc_mode_values).
!>
!> Load and stress rates of the control component (step 8.2c, addendum §3.1-§3.4, iteration
!> (B); decisions Q-T1, Q-T5, Q-T7, Q-T8 (i) and (iii), Q-T9): the control component is then
!> load-controlled too, its value to reach at the attempt's end time linear in time to its load
!> or stress target (a step of duration |v* - v_start|/|rate| from the value measured at its
!> start, Q-T4, D44's increments, Q-T5), or the start value plus the rate times the time (a
!> time, a target_duration, a strain target; a rate 0 holds it: a dwell). Every attempt (and
!> sub-increment) runs v2's triaxial CLR iteration (bc_iterate_rate; the pure parts in
!> bc_iteration_mod): a scaling loop (every load-controlled face velocity scaled, 1.1 or 0.9 in
!> a dwell, then the secant on the control, done when the control is within max(tol, 0.1 of the
!> increment's change) of its value; not iterated at the start of a dwell or a reload, in an
!> unloading step, or with a zero control velocity), then a finite-difference Newton loop over
!> the load-controlled face velocities until every value is within its tolerance (v2's
!> perturbation sizes with v2-fixed C60's lower bound, Q-T1, and its test before the first
!> perturbation), each loop at most max_bc_iter iterations (then a fatal error). A step is
!> loading, unloading or a dwell by the control value's magnitude (Q-T9; bc_step_action; a
!> strain-rate step by its measured control-face force); the first attempt of an episode starts
!> as v2's: at the run's first increment and at a reload (loading after a step that was not),
!> v2's elastic estimate of the face velocities for the increment's change (after a solve at
!> zero velocity at a reload); at an unload start the velocities over 10; at a dwell start the
!> velocities stored at the first increment of the last dwell after the same kind of step (0
!> the first time). A dwell's increments: 0.01 dtime at its start, then
!> min(max_strain_incr L0 / |V|, dtime), ending at its time exactly (step_dwell_dt). A strain
!> target at a load or stress rate (Q-T8 (i)) ends at the first increment that reaches it, the
!> increments chosen by v2's extrapolation (as a load target, on the bounding-box strain, with
!> tolerance 0). Necking is not tested under a load or stress rate (addendum §3.5); cyclic
!> hardening is reset at a step whose control rate (strain, load or stress) reverses. A run whose
!> step 1 has a load or stress rate has no initialization step (D49; its first increment starts
!> from v2's elastic estimate).
!>
!> Order of operations of an increment: mid-node linearization of the coordinates (D26, when
!> asked; v2 resets the stored coordinates at every iteration, the same as once per increment);
!> increment_solve from the current coordinates; coordinates updated with the converged
!> velocities; time = the step's end time at the last increment of a strain or time target (no
!> accumulated round-off), otherwise time + dt; face forces written; completion test (load and
!> stress targets: the next increment), the step results when complete, the necking test, then
!> the run limits (max_incr,
!> max_total_time, max_strain and max_eqstrain on the bounding-box lengths, v2
!> mesh_lengths): a tripped limit writes the results of that increment as the current step and
!> stops with an error. The step files are written at the end of each step with
!> print_results 1.
!>
!> Sub-incrementation (subincrementation on, Q35, D63, addendum_subincrementation.md §2-§5): an
!> increment attempt is rejected when a material point solve did not converge at the converged
!> velocities or the global iteration failed (global_max_iters, diverging; a CG failure stays
!> fatal); the attempt is undone (coordinates, state, velocities) and its interval is solved as
!> two halves, each subdivided again only if it is rejected itself (recursive bisection,
!> solve_interval), at most subincrement_max_depth times (limit L1); at most
!> subincrement_max_successive successive increments may be subdivided (limit L2); beyond
!> either limit the run stops with an error. Each accepted sub-increment is a full increment
!> (mid-node reset, increment_solve, geometry update); the face forces, the output and the run
!> limits follow the whole increment. The initialization step is not subdivided.
!> Source: Baeker (2018), arXiv:1811.05753, §2.1-1 and §2.2-2 to §2.2-4 (cutback of a
!> non-converged increment, a minimum increment and a limited number of cutbacks, each ending
!> the run with an error).
!>
!> Restart (step 8.4, dev/theory/addendum_restart.md): with print restart, at the end of every
!> step every rank writes its restart file simulation.sim/restart/restart.rank<r> (the run's
!> counters and the step just completed, the triaxial carried values, the output file counts,
!> the configuration's fingerprint; the rank's mesh part as a check, the coordinates, the
!> velocities and the whole material state; restart_file_mod), as <file>.new, renamed over the
!> previous one when every rank has written (the latest restart point only). With restart on,
!> the run continues simulation.sim from those files: its configuration lists the steps after
!> the restart point (its step 1 follows the restart point's step, whose record is the previous
!> step of step 1: planning from its end, the velocity-guess ratio, the cyclic-hardening reset,
!> the episode starts); its steps are written as step_offset + 1, ...; the material, the
!> boundary conditions and the input files must be unchanged (fingerprint; other changes a
!> warning); the same number of ranks; no step 0, initialization step or
!> isotropic viscoplastic start; the output files continued (sim_reopen).
!>
!> Locals:
!>   config       type(config_t)            the configuration
!>   mesh         type(mesh_t)              rank-local mesh
!>   elt          type(element_t)           element descriptor
!>   phases       type(phase_model_t)(:)    one model per phase, with its elements
!>   model        type(elastic_model_t)     an elastic phase's model (set up, then copied)
!>   cp_model     type(cp_model_t)          a crystal-plasticity phase's model (likewise)
!>   crys         type(crystal_t)           slip geometry of a crystal-plasticity phase
!>   hard         type(hardening_params_t)  hardening of a crystal-plasticity phase
!>   irr          type(hardening_irr_t)     its irradiation softening (step 7.3)
!>   vertices     real(rk)(:, :)            (6, :) yield-surface vertices of a crystal-
!>                                          plasticity phase, one of each +- pair, deviatoric
!>                                          Mandel, crystal frame: computed (initialization
!>                                          isotropic_viscoplastic) [MPa]; else size 0
!>   num_independent  integer               independent combinations of 5 slip systems of the
!>                                          phase (yield_vertices_compute)
!>   tau_0_phase  real(rk)(:)               the phase's initial strengths with the precipitate
!>                                          strength (step 7.1) [MPa]
!>   tau_vertex_phase real(rk)(:)           the strengths the vertices are computed from
!>                                          (tau_0_phase; Kocks-Mecking: with b mu_0
!>                                          sqrt(h_n rho_0) added, step 7.2) [MPa]
!>   clock_start, clock_end, clock_rate  integer(int64)  wall clock of the vertex computation
!>   num_slip_max integer(ik)               the largest number of slip systems of the phases
!>   num_slip_elt integer(ik)(:)            (num_elts) slip systems of each element's phase
!>   dt_init      real(rk)                  duration of the initialization step (0: none) [s]
!>   init_done    logical                   the initialization (elastic step or isotropic
!>                                          viscoplastic solve) was run
!>   state_prev, state  type(material_state_t)  states at the start and end of an increment
!>   specs        type(bc_spec_t)(:)        velocity conditions of the step
!>   cons         type(constraint_t)        multi-point constraints (inactive without them)
!>   mpc_specs    type(mpc_spec_t)(:)       the configuration's MPC relations (constraint_mod
!>                                          form)
!>   bc_mask      logical(:)                prescribed dofs
!>   bc_vel, vel  real(rk)(:)               prescribed and solved velocities [mm/s]
!>   coo          real(rk)(:, :)            (3, num_nodes_loc) current coordinates [mm]
!>   opts         type(increment_opts_t)    iteration parameters
!>   info         type(increment_info_t)    iteration report
!>   sim          type(sim_writer_t)        output writer
!>   stress_c, strain_el_c  real(rk)(:, :)  (6, num_elts) centroid-point stress [MPa] and
!>                                          elastic strain [-], Mandel, sample frame
!>   rot_c      real(rk)(:, :, :)           (3, 3, num_elts) centroid-point rot_c2s [-]
!>   out6       real(rk)(:, :)              (6, num_elts) a centroid-point Mandel output
!>   out3       real(rk)(:, :)              (3, num_elts) a centroid-point 3-vector output
!>   out1       real(rk)(:, :)              (1, num_elts) a centroid-point scalar output
!>   out9       real(rk)(:, :)              (9, num_elts) a centroid-point 3x3 output (velgrad)
!>   loads, areas real(rk)(:, :), (:)       face-set forces [N] and areas [mm^2]
!>   aux, inputs  character(64)(:)          auxiliary files of the mesh; input files to copy
!>   tess_exists  logical                   simulation.tess is in the working directory (copied
!>                                          with the inputs, not read)
!>   sched        type(step_schedule_t)     the steps as a time schedule
!>   dt           real(rk)                  time increment [s]
!>   time         real(rk)                  time at the end of the current increment [s]
!>   rate_general real(rk)                  general: max |set_bc velocity| / gage length of the
!>                                          loading axis (0 without one) [1/s]
!>   lengths0     real(rk)(3)               original bounding-box lengths [mm]
!>   axis         integer                   control component (0: general without strain
!>                                          targets)
!>   i_face       integer                   control face set (x1, y1 or z1 of the control
!>                                          component; 0 when not needed)
!>   is_load      logical                   the current step ends at a load or stress target
!>   step_done    logical                   the current step is complete
!>   written      logical                   the increment's results were written as the step's
!>   dt_next      real(rk)                  next increment of a load or stress target [s]
!>   coo_land, vel_land, state_land         coordinates [mm], velocities [mm/s] and state at
!>                                          the start of an increment of a load or stress
!>                                          target (to discard an attempt that passes it)
!>   value, value_prev  real(rk)            load [N] or stress [MPa] of the control face after
!>                                          and before the increment (load and stress targets)
!>   force_prev   real(rk)                  control-face force before the increment (necking
!>                                          test) [N]
!>   strain_now   real(rk)(3)               engineering strains of the prescribed components
!>                                          reached [-]
!>   axis_letter  character(1)(3)           x, y, z
!>   i_step, i_incr  integer                step, increment within the step
!>   incr         integer                   increments done (cumulative)
!>   num_successive  integer                successive increments subdivided so far (L2)
!>   num_attempts integer                   attempts of the current increment
!>   num_subincrs integer                   accepted sub-increments of the current increment
!>   dt_smallest  real(rk)                  smallest accepted sub-increment of the increment [s]
!>   i_phase, i   integer                   indices
!>   i_elt, i_qpt integer(ik)               element and point
!>   text         character(200)            message buffer
!>   vel_face     real(rk)(3)               normal velocities of the faces x1, y1, z1 (triaxial;
!>                                          load-controlled ones found by the iteration) [mm/s]
!>   vel_face_land  real(rk)(3)             vel_face at the start of a landed increment [mm/s]
!>   face_set     integer(3)                face-set index of x1, y1, z1 of the load-controlled
!>                                          components (0: not needed)
!>   face_dof     integer(:)                (3 num_nodes_loc) component c of a dof that is the
!>                                          c-velocity of a node on face c1 of a load-controlled
!>                                          component, else 0
!>   bc_any       logical                   some step has load-controlled components (8.2b)
!>   bc_on        logical                   the current step has load-controlled components
!>   bc_path      logical                   their values follow the load path (control load or
!>                                          stress target), else they are interpolated in time
!>   bc_comps     integer(:)                the load-controlled components of the step
!>   bc_q_start   real(rk)(3)               their values at the step's start [N] or [MPa]
!>   bc_ctrl_start  real(rk)                the control value at the step's start (path) [N] or
!>                                          [MPa]
!>   bc_delta     real(rk)(3)               perturbation sizes of the face velocities (0: not set
!>                                          yet) [mm/s]
!>   box_lo, box_len  real(rk)(3)           bounding box of the attempt's start coordinates:
!>                                          lowest coordinates and lengths [mm]
!>   time_solved  real(rk)                  time at the end of the last accepted (sub-)increment
!>                                          [s]
!>   strain_box   real(rk)(3)               engineering strains of the bounding box [-]
!>   areas0       real(rk)(3)               original areas of the faces x1, y1, z1 of the
!>                                          load-controlled components (0 for the others) [mm^2]
!>   e_avg, nu_avg  real(rk)                Voigt-averaged E [MPa] and nu [-] of the phases
!>                                          (by volume fraction, Q-T12)
!>   stress_max   real(rk)(3)               largest target magnitude of each component over the
!>                                          run, as a stress (loads over areas0) [MPa]
!>                                          (iteration (B)'s perturbation floor, C60)
!>   action, prev_action  integer           kind of the current and the previous step
!>                                          (bc_action_*; 0: none yet)
!>   start_run, start_reload, start_unload, start_dwell  logical  the current step's first
!>                                          attempt starts the run, a reload, an unload, a dwell
!>                                          episode (iteration (B), v2 CLR)
!>   bc_at_start  logical                   the next attempt starts at the step's start (the
!>                                          episode starts apply)
!>   vel_dwell_load, vel_dwell_unload  real(rk)(3)  face velocities stored at the first
!>                                          increment of the last dwell after loading,
!>                                          unloading [mm/s]
!>   force_start  real(rk)                  control-face force at the step's start [N]
!>   last_incr    logical                   the increment ends a step of known duration
!>   is_open      logical                   the current step's duration is not known (a load or
!>                                          stress target at a strain rate, a strain target at
!>                                          a load or stress rate)
!>   value_start  real(rk)                  control value at the start of a load- or stress-rate
!>                                          step (planning) [N] or [MPa]
!>   restarted    logical                   the run continues from restart files (restart on)
!>   step_offset  integer                   steps of the run before this one's step 1 (the
!>                                          restart point's step number; 0 without restart)
!>   prev_rec     type(step_record_t)       the restart point's step (the previous step of
!>                                          step 1)
!>   rst_in       type(restart_file_t)      the restart file being read
!>   rst_incr, rst_num_successive  integer  increments and successive subdivided increments at
!>                                          the restart point
!>   rst_time, rst_force_prev  real(rk)     time [s] and control-face force [N] there
!>   rst_strain_now, rst_stress_max  real(rk)(3)  strains reached [-] and the perturbation
!>                                          floors [MPa] there
!>   rst_lines    integer(3)                lines of the force, convergence and macroscopic
!>                                          files there
!>   rst_rows     integer(3)                rows of the HDF5 force, convergence and
!>                                          macroscopic tables there (step 9.1)
!>   rst_conv     logical(2)                the convergence file's columns there (attempts,
!>                                          boundary iterations)
!>   rst_num_steps integer                  largest step written there
program fepx

  use, intrinsic :: iso_fortran_env, only: int64, compiler_version, output_unit
  use hdf5_file_mod, only: hdf5_available, h5_library_version
  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error, set_abort_handler
  use log_mod, only: log_info, log_warning, log_file, log_text, log_warning_count, &
    & log_hold_start, log_hold_take, log_lines_write, log_real_text, log_int_text, log_count_text
  use build_info_mod, only: build_version_text
  use parallel_mod, only: par_init, par_finalize, par_library_version, num_ranks, my_rank, &
    & par_sum, par_min, par_max, par_bcast, par_barrier, par_any, par_wtime, par_log_flush
  use config_mod, only: config_t, config_read, cfg_phase_t, config_log_material, &
    & config_log_phase, config_log_options
  use mesh_file_mod, only: elt_type_name
  use mesh_mod, only: mesh_t
  use mesh_distribute_mod, only: mesh_load
  use element_mod, only: element_t, element_init, element_volume
  use tensor_mod, only: equiv_strain, equiv_stress, matrix_from_mandel, matrix_from_skew
  use rotation_mod, only: transform_crystal_to_sample_mandel
  use material_state_mod, only: material_state_t, material_state_alloc
  use elastic_model_mod, only: elastic_model_t, elastic_model_init
  use crystal_mod, only: crystal_t, crystal_init, crystal_structure_from_name, &
    & crystal_family_size, shear_modulus_hill
  use crystal_plasticity_model_mod, only: cp_model_t, cp_model_init, cp_model_set_element_values
  use yield_vertices_mod, only: yield_vertices_compute
  use hardening_mod, only: hardening_params_t, hardening_none, hardening_isotropic, &
    & hardening_anisotropic, hardening_kocks_mecking, hardening_km_strengths, hardening_irr_t, &
    & hardening_irr_loop_average, hardening_irr_planes, hardening_irr_dmg_0, &
    & hardening_irr_strengths
  use bc_mod, only: bc_spec_t, bc_preset_specs, bc_triaxial_specs, bc_resolve, bc_gage_length, &
    & bc_lengths, bc_preset_grip, bc_preset_symmetry, bc_preset_minimal
  use increment_mod, only: phase_model_t, increment_opts_t, increment_info_t, increment_solve, &
    & midnode_linearize, increment_converged, increment_max_iters, increment_diverging
  use face_force_mod, only: face_forces
  use constraint_mod, only: constraint_t, mpc_spec_t, constraint_build_mpc, &
    & constraint_propagate_bcs, constraint_project, constraint_check_bcs, constraint_set_rates, &
    & constraint_zero_mean, constraint_macro_rate
  use periodic_mod, only: periodic_build, periodic_stress_average
  use isotropic_start_mod, only: isotropic_start_solve
  use pcg_mod, only: pcg_precond_none, pcg_precond_jacobi, pcg_precond_block_jacobi
  use sim_writer_mod, only: sim_writer_t, sim_open, sim_write_nodes, sim_write_elts, &
    & sim_end_step, sim_write_elts_values, sim_write_elts_slip, sim_write_ori, sim_write_forces, &
    & sim_write_convergence, sim_write_macro, sim_close, sim_reopen, sim_restart_dir
  use step_schedule_mod, only: step_record_t, step_schedule_record
  use restart_file_mod, only: restart_file_t, restart_create, restart_open, restart_close, &
    & restart_put, restart_get, restart_put_state, restart_get_state, restart_put_step, &
    & restart_get_step, restart_rename
  use step_schedule_mod, only: step_schedule_t, step_schedule_build, step_schedule_plan, &
    & step_schedule_load_start, step_load_next_dt, step_load_landing, step_load_false_position, &
    & landing_within, landing_passed, step_load_complete, step_dwell_dt
  use config_mod, only: target_by_stress, rate_by_stress
  use bc_iteration_mod, only: bc_converged, bc_pert_signs, bc_correction, bc_pert_update, &
    & bc_elastic_moduli, bc_elastic_velocities, bc_initial_velocities, bc_mode_values, &
    & bc_step_action, bc_rate_first_converged, bc_rate_scale, bc_rate_pert_floor, &
    & bc_rate_pert_initial, bc_rate_pert_update, bc_mode_path, bc_mode_time, bc_mode_rate, &
    & bc_action_loading, bc_action_unloading, bc_action_dwell

  implicit none

  type(config_t) :: config
  type(mesh_t) :: mesh
  type(element_t) :: elt
  type(phase_model_t), allocatable :: phases(:)
  type(elastic_model_t) :: model
  type(cp_model_t) :: cp_model
  real(rk), allocatable :: vertices(:, :), tau_0_phase(:), tau_vertex_phase(:)
  type(crystal_t) :: crys
  type(hardening_params_t) :: hard
  type(hardening_irr_t) :: irr
  integer(ik) :: num_slip_max
  integer(ik), allocatable :: num_slip_elt(:)
  real(rk) :: dt_init
  logical :: init_done
  type(material_state_t) :: state_prev, state, state_land
  type(bc_spec_t), allocatable :: specs(:)
  type(constraint_t) :: cons
  type(mpc_spec_t), allocatable :: mpc_specs(:)
  logical, allocatable :: bc_mask(:)
  real(rk), allocatable :: bc_vel(:), vel(:), coo(:, :), coo_land(:, :), vel_land(:)
  type(increment_opts_t) :: opts
  type(increment_info_t) :: info
  character(len=64), allocatable :: aux(:), inputs(:)
  type(sim_writer_t) :: sim
  real(rk), allocatable :: stress_c(:, :), strain_el_c(:, :), rot_c(:, :, :), loads(:, :), &
    & areas(:), out6(:, :), out3(:, :), out1(:, :), out9(:, :)
  type(step_schedule_t) :: sched
  real(rk) :: dt, time, rate_general, lengths0(3)
  integer :: axis, i_phase, i, i_step, i_incr, incr, i_face
  logical :: is_load, step_done, written
  real(rk) :: dt_next, value, value_prev, force_prev, strain_now(3)
  character(len=1), parameter :: axis_letter(3) = ['x', 'y', 'z']
  integer :: num_successive, num_attempts, num_subincrs
  real(rk) :: dt_smallest
  integer :: num_independent
  integer(int64) :: clock_start, clock_end, clock_rate
  integer(ik) :: i_elt, i_qpt
  character(len=200) :: text
  real(rk) :: vel_face(3), vel_face_land(3), bc_q_start(3), bc_ctrl_start, bc_delta(3)
  real(rk) :: box_lo(3), box_len(3), time_solved, strain_box(3)
  integer :: face_set(3)
  integer, allocatable :: face_dof(:), bc_comps(:)
  logical :: bc_any, bc_on, bc_path
  real(rk) :: areas0(3), e_avg, nu_avg, stress_max(3), vel_dwell_load(3), vel_dwell_unload(3)
  real(rk) :: force_start, value_start
  integer :: action, prev_action
  logical :: start_run, start_reload, start_unload, start_dwell, bc_at_start, last_incr, is_open
  logical :: restarted
  logical :: tess_exists
  integer :: step_offset
  type(step_record_t) :: prev_rec
  type(restart_file_t) :: rst_in
  integer :: rst_incr, rst_num_successive, rst_lines(3), rst_rows(3), rst_num_steps
  real(rk) :: rst_time, rst_force_prev, rst_strain_now(3), rst_stress_max(3)
  logical :: rst_conv(2)

  !> Log lines kept by the hold mode for the setup summary (step 9.7).
  type :: held_text_t
    character(len=:), allocatable :: text
  end type held_text_t
  type(held_text_t), allocatable :: phase_held(:)
  character(len=:), allocatable :: loading_held, restart_held, held_part
  real(rk) :: time_start, time_setup_end, time_init_end, time_loop_end, time_mark
  real(rk) :: time_output, time_output_setup, time_output_run, time_elements_run, time_cg_run
  integer :: num_iters_run, num_nr_run, num_cg_run, incr_start
  integer(ik) :: num_warnings
  logical :: in_run

  ! Command line (step 9.7, Q64): --help and --version before MPI starts; any other argument
  ! is an error.
  call command_line()

  ! MPI, the wall clock of the run, the banner (step 9.7, Q58, Q59).
  call par_init()
  time_start = par_wtime()
  time_output = 0.0_rk
  time_elements_run = 0.0_rk
  time_cg_run = 0.0_rk
  num_iters_run = 0
  num_nr_run = 0
  num_cg_run = 0
  in_run = .false.
  call write_banner(.true.)

  ! --- Input ---

  ! Configuration (read by every rank), then the mesh and its auxiliary files; the file lines
  ! are written as the files are read, the summary after the setup (log_setup_summary).
  call log_info('Loading simulation...')
  call log_file('i', "Parsing file `simulation.cfg'...")
  call config_read(config, 'simulation.cfg')
  call log_file('i', "Parsed file `simulation.cfg'.")

  allocate (aux(0))
  if (config%read_ori_from_file) aux = [character(len=64) :: aux, 'simulation.ori']
  if (config%read_phase_from_file) aux = [character(len=64) :: aux, 'simulation.phase']
  if (config%read_from_opt_file) aux = [character(len=64) :: aux, 'simulation.opt']
  call mesh_load(mesh, 'simulation.msh', aux)
  call element_init(elt, mesh%elt_type)

  ! The mesh must match the configuration.
  if (mesh%num_phases /= int(config%num_phases, ik)) then
    write (text, '(a,i0,a,i0,a)') 'the mesh has ', mesh%num_phases, ' phase(s), the ' // &
      & 'configuration ', config%num_phases, ' (number_of_phases)'
    call fatal_error('fepx', trim(text))
  end if
  if (.not. mesh%has_ori) call fatal_error('fepx', 'the mesh has no orientations ' // &
    & '($ElsetOrientations, $ElementOrientations or simulation.ori)')

  ! Multi-point constraints (step 8.3a): the replicated table and its local rows. (Their log
  ! lines, those of periodic conditions and of the step schedule are kept for the summary's
  ! "- Loading" part.)
  call log_hold_start()
  if (size(config%mpcs) > 0) then
    allocate (mpc_specs(size(config%mpcs)))
    do i = 1, size(config%mpcs)
      mpc_specs(i)%nset = config%mpcs(i)%nset
      mpc_specs(i)%comp = config%mpcs(i)%comp
      mpc_specs(i)%coeff = config%mpcs(i)%coeff
      mpc_specs(i)%offset = config%mpcs(i)%offset
      mpc_specs(i)%line = config%mpcs(i)%line
    end do
    call constraint_build_mpc(cons, mesh, mpc_specs)
  end if

  ! Periodic conditions (step 8.3b): the drive pairs, the replicated table of the periodicity
  ! relations, its local rows and the zero-mean weights.
  if (config%bc_type == 'periodic') call periodic_build(cons, mesh, elt, &
    & config%periodic_kinematics == 'small_strain', config%per_imposed, &
    & config%periodic_preconditioner == 'exact_diagonal')

  call log_hold_take(loading_held)

  ! Restart (step 8.4): the restart point's counters, step and checks, from this rank's file
  ! (its log lines kept for the summary's "- Restart" part).
  restarted = .false.
  step_offset = 0
  call log_hold_start()
  if (config%restart) call restart_read_globals()
  call log_hold_take(restart_held)

  ! The steps as a time schedule (with the checks that need the start-of-step values), from the
  ! strain rate of general conditions: max |set_bc velocity| over the gage length (v2); after
  ! a restart point, from its step.
  axis = config%control_axis
  rate_general = 0.0_rk
  if (config%bc_type == 'general' .and. axis > 0) rate_general = &
    & maxval(abs(config%vel_bcs(:)%value))/bc_gage_length(mesh, axis)
  call log_hold_start()
  if (restarted) then
    call step_schedule_build(sched, config, rate_general, prev_rec, rst_incr)
  else
    call step_schedule_build(sched, config, rate_general)
  end if
  call log_hold_take(held_part)
  loading_held = loading_held // held_part

  ! The control face x1, y1 or z1 of the control component, whose force is the load of load
  ! and stress targets and of the necking test (addendum §3.1; v2's loading_face).
  i_face = 0
  if (any(sched%is_load) .or. config%check_necking .or. any(sched%is_rate)) then
    do i = 1, size(mesh%fasets)
      if (mesh%fasets(i)%label == axis_letter(axis) // '1') i_face = i
    end do
    if (i_face == 0) call fatal_error('fepx', 'the mesh has no face set ' // &
      & axis_letter(axis) // '1, whose force is the load of the control component ' // &
      & int_str(int(11*axis, ik)) // ' (target_load, target_stress, check_necking)')
  end if

  ! The faces x1, y1, z1 of the load-controlled components (step 8.2b): their face sets (whose
  ! forces are controlled) and the dofs of their normal velocities.
  bc_any = any(sched%bc_kind /= 0)
  bc_on = .false.
  bc_path = .false.
  allocate (bc_comps(0))
  face_set = 0
  vel_face = 0.0_rk
  bc_delta = 0.0_rk
  bc_q_start = 0.0_rk
  bc_ctrl_start = 0.0_rk
  areas0 = 0.0_rk
  e_avg = 0.0_rk
  nu_avg = 0.0_rk
  stress_max = 0.0_rk
  vel_dwell_load = 0.0_rk
  vel_dwell_unload = 0.0_rk
  action = 0
  prev_action = 0
  start_run = .false.
  start_reload = .false.
  start_unload = .false.
  start_dwell = .false.
  bc_at_start = .false.
  force_start = 0.0_rk
  if (bc_any) call bc_faces()

  ! --- Material: one model per phase, the state with the mesh orientations ---

  allocate (phases(config%num_phases))
  allocate (phase_held(config%num_phases))
  num_slip_max = 0_ik
  allocate (num_slip_elt(mesh%num_elts))
  do i_phase = 1, config%num_phases
    ! (The phase's derived quantities on the log are kept for its part of the summary.)
    call log_hold_start()
    associate (ph => config%phases(i_phase))
      if (ph%model == 'crystal_plasticity') then
        ! Slip geometry of the active families (D48), then the model (constant strengths).
        call crystal_init(crys, crystal_structure_from_name(ph%crystal_type), &
          & c_over_a=ph%c_over_a, families=ph%families)
        ! Hardening law and its per-system parameters (D51).
        hard%law = hardening_none
        if (ph%hardening == 'saturation_isotropic') hard%law = hardening_isotropic
        if (ph%hardening == 'saturation_anisotropic') hard%law = hardening_anisotropic
        if (ph%hardening == 'kocks_mecking') hard%law = hardening_kocks_mecking
        hard%h_0 = ph%h_0
        hard%tau_s = ph%tau_s
        hard%n = ph%n_hard
        hard%interaction = ph%interaction
        hard%max_iters = config%point_max_iters_hardening
        hard%tol_hardening = config%point_tol_hardening
        ! Modifiers (step 7.1, D59); the precipitate strength raises the phase's tau_0.
        hard%saturation_evolution = ph%saturation_evolution
        hard%gammadot_s0 = ph%gammadot_s0
        hard%m_prime = ph%m_prime
        hard%cyclic = ph%cyclic
        hard%cyclic_a = ph%cyclic_a
        hard%cyclic_c = ph%cyclic_c
        ! Kocks-Mecking (step 7.2, D60): its parameters; mu_0 is Hill's average shear modulus
        ! of the phase's elastic constants (addendum_kocks_mecking.md §4).
        hard%km_b = ph%km_b
        hard%km_h_n = ph%km_h_n
        hard%km_rho_0 = ph%km_rho_0
        hard%km_k1 = ph%km_k1
        hard%km_k20 = ph%km_k20
        hard%km_gammadot_k0 = ph%km_gammadot_k0
        hard%km_n = ph%km_n
        hard%km_mu_0 = 0.0_rk
        if (hard%law == hardening_kocks_mecking) then
          hard%km_mu_0 = shear_modulus_hill(ph%stiff_el_c)
          call log_info('Kocks-Mecking shear modulus mu_0 (Hill): ' // &
            & log_real_text(hard%km_mu_0, 4) // ' MPa; initial strength increment b mu_0 ' // &
            & 'sqrt(h_n rho_0): ' // &
            & log_real_text(hard%km_b*hard%km_mu_0*sqrt(hard%km_h_n*hard%km_rho_0), 4) // &
            & ' MPa', 3)
        end if
        ! Irradiation softening (step 7.3, D61): its parameters, the average loop tensor of the
        ! loop-normal families and the slip planes N^a (addendum_irradiation.md §5-§6).
        irr = hardening_irr_t()
        if (ph%irradiation) then
          irr%h_d = ph%irr_h_d
          irr%eta = ph%irr_eta
          irr%loop_density = ph%irr_loop_density
          irr%loop_diameter = ph%irr_loop_diameter
          irr%loop_avg_c = hardening_irr_loop_average(ph%irr_loop_normals_c)
          irr%plane_c = hardening_irr_planes(crys%slip_normal_c)
          call log_info('Irradiation, initial defect strength term (phase loops, system 1): ' &
            & // log_real_text(hard%km_b*hard%km_mu_0*sqrt(irr%h_d* &
            & max(dot_product(irr%plane_c(:, 1), hardening_irr_dmg_0(irr, irr%loop_density, &
            & irr%loop_diameter)), 0.0_rk)), 4) // ' MPa', 3)
        end if
        tau_0_phase = ph%tau_0 + ph%tau_p
        if (ph%precipitation) then
          call log_info('Precipitate strength: ' // log_real_text(ph%tau_p, 4) // ' MPa (' // &
            & trim(ph%precip_regime) // ') added to tau_0', 3)
        end if
        ! Yield-surface vertices of the increment-1 estimate: computed from the phase's active
        ! systems and tau_0 (initialization isotropic_viscoplastic, D57; at least five
        ! independent systems).
        if (config%initialization == 'isotropic_viscoplastic') then
          call system_clock(clock_start, clock_rate)
          ! With Kocks-Mecking, from the initial strengths tau_0 + b mu_0 sqrt(h_n rho_0) of the
          ! phase (addendum_kocks_mecking.md §5; per-element rho_0: an approximate start, D58).
          tau_vertex_phase = tau_0_phase
          if (hard%law == hardening_kocks_mecking) tau_vertex_phase = &
            & hardening_km_strengths(hard, tau_0_phase, hard%km_rho_0)
          ! With irradiation, also the defect term of the phase's loops (addendum_irradiation.md
          ! §7).
          if (ph%irradiation) tau_vertex_phase = hardening_irr_strengths(hard, irr, &
            & tau_0_phase, hard%km_rho_0, hardening_irr_dmg_0(irr, irr%loop_density, &
            & irr%loop_diameter))
          call yield_vertices_compute(crys%schmid_sym_c, tau_vertex_phase, vertices, &
            & num_independent)
          call system_clock(clock_end)
          if (num_independent == 0) then
            write (text, '(a,i0,a)') 'initialization isotropic_viscoplastic: the active ' // &
              & 'slip systems of ' // &
              & 'phase ', i_phase, ' span fewer than 5 deviatoric dimensions (at least ' // &
              & 'five independent systems are needed)'
            call fatal_error('fepx', trim(text))
          end if
          write (text, '(a,i0,a,i0,a,f0.3,a)') 'Yield-surface vertices: ', size(vertices, 2), &
            & ' +- pairs (', num_independent, ' independent combinations of 5 systems), ', &
            & real(clock_end - clock_start, rk)/real(clock_rate, rk), ' s'
          call log_info(trim(text), 3)
        else
          if (allocated(vertices)) deallocate (vertices)
          allocate (vertices(6, 0))
        end if
        call cp_model_init(cp_model, ph%stiff_el_c, crys%schmid_sym_c, crys%schmid_skew_c, &
          & ph%m, ph%gammadot_0, tau_0_phase, tol_newton=config%point_tol, &
          & max_iters_newton=config%point_max_iters_newton, &
          & max_iters_state=config%point_max_iters_state, tol_state=config%point_tol_state, &
          & hard=hard, vertices_c=vertices)
        if (ph%irradiation) then
          cp_model%irradiation = .true.
          cp_model%irr = irr
        end if
        call element_values(i_phase, cp_model)
        allocate (phases(i_phase)%model, source=cp_model)
        num_slip_max = max(num_slip_max, cp_model%num_slip)
      else
        call elastic_model_init(model, ph%stiff_el_c, tol_state=config%point_tol_state, &
          & max_iters_state=config%point_max_iters_state)
        allocate (phases(i_phase)%model, source=model)
      end if
    end associate
    phases(i_phase)%elts = pack([(i_elt, i_elt = 1_ik, mesh%num_elts)], &
      & mesh%elt_phase == int(i_phase, ik))
    num_slip_elt(phases(i_phase)%elts) = phases(i_phase)%model%num_slip
    call log_hold_take(phase_held(i_phase)%text)
  end do
  call opt_warnings()
  call material_state_alloc(state, elt%num_qpts, mesh%num_elts, num_slip_max)
  do i_elt = 1_ik, mesh%num_elts
    do i_qpt = 1_ik, elt%num_qpts
      state%rot_c2s(:, :, i_qpt, i_elt) = mesh%rot_c2s(:, :, i_elt)
    end do
  end do
  do i_phase = 1, config%num_phases
    call phases(i_phase)%model%init_state(phases(i_phase)%elts, state)
  end do

  ! --- Output of the initial state (step 0) ---

  ! Original coordinates, zero displacement, stress and elastic strain, the mesh orientations,
  ! face areas with zero forces, and as velocities the prescribed velocities of step 1, zero
  ! elsewhere (as v2, which initializes the velocities with them before writing step 0).
  coo = mesh%coo_s
  allocate (vel(3*mesh%num_nodes_loc), bc_vel(3*mesh%num_nodes_loc))
  allocate (bc_mask(3*mesh%num_nodes_loc))
  vel = 0.0_rk
  init_done = .false.
  allocate (stress_c(6, mesh%num_elts), strain_el_c(6, mesh%num_elts))
  allocate (rot_c(3, 3, mesh%num_elts))
  allocate (out6(6, mesh%num_elts), out3(3, mesh%num_elts), out1(1, mesh%num_elts))
  allocate (out9(9, mesh%num_elts))
  allocate (loads(3, size(mesh%fasets)), areas(size(mesh%fasets)))
  lengths0 = bc_lengths(mesh, mesh%coo_s)
  ! (Also with print restart: the restart point carries the perturbation floors, step 8.4.)
  ! The triaxial estimates' log lines are kept for the summary's "- Loading" part.
  call log_hold_start()
  if (bc_any .or. config%print_restart .or. restarted) call bc_setup()
  if (.not. restarted) then
    if (any(sched%bc_kind(:, 1) /= 0) .and. .not. sched%is_rate(1)) call bc_initial_estimate()
    call step_conditions(1)
  end if
  call log_hold_take(held_part)
  loading_held = loading_held // held_part
  ! The input files copied to simulation.sim/inputs: those read, and the tessellation
  ! simulation.tess when it is in the working directory (not read, but kept with the results
  ! for Neper's post-processing and visualization).
  inputs = [character(len=64) :: 'simulation.cfg', 'simulation.msh', aux]
  inquire (file='simulation.tess', exist=tess_exists)
  if (tess_exists) inputs = [character(len=64) :: inputs, 'simulation.tess']

  ! The setup summary (step 9.7, dev/plan/terminal_output.md §3.2), then the output files.
  call log_setup_summary()
  call log_file('o', "Writing results to `simulation.sim'...")
  time_mark = par_wtime()
  if (restarted) then
    ! A restarted run (step 8.4): the output directory continued from the restart point (the
    ! convergence columns of both runs).
    call sim_reopen(sim, 'simulation.sim', config%print_vars, mesh, inputs, &
      & config%subincrementation .or. any(sched%is_load) .or. rst_conv(1), bc_any .or. &
      & rst_conv(2), config%bc_type == 'periodic', rst_num_steps, rst_lines(1), rst_lines(2), &
      & rst_lines(3), rst_conv(1), rst_conv(2), config%output_ascii, config%output_hdf5, rst_rows)
  else
    ! (The convergence file numbers the attempts of an increment with sub-incrementation and
    ! with load or stress targets, whose increments may be redone to land on the target.)
    call sim_open(sim, 'simulation.sim', config%print_vars, mesh, inputs, &
      & attempts=config%subincrementation .or. any(sched%is_load), bc_iterations=bc_any, &
      & macro=config%bc_type == 'periodic', ascii=config%output_ascii, hdf5=config%output_hdf5)
    time = 0.0_rk
    call write_step(0)
    call face_forces(mesh, elt, coo, stress_c, loads, areas)
    call sim_write_forces(sim, 0, 0, loads, areas, 0.0_rk)
    call write_macro(0, 0)
  end if
  time_output = time_output + (par_wtime() - time_mark)
  time_output_setup = time_output
  call par_log_flush()
  time_setup_end = par_wtime()

  ! --- Time stepping ---

  opts%global_max_iters = config%global_max_iters
  opts%global_tol = config%global_tol
  opts%global_tol_abs = config%global_tol_abs
  opts%global_tol_slow = config%global_tol_slow
  opts%global_nr_switch = config%global_nr_switch
  opts%global_nr_slow_ratio = config%global_nr_slow_ratio
  opts%cg_tol = config%cg_tol
  opts%cg_max_iters = config%cg_max_iters
  ! The CG preconditioner (step 9.2): the keyword's word as pcg_mod's code.
  select case (config%cg_preconditioner)
  case ('block_jacobi')
    opts%cg_precond = pcg_precond_block_jacobi
  case ('none')
    opts%cg_precond = pcg_precond_none
  case default
    opts%cg_precond = pcg_precond_jacobi
  end select
  time = 0.0_rk
  strain_now = 0.0_rk
  force_prev = 0.0_rk
  dt_next = 0.0_rk
  incr = 0
  num_successive = 0

  ! Initialization step (D49): part of step 1, not written, not counted as an increment (none
  ! after a restart point).
  init_done = .false.
  dt_init = 0.0_rk
  if (.not. restarted .and. any(config%phases(:)%model == 'crystal_plasticity')) &
    & call log_info('Initializing simulation...')
  if (.not. restarted) dt_init = initialization_dt()
  if (dt_init > 0.0_rk) then
    call step_conditions(1)
    call log_info('Elastic step of ' // log_real_text(dt_init, 4) // ' s (step 1)', 1)
    if (config%midnode_linearization) call midnode_linearize(mesh, elt, coo)
    state_prev = state
    call increment_solve(mesh, elt, phases, coo, bc_mask, dt_init, opts, state_prev, state, &
      & vel, info, cons)
    call report_increment(1, .true.)
    coo = coo + dt_init*reshape(vel, [3_ik, mesh%num_nodes_loc])
    time = dt_init
    init_done = .true.
  end if

  ! Isotropic viscoplastic start (D78, v2 solveit_vp): with a crystal-plasticity phase, the
  ! velocity guess of increment 1 from the step-1 conditions (coordinates and time unchanged,
  ! nothing written; step 0 holds the conditions' guess, as v2).
  if (config%initialization == 'isotropic_viscoplastic' .and. &
    & any(config%phases(:)%model == 'crystal_plasticity') .and. .not. restarted) then
    call log_info('Isotropic viscoplastic solve (velocity guess of increment 1)', 1)
    call isotropic_start_solve(mesh, elt, coo, bc_mask, config%cg_tol, config%cg_max_iters, &
      & vel, cons, opts%cg_precond)
    init_done = .true.
  end if

  ! A restarted run (step 8.4): the coordinates, velocities, state and carried values of the
  ! restart point.
  if (restarted) call restart_read_fields()

  ! The control face's force after the initialization step (start of the necking test; after a
  ! restart point, its value there).
  if (i_face > 0 .and. .not. restarted) then
    stress_c = state%stress_s(:, elt%qpt_centroid, :)
    call face_forces(mesh, elt, coo, stress_c, loads, areas)
    force_prev = loads(axis, i_face)
  end if
  call par_log_flush()
  time_init_end = par_wtime()
  in_run = .true.
  incr_start = incr

  do i_step = 1, sched%num_steps
    call log_info('Running step ' // log_int_text(step_offset + i_step) // '...')
    ! Plan a step whose start was not known before the run (after a load or stress target,
    ! addendum §3.5) from the time and strains reached.
    if (i_step > sched%num_planned) then
      ! (A load or stress rate: from the control value measured now, step 8.2c.)
      value_start = 0.0_rk
      if (sched%is_rate(i_step)) value_start = rate_control_value(i_step)
      call step_schedule_plan(sched, config, i_step, time, strain_now, rate_general, &
        & value_start)
      call log_step_plan(i_step, 1)
    end if

    ! A load or stress target: its direction from the current value (addendum §3.3).
    is_load = sched%is_load(i_step)
    if (is_load) then
      stress_c = state%stress_s(:, elt%qpt_centroid, :)
      call face_forces(mesh, elt, coo, stress_c, loads, areas)
      value_prev = control_value(i_step)
      call step_schedule_load_start(sched, config, i_step, value_prev, rate_general)
      dt_next = sched%dt(i_step)
    end if
    ! A strain target at a load or stress rate (step 8.2c): the strain reached, increments as a
    ! load target's.
    is_open = is_load .or. sched%is_open(i_step)
    if (sched%is_open(i_step)) then
      value_prev = strain_now(axis)
      dt_next = sched%dt(i_step)
    end if

    ! Load-controlled components of the step and their start values (step 8.2b); at a load or
    ! stress rate the kind of step and its episode start (step 8.2c); the control-face force at
    ! the start (the kind of a strain-rate step, measured).
    call bc_step_start(i_step)
    start_run = .false.
    start_reload = .false.
    start_unload = .false.
    start_dwell = .false.
    bc_at_start = .false.
    if (sched%is_rate(i_step)) call bc_rate_step_start(i_step)
    if (any(sched%is_rate)) force_start = control_force()

    ! Velocity conditions of the step, and the initial velocity guess.
    call step_conditions(i_step)

    i_incr = 0
    step_done = .false.
    do while (.not. step_done)
      i_incr = i_incr + 1
      ! Time increment: strain and time targets, full increments, then the last one (D44);
      ! load and stress targets, the increment chosen after the previous one (v2).
      ! A dwell (step 8.2c): v2's increments, the last one ending at the step's time.
      last_incr = .false.
      if (is_open) then
        dt = dt_next
      else if (sched%is_dwell(i_step)) then
        call step_dwell_dt(i_incr == 1 .and. start_dwell, vel_face(axis), lengths0(axis), &
          & sched%max_strain_incr, sched%dt(i_step), sched%time_end(i_step) - time, dt, &
          & last_incr)
      else
        dt = sched%dt(i_step)
        last_incr = i_incr == sched%num_incrs(i_step)
        if (last_incr) dt = sched%dt_last(i_step)
      end if
      if (i_step == 1 .and. i_incr == 1) dt = dt - dt_init
      incr = incr + 1
      ! The increment's end time (the step's end at its last increment; a landed increment of
      ! a load or stress target may end earlier, land_increment) and length.
      if (last_incr) then
        call log_info('Increment ' // log_int_text(incr) // ': t = ' // &
          & log_real_text(sched%time_end(i_step), 4) // ' s, dt = ' // log_real_text(dt, 4) // &
          & ' s', 1)
      else
        call log_info('Increment ' // log_int_text(incr) // ': t = ' // &
          & log_real_text(time + dt, 4) // ' s, dt = ' // log_real_text(dt, 4) // ' s', 1)
      end if

      ! Increment 1 with initialization isotropic_viscoplastic: the point solves start from
      ! v2's vertex estimate (D54 (3), D57, D78); later increments from the previous stress.
      state%start_from_vertex = incr == 1 .and. &
        & config%initialization == 'isotropic_viscoplastic'
      ! A step that reverses the control rate (strain, load or stress rate, step 8.2c, addendum
      ! §3.5) resets the accumulated slip of cyclic hardening (D59; v2 only for a reversal at
      ! the same rate, v2-fixed C47).
      ! (Step 1 of a restarted run: against the restart point's step, step 8.4.)
      if (i_incr == 1 .and. (i_step > 1 .or. restarted) .and. axis > 0) then
        if (control_rate_sign(step_schedule_record(sched, config, i_step))* &
          & control_rate_sign(previous_step(i_step)) < 0.0_rk) then
          state%slip_accum = 0.0_rk
          call log_info('Control rate reversed: accumulated slip (cyclic hardening) reset', 2)
        end if
      end if
      ! The increment (with sub-incrementation: subdivided where rejected), its geometry
      ! update included; for a load or stress target, redone shorter when it passes the target
      ! by more than its tolerance (landing rule, step 8.2a), from the start of the increment
      ! kept here. Then the time (the step's end time exactly at the last increment of a
      ! strain or time target).
      num_attempts = 0
      time_solved = time
      if (is_load) then
        coo_land = coo
        vel_land = vel
        vel_face_land = vel_face
        state_land = state
      end if
      call solve_increment(dt)
      if (is_load) call land_increment(i_step, dt)
      if (last_incr) then
        time = sched%time_end(i_step)
      else
        time = time + dt
      end if

      ! The first increment of a dwell episode stores its face velocities (v2's
      ! initial_load_dwell_vel, initial_unload_dwell_vel; step 8.2c).
      if (i_incr == 1 .and. start_dwell) then
        if (prev_action == bc_action_loading) vel_dwell_load = vel_face
        if (prev_action == bc_action_unloading) vel_dwell_unload = vel_face
      end if

      ! Centroid-point results and the face forces of the increment.
      stress_c = state%stress_s(:, elt%qpt_centroid, :)
      strain_el_c = state%strain_el_s(:, elt%qpt_centroid, :)
      rot_c = state%rot_c2s(:, :, elt%qpt_centroid, :)
      call face_forces(mesh, elt, coo, stress_c, loads, areas)
      time_mark = par_wtime()
      call sim_write_forces(sim, step_offset + i_step, incr, loads, areas, time)
      call write_macro(step_offset + i_step, incr)
      time_output = time_output + (par_wtime() - time_mark)

      ! Completion: a strain or time target at the step's last increment; a load or stress
      ! target when within its tolerance (landing rule; an increment past it by more was
      ! redone by land_increment), otherwise the next increment from the extrapolated value
      ! (v2 loading_timeincr).
      ! A strain target at a load or stress rate (step 8.2c): when reached (tolerance 0), else
      ! the next increment as a load target's, on the bounding-box strain.
      if (bc_on) strain_box = bc_lengths(mesh, coo)/lengths0 - 1.0_rk
      if (is_load) then
        value = control_value(i_step)
        step_done = step_load_landing(value, sched%target_value(i_step), &
          & sched%load_sign(i_step), sched%target_tol(i_step)) == landing_within
        if (.not. step_done) dt_next = step_load_next_dt(value, value_prev, dt, &
          & sched%target_value(i_step), sched%load_sign(i_step), sched%target_tol(i_step), &
          & sched%dt(i_step), sched%dt_min(i_step), sched%dtime_factor)
        value_prev = value
      else if (is_open) then
        value = strain_box(axis)
        step_done = step_load_complete(value, sched%target_value(i_step), &
          & sched%load_sign(i_step), 0.0_rk)
        if (.not. step_done) dt_next = step_load_next_dt(value, value_prev, dt, &
          & sched%target_value(i_step), sched%load_sign(i_step), 0.0_rk, sched%dt(i_step), &
          & sched%dt_min(i_step), sched%dtime_factor)
        value_prev = value
      else
        step_done = last_incr
      end if

      ! Strains of the prescribed components reached (the planned end strains at the end of a
      ! strain or time target); a completed load or stress target records its end.
      if (step_done .and. .not. is_open) then
        strain_now = sched%strain_end(:, i_step)
      else
        strain_now = sched%strain_start(:, i_step) + sched%strain_rate(:, i_step)* &
          & (time - sched%time_start(i_step))
      end if
      ! Load-controlled components: the strain of the bounding box (step 8.2b).
      if (bc_on) then
        where (sched%bc_kind(:, i_step) /= 0) strain_now = strain_box
      end if
      if (step_done .and. is_open) then
        sched%time_end(i_step) = time
        sched%strain_end(:, i_step) = strain_now
        sched%num_incrs(i_step) = i_incr
        call log_load_end(i_step, value)
      end if
      ! A dwell records its increments; a step at a load or stress rate logs its end (step 8.2c).
      if (step_done .and. sched%is_dwell(i_step)) sched%num_incrs(i_step) = i_incr
      if (step_done .and. sched%is_rate(i_step) .and. .not. is_open) call log_rate_end(i_step)

      ! The kind of a completed step, for the episode starts of iteration (B) (step 8.2c): its
      ! own at a load or stress rate, else by its measured control-face force.
      if (step_done .and. any(sched%is_rate)) then
        if (sched%is_rate(i_step)) then
          prev_action = action
        else
          prev_action = bc_step_action(loads(axis, i_face) - force_start, force_start, &
            & loads(axis, i_face), .true., config%load_tol)
        end if
      end if

      ! End of the step: its results when printed.
      written = step_done .and. sched%print_results(i_step)
      time_mark = par_wtime()
      if (written) call write_step(i_step)

      ! The restart point of a completed step (print restart, step 8.4).
      if (step_done .and. config%print_restart) call restart_write(i_step)
      time_output = time_output + (par_wtime() - time_mark)

      ! Necking (check_necking on, addendum §3.6, v2-fixed C54): after every increment.
      if (config%check_necking) call check_necking(i_step, written)

      ! Every rank's buffered log lines (all_ranks), then completion, then the run limits
      ! (D44, v2-fixed C32).
      call par_log_flush()
      if (i_step == sched%num_steps .and. step_done) exit
      call check_limits(i_step, written)
    end do
  end do
  time_loop_end = par_wtime()
  time_output_run = time_output - time_output_setup

  time_mark = par_wtime()
  call sim_close(sim)
  time_output = time_output + (par_wtime() - time_mark)

  ! The closing summary (step 9.7, dev/plan/terminal_output.md §3.4, Q63).
  call par_log_flush()
  call par_sum(int(log_warning_count(), ik), num_warnings)
  call log_closing_summary()
  call par_finalize()

contains

  !> Per-elset or per-element values of a crystal-plasticity phase from simulation.opt
  !> (step 7.1, D59, addendum_hardening_modifiers.md §6): for each variable (Crss: tau_0,
  !> CrssSat: tau_s, HardCoeff: h_0; with Kocks-Mecking hardening only Rho0: element_values_km)
  !> the per-elset or the per-element section, not both, each
  !> covering every elset or element (mesh_file_mod); the values of the phase's local elements
  !> replace the phase's values. Each entity gives one value, one per active family or one per
  !> slip system (tau_s and h_0 one value with saturation_isotropic, and tau_s one value with
  !> saturation_evolution on). Strengths from the file are used as given (no precipitate
  !> strength, a warning when precipitate_strengthening is on); values for other phases are
  !> ignored, with a warning (opt_warnings). Checks per element: tau_0 > 0, h_0 >= 0, tau_s > tau_0
  !> (system 1 with saturation_isotropic, every system otherwise).
  !>
  !> Arguments:
  !>   i_phase  (in)    integer           phase (crystal plasticity)
  !>   cp       (inout) type(cp_model_t)  its model (per-element values set)
  !> Locals:
  !>   names    character(10)(3)          the variables (Crss, CrssSat, HardCoeff)
  !>   var, i_sec, k, nv  integer         variable, section, entity index, values given
  !>   i_e      integer(ik)               local element
  !>   arrays   real(rk)(num_slip, num_elts, 3)  per-element values of the three variables
  !>   given    logical(3)                a section is given for the variable
  !>   single   logical                   only one value per entity is accepted
  !>   ok       logical                   the count matched
  !>   tau_0, tau_s, h  real(rk)(num_slip)  an element's initial and saturation strengths and h_0
  !>                                        [MPa]
  !>   entity   character(:), allocatable "elset <id>" or "element <id>", for messages
  !>   counts   character(:), allocatable the accepted counts, for messages
  subroutine element_values(i_phase, cp)

    integer, intent(in) :: i_phase
    type(cp_model_t), intent(inout) :: cp

    character(len=10), parameter :: names(3) = [character(len=10) :: 'Crss', 'CrssSat', &
      & 'HardCoeff']
    integer :: var, i_sec, nv
    integer(ik) :: k, i_e
    real(rk), allocatable :: arrays(:, :, :)
    logical :: given(3), single, ok
    real(rk) :: tau_0(cp%num_slip), tau_s(cp%num_slip), h(cp%num_slip)
    character(len=:), allocatable :: entity, counts

    ! Kocks-Mecking: only rho_0 is read from the file (step 7.2, D60).
    if (config%phases(i_phase)%hardening == 'kocks_mecking') then
      call element_values_km(i_phase, cp)
      return
    end if

    associate (ph => config%phases(i_phase))
      allocate (arrays(cp%num_slip, mesh%num_elts, 3))
      arrays = 0.0_rk
      given = .false.
      do var = 1, 3
        ! One form per variable.
        if (mesh%values(var)%defined .and. mesh%values(var + 3)%defined) call fatal_error( &
          & 'fepx', 'simulation.opt: give $Elset' // trim(names(var)) // ' or $Element' // &
          & trim(names(var)) // ', not both')
        i_sec = 0
        if (mesh%values(var)%defined) i_sec = var
        if (mesh%values(var + 3)%defined) i_sec = var + 3
        if (i_sec == 0) cycle
        given(var) = .true.
        single = var > 1 .and. (ph%hardening == 'saturation_isotropic' .or. &
          & (var == 2 .and. ph%saturation_evolution))

        ! The values of each local element of the phase, expanded to one per system.
        do i_e = 1_ik, mesh%num_elts
          if (mesh%elt_phase(i_e) /= int(i_phase, ik)) cycle
          if (i_sec <= 3) then
            k = mesh%elt_elset(i_e)
            entity = 'elset ' // int_str(mesh%elset_id(k))
          else
            k = i_e
            entity = 'element ' // int_str(mesh%elt_id_glob(i_e))
          end if
          nv = int(mesh%values(i_sec)%num_values(k))
          call expand_values(ph, mesh%values(i_sec)%values(1:nv, k), single .or. var == 1, &
            & single, arrays(:, i_e, var), ok)
          if (.not. ok) then
            counts = 'one value'
            if (.not. single) counts = 'one value, one per active slip family or one per ' // &
              & 'slip system'
            call fatal_error('fepx', 'simulation.opt: ' // trim(names(var)) // ' of ' // &
              & entity // ' (phase ' // int_str(int(i_phase, ik)) // '): ' // &
              & int_str(int(nv, ik)) // ' value(s); give ' // counts)
          end if
        end do
      end do

      ! Warning: strengths from the file without the precipitate strength (the values of other
      ! phases: opt_warnings).
      if (given(1) .and. ph%precipitation) call log_warning('Initial strengths from ' // &
        & 'simulation.opt are used as given; precipitate_strengthening is ignored for them.', 3)

      ! Checks per element (when the file gives values), then the values in the model.
      if (.not. any(given)) return
      do i_e = 1_ik, mesh%num_elts
        if (mesh%elt_phase(i_e) /= int(i_phase, ik)) cycle
        tau_0 = cp%tau_0
        if (given(1)) tau_0 = arrays(:, i_e, 1)
        entity = 'element ' // int_str(mesh%elt_id_glob(i_e))
        if (.not. all(tau_0 > 0.0_rk)) call fatal_error('fepx', 'simulation.opt: ' // entity // &
          & ': initial strengths must be positive')
        if (ph%hardening == 'none') cycle
        tau_s = ph%tau_s
        h = ph%h_0
        if (given(2)) tau_s = arrays(:, i_e, 2)
        if (given(3)) h = arrays(:, i_e, 3)
        if (any(h < 0.0_rk)) call fatal_error('fepx', 'simulation.opt: ' // entity // &
          & ': hardening coefficients must not be negative')
        if (ph%hardening == 'saturation_isotropic') then
          ok = tau_s(1) > tau_0(1)
        else
          ok = all(tau_s > tau_0)
        end if
        if (.not. ok) call fatal_error('fepx', 'simulation.opt: ' // entity // &
          & ': the saturation strength must exceed the initial strength')
      end do
      if (given(1)) call cp_model_set_element_values(cp, tau_0=arrays(:, :, 1))
      if (given(2)) call cp_model_set_element_values(cp, tau_s=arrays(:, :, 2))
      if (given(3)) call cp_model_set_element_values(cp, h_0=arrays(:, :, 3))
    end associate

  end subroutine element_values

  !> Per-elset or per-element values of a Kocks-Mecking phase from simulation.opt: the initial
  !> network dislocation density (step 7.2, D60, addendum_kocks_mecking.md §6: $ElsetRho0 or
  !> $ElementRho0, > 0 [1/mm^2], replacing km_rho_0) and, with irradiation on (step 7.3, D61,
  !> addendum_irradiation.md §6), the loop number density ($ElsetIrrLoopDensity or
  !> $ElementIrrLoopDensity, >= 0 [1/mm^3]) and the loop diameter ($ElsetIrrLoopDiameter or
  !> $ElementIrrLoopDiameter, > 0 [mm]), replacing irr_loop_density and irr_loop_diameter; a
  !> variable not in the file keeps the cfg value for every element. The other sections are
  !> ignored for this phase (opt_warnings).
  !>
  !> Arguments:
  !>   i_phase  (in)    integer           phase (crystal plasticity, hardening kocks_mecking)
  !>   cp       (inout) type(cp_model_t)  its model (per-element values set)
  !> Locals:
  !>   values   real(rk)(num_elts)        per-element values of one variable
  !>   given    logical                   the variable's section is given
  subroutine element_values_km(i_phase, cp)

    integer, intent(in) :: i_phase
    type(cp_model_t), intent(inout) :: cp

    real(rk), allocatable :: values(:)
    logical :: given

    ! Initial density.
    call opt_scalar(i_phase, 7, 'Rho0', 'the initial dislocation density', .true., values, given)
    if (given) call cp_model_set_element_values(cp, rho_0=values)

    ! Irradiation: loop densities and diameters.
    if (.not. config%phases(i_phase)%irradiation) return
    call opt_scalar(i_phase, 9, 'IrrLoopDensity', 'the loop density', .false., values, given)
    if (given) call cp_model_set_element_values(cp, loop_density=values)
    call opt_scalar(i_phase, 11, 'IrrLoopDiameter', 'the loop diameter', .true., values, given)
    if (given) call cp_model_set_element_values(cp, loop_diameter=values)

  end subroutine element_values_km

  !> One scalar variable of simulation.opt for the local elements of a phase (steps 7.2, 7.3):
  !> the per-elset section i_sec or the per-element section i_sec + 1, not both, covering every
  !> elset or element (mesh_file_mod), one value per entity, > 0 (positive) or >= 0.
  !>
  !> Arguments:
  !>   i_phase   (in)  integer            phase
  !>   i_sec     (in)  integer            the per-elset section (the per-element one follows)
  !>   name      (in)  character          variable name in the section names (e.g. Rho0)
  !>   what      (in)  character          its description, for messages
  !>   positive  (in)  logical            values must be > 0 (otherwise >= 0)
  !>   values    (out) real(rk)(num_elts), allocatable  the values (0 for other phases)
  !>   given     (out) logical            a section is given
  !> Locals:
  !>   sec      integer                   the section used
  !>   k, i_e   integer(ik)               entity index, local element
  !>   nv       integer                   values given for the entity
  !>   entity   character(:), allocatable "elset <id>" or "element <id>", for messages
  subroutine opt_scalar(i_phase, i_sec, name, what, positive, values, given)

    integer, intent(in) :: i_phase, i_sec
    character(len=*), intent(in) :: name, what
    logical, intent(in) :: positive
    real(rk), allocatable, intent(out) :: values(:)
    logical, intent(out) :: given

    integer :: sec, nv
    integer(ik) :: k, i_e
    character(len=:), allocatable :: entity

    ! One form; none given: the cfg value.
    if (mesh%values(i_sec)%defined .and. mesh%values(i_sec + 1)%defined) call fatal_error( &
      & 'fepx', 'simulation.opt: give $Elset' // name // ' or $Element' // name // ', not both')
    sec = 0
    if (mesh%values(i_sec)%defined) sec = i_sec
    if (mesh%values(i_sec + 1)%defined) sec = i_sec + 1
    given = sec > 0
    allocate (values(mesh%num_elts), source=0.0_rk)
    if (.not. given) return

    ! One value in range for each local element of the phase.
    do i_e = 1_ik, mesh%num_elts
      if (mesh%elt_phase(i_e) /= int(i_phase, ik)) cycle
      if (sec == i_sec) then
        k = mesh%elt_elset(i_e)
        entity = 'elset ' // int_str(mesh%elset_id(k))
      else
        k = i_e
        entity = 'element ' // int_str(mesh%elt_id_glob(i_e))
      end if
      nv = int(mesh%values(sec)%num_values(k))
      if (nv /= 1) call fatal_error('fepx', 'simulation.opt: ' // name // ' of ' // entity // &
        & ' (phase ' // int_str(int(i_phase, ik)) // '): ' // int_str(int(nv, ik)) // &
        & ' value(s); give one value')
      values(i_e) = mesh%values(sec)%values(1, k)
      if (positive .and. .not. values(i_e) > 0.0_rk) call fatal_error('fepx', &
        & 'simulation.opt: ' // entity // ': ' // what // ' ' // name // ' must be positive')
      if (.not. positive .and. values(i_e) < 0.0_rk) call fatal_error('fepx', &
        & 'simulation.opt: ' // entity // ': ' // what // ' ' // name // &
        & ' must not be negative')
    end do

  end subroutine opt_scalar

  !> Warnings for simulation.opt values a phase does not use (steps 7.1-7.3, D59 (4), D60,
  !> D61): any value section for an elastic phase; Crss, CrssSat or HardCoeff for a
  !> Kocks-Mecking phase; Rho0 for a crystal-plasticity phase with another law; IrrLoopDensity
  !> or IrrLoopDiameter for a phase without irradiation. Only phases with elements are named;
  !> one warning per phase and kind.
  !>
  !> Locals:
  !>   i_p        integer   phase
  !>   strengths  logical   a Crss, CrssSat or HardCoeff section is given
  !>   densities  logical   a Rho0 section is given
  !>   loops      logical   an IrrLoopDensity or IrrLoopDiameter section is given
  subroutine opt_warnings()

    integer :: i_p
    logical :: strengths, densities, loops

    strengths = any(mesh%values(1:6)%defined)
    densities = any(mesh%values(7:8)%defined)
    loops = any(mesh%values(9:12)%defined)
    do i_p = 1, config%num_phases
      if (.not. any(mesh%elset_phase == int(i_p, ik))) cycle
      associate (ph => config%phases(i_p))
        if (ph%model /= 'crystal_plasticity') then
          if (strengths .or. densities .or. loops) call log_warning('simulation.opt gives ' // &
            & 'plasticity parameters for phase ' // int_str(int(i_p, ik)) // &
            & ', which is elastic: ignored.', 1)
          cycle
        else if (ph%hardening == 'kocks_mecking') then
          if (strengths) call log_warning('simulation.opt gives Crss, CrssSat or HardCoeff ' // &
            & 'values for phase ' // int_str(int(i_p, ik)) // ', whose hardening is ' // &
            & 'kocks_mecking (only Rho0 and the irradiation loop sections are read): ignored.', 1)
        else
          if (densities) call log_warning('simulation.opt gives Rho0 values for phase ' // &
            & int_str(int(i_p, ik)) // ', whose hardening is not kocks_mecking: ignored.', 1)
        end if
        if (loops .and. .not. ph%irradiation) call log_warning('simulation.opt gives ' // &
          & 'IrrLoopDensity or IrrLoopDiameter values for phase ' // int_str(int(i_p, ik)) // &
          & ', which has no irradiation: ignored.', 1)
      end associate
    end do

  end subroutine opt_warnings

  !> Values of one entity of simulation.opt expanded to one per slip system (the counting rule
  !> of dev/cfg_keywords.md §2.1): one value (when allowed), one per active family (in the
  !> slip_families order) or one per slip system.
  !>
  !> Arguments:
  !>   ph        (in)  type(cfg_phase_t)  the phase (families, num_slip)
  !>   vals      (in)  real(rk)(:)        the values given
  !>   one_ok    (in)  logical            one value for all systems is accepted
  !>   only_one  (in)  logical            only one value is accepted
  !>   out       (out) real(rk)(num_slip) one value per system
  !>   ok        (out) logical            the count is accepted
  !> Locals:
  !>   i, a, n   integer                  family, system of the family, system counter
  subroutine expand_values(ph, vals, one_ok, only_one, out, ok)

    type(cfg_phase_t), intent(in) :: ph
    real(rk), intent(in) :: vals(:)
    logical, intent(in) :: one_ok, only_one
    real(rk), intent(out) :: out(:)
    logical, intent(out) :: ok

    integer :: i, a, n

    ok = .true.
    out = 0.0_rk
    if (size(vals) == 1 .and. one_ok) then
      out = vals(1)
    else if (only_one) then
      ok = .false.
    else if (size(vals) == size(ph%families)) then
      n = 0
      do i = 1, size(ph%families)
        do a = 1, crystal_family_size(crystal_structure_from_name(ph%crystal_type), &
          & ph%families(i))
          n = n + 1
          out(n) = vals(i)
        end do
      end do
    else if (size(vals) == ph%num_slip) then
      out = vals
    else
      ok = .false.
    end if

  end subroutine expand_values

  !> Velocity conditions of a step (bc_mask, bc_vel) and the initial velocity guess of its
  !> first increment: the prescribed values at the first increment of the run; otherwise the
  !> previous velocities scaled by the ratio of the step strain rates (v2 vel_factor; restarted
  !> from the prescribed values after a step with zero rate; when either step has
  !> load-controlled components, as they are, or with triaxial conditions and the control
  !> component at a strain rate in both steps times the ratio of its strain rates, the
  !> load-controlled faces' velocities included: v2 driver_triaxcsr, addendum §3.4 revised,
  !> D86), with the prescribed values on the prescribed dofs;
  !> a run whose step 1 has load-controlled components starts from v2's homogeneous field of the
  !> face velocities (step 8.2b).
  !>
  !> Arguments:
  !>   i_step  (in)  integer  step
  !> Locals:
  !>   preset  integer        bc_mod preset code
  !>   ratio   real(rk)       scale of the previous velocities [-]
  !>   ratio_control real(rk) triaxial: ratio of the control component's strain rates, 1 unless
  !>                          it is at a strain rate in both steps (v2 vel_factor) [-]
  !>   ratios  real(rk)(3)    ratio of the step strain rates of each prescribed component [-]
  !>   ratios9 real(rk)(9)    periodic: ratio of the imposed rates of Fbar-dot [-]
  !>   lengths real(rk)(3)    original lengths along the axes (gage lengths) [mm]
  !>   j       integer        condition index, component
  !>   i_node  integer(ik)    local node
  !>   prev    type(step_record_t)  the previous step (of step 1 after a restart point: the
  !>                          restart point's step, step 8.4)
  subroutine step_conditions(i_step)

    integer, intent(in) :: i_step

    integer :: preset, j
    integer(ik) :: i_node
    real(rk) :: ratio, ratios(3), lengths(3), ratios9(9), ratio_control
    type(step_record_t) :: prev

    ! The conditions of the step.
    if (allocated(specs)) deallocate (specs)
    select case (config%bc_type)
    case ('uniaxial_grip', 'uniaxial_symmetry', 'uniaxial_minimal')
      ! Preset along the loading axis, moving face at the step's strain rate times the original
      ! gage length (v2).
      select case (config%bc_type)
      case ('uniaxial_grip')
        preset = bc_preset_grip
      case ('uniaxial_symmetry')
        preset = bc_preset_symmetry
      case default
        preset = bc_preset_minimal
      end select
      ! (A load or stress rate, step 8.2c: the face at its current velocity, iterated.)
      if (sched%is_rate(i_step)) then
        call bc_preset_specs(preset, axis, vel_face(axis), specs)
      else
        call bc_preset_specs(preset, axis, sched%strain_rate(axis, i_step)* &
          & bc_gage_length(mesh, axis), specs)
      end if
    case ('triaxial')
      ! Symmetry faces x0, y0, z0; each prescribed face x1, y1, z1 at its component's strain
      ! rate times the original length; free faces without a condition (D45); a
      ! load-controlled face at its current velocity (step 8.2b: carried from the previous
      ! increment).
      lengths = bc_lengths(mesh, mesh%coo_s)
      ! At a step change with the control component at a strain rate in both steps, the
      ! load-controlled faces' velocities times the ratio of its strain rates (v2 vel_factor,
      ! driver_triaxcsr; addendum §3.4 revised, D86).
      ratio_control = triaxial_control_ratio(i_step)
      where (sched%bc_kind(:, i_step) /= 0) vel_face = ratio_control*vel_face
      where (sched%bc_kind(:, i_step) == 0) vel_face = sched%strain_rate(:, i_step)*lengths
      call bc_triaxial_specs(vel_face, sched%prescribed, specs)
    case ('periodic')
      ! Periodic conditions (step 8.3b, Q-P6): the set_bc vel conditions as given (v2).
      allocate (specs(size(config%vel_bcs)))
      do j = 1, size(config%vel_bcs)
        specs(j)%nset = config%vel_bcs(j)%nset
        specs(j)%comp = config%vel_bcs(j)%comp
        specs(j)%value = config%vel_bcs(j)%value
      end do
    case default
      ! general: the set_bc vel conditions times the step's velocity factor (v2).
      allocate (specs(size(config%vel_bcs)))
      do j = 1, size(config%vel_bcs)
        specs(j)%nset = config%vel_bcs(j)%nset
        specs(j)%comp = config%vel_bcs(j)%comp
        specs(j)%value = sched%vel_factor(i_step)*config%vel_bcs(j)%value
      end do
    end select
    call bc_resolve(mesh, specs, bc_mask, bc_vel)
    if (config%bc_type == 'periodic') then
      ! Periodic conditions: no condition on a constrained dof (Q-P6), the offsets of the
      ! step's imposed rates (Q-P5).
      call constraint_check_bcs(cons, mesh, bc_mask)
      call constraint_set_rates(cons, config%per_rate(:, i_step))
    else if (cons%active) then
      ! Multi-point constraints: the conditions carried through their groups (C66).
      call constraint_propagate_bcs(cons, mesh, bc_mask, bc_vel)
    end if

    ! Initial guess (v2): the prescribed values at the start; then the previous velocities
    ! scaled by the ratio of the step rates when every prescribed component changes by the
    ! same ratio, the previous velocities as they are when the ratios differ (triaxial), and
    ! the prescribed values again after a step with a zero rate.
    ! (Step 1 of a restarted run: against the restart point's step, step 8.4.)
    ratio = 0.0_rk
    if (i_step == 1 .and. init_done .and. .not. restarted) ratio = 1.0_rk
    if (i_step > 1 .or. restarted) then
      prev = previous_step(i_step)
      if (any(sched%bc_kind(:, i_step) /= 0) .or. any(prev%bc_kind /= 0)) then
        ! Load-controlled components in either step: the previous velocities (addendum §3.4),
        ! triaxial: times the ratio of the control component's strain rates (revised, D86).
        ratio = 1.0_rk
        if (config%bc_type == 'triaxial') ratio = triaxial_control_ratio(i_step)
      else if (config%bc_type == 'periodic') then
        ! Periodic conditions: the ratio of the imposed rates when they all change by it (and
        ! none was 0), else the previous velocities as they are.
        ratio = 1.0_rk
        if (all(abs(prev%per_rate) > 0.0_rk .or. .not. config%per_imposed)) then
          ratios9 = 1.0_rk
          where (config%per_imposed) ratios9 = config%per_rate(:, i_step)/prev%per_rate
          ratio = ratios9(config%per_control)
          if (any(abs(ratios9 - ratio) > 1.0e-12_rk*abs(ratio) .and. config%per_imposed)) &
            & ratio = 1.0_rk
        end if
      else if (axis == 0) then
        ratio = sched%vel_factor(i_step)/prev%vel_factor
      else if (all(abs(prev%strain_rate) > 0.0_rk .or. .not. sched%prescribed)) then
        ratios = 1.0_rk
        where (sched%prescribed) ratios = sched%strain_rate(:, i_step)/prev%strain_rate
        ratio = ratios(axis)
        if (any(abs(ratios - ratio) > 1.0e-12_rk*abs(ratio) .and. sched%prescribed)) &
          & ratio = 1.0_rk
      end if
    end if
    vel = ratio*vel

    ! Triaxial load control from the start of the run (step 8.2b): v2's homogeneous field of the
    ! face velocities, V_c (x_c - x_min,c)/L_c on every prescribed component (driver_triaxcsr).
    if (i_step == 1 .and. .not. init_done .and. .not. restarted .and. &
      & any(sched%bc_kind(:, 1) /= 0)) then
      call box_bounds(coo, box_lo, box_len)
      do j = 1, 3
        if (.not. sched%prescribed(j)) cycle
        do i_node = 1_ik, mesh%num_nodes_loc
          vel(3*(i_node - 1) + j) = vel_face(j)*(coo(j, i_node) - box_lo(j))/box_len(j)
        end do
      end do
    end if
    where (bc_mask) vel = bc_vel
    ! Multi-point constraints and periodicity: the guess projected onto them (v2's uniaxial
    ! driver; periodic: v2's offset field at the start); periodic: the zero-mean fluctuation.
    if (cons%active) call constraint_project(cons, vel)
    if (cons%zero_mean) call constraint_zero_mean(cons, mesh, vel)

  end subroutine step_conditions

  !> Triaxial conditions at a step change: the ratio of the control component's strain rates
  !> when it is at a strain rate in this step and the previous one (v2 vel_factor =
  !> step_velocity(k)/step_velocity(k-1), driver_triaxcsr; addendum §3.4 revised, D86); 1 at
  !> the first step of a run, when either step has the control component at a load or stress
  !> rate (v2 CLR keeps its velocities), or after a step with a zero strain rate.
  !>
  !> Arguments:
  !>   i_step  (in)  integer   step
  !>   ratio   (out, result) real(rk)  the factor on the previous velocities [-]
  !> Locals:
  !>   prev    type(step_record_t)  the previous step (of step 1 after a restart point: the
  !>                          restart point's step, step 8.4)
  function triaxial_control_ratio(i_step) result(ratio)

    integer, intent(in) :: i_step
    real(rk) :: ratio

    type(step_record_t) :: prev

    ! No previous step: no change.
    ratio = 1.0_rk
    if (i_step == 1 .and. .not. restarted) return

    ! The control component at a strain rate in both steps, the previous one nonzero.
    prev = previous_step(i_step)
    if (sched%bc_kind(axis, i_step) /= 0 .or. prev%bc_kind(axis) /= 0) return
    if (.not. abs(prev%strain_rate(axis)) > 0.0_rk) return
    ratio = sched%strain_rate(axis, i_step)/prev%strain_rate(axis)

  end function triaxial_control_ratio

  !> Duration of the initialization step (D49; program header), 0 when there is none: no
  !> crystal-plasticity phase, initialization not elastic_step, or no loading in step 1.
  !>
  !> Locals:
  !>   tau_0_min    real(rk)  smallest tau_0 of the crystal-plasticity phases [MPa]
  !>   stiff_bound  real(rk)  largest absolute row sum of the phases' [C] [MPa]
  !>   rate         real(rk)  fastest strain rate of step 1 [1/s]
  !>   dt_first     real(rk)  first increment of step 1 [s]
  !>   j            integer   phase
  function initialization_dt() result(dt_out)

    real(rk) :: dt_out

    real(rk) :: tau_0_min, stiff_bound, rate, dt_first
    integer :: j

    dt_out = 0.0_rk
    if (config%initialization /= 'elastic_step') return
    if (.not. any(config%phases(:)%model == 'crystal_plasticity')) return
    ! A load or stress rate in step 1 (step 8.2c): its first increment starts from v2's elastic
    ! estimate of the face velocities.
    if (sched%is_rate(1)) then
      call log_info('No elastic step (step 1 has a load or stress rate)', 1)
      return
    end if

    ! Expected yield 3 tau_0 (smallest tau_0) and the stiffness bound over the phases.
    tau_0_min = huge(1.0_rk)
    stiff_bound = 0.0_rk
    do j = 1, config%num_phases
      if (config%phases(j)%model == 'crystal_plasticity') &
        & tau_0_min = min(tau_0_min, minval(config%phases(j)%tau_0))
      stiff_bound = max(stiff_bound, maxval(sum(abs(config%phases(j)%stiff_el_c), dim=2)))
    end do

    ! Fastest strain rate of step 1.
    if (config%bc_type == 'general') then
      rate = maxval(abs(config%vel_bcs(:)%value))*abs(sched%vel_factor(1))/minval(lengths0)
    else if (config%bc_type == 'periodic') then
      ! Periodic: the largest imposed rate of Fbar-dot.
      rate = maxval(abs(config%per_rate(:, 1)))
    else
      rate = maxval(abs(sched%strain_rate(:, 1)), mask=sched%prescribed)
    end if
    if (.not. (rate > 0.0_rk)) then
      call log_info('No elastic step (step 1 has no loading rate)', 1)
      return
    end if

    ! 1/100 of 3 tau_0 over the stiffness bound, at the rate; at most half the first increment.
    dt_out = 0.01_rk*3.0_rk*tau_0_min/stiff_bound/rate
    dt_first = sched%dt(1)
    if (sched%num_incrs(1) == 1) dt_first = sched%dt_last(1)
    if (.not. (dt_out < dt_first)) dt_out = 0.5_rk*dt_first

  end function initialization_dt

  !> Solve the interval dt_sub of the current increment from the current coordinates, state
  !> and velocities: mid-nodes reset when asked (D26), increment_solve, geometry update
  !> x = x + dt_sub v. With sub-incrementation (Q35, D63), a rejected attempt (a point solve
  !> not converged at the converged velocities, or the global iteration failed with
  !> global_max_iters or diverging) is undone and the interval solved as two halves, each by
  !> this routine (depth + 1); a rejection at depth subincrement_max_depth (L1), or the first
  !> rejection of an increment after subincrement_max_successive successive subdivided
  !> increments (L2), ends the run with an error. Without sub-incrementation: one attempt, its
  !> failures handled by report_increment (D49 (2): failed points keep their last iterate).
  !> With load-controlled components (step 8.2b) the attempt is the boundary-velocity iteration
  !> (bc_iterate), any of whose solves can reject it; the face velocities are undone with it.
  !>
  !> Source : addendum_subincrementation.md §2 (trigger), §3 (procedure solve_interval), §4
  !>          (limits L1, L2 and their messages), §5 (log); Baeker (2018), arXiv:1811.05753,
  !>          §2.1-1 and §2.2-2 to §2.2-4 (cutback of a non-converged increment, a minimum
  !>          increment and a limited number of cutbacks ending the run with an error).
  !>
  !> Arguments:
  !>   dt_sub  (in)  real(rk)  the interval: dt / 2^depth [s]
  !>   depth   (in)  integer   subdivisions so far
  !> Locals:
  !>   coo_save    real(rk)(:, :)          (3, num_nodes_loc) coordinates at the start of the
  !>                                       interval [mm]
  !>   vel_save    real(rk)(:)             (3 num_nodes_loc) velocities (initial guess) at the
  !>                                       start of the interval [mm/s]
  !>   vel_face_save  real(rk)(3)          face velocities at the start of the interval [mm/s]
  !>   state_save  type(material_state_t)  state at the start of the interval
  !>   rejected    logical                 the attempt is rejected (trigger T1 or T2)
  !>   reason      character(120)          why, for the messages
  !>   message     character(400)          message buffer
  recursive subroutine solve_interval(dt_sub, depth)

    real(rk), intent(in) :: dt_sub
    integer, intent(in) :: depth

    real(rk), allocatable :: coo_save(:, :), vel_save(:)
    real(rk) :: vel_face_save(3)
    type(material_state_t) :: state_save
    logical :: rejected
    character(len=120) :: reason
    character(len=400) :: message

    ! The start of the interval, to undo a rejected attempt (the state only with
    ! sub-incrementation: the only case that undoes an attempt).
    coo_save = coo
    vel_save = vel
    vel_face_save = vel_face
    if (config%subincrementation) state_save = state

    ! The attempt: mid-nodes reset when asked (D26), then the increment over dt_sub: one solve,
    ! or with load-controlled components the boundary-velocity iteration (step 8.2b).
    if (config%midnode_linearization) call midnode_linearize(mesh, elt, coo)
    state_prev = state
    num_attempts = num_attempts + 1
    if (bc_on .and. sched%is_rate(i_step)) then
      call bc_iterate_rate(dt_sub, rejected, reason)
    else if (bc_on) then
      call bc_iterate(dt_sub, rejected, reason)
    else
      call increment_solve(mesh, elt, phases, coo, bc_mask, dt_sub, opts, state_prev, state, &
        & vel, info, cons)
      rejected = attempt_rejected(reason)
      call report_increment(num_attempts, .not. rejected)
    end if

    ! Accepted: the geometry update; later sub-increments start from the previous stress.
    if (.not. rejected) then
      coo = coo + dt_sub*reshape(vel, [3_ik, mesh%num_nodes_loc])
      time_solved = time_solved + dt_sub
      state%start_from_vertex = .false.
      bc_at_start = .false.
      num_subincrs = num_subincrs + 1
      dt_smallest = min(dt_smallest, dt_sub)
      return
    end if

    ! Rejected: the limits (addendum §4), then the two halves from the undone attempt.
    if (depth >= config%subincrement_max_depth) then
      write (message, '(a,i0,a,i0,a,i0,a,es12.5,a)') 'increment ', incr, ': ' // trim(reason) // &
        & ' after ', config%subincrement_max_depth, ' subdivisions (sub-increment dt/2^', &
        & config%subincrement_max_depth, ' = ', dt_sub, ' s); refine the time stepping ' // &
        & '(number_of_incrs, dtime or dstrain) or raise subincrement_max_depth'
      call fatal_error('fepx', trim(message))
    end if
    if (num_attempts == 1) then
      if (num_successive >= config%subincrement_max_successive) then
        write (message, '(i0,a,i0,a)') num_successive + 1, ' successive increments needed ' // &
          & 'sub-incrementation (increment ', incr, '; ' // trim(reason) // &
          & '); refine the time stepping (number_of_incrs, dtime or dstrain) or raise ' // &
          & 'subincrement_max_successive'
        call fatal_error('fepx', trim(message))
      end if
      num_successive = num_successive + 1
    end if
    call log_info('Attempt over ' // log_real_text(dt_sub, 4) // ' s rejected (' // &
      & trim(reason) // '); subdividing (depth ' // log_int_text(depth + 1) // ')', 2)
    coo = coo_save
    vel = vel_save
    vel_face = vel_face_save
    state = state_save
    call solve_interval(0.5_rk*dt_sub, depth + 1)
    call solve_interval(0.5_rk*dt_sub, depth + 1)

  end subroutine solve_interval

  !> Whether the solve just made rejects the attempt (sub-incrementation trigger,
  !> addendum_subincrementation.md §2): T1, a point solve not converged at the converged
  !> velocities (num_point_failed is summed over the ranks: the same decision on every rank);
  !> T2, a global failure other than the linear solver's. Never without sub-incrementation.
  !>
  !> Arguments:
  !>   reason    (out) character(120)  why, for the messages (empty when not rejected)
  !>   rejected  (out, result) logical the attempt is rejected
  function attempt_rejected(reason) result(rejected)

    character(len=120), intent(out) :: reason
    logical :: rejected

    rejected = .false.
    reason = ''
    if (.not. config%subincrementation) return
    if (info%status == increment_converged .and. info%num_point_failed > 0_ik) then
      rejected = .true.
      write (reason, '(i0,a)') info%num_point_failed, ' point solve(s) did not converge ' // &
        & 'at the converged velocities'
    else if (info%status == increment_max_iters) then
      rejected = .true.
      reason = 'the global iteration did not converge in global_max_iters iterations'
    else if (info%status == increment_diverging) then
      rejected = .true.
      reason = 'the global iteration is diverging'
    end if

  end function attempt_rejected

  !> Boundary-velocity iteration (A) of an increment attempt over dt_sub (step 8.2b; addendum
  !> §3.4; v2 driver_triaxcsr's loop bc_iteration): at most max_bc_iter iterations, each a solve
  !> at the current face velocities (trial_solve), its values compared with the values to reach
  !> (bc_mode_values: on the load path, in time at the attempt's end time, or by a rate, step
  !> 8.2c; path and other rows together: the path rows on the path, the others at their
  !> values);
  !> within the tolerances the solve is the attempt (its state and velocities kept); otherwise one
  !> solve with each load-controlled face velocity perturbed by s_k delta_k (bc_pert_signs), the
  !> correction from v2's system (bc_correction), and the new perturbation sizes
  !> (bc_pert_update). Every solve starts from the start of the attempt. A rejected solve
  !> (sub-incrementation) ends the iteration with the attempt rejected; no convergence in
  !> max_bc_iter iterations, or a singular system, ends the run with an error.
  !>
  !> Arguments:
  !>   dt_sub    (in)  real(rk)        length of the attempt [s]
  !>   rejected  (out) logical         a solve rejected the attempt
  !>   reason    (out) character(120)  why, for the messages
  !> Locals:
  !>   n, m      integer               load-controlled components, values compared (n, + 1 with
  !>                                   the path's control value)
  !>   iter, k, c  integer             iteration, index in bc_comps, component
  !>   status    integer               solve status of the correction
  !>   q, q_pert real(rk)(4)           values at the iteration's velocities and with a
  !>                                   perturbation: the load-controlled components, then the
  !>                                   control (path) [N] or [MPa]
  !>   q_ideal   real(rk)(3)           values to reach [N] or [MPa]
  !>   q_aim, ratio  real(rk)(4)       right-hand side values and path directions of the system
  !>   signs, pert, dvel, tol  real(rk)(3)  perturbation signs, signed perturbations [mm/s],
  !>                                   velocity corrections [mm/s], tolerances [N] or [MPa]
  !>   frac      real(rk)              fraction of the step at the attempt's end (time-based)
  !>   elapsed   real(rk)              time since the step's start at the attempt's end [s]
  !>   lambda    real(rk)              progress of the control along its step (path) [-]
  !>   modes     integer(3)            how each row's value to reach is given (bc_mode_*)
  !>   dq        real(rk)(:, :)        (m, n) value changes of the perturbed solves
  !>   vel_base  real(rk)(:)           velocities of the iteration's solve [mm/s]
  !>   vel_face_base  real(rk)(3)      its face velocities [mm/s]
  !>   message   character(600)        error message
  subroutine bc_iterate(dt_sub, rejected, reason)

    real(rk), intent(in) :: dt_sub
    logical, intent(out) :: rejected
    character(len=120), intent(out) :: reason

    integer :: n, m, iter, k, c, status, modes(3)
    real(rk) :: q(4), q_pert(4), q_ideal(3), q_aim(4), ratio(4), signs(3), pert(3), dvel(3), &
      & tol(3), frac, elapsed, lambda
    real(rk), allocatable :: dq(:, :), vel_base(:)
    real(rk) :: vel_face_base(3)
    character(len=600) :: message

    ! Sizes, tolerances, and the box of the attempt's start coordinates (v2: x / L).
    n = size(bc_comps)
    m = n
    if (bc_path) m = n + 1
    allocate (dq(m, n))
    q = 0.0_rk
    q_ideal = 0.0_rk
    tol(1:n) = sched%bc_tol(bc_comps, i_step)
    modes(1:n) = sched%bc_mode(bc_comps, i_step)
    call box_bounds(coo, box_lo, box_len)

    ! The attempt's end in the step: the time since its start, and the fraction of a step of
    ! known duration.
    elapsed = time_solved + dt_sub - sched%time_start(i_step)
    frac = 0.0_rk
    if (sched%time_end(i_step) > sched%time_start(i_step)) frac = elapsed/ &
      & (sched%time_end(i_step) - sched%time_start(i_step))

    ! Perturbation sizes not yet set: min_pert_frac times the control face velocity (v2).
    do k = 1, n
      c = bc_comps(k)
      if (.not. bc_delta(c) > 0.0_rk) bc_delta(c) = config%min_pert_frac*abs(vel_face(axis))
    end do

    do iter = 1, config%max_bc_iter
      ! The solve at the current face velocities, and the values to reach.
      call trial_solve(dt_sub, iter, 0, rejected, reason, q(1:m))
      if (rejected) return
      lambda = 0.0_rk
      if (bc_path) lambda = (q(m) - bc_ctrl_start)/(sched%target_value(i_step) - bc_ctrl_start)
      q_ideal(1:n) = bc_mode_values(modes(1:n), bc_q_start(bc_comps), &
        & sched%bc_target(bc_comps, i_step), sched%bc_rate(bc_comps, i_step), lambda, frac, &
        & elapsed)
      call log_bc_iteration(iter, q(1:n), q_ideal(1:n))

      ! Within the tolerances: this solve is the attempt.
      if (bc_converged(q(1:n), q_ideal(1:n), tol(1:n))) return
      if (iter == config%max_bc_iter) exit

      ! One solve with each load-controlled face velocity perturbed towards its value, from the
      ! iteration's velocities.
      signs(1:n) = bc_pert_signs(q(1:n), q_ideal(1:n))
      vel_base = vel
      vel_face_base = vel_face
      do k = 1, n
        c = bc_comps(k)
        pert(k) = signs(k)*bc_delta(c)
        call shift_face(c, pert(k))
        call trial_solve(dt_sub, iter, 11*c, rejected, reason, q_pert(1:m))
        if (rejected) return
        dq(:, k) = q_pert(1:m) - q(1:m)
        vel = vel_base
        vel_face = vel_face_base
      end do

      ! The correction (path rows: the values back on the load path; the others: their values
      ! to reach), applied as homogeneous fields; the next perturbation sizes.
      ratio = 0.0_rk
      q_aim(1:n) = q_ideal(1:n)
      if (bc_path) then
        where (modes(1:n) == bc_mode_path)
          q_aim(1:n) = bc_q_start(bc_comps)
          ratio(1:n) = (sched%bc_target(bc_comps, i_step) - bc_q_start(bc_comps))/ &
            & (sched%target_value(i_step) - bc_ctrl_start)
        end where
        q_aim(m) = bc_ctrl_start
        ratio(m) = 1.0_rk
      end if
      call bc_correction(dq, q(1:m), q_aim(1:m), ratio(1:m), pert(1:n), bc_path, dvel(1:n), &
        & status)
      if (status /= 0) then
        call sim_close(sim)
        write (message, '(a,i0,a,i0,a,i0,a)') 'step ', i_step, ', increment ', incr, &
          & ', boundary-velocity iteration ', iter, ': the perturbed solves give a singular ' // &
          & 'system (the face forces do not respond to the face velocities); check the ' // &
          & 'load-controlled components and min_pert_frac'
        call fatal_error('fepx', trim(message))
      end if
      do k = 1, n
        call shift_face(bc_comps(k), dvel(k))
      end do
      bc_delta(bc_comps) = bc_pert_update(dvel(1:n), vel_face(axis), config%min_pert_frac)
    end do

    ! Not converged in max_bc_iter iterations.
    call sim_close(sim)
    write (message, '(a,i0,a,i0,a,i0,a)') 'step ', i_step, ', increment ', incr, &
      & ': the boundary-velocity iteration did not converge in max_bc_iter ', &
      & config%max_bc_iter, ' iterations ('
    do k = 1, n
      c = bc_comps(k)
      if (k == 1) then
        write (message, '(a,i0,a,es13.6,a,es13.6,a)') trim(message) // 'component ', 11*c, &
          & ': ', q(k), ', to reach ', q_ideal(k), ' ' // trim(bc_unit(c))
      else
        write (message, '(a,i0,a,es13.6,a,es13.6,a)') trim(message) // '; component ', 11*c, &
          & ': ', q(k), ', to reach ', q_ideal(k), ' ' // trim(bc_unit(c))
      end if
    end do
    message = trim(message) // '); raise max_bc_iter or the tolerances (load_tol, ' // &
      & 'stress_tol), or shorten the increments'
    call fatal_error('fepx', trim(message))

  end subroutine bc_iterate

  !> One solve of the boundary-velocity iteration (step 8.2b): from the start of the attempt
  !> (state_prev, the attempt's start coordinates) with the current velocities as the initial
  !> guess (the face velocities prescribed), its convergence record and log (report_increment),
  !> the sub-incrementation trigger, then the controlled values of the load-controlled
  !> components and, on the load path, of the control component, from the face forces of the
  !> end-of-increment stresses on the end-of-increment geometry x + dt v (Q-T10).
  !>
  !> Arguments:
  !>   dt_sub    (in)  real(rk)        length of the attempt [s]
  !>   iter      (in)  integer         boundary-velocity iteration
  !>   solve     (in)  integer         0, or the component (11, 22, 33) whose face velocity is
  !>                                   perturbed
  !>   rejected  (out) logical         the solve rejects the attempt
  !>   reason    (out) character(120)  why, for the messages
  !>   q         (out) real(rk)(:)     values of bc_comps, then the control's (path) [N] or [MPa]
  !> Locals:
  !>   coo_end   real(rk)(:, :)        (3, num_nodes_loc) end-of-increment coordinates [mm]
  !>   k, c      integer               index in bc_comps, component
  subroutine trial_solve(dt_sub, iter, solve, rejected, reason, q)

    real(rk), intent(in) :: dt_sub
    integer, intent(in) :: iter, solve
    logical, intent(out) :: rejected
    character(len=120), intent(out) :: reason
    real(rk), intent(out) :: q(:)

    real(rk), allocatable :: coo_end(:, :)
    integer :: k, c

    ! The solve from the start of the attempt.
    q = 0.0_rk
    state = state_prev
    call increment_solve(mesh, elt, phases, coo, bc_mask, dt_sub, opts, state_prev, state, vel, &
      & info, cons)
    rejected = attempt_rejected(reason)
    call report_increment(num_attempts, .not. rejected, iter, solve)
    if (rejected) return

    ! Face forces on the end-of-increment geometry, then the controlled values.
    coo_end = coo + dt_sub*reshape(vel, [3_ik, mesh%num_nodes_loc])
    stress_c = state%stress_s(:, elt%qpt_centroid, :)
    call face_forces(mesh, elt, coo_end, stress_c, loads, areas)
    do k = 1, size(bc_comps)
      c = bc_comps(k)
      q(k) = loads(c, face_set(c))
      if (sched%bc_kind(c, i_step) == target_by_stress) q(k) = q(k)/areas(face_set(c))
    end do
    if (bc_path) q(size(bc_comps) + 1) = control_value(i_step)

  end subroutine trial_solve

  !> Change the velocity of the face c1 of component c by dv (step 8.2b): the c-velocities of
  !> every node by dv (x_c - x_min,c)/L_c (v2's homogeneous field: the faces c0 keep theirs, c1
  !> moves by dv), the dofs of face c1 set to the new face velocity exactly.
  !>
  !> Arguments:
  !>   c   (in)  integer   component (1..3)
  !>   dv  (in)  real(rk)  velocity change [mm/s]
  !> Locals:
  !>   i_node  integer(ik)  local node
  subroutine shift_face(c, dv)

    integer, intent(in) :: c
    real(rk), intent(in) :: dv

    integer(ik) :: i_node

    vel_face(c) = vel_face(c) + dv
    do i_node = 1_ik, mesh%num_nodes_loc
      vel(3*(i_node - 1) + c) = vel(3*(i_node - 1) + c) + dv*(coo(c, i_node) - box_lo(c))/ &
        & box_len(c)
    end do
    where (face_dof == c) vel = vel_face(c)

  end subroutine shift_face

  !> Bounding box of coordinates over all ranks: lowest coordinates and lengths. Collective.
  !>
  !> Arguments:
  !>   coo_box  (in)  real(rk)(:, :)  (3, num_nodes_loc) coordinates [mm]
  !>   lo       (out) real(rk)(3)     lowest coordinate along x, y, z [mm]
  !>   len      (out) real(rk)(3)     lengths along x, y, z [mm]
  !> Locals:
  !>   lo_loc, hi_loc, hi  real(rk)(3)  rank-local extents, highest coordinates [mm]
  !>   j        integer               axis
  subroutine box_bounds(coo_box, lo, len)

    real(rk), intent(in) :: coo_box(:, :)
    real(rk), intent(out) :: lo(3), len(3)

    real(rk) :: lo_loc(3), hi_loc(3), hi(3)
    integer :: j

    lo_loc = huge(1.0_rk)
    hi_loc = -huge(1.0_rk)
    if (mesh%num_nodes_loc > 0_ik) then
      lo_loc = minval(coo_box(:, 1:mesh%num_nodes_loc), dim=2)
      hi_loc = maxval(coo_box(:, 1:mesh%num_nodes_loc), dim=2)
    end if
    do j = 1, 3
      call par_min(lo_loc(j), lo(j))
      call par_max(hi_loc(j), hi(j))
    end do
    len = hi - lo

  end subroutine box_bounds

  !> Log one boundary-velocity iteration (v2: velocities, loads, targets): the face velocities,
  !> then each load-controlled component's value and the value to reach.
  !>
  !> Arguments:
  !>   iter     (in)  integer       iteration
  !>   q        (in)  real(rk)(:)   values reached [N] or [MPa]
  !>   q_ideal  (in)  real(rk)(:)   values to reach [N] or [MPa]
  !> Locals:
  !>   line     character(400)      the message
  !>   k, c     integer             index in bc_comps, component
  subroutine log_bc_iteration(iter, q, q_ideal)

    integer, intent(in) :: iter
    real(rk), intent(in) :: q(:), q_ideal(:)

    character(len=400) :: line
    integer :: k, c

    ! (Six significant digits: the verification checks read the face velocities back.)
    line = 'Load iteration ' // log_int_text(iter) // ': face velocities ' // &
      & log_real_text(vel_face(1), 5) // ' ' // log_real_text(vel_face(2), 5) // ' ' // &
      & log_real_text(vel_face(3), 5) // ' mm/s;'
    do k = 1, size(bc_comps)
      c = bc_comps(k)
      line = trim(line) // ' ' // log_int_text(11*c) // ': ' // log_real_text(q(k), 6) // &
        & ' (to reach ' // log_real_text(q_ideal(k), 6) // ') ' // trim(bc_unit(c)) // ';'
    end do
    call log_info(trim(line), 2)

  end subroutine log_bc_iteration

  !> Unit of a load-controlled component's value in the current step: N or MPa.
  !>
  !> Arguments:
  !>   c     (in)  integer  component
  !>   unit  (out, result) character(3)
  function bc_unit(c) result(unit)

    integer, intent(in) :: c
    character(len=3) :: unit

    unit = 'N'
    if (sched%bc_kind(c, i_step) == target_by_stress) unit = 'MPa'

  end function bc_unit

  !> Start of a step for triaxial load control (step 8.2b, addendum §3.1, Q-T4): its
  !> load-controlled components (with the control component at a load or stress rate, step
  !> 8.2c), whether some follow the load path (a load or stress target of the control component
  !> at a strain rate), and their start values measured now (step 1: the initial values 0), with
  !> the control's value on the path.
  !>
  !> Arguments:
  !>   i_bc  (in)  integer  the step
  !> Locals:
  !>   all_comps  integer(3)  1, 2, 3
  !>   c          integer     component
  subroutine bc_step_start(i_bc)

    integer, intent(in) :: i_bc

    integer, parameter :: all_comps(3) = [1, 2, 3]
    integer :: c

    bc_comps = pack(all_comps, sched%bc_kind(:, i_bc) /= 0)
    bc_on = size(bc_comps) > 0
    bc_path = any(sched%bc_mode(bc_comps, i_bc) == bc_mode_path)
    bc_q_start = 0.0_rk
    bc_ctrl_start = 0.0_rk
    if (.not. bc_on .or. (i_bc == 1 .and. .not. restarted)) return

    ! The values at the step's start, on the current geometry.
    stress_c = state%stress_s(:, elt%qpt_centroid, :)
    call face_forces(mesh, elt, coo, stress_c, loads, areas)
    do c = 1, 3
      if (sched%bc_kind(c, i_bc) == 0) cycle
      bc_q_start(c) = loads(c, face_set(c))
      if (sched%bc_kind(c, i_bc) == target_by_stress) bc_q_start(c) = bc_q_start(c)/ &
        & areas(face_set(c))
    end do
    if (bc_path) bc_ctrl_start = control_value(i_bc)

  end subroutine bc_step_start

  !> The faces of the load-controlled components (step 8.2b): the face sets x1, y1, z1 of every
  !> component load-controlled in some step (an error when the mesh lacks one), and face_dof, the
  !> normal-velocity dofs of their nodes.
  !>
  !> Locals:
  !>   c, j, n   integer               component, face set, conditions set
  !>   specs_f   type(bc_spec_t)(:)    one condition per face, valued with its component (set
  !>                                   one by one: gfortran 13 loses an appended allocatable
  !>                                   character component, ROADMAP §G)
  !>   mask_f    logical(:)            (3 num_nodes_loc) dofs of the faces
  !>   values_f  real(rk)(:)           (3 num_nodes_loc) the component of each face dof
  subroutine bc_faces()

    integer :: c, j, n
    type(bc_spec_t), allocatable :: specs_f(:)
    logical, allocatable :: mask_f(:)
    real(rk), allocatable :: values_f(:)

    ! The face sets.
    allocate (specs_f(count([(any(sched%bc_kind(c, :) /= 0), c = 1, 3)])))
    n = 0
    do c = 1, 3
      if (.not. any(sched%bc_kind(c, :) /= 0)) cycle
      do j = 1, size(mesh%fasets)
        if (mesh%fasets(j)%label == axis_letter(c) // '1') face_set(c) = j
      end do
      if (face_set(c) == 0) call fatal_error('fepx', 'the mesh has no face set ' // &
        & axis_letter(c) // '1, whose force is the load of the load-controlled component ' // &
        & int_str(int(11*c, ik)) // ' (target_load, target_stress)')
      n = n + 1
      specs_f(n)%nset = axis_letter(c) // '1'
      specs_f(n)%comp = c
      specs_f(n)%value = real(c, rk)
    end do

    ! The normal-velocity dofs of their nodes.
    allocate (mask_f(3*mesh%num_nodes_loc), values_f(3*mesh%num_nodes_loc))
    call bc_resolve(mesh, specs_f, mask_f, values_f)
    allocate (face_dof(3*mesh%num_nodes_loc))
    face_dof = merge(nint(values_f), 0, mask_f)

  end subroutine bc_faces

  !> Setup of triaxial load control (steps 8.2b, 8.2c): the Voigt-averaged E and nu of the
  !> phases weighted by their volume fractions (Q-T12; v2 est_avg_mod by element counts), the
  !> original areas of the faces x1, y1, z1 of the load-controlled components (v2's area0), and
  !> each component's largest target magnitude over the run as a stress (loads over the
  !> original area of its face x1, y1 or z1; iteration (B)'s perturbation floor, v2-fixed C60;
  !> for every component with that face set, load-controlled in this run or not, so that a
  !> restart point carries the floor of the steps before it, step 8.4).
  !>
  !> Locals:
  !>   stiffs    real(rk)(:, :, :)  (6, 6, num_phases) the phases' stiffnesses [MPa]
  !>   vol_loc, vol  real(rk)(:)    (num_phases) phase volumes, rank-local and total [mm^3]
  !>   i_e       integer(ik)        local element
  !>   j, c      integer            phase, component
  !>   k         integer            face set x1, y1 or z1 of a component (0: none; one of zero
  !>                                  area is skipped)
  subroutine bc_setup()

    real(rk), allocatable :: stiffs(:, :, :), vol_loc(:), vol(:)
    integer(ik) :: i_e
    integer :: j, c, k

    ! Phase volume fractions and the Voigt moduli.
    allocate (stiffs(6, 6, config%num_phases), vol_loc(config%num_phases))
    allocate (vol(config%num_phases))
    vol_loc = 0.0_rk
    do i_e = 1_ik, mesh%num_elts
      j = int(mesh%elt_phase(i_e))
      vol_loc(j) = vol_loc(j) + element_volume(elt, coo(:, mesh%elt_nodes_loc(:, i_e)))
    end do
    call par_sum(vol_loc, vol)
    do j = 1, config%num_phases
      stiffs(:, :, j) = config%phases(j)%stiff_el_c
    end do
    call bc_elastic_moduli(stiffs, vol/sum(vol), e_avg, nu_avg)

    ! The original face areas, and the largest target magnitudes as stresses.
    stress_c = 0.0_rk
    call face_forces(mesh, elt, coo, stress_c, loads, areas)
    do c = 1, 3
      k = 0
      do j = 1, size(mesh%fasets)
        if (mesh%fasets(j)%label == axis_letter(c) // '1') k = j
      end do
      if (face_set(c) /= 0) areas0(c) = areas(face_set(c))
      ! (No such face set, or one without area: periodic meshes, no load control.)
      if (k == 0) cycle
      if (.not. areas(k) > 0.0_rk) cycle
      stress_max(c) = max(maxval(abs(config%target_load(c, :)))/areas(k), &
        & maxval(abs(config%target_stress(c, :))))
    end do
    if (bc_any) call log_info('Triaxial load control: elastic estimate E = ' // &
      & log_real_text(e_avg, 4) // ' MPa, nu = ' // log_real_text(nu_avg, 4), 2)

  end subroutine bc_setup

  !> Initial face velocities of the load-controlled components of step 1 at a strain rate (step
  !> 8.2b, addendum §3.4 "Initial velocities"; v2 driver_triaxcsr's initialization):
  !> bc_initial_velocities with the moduli of bc_setup: with a load or stress target of the
  !> control component the estimate of step 1's target stresses (loads over the original face
  !> areas) scaled to the control strain rate (v2; rate rows, step 8.2c, count 0 there);
  !> otherwise the control's stress rate from Hooke's law and the others' rates, of their
  !> targets over the step's duration or their given rates (step 8.2c).
  !>
  !> Locals:
  !>   stress    real(rk)(3)        target stresses [MPa] or their rates [MPa/s]
  !>   vel_est   real(rk)(3)        estimated face velocities [mm/s]
  !>   duration  real(rk)           duration of step 1 (time-based rows) [s]
  !>   c         integer            component
  subroutine bc_initial_estimate()

    real(rk) :: stress(3), vel_est(3), duration
    integer :: c

    ! Step 1's target stresses (original face areas), as rates over a known duration; rate
    ! rows as their stress rates.
    duration = sched%time_end(1) - sched%time_start(1)
    stress = 0.0_rk
    do c = 1, 3
      if (sched%bc_kind(c, 1) == 0) cycle
      if (sched%bc_mode(c, 1) == bc_mode_rate) then
        if (.not. sched%is_load(1)) stress(c) = sched%bc_rate(c, 1)
      else if (sched%is_load(1)) then
        stress(c) = sched%bc_target(c, 1)
      else
        stress(c) = sched%bc_target(c, 1)/duration
      end if
      if (sched%bc_kind(c, 1) /= target_by_stress) stress(c) = stress(c)/areas0(c)
    end do
    if (sched%is_load(1)) then
      stress(axis) = sched%target_value(1)
      if (sched%target_kind(1) /= target_by_stress) stress(axis) = stress(axis)/areas(i_face)
    end if

    ! The face velocities.
    vel_est = bc_initial_velocities(stress, axis, sched%strain_rate(axis, 1)*lengths0(axis), &
      & lengths0, e_avg, nu_avg, sched%is_load(1))
    where (sched%bc_kind(:, 1) /= 0) vel_face = vel_est
    call log_info('Triaxial load control: initial face velocities ' // &
      & log_real_text(vel_est(1), 4) // ' ' // log_real_text(vel_est(2), 4) // ' ' // &
      & log_real_text(vel_est(3), 4) // ' mm/s', 2)

  end subroutine bc_initial_estimate

  !> Start of a step at a load or stress rate of the control component (step 8.2c, addendum
  !> §3.4; v2 driver_triaxclr's curr_action and episode starts): its kind by the control
  !> value's magnitude (bc_step_action, Q-T9; the end value is the target, the start value plus
  !> the rate times a known duration, or not known for a strain target), and the episode its
  !> first attempt starts: the run (step 1: v2's estimate at increment 1), a reload (loading
  !> after a step that was not loading), an unload (unloading after one that was not), a dwell
  !> (a dwell after a step that was not a dwell).
  !>
  !> Arguments:
  !>   i_bc  (in)  integer  the step
  !> Locals:
  !>   q_end      real(rk)  control value at the step's end (when known) [N] or [MPa]
  !>   end_known  logical   it is known
  subroutine bc_rate_step_start(i_bc)

    integer, intent(in) :: i_bc

    real(rk) :: q_end
    logical :: end_known

    ! The end value: the target, or the rate over a known duration.
    end_known = .not. sched%is_open(i_bc)
    if (sched%bc_mode(axis, i_bc) == bc_mode_time) then
      q_end = sched%bc_target(axis, i_bc)
    else
      q_end = bc_q_start(axis) + sched%ctrl_rate(i_bc)*(sched%time_end(i_bc) - &
        & sched%time_start(i_bc))
    end if
    action = bc_step_action(sched%ctrl_rate(i_bc), bc_q_start(axis), q_end, end_known, &
      & sched%bc_tol(axis, i_bc))

    ! The episode starts (v2: reload not at step 1; the first increment of the run takes the
    ! estimate).
    ! (Step 1 of a restarted run follows the restart point's step, step 8.4.)
    start_run = i_bc == 1 .and. .not. restarted
    start_reload = action == bc_action_loading .and. prev_action /= bc_action_loading .and. &
      & (i_bc > 1 .or. restarted)
    start_unload = action == bc_action_unloading .and. prev_action /= bc_action_unloading
    start_dwell = action == bc_action_dwell .and. prev_action /= bc_action_dwell
    bc_at_start = .true.
    select case (action)
    case (bc_action_loading)
      text = 'loading'
    case (bc_action_unloading)
      text = 'unloading'
    case default
      text = 'dwell'
    end select
    if (start_reload) text = trim(text) // ', starting a reload episode'
    if (start_unload) text = trim(text) // ', starting an unload episode'
    if (start_dwell) text = trim(text) // ', starting a dwell episode'
    call log_info('Step ' // int_str(int(step_offset + i_bc, ik)) // ': ' // trim(text), 1)

  end subroutine bc_rate_step_start

  !> Boundary-velocity iteration (B) of an increment attempt over dt_sub at a load or stress
  !> rate of the control component (step 8.2c; addendum §3.4; v2 driver_triaxclr's
  !> bc_iteration_1 and bc_iteration_2, with v2-fixed C60): the values to reach of the
  !> load-controlled components (the control included) at the attempt's end time
  !> (bc_mode_values); at the step's first attempt v2's episode starts (the run's first
  !> increment and a reload: the elastic estimate of the face velocities for the increment's
  !> change, after a solve at zero velocity at a reload; an unload: the velocities over 10; a
  !> dwell: the velocities stored at the last dwell after the same kind of step); the scaling
  !> loop (every load-controlled face velocity scaled by bc_rate_scale until the control is
  !> within bc_rate_first_converged; not iterated at a dwell or reload start, in an unloading
  !> step, with a zero control velocity); the Newton loop (each load-controlled face velocity
  !> perturbed by +- delta, the correction of bc_correction from the finite-difference
  !> Jacobian, a solve, until every value is within its tolerance; the test also before the
  !> first perturbation, C60; delta from bc_rate_pert_initial, then bc_rate_pert_update, at
  !> least bc_rate_pert_floor). Every solve starts from the start of the attempt (trial_solve).
  !> A rejected solve ends the iteration with the attempt rejected; a loop not converged in
  !> max_bc_iter iterations, or a singular system, ends the run with an error. The iterations
  !> are numbered on across the two loops in the convergence file (the zero-velocity solve 0).
  !>
  !> Arguments:
  !>   dt_sub    (in)  real(rk)        length of the attempt [s]
  !>   rejected  (out) logical         a solve rejected the attempt
  !>   reason    (out) character(120)  why, for the messages
  !> Locals:
  !>   n, kc     integer               load-controlled components, the control's index in them
  !>   iter, iter2, k, c, status  integer  scaling and Newton iterations, index, component,
  !>                                   solve status of the correction
  !>   q, q_pert real(rk)(3)           values at the iteration's velocities and with a
  !>                                   perturbation [N] or [MPa]
  !>   q_init    real(rk)(3)           values at the attempt's start [N] or [MPa]
  !>   q_target  real(rk)(3)           values to reach at the attempt's end [N] or [MPa]
  !>   q_prev, vel_prev  real(rk)      control value and velocity of the solve before (scaling)
  !>   vel_init  real(rk)              control face velocity at the attempt's start (after the
  !>                                   episode start) [mm/s]
  !>   ds, vel_est  real(rk)(3)        stress changes [MPa] and face velocities [mm/s] of the
  !>                                   elastic estimate
  !>   signs, pert, dvel, tol  real(rk)(3)  perturbation signs, signed perturbations [mm/s],
  !>                                   velocity corrections [mm/s], tolerances [N] or [MPa]
  !>   delta, floor  real(rk)          perturbation size and its lower bound [mm/s]
  !>   frac, elapsed  real(rk)         fraction of the step [-] and time since its start [s] at
  !>                                   the attempt's end
  !>   scale     real(rk)              scale of the face velocities [-]
  !>   dq        real(rk)(3, 3)        value changes of the perturbed solves
  !>   zeros     real(rk)(3)           path directions (none) [-]
  !>   vel_base  real(rk)(:)           velocities of the iteration's solve [mm/s]
  !>   vel_face_base  real(rk)(3)      its face velocities [mm/s]
  subroutine bc_iterate_rate(dt_sub, rejected, reason)

    real(rk), intent(in) :: dt_sub
    logical, intent(out) :: rejected
    character(len=120), intent(out) :: reason

    integer :: n, kc, iter, iter2, k, c, status
    real(rk) :: q(3), q_pert(3), q_init(3), q_target(3), q_prev, vel_prev, vel_init, ds(3), &
      & vel_est(3), signs(3), pert(3), dvel(3), tol(3), delta, floor, frac, elapsed, scale, &
      & dq(3, 3), zeros(3), vel_face_base(3)
    real(rk), allocatable :: vel_base(:)

    ! Sizes, tolerances, the box of the attempt's start coordinates, and the values at the
    ! attempt's start (v2's initial_load).
    n = size(bc_comps)
    kc = findloc(bc_comps, axis, dim=1)
    tol(1:n) = sched%bc_tol(bc_comps, i_step)
    call box_bounds(coo, box_lo, box_len)
    call bc_values_now(q_init(1:n))
    q = q_init
    zeros = 0.0_rk

    ! The values to reach at the attempt's end.
    elapsed = time_solved + dt_sub - sched%time_start(i_step)
    frac = 0.0_rk
    if (sched%time_end(i_step) > sched%time_start(i_step)) frac = elapsed/ &
      & (sched%time_end(i_step) - sched%time_start(i_step))
    q_target(1:n) = bc_mode_values(sched%bc_mode(bc_comps, i_step), bc_q_start(bc_comps), &
      & sched%bc_target(bc_comps, i_step), sched%bc_rate(bc_comps, i_step), 0.0_rk, frac, &
      & elapsed)

    ! Episode starts (v2), at the step's first attempt: a solve at zero velocity at a reload,
    ! then the elastic estimate (run start, reload); an unload start backs off by 10; a dwell
    ! start takes the stored velocities.
    if (bc_at_start) then
      if (start_reload) then
        call log_info('Reload: solve at zero velocity of the load-controlled faces', 2)
        call set_faces(bc_comps, spread(0.0_rk, 1, n))
        call trial_solve(dt_sub, 0, 0, rejected, reason, q(1:n))
        if (rejected) return
      end if
      if (start_run .or. start_reload) then
        ds = 0.0_rk
        do k = 1, n
          c = bc_comps(k)
          ds(c) = q_target(k) - q(k)
          if (sched%bc_kind(c, i_step) /= target_by_stress) ds(c) = ds(c)/areas0(c)
        end do
        vel_est = bc_elastic_velocities(ds/dt_sub, box_len, e_avg, nu_avg)
        call set_faces(bc_comps, vel_est(bc_comps))
      end if
      if (start_unload) call scale_faces(bc_comps, 0.1_rk)
      if (start_dwell) then
        if (prev_action == bc_action_unloading) then
          call set_faces(bc_comps, vel_dwell_unload(bc_comps))
        else
          call set_faces(bc_comps, vel_dwell_load(bc_comps))
        end if
      end if
    end if
    vel_init = vel_face(axis)

    ! Scaling loop (v2 bc_iteration_1).
    q_prev = 0.0_rk
    vel_prev = 0.0_rk
    do iter = 1, config%max_bc_iter
      call trial_solve(dt_sub, iter, 0, rejected, reason, q(1:n))
      if (rejected) return
      call log_bc_iteration(iter, q(1:n), q_target(1:n))
      if (bc_rate_first_converged(q(kc), q_target(kc), q_init(kc), tol(kc))) exit
      if (bc_at_start .and. (start_dwell .or. start_reload)) exit
      if (action == bc_action_unloading .or. .not. abs(vel_face(axis)) > 0.0_rk) exit
      scale = bc_rate_scale(iter, action == bc_action_dwell, vel_face(axis), vel_prev, &
        & q_target(kc), q(kc), q_prev)
      q_prev = q(kc)
      vel_prev = vel_face(axis)
      call scale_faces(bc_comps, scale)
      if (iter == config%max_bc_iter) call bc_not_converged('scaling loop', q(1:n), &
        & q_target(1:n))
    end do

    ! Newton loop (v2 bc_iteration_2, v2-fixed C60: the lower bound of the perturbations and
    ! the test before the first one).
    floor = bc_rate_pert_floor(stress_max, box_len, e_avg, nu_avg, dt_sub, config%min_pert_frac)
    delta = bc_rate_pert_initial(vel_face(axis), vel_init, bc_at_start .and. start_dwell, floor)
    iter2 = 0
    do while (.not. bc_converged(q(1:n), q_target(1:n), tol(1:n)))
      iter2 = iter2 + 1
      ! One solve with each load-controlled face velocity perturbed towards its value.
      signs(1:n) = bc_pert_signs(q(1:n), q_target(1:n))
      vel_base = vel
      vel_face_base = vel_face
      do k = 1, n
        c = bc_comps(k)
        pert(k) = signs(k)*delta
        call shift_face(c, pert(k))
        call trial_solve(dt_sub, iter + iter2, 11*c, rejected, reason, q_pert(1:n))
        if (rejected) return
        dq(1:n, k) = q_pert(1:n) - q(1:n)
        vel = vel_base
        vel_face = vel_face_base
      end do

      ! The correction towards the values to reach (the finite-difference Jacobian), applied
      ! as homogeneous fields, then the solve at the corrected velocities.
      call bc_correction(dq(1:n, 1:n), q(1:n), q_target(1:n), zeros(1:n), pert(1:n), .false., &
        & dvel(1:n), status)
      if (status /= 0) then
        call sim_close(sim)
        write (text, '(a,i0,a,i0,a)') 'step ', i_step, ', increment ', incr, &
          & ': the perturbed solves of the boundary-velocity iteration give a singular ' // &
          & 'system (the face forces do not respond to the face velocities); check the ' // &
          & 'load-controlled components and min_pert_frac'
        call fatal_error('fepx', trim(text))
      end if
      do k = 1, n
        call shift_face(bc_comps(k), dvel(k))
      end do
      call trial_solve(dt_sub, iter + iter2, 0, rejected, reason, q(1:n))
      if (rejected) return
      call log_bc_iteration(iter + iter2, q(1:n), q_target(1:n))
      if (bc_converged(q(1:n), q_target(1:n), tol(1:n))) exit
      if (iter2 == config%max_bc_iter) call bc_not_converged('Newton loop', q(1:n), &
        & q_target(1:n))
      delta = bc_rate_pert_update(dvel(1:n), floor)
    end do

  end subroutine bc_iterate_rate

  !> The run stops: a loop of iteration (B) did not converge in max_bc_iter iterations (step
  !> 8.2c); the message names the step, the increment, the loop, the values and the values to
  !> reach.
  !>
  !> Arguments:
  !>   loop     (in)  character    the loop ("scaling loop", "Newton loop")
  !>   q        (in)  real(rk)(:)  values reached [N] or [MPa]
  !>   q_ideal  (in)  real(rk)(:)  values to reach [N] or [MPa]
  !> Locals:
  !>   message  character(600)     the message
  !>   k, c     integer            index in bc_comps, component
  subroutine bc_not_converged(loop, q, q_ideal)

    character(len=*), intent(in) :: loop
    real(rk), intent(in) :: q(:), q_ideal(:)

    character(len=600) :: message
    integer :: k, c

    call sim_close(sim)
    write (message, '(a,i0,a,i0,a,i0,a)') 'step ', i_step, ', increment ', incr, &
      & ': the boundary-velocity iteration (' // loop // ') did not converge in max_bc_iter ', &
      & config%max_bc_iter, ' iterations ('
    do k = 1, size(bc_comps)
      c = bc_comps(k)
      if (k > 1) message = trim(message) // '; '
      write (message, '(a,i0,a,es13.6,a,es13.6,a)') trim(message) // 'component ', 11*c, &
        & ': ', q(k), ', to reach ', q_ideal(k), ' ' // trim(bc_unit(c))
    end do
    message = trim(message) // '); raise max_bc_iter or the tolerances (load_tol, ' // &
      & 'stress_tol), or shorten the increments'
    call fatal_error('fepx', trim(message))

  end subroutine bc_not_converged

  !> The values of the load-controlled components now: the face forces of the current state on
  !> the current coordinates (collective), as loads [N] or stresses [MPa] (step 8.2c).
  !>
  !> Arguments:
  !>   q  (out) real(rk)(:)  values of bc_comps [N] or [MPa]
  !> Locals:
  !>   k, c  integer         index in bc_comps, component
  subroutine bc_values_now(q)

    real(rk), intent(out) :: q(:)

    integer :: k, c

    stress_c = state%stress_s(:, elt%qpt_centroid, :)
    call face_forces(mesh, elt, coo, stress_c, loads, areas)
    do k = 1, size(bc_comps)
      c = bc_comps(k)
      q(k) = loads(c, face_set(c))
      if (sched%bc_kind(c, i_step) == target_by_stress) q(k) = q(k)/areas(face_set(c))
    end do

  end subroutine bc_values_now

  !> Set the velocities of the faces c1 of components (step 8.2c): v2's homogeneous field
  !> V_c (x_c - x_min,c)/L_c on every node's c-velocity, the dofs of face c1 exactly V_c.
  !>
  !> Arguments:
  !>   comps   (in)  integer(:)   components (1..3)
  !>   values  (in)  real(rk)(:)  their face velocities [mm/s]
  !> Locals:
  !>   k, c    integer            index, component
  !>   i_node  integer(ik)        local node
  subroutine set_faces(comps, values)

    integer, intent(in) :: comps(:)
    real(rk), intent(in) :: values(:)

    integer :: k, c
    integer(ik) :: i_node

    do k = 1, size(comps)
      c = comps(k)
      vel_face(c) = values(k)
      do i_node = 1_ik, mesh%num_nodes_loc
        vel(3*(i_node - 1) + c) = values(k)*(coo(c, i_node) - box_lo(c))/box_len(c)
      end do
      where (face_dof == c) vel = vel_face(c)
    end do

  end subroutine set_faces

  !> Scale the velocities of components (step 8.2c; v2 bc_iteration_1 and the unload start):
  !> every node's c-velocity and the face velocity V_c by a factor.
  !>
  !> Arguments:
  !>   comps   (in)  integer(:)  components (1..3)
  !>   factor  (in)  real(rk)    the factor [-]
  !> Locals:
  !>   k, c    integer           index, component
  !>   i_node  integer(ik)       local node
  subroutine scale_faces(comps, factor)

    integer, intent(in) :: comps(:)
    real(rk), intent(in) :: factor

    integer :: k, c
    integer(ik) :: i_node

    do k = 1, size(comps)
      c = comps(k)
      vel_face(c) = factor*vel_face(c)
      do i_node = 1_ik, mesh%num_nodes_loc
        vel(3*(i_node - 1) + c) = factor*vel(3*(i_node - 1) + c)
      end do
      where (face_dof == c) vel = vel_face(c)
    end do

  end subroutine scale_faces

  !> The control value of a step at a load or stress rate now (step 8.2c, for its plan): the
  !> control-face force [N], or over the face's current area [MPa] (collective).
  !>
  !> Arguments:
  !>   i_val  (in)  integer  the step
  !>   v      (out, result) real(rk)  the value [N] or [MPa]
  function rate_control_value(i_val) result(v)

    integer, intent(in) :: i_val
    real(rk) :: v

    v = control_force()
    if (sched%rate_kind(i_val) == rate_by_stress) v = v/areas(i_face)

  end function rate_control_value

  !> The control-face force now: the face forces of the current state on the current
  !> coordinates (collective), component axis of face i_face [N] (step 8.2c).
  !>
  !> Arguments:
  !>   f  (out, result) real(rk)  the force [N]
  function control_force() result(f)

    real(rk) :: f

    stress_c = state%stress_s(:, elt%qpt_centroid, :)
    call face_forces(mesh, elt, coo, stress_c, loads, areas)
    f = loads(axis, i_face)

  end function control_force

  !> Sign of a step's control rate (step 8.2c, addendum §3.5, the cyclic-hardening reset): its
  !> load or stress rate, or its strain rate (0 for a dwell or a zero strain rate).
  !>
  !> Arguments:
  !>   rec  (in)  type(step_record_t)  the step
  !>   r    (out, result) real(rk)     +1, -1 or 0
  function control_rate_sign(rec) result(r)

    type(step_record_t), intent(in) :: rec
    real(rk) :: r

    if (rec%is_rate) then
      r = rec%ctrl_rate
    else
      r = rec%strain_rate(axis)
    end if
    if (abs(r) > 0.0_rk) then
      r = sign(1.0_rk, r)
    else
      r = 0.0_rk
    end if

  end function control_rate_sign

  !> Log the end of a step at a load or stress rate (step 8.2c): the control value reached, the
  !> increments and the time.
  !>
  !> Arguments:
  !>   i_end  (in)  integer  the step
  !> Locals:
  !>   v  real(rk)  the control value [N] or [MPa]
  subroutine log_rate_end(i_end)

    integer, intent(in) :: i_end

    real(rk) :: v

    v = rate_control_value(i_end)
    call log_info('Step ' // log_int_text(step_offset + i_end) // ': control ' // &
      & log_int_text(11*axis) // ' at ' // log_real_text(v, 6) // ' ' // trim(bc_unit(axis)) // &
      & ' after ' // log_count_text(sched%num_incrs(i_end), 'increment') // ', t = ' // &
      & log_real_text(time, 4) // ' s', 1)

  end subroutine log_rate_end

  !> Solve one increment (or landing attempt) of length dt_incr from the current coordinates,
  !> state and velocities (solve_interval: with sub-incrementation, subdivided where
  !> rejected), and log its sub-increments; an increment solved without subdivision resets
  !> the count of successive subdivided increments (limit L2, D63).
  !>
  !> Arguments:
  !>   dt_incr  (in)  real(rk)  length of the increment [s]
  subroutine solve_increment(dt_incr)

    real(rk), intent(in) :: dt_incr

    ! The increment, subdivided where rejected.
    num_subincrs = 0
    dt_smallest = dt_incr
    call solve_interval(dt_incr, 0)

    ! Its sub-increments (sub-incrementation, D63).
    if (config%subincrementation) then
      if (num_subincrs > 1) then
        call log_info(log_count_text(num_subincrs, 'sub-increment') // ' (' // &
          & log_count_text(num_attempts, 'attempt') // '), smallest ' // &
          & log_real_text(dt_smallest, 4) // ' s', 2)
      else
        num_successive = 0
      end if
    end if

  end subroutine solve_increment

  !> Landing on a load or stress target (step 8.2a; addendum_triaxial_load_control.md §3.6,
  !> decisions Q-T13, Q-T14): after an increment of a load- or stress-target step solved from
  !> the start state kept in coo_land, vel_land, state_land (value value_prev there), an
  !> increment whose value passes the target by more than the tolerance delta is discarded (the
  !> start restored, as a rejected sub-incrementation attempt) and redone with the length
  !> given by false position on the bracket (dt_lo, v_lo) = (0, value_prev), (dt_hi, v_hi) =
  !> (dt, value) (step_load_false_position; Press et al. 2007, Numerical Recipes 3rd ed.,
  !> §9.2): an attempt within delta of the target is accepted, one past it by more replaces the
  !> passed end of the bracket, one short of it by more the short end; at most max_bc_iter
  !> attempts, then the start state (the last accepted increment) is written as the current step
  !> and the run stops with an error. The strain rate is unchanged; dtime_min and dtime_factor
  !> do not apply to the attempts. On return the state is that of the accepted attempt and
  !> dt_incr its length.
  !>
  !> Arguments:
  !>   i_land   (in)    integer   the step (a load or stress target)
  !>   dt_incr  (inout) real(rk)  length of the increment solved; on return, of the accepted
  !>                              attempt [s]
  !> Locals:
  !>   v_now          real(rk)   value of the last attempt [N] or [MPa]
  !>   dt_lo, v_lo    real(rk)   short end of the bracket: length [s] and value
  !>   dt_hi, v_hi    real(rk)   passed end of the bracket: length [s] and value
  !>   dt_try         real(rk)   length of an attempt [s]
  !>   i_attempt      integer    landing attempt
  !>   unit           character(4)  unit of the values (N or MPa), for messages
  !>   message        character(400)  the error message
  subroutine land_increment(i_land, dt_incr)

    integer, intent(in) :: i_land
    real(rk), intent(inout) :: dt_incr

    real(rk) :: v_now, dt_lo, v_lo, dt_hi, v_hi, dt_try
    integer :: i_attempt
    character(len=4) :: unit
    character(len=400) :: message

    ! The value reached; nothing to do unless the target is passed by more than delta.
    v_now = attempt_value(i_land)
    if (step_load_landing(v_now, sched%target_value(i_land), sched%load_sign(i_land), &
      & sched%target_tol(i_land)) /= landing_passed) return
    unit = merge('MPa ', 'N   ', sched%target_kind(i_land) == target_by_stress)

    ! The bracket: the start of the increment (short of the target by more than delta: the step
    ! was not complete) and the increment just solved.
    dt_lo = 0.0_rk
    v_lo = value_prev
    dt_hi = dt_incr
    v_hi = v_now
    do i_attempt = 1, config%max_bc_iter
      ! The next attempt by false position, from the start of the increment.
      dt_try = step_load_false_position(dt_lo, v_lo, dt_hi, v_hi, sched%target_value(i_land))
      call log_info('Value ' // log_real_text(v_hi, 6) // ' ' // trim(unit) // ' after ' // &
        & log_real_text(dt_hi, 4) // ' s passes the target by more than the tolerance; ' // &
        & 'landing attempt ' // log_int_text(i_attempt) // ' over ' // &
        & log_real_text(dt_try, 4) // ' s', 2)
      coo = coo_land
      vel = vel_land
      vel_face = vel_face_land
      state = state_land
      time_solved = time
      call solve_increment(dt_try)
      v_now = attempt_value(i_land)

      ! Within delta: accepted; otherwise the bracket's end on the same side is replaced.
      select case (step_load_landing(v_now, sched%target_value(i_land), &
        & sched%load_sign(i_land), sched%target_tol(i_land)))
      case (landing_within)
        dt_incr = dt_try
        call log_info('Landed at ' // log_real_text(v_now, 6) // ' ' // trim(unit) // ' in ' // &
          & log_count_text(i_attempt, 'attempt') // ' (dt = ' // log_real_text(dt_try, 4) // &
          & ' s)', 2)
        return
      case (landing_passed)
        dt_hi = dt_try
        v_hi = v_now
      case default
        dt_lo = dt_try
        v_lo = v_now
      end select
    end do

    ! Not landed: the last accepted increment (the start state) as the current step, then stop.
    coo = coo_land
    vel = vel_land
    vel_face = vel_face_land
    state = state_land
    call write_step(i_land)
    call sim_close(sim)
    write (message, '(a,i0,a,i0,a,es13.6,a,es10.3,a,i0,a,es12.5,a,es13.6,a,es12.5,a,es13.6,a)') &
      & 'step ', step_offset + i_land, ', increment ', incr, ': no landing on the target ', &
      & sched%target_value(i_land), ' ' // trim(unit) // ' within the tolerance ', &
      & sched%target_tol(i_land), ' in max_bc_iter ', config%max_bc_iter, &
      & ' attempts (bracket: ', dt_lo, ' s -> ', v_lo, ', ', dt_hi, ' s -> ', v_hi, &
      & '); the last accepted increment is written; raise max_bc_iter or the tolerance'
    call fatal_error('fepx', trim(message))

  end subroutine land_increment

  !> The controlled value of a load or stress target after an attempt: the face forces of the
  !> current state on the current coordinates (collective), then control_value.
  !>
  !> Arguments:
  !>   i_val  (in)  integer  the step
  !>   v      (out, result) real(rk)  the value [N] or [MPa]
  function attempt_value(i_val) result(v)

    integer, intent(in) :: i_val
    real(rk) :: v

    stress_c = state%stress_s(:, elt%qpt_centroid, :)
    call face_forces(mesh, elt, coo, stress_c, loads, areas)
    v = control_value(i_val)

  end function attempt_value

  !> Report the global iteration of an increment attempt (or of one solve of its
  !> boundary-velocity iteration): its convergence record (when printed), then the log; a
  !> failure ends the run. An attempt rejected by sub-incrementation (Q35, D63) is only recorded
  !> (solve_interval logs it).
  !>
  !> Arguments:
  !>   attempt   (in)  integer  attempt within the increment (1 without sub-incrementation)
  !>   accepted  (in)  logical  the attempt was accepted
  !>   bc_iter   (in, optional) integer  boundary-velocity iteration (step 8.2b)
  !>   bc_solve  (in, optional) integer  0, or the component whose face velocity is perturbed
  subroutine report_increment(attempt, accepted, bc_iter, bc_solve)

    integer, intent(in) :: attempt
    logical, intent(in) :: accepted
    integer, intent(in), optional :: bc_iter, bc_solve

    ! The run's totals (after the initialization; the closing summary, step 9.7).
    if (in_run) then
      num_iters_run = num_iters_run + info%num_iters
      num_nr_run = num_nr_run + info%num_nr_iters
      num_cg_run = num_cg_run + info%cg_iters
      time_elements_run = time_elements_run + info%time_elements
      time_cg_run = time_cg_run + info%time_cg
    end if

    ! Convergence record of the increment's global iterations (incr 0: the initialization).
    call sim_write_convergence(sim, incr, info%iter_method(1:info%num_iters), &
      & info%iter_norms(:, 1:info%num_iters), info%iter_cg(1:info%num_iters), attempt, accepted, &
      & bc_iter, bc_solve)
    if (.not. accepted) return

    select case (info%status)
    case (increment_converged)
      if (info%max_point_failed > 0_ik) then
        write (text, '(a,i0,a,i0,a)') 'Material point solves did not converge (', &
          & info%num_point_failed, ' at the converged velocities, at most ', &
          & info%max_point_failed, ' in an iteration); they keep their last iterate.'
        call log_warning(trim(text), 2)
      end if
      ! (Q60 (b1): the solve's wall time; its split goes to the closing summary.)
      call log_info('Converged in ' // log_count_text(info%num_iters, 'iteration') // ' (' // &
        & seconds_text(info%time_total) // ' s)', 2)
    case (increment_max_iters)
      call fatal_error('fepx', 'the global iteration did not converge in global_max_iters ' // &
        & 'iterations')
    case (increment_diverging)
      call fatal_error('fepx', 'the global iteration is diverging')
    case default
      call fatal_error('fepx', 'the linear solver (CG) failed: it reached cg_max_iters, or ' // &
        & 'a nodal block of the block-Jacobi preconditioner is not positive definite')
    end select

  end subroutine report_increment

  !> Log the plan of a step: times and increments, then, one deeper and one per line, the strain
  !> of the control component and the other prescribed components (triaxial); for a load or
  !> stress target, the target and the bounds of its increments; at a load or stress rate (step 8.2c) the rate, then the target, a dwell,
  !> or the strain target with the bounds of its increments. The step's line at the given depth,
  !> the other components one deeper (step 9.7: depth 2 in the setup summary, 1 under
  !> "Running step" for a step planned during the run).
  !>
  !> Arguments:
  !>   i_plan  (in)  integer  the step
  !>   depth   (in)  integer  depth of the step's line (log_mod)
  !> Locals:
  !>   j_comp  integer        a normal component
  !>   head    character(:), allocatable  "Step <n>"
  !>   line    character(:), allocatable  the message
  !>   rate_row, unit  character(16), character(6)  the control's rate row and its unit
  subroutine log_step_plan(i_plan, depth)

    integer, intent(in) :: i_plan, depth

    integer :: j_comp
    character(len=:), allocatable :: head, line
    character(len=16) :: rate_row
    character(len=6) :: unit

    head = 'Step ' // log_int_text(step_offset + i_plan)
    if (sched%is_rate(i_plan)) then
      ! A load or stress rate of the control component (step 8.2c).
      rate_row = 'load_rate'
      unit = 'N/s'
      if (sched%rate_kind(i_plan) == rate_by_stress) then
        rate_row = 'stress_rate'
        unit = 'MPa/s'
      end if
      line = head // ': ' // trim(rate_row) // ' ' // log_int_text(11*axis) // ' ' // &
        & log_real_text(sched%ctrl_rate(i_plan), 4) // ' ' // trim(unit)
      if (sched%is_open(i_plan)) then
        line = line // ' to target_strain ' // log_int_text(11*axis) // ' ' // &
          & log_real_text(sched%target_value(i_plan), 4) // ' from t = ' // &
          & log_real_text(sched%time_start(i_plan), 4) // ' s, increments ' // &
          & log_real_text(sched%dt_min(i_plan), 4) // ' to ' // &
          & log_real_text(sched%dt(i_plan), 4) // ' s'
      else if (sched%is_dwell(i_plan)) then
        line = line // ' (dwell), t = ' // log_real_text(sched%time_start(i_plan), 4) // &
          & ' -> ' // log_real_text(sched%time_end(i_plan), 4) // ' s, increments of at ' // &
          & 'most ' // log_real_text(sched%dt(i_plan), 4) // ' s'
      else
        line = line // ', t = ' // log_real_text(sched%time_start(i_plan), 4) // ' -> ' // &
          & log_real_text(sched%time_end(i_plan), 4) // ' s, ' // &
          & log_count_text(sched%num_incrs(i_plan), 'increment')
        if (sched%bc_mode(axis, i_plan) == bc_mode_time) line = line // ' to ' // &
          & log_real_text(sched%bc_target(axis, i_plan), 4)
      end if
    else if (sched%is_load(i_plan)) then
      if (sched%target_kind(i_plan) == target_by_stress) then
        line = head // ': target_stress ' // log_int_text(11*axis) // ' ' // &
          & log_real_text(sched%target_value(i_plan), 4) // ' MPa'
      else
        line = head // ': target_load ' // log_int_text(11*axis) // ' ' // &
          & log_real_text(sched%target_value(i_plan), 4) // ' N'
      end if
      line = line // ' from t = ' // log_real_text(sched%time_start(i_plan), 4) // &
        & ' s, increments ' // log_real_text(sched%dt_min(i_plan), 4) // ' to ' // &
        & log_real_text(sched%dt(i_plan), 4) // ' s'
    else
      line = head // ': t = ' // log_real_text(sched%time_start(i_plan), 4) // ' -> ' // &
        & log_real_text(sched%time_end(i_plan), 4) // ' s, ' // &
        & log_count_text(sched%num_incrs(i_plan), 'increment')
    end if
    call log_info(line, depth)
    ! A strain step: the control component's strain on its own line, one deeper, as the other
    ! components' lines below.
    if (.not. sched%is_rate(i_plan) .and. .not. sched%is_load(i_plan) .and. axis > 0) &
      & call log_info(strain_line(trim(sched%comp_label(axis)), axis, i_plan), depth + 1)
    do j_comp = 1, 3
      if (j_comp == axis .or. .not. sched%prescribed(j_comp)) cycle
      if (sched%bc_mode(j_comp, i_plan) == bc_mode_rate) then
        ! A rate row (step 8.2c).
        if (sched%bc_kind(j_comp, i_plan) == target_by_stress) then
          line = 'Component ' // log_int_text(11*j_comp) // ': stress_rate ' // &
            & log_real_text(sched%bc_rate(j_comp, i_plan), 4) // ' MPa/s (load-controlled)'
        else
          line = 'Component ' // log_int_text(11*j_comp) // ': load_rate ' // &
            & log_real_text(sched%bc_rate(j_comp, i_plan), 4) // ' N/s (load-controlled)'
        end if
      else if (sched%bc_kind(j_comp, i_plan) == target_by_stress) then
        line = 'Component ' // log_int_text(11*j_comp) // ': target_stress ' // &
          & log_real_text(sched%bc_target(j_comp, i_plan), 4) // ' MPa (load-controlled)'
      else if (sched%bc_kind(j_comp, i_plan) /= 0) then
        line = 'Component ' // log_int_text(11*j_comp) // ': target_load ' // &
          & log_real_text(sched%bc_target(j_comp, i_plan), 4) // ' N (load-controlled)'
      else if (sched%is_load(i_plan) .or. sched%is_open(i_plan)) then
        line = 'Strain ' // log_int_text(11*j_comp) // ': rate ' // &
          & log_real_text(sched%strain_rate(j_comp, i_plan), 4) // ' 1/s'
      else
        line = strain_line(log_int_text(11*j_comp), j_comp, i_plan)
      end if
      call log_info(line, depth + 1)
    end do

  end subroutine log_step_plan

  !> The plan line of a strain-controlled component of a step:
  !> "Strain <ij>: <start> -> <end> at <rate>/s" (log_step_plan; at most 71 columns at depth 3,
  !> all three values negative included, so that it never wraps at log_width).
  !>
  !> Arguments:
  !>   label   (in)  character  the component as text, e.g. "22"
  !>   i_comp  (in)  integer    the normal component (1 to 3)
  !>   i_plan  (in)  integer    the step
  !> Result:
  !>   text    character(:), allocatable  the line
  function strain_line(label, i_comp, i_plan) result(text)

    character(len=*), intent(in) :: label
    integer, intent(in) :: i_comp, i_plan
    character(len=:), allocatable :: text

    ! Start and end strains of the step and the strain rate.
    text = 'Strain ' // label // ': ' // &
      & log_real_text(sched%strain_start(i_comp, i_plan), 4) // ' -> ' // &
      & log_real_text(sched%strain_end(i_comp, i_plan), 4) // ' at ' // &
      & log_real_text(sched%strain_rate(i_comp, i_plan), 4) // '/s'

  end function strain_line

  !> Log the end of a load or stress target step: the target, the value reached, the
  !> difference in the step's direction (past the target when positive; within the tolerance,
  !> landing rule) and the increments.
  !>
  !> Arguments:
  !>   i_end  (in)  integer   the step
  !>   reached (in) real(rk)  the value reached [N] or [MPa]
  subroutine log_load_end(i_end, reached)

    integer, intent(in) :: i_end
    real(rk), intent(in) :: reached

    call log_info('Step ' // log_int_text(step_offset + i_end) // ': target ' // &
      & log_real_text(sched%target_value(i_end), 6) // ', reached ' // &
      & log_real_text(reached, 6) // ' (past it by ' // &
      & log_real_text((reached - sched%target_value(i_end))*sched%load_sign(i_end), 3) // &
      & ') in ' // log_count_text(sched%num_incrs(i_end), 'increment') // ', t = ' // &
      & log_real_text(time, 4) // ' s', 1)

  end subroutine log_load_end

  !> The controlled value of a load or stress target step from the current face forces: the
  !> control component of the force on the control face [N], or that force over the face's
  !> current area [MPa] (addendum §3.1).
  !>
  !> Arguments:
  !>   i_val  (in)  integer  the step
  !>   v      (out, result) real(rk)  the value [N] or [MPa]
  function control_value(i_val) result(v)

    integer, intent(in) :: i_val
    real(rk) :: v

    v = loads(axis, i_face)
    if (sched%target_kind(i_val) == target_by_stress) v = v/areas(i_face)

  end function control_value

  !> Necking test after an increment (check_necking on; addendum §3.6, v2-fixed C54): the
  !> force on the control face changed against the step's direction (load and stress targets:
  !> the target's sign; strain and time targets: the sign of the control strain rate; not at a
  !> load or stress rate of the control component, step 8.2c, addendum §3.5). The
  !> results of this increment are written as the current step (unless just written), then the
  !> run stops with an error naming the step, the increment and the forces.
  !>
  !> Arguments:
  !>   i_neck   (in)  integer  current step
  !>   written  (in)  logical  the results of this increment were just written as the step's
  !> Locals:
  !>   force    real(rk)       control-face force after the increment [N]
  !>   dir      real(rk)       direction of the step, +1, -1 (0: no direction, no test)
  subroutine check_necking(i_neck, written)

    integer, intent(in) :: i_neck
    logical, intent(in) :: written

    real(rk) :: force, dir

    force = loads(axis, i_face)
    if (sched%is_rate(i_neck)) then
      ! Not under a load or stress rate (addendum §3.5: a load that cannot be carried shows as
      ! an iteration failure).
      dir = 0.0_rk
    else if (sched%is_load(i_neck)) then
      dir = sched%load_sign(i_neck)
    else if (abs(sched%strain_rate(axis, i_neck)) > 0.0_rk) then
      dir = sign(1.0_rk, sched%strain_rate(axis, i_neck))
    else
      dir = 0.0_rk
    end if
    if ((force - force_prev)*dir < 0.0_rk) then
      if (.not. written) call write_step(i_neck)
      call sim_close(sim)
      write (text, '(a,i0,a,i0,a,es11.4,a,a,a,es13.6,a,es13.6,a)') 'stopped in step ', &
        & step_offset + i_neck, ' at increment ', incr, ' (t = ', time, ' s): necking (check_necking on): ' &
        & // 'the load on face ', axis_letter(axis) // '1', ' fell from ', force_prev, &
        & ' N to ', force, ' N'
      call fatal_error('fepx', trim(text))
    end if
    force_prev = force

  end subroutine check_necking

  !> Write the results of a step: current coordinates, displacement from the original
  !> coordinates, velocities, and the centroid-point element results of the state (stress,
  !> elastic strain, orientation; step 6.5: slip-system strengths, slip, slip rates, resolved
  !> shear stresses, plastic strain and its equivalent, plastic deformation rate, plastic spin
  !> rate, lattice rotation rate, work; step 6.5b: total strain, the equivalent values, the
  !> deformation rate and velocity gradient, the spin and slip parts of the rotation rate,
  !> plastic work, work rates, the sample-frame plastic strain; step 7.2: the Kocks-Mecking
  !> density km_rho [1/mm^2] and recovery coefficient km_k2 [-], rate km_rhorate [1/(mm^2 s)], zero for other phases).
  !>
  !> Output forms (v2 write_res; manual, output page): symmetric tensors as tensor components
  !> 11 22 33 23 13 12 in the sample frame (sim_write_elts); strain_pl accumulated in the crystal
  !> frame and written in the sample frame of the current orientation (v2); spinrate as the
  !> components 12 13 23 of the plastic spin tensor sum gammadot q (v2's order); rotrate as the
  !> axial vector omega (W v = omega x v) of the lattice rotation rate, i.e. (-w23, w13, -w12)
  !> (v2's rotation vector of d_rstar over dtime); rotrate_spin and rotrate_slip in the same
  !> form, of the w^p (the spin of the stress solve) and the slip spin sum gammadot q used in
  !> the orientation update of the increment (v2's d_rstar_spin, d_rstar_slip: exp(wp_hat dt)
  !> and exp(wp_ss dt) in the crystal frame, whose rotation vectors over dt are the spins,
  !> written in the sample frame; rotrate = rotrate_spin - rotrate_slip); velgrad as the nine
  !> components L11 L12 L13 L21 ... L33, L_ij = dv_i/dx_j (v2 write_res_tensor); the equivalent
  !> values sqrt(2/3 e':e') of strains and rates, sqrt(3/2 s':s') of the stress (v2
  !> strain_equiv_3x3, stress_equiv_3x3; defrate_eq of the deviatoric d, as v2's defrate_eq of
  !> results%d); slip-system values of each element's phase.
  !>
  !> Arguments:
  !>   step_in  (in)  integer  step number of the configuration (0: initial state)
  !> Locals:
  !>   step   integer        the step number written (step_offset + step_in after a restart
  !>                         point)
  !>   q      integer(ik)    centroid point
  !>   i_elt  integer(ik)    element
  !>   velgrad  real(rk)(3, 3)  velocity gradient d + w of an element [1/s]
  subroutine write_step(step_in)

    integer, intent(in) :: step_in

    integer(ik) :: q, i_elt
    real(rk) :: velgrad(3, 3)
    integer :: step

    ! The step's number in the output (after a restart point: step_offset + step_in, step 8.4).
    step = step_in
    if (step_in > 0) step = step_offset + step_in

    ! Centroid-point values of the current state.
    q = elt%qpt_centroid
    stress_c = state%stress_s(:, q, :)
    strain_el_c = state%strain_el_s(:, q, :)
    rot_c = state%rot_c2s(:, :, q, :)

    ! Node and element variables (each written only when printed).
    call sim_write_nodes(sim, mesh, 'coo', coo, step)
    call sim_write_nodes(sim, mesh, 'disp', coo - mesh%coo_s, step)
    call sim_write_nodes(sim, mesh, 'vel', reshape(vel, [3_ik, mesh%num_nodes_loc]), step)
    call sim_write_elts(sim, mesh, 'stress', stress_c, step)
    call sim_write_elts(sim, mesh, 'strain_el', strain_el_c, step)
    call sim_write_ori(sim, mesh, rot_c, step)

    ! Slip-system values (each element's phase's systems).
    call sim_write_elts_slip(sim, mesh, 'crss', state%crss(:, q, :), num_slip_elt, step)
    call sim_write_elts_slip(sim, mesh, 'slip', state%slip(:, q, :), num_slip_elt, step)
    call sim_write_elts_slip(sim, mesh, 'sliprate', state%sliprate(:, q, :), num_slip_elt, step)
    call sim_write_elts_slip(sim, mesh, 'rss', state%rss(:, q, :), num_slip_elt, step)

    ! Plastic strain (crystal-frame sum, in the sample frame of the current orientation) and
    ! its equivalent value; plastic deformation rate.
    do i_elt = 1_ik, mesh%num_elts
      out6(:, i_elt) = transform_crystal_to_sample_mandel(rot_c(:, :, i_elt), &
        & state%strain_pl_c(:, q, i_elt))
      out1(1, i_elt) = equiv_strain(out6(:, i_elt))
    end do
    call sim_write_elts(sim, mesh, 'strain_pl', out6, step)
    call sim_write_elts_values(sim, mesh, 'strain_pl_eq', out1, step)
    call sim_write_elts(sim, mesh, 'defrate_pl', state%defrate_pl_s(:, q, :), step)

    ! Plastic spin rate, components 12 13 23 of the skew tensor (stored w23 w13 w12).
    out3(1, :) = state%spin_pl_s(3, q, :)
    out3(2, :) = state%spin_pl_s(2, q, :)
    out3(3, :) = state%spin_pl_s(1, q, :)
    call sim_write_elts_values(sim, mesh, 'spinrate', out3, step)

    ! Lattice rotation rate, axial vector (-w23, w13, -w12).
    out3(1, :) = -state%rotrate_s(1, q, :)
    out3(2, :) = state%rotrate_s(2, q, :)
    out3(3, :) = -state%rotrate_s(3, q, :)
    call sim_write_elts_values(sim, mesh, 'rotrate', out3, step)

    ! Work.
    out1(1, :) = state%work(q, :)
    call sim_write_elts_values(sim, mesh, 'work', out1, step)

    ! Equivalent stress and elastic strain.
    do i_elt = 1_ik, mesh%num_elts
      out1(1, i_elt) = equiv_stress(stress_c(:, i_elt))
    end do
    call sim_write_elts_values(sim, mesh, 'stress_eq', out1, step)
    do i_elt = 1_ik, mesh%num_elts
      out1(1, i_elt) = equiv_strain(strain_el_c(:, i_elt))
    end do
    call sim_write_elts_values(sim, mesh, 'strain_el_eq', out1, step)

    ! Total strain (sample-frame sum of d dt) and its equivalent value.
    call sim_write_elts(sim, mesh, 'strain', state%strain_s(:, q, :), step)
    do i_elt = 1_ik, mesh%num_elts
      out1(1, i_elt) = equiv_strain(state%strain_s(:, q, i_elt))
    end do
    call sim_write_elts_values(sim, mesh, 'strain_eq', out1, step)

    ! Plastic strain accumulated in the sample frame (D53) and its equivalent value.
    call sim_write_elts(sim, mesh, 'strain_pl_sample', state%strain_pl_s(:, q, :), step)
    do i_elt = 1_ik, mesh%num_elts
      out1(1, i_elt) = equiv_strain(state%strain_pl_s(:, q, i_elt))
    end do
    call sim_write_elts_values(sim, mesh, 'strain_pl_sample_eq', out1, step)

    ! Deformation rate, its equivalent value (deviatoric part), and that of d^p.
    call sim_write_elts(sim, mesh, 'defrate', state%defrate_s(:, q, :), step)
    do i_elt = 1_ik, mesh%num_elts
      out1(1, i_elt) = equiv_strain(state%defrate_s(:, q, i_elt))
    end do
    call sim_write_elts_values(sim, mesh, 'defrate_eq', out1, step)
    do i_elt = 1_ik, mesh%num_elts
      out1(1, i_elt) = equiv_strain(state%defrate_pl_s(:, q, i_elt))
    end do
    call sim_write_elts_values(sim, mesh, 'defrate_pl_eq', out1, step)

    ! Velocity gradient L = d + w, row by row (11 12 13 21 22 23 31 32 33).
    do i_elt = 1_ik, mesh%num_elts
      velgrad = matrix_from_mandel(state%defrate_s(:, q, i_elt)) + &
        & matrix_from_skew(state%spin_s(:, q, i_elt))
      out9(:, i_elt) = reshape(transpose(velgrad), [9])
    end do
    call sim_write_elts_values(sim, mesh, 'velgrad', out9, step)

    ! Spin and slip parts of the rotation rate, axial vectors (-w23, w13, -w12) of the w^p and
    ! the slip spin of the orientation update (v2 d_rstar_spin, d_rstar_slip).
    out3(1, :) = -state%spin_hat_update_s(1, q, :)
    out3(2, :) = state%spin_hat_update_s(2, q, :)
    out3(3, :) = -state%spin_hat_update_s(3, q, :)
    call sim_write_elts_values(sim, mesh, 'rotrate_spin', out3, step)
    out3(1, :) = -state%spin_slip_update_s(1, q, :)
    out3(2, :) = state%spin_slip_update_s(2, q, :)
    out3(3, :) = -state%spin_slip_update_s(3, q, :)
    call sim_write_elts_values(sim, mesh, 'rotrate_slip', out3, step)

    ! Plastic work, work rate and plastic work rate.
    out1(1, :) = state%work_pl(q, :)
    call sim_write_elts_values(sim, mesh, 'work_pl', out1, step)
    out1(1, :) = state%workrate(q, :)
    call sim_write_elts_values(sim, mesh, 'workrate', out1, step)
    out1(1, :) = state%workrate_pl(q, :)
    call sim_write_elts_values(sim, mesh, 'workrate_pl', out1, step)

    ! Kocks-Mecking (step 7.2, D60; names and the rate step 7.3, D61): network dislocation
    ! density, recovery coefficient k_2 and density rate (zero for elements of other phases,
    ! and k_2 and the rate where Gammadot = 0).
    out1(1, :) = state%km_rho(q, :)
    call sim_write_elts_values(sim, mesh, 'km_rho', out1, step)
    out1(1, :) = state%km_k2(q, :)
    call sim_write_elts_values(sim, mesh, 'km_k2', out1, step)
    out1(1, :) = state%km_rhorate(q, :)
    call sim_write_elts_values(sim, mesh, 'km_rhorate', out1, step)

    ! Irradiation (step 7.3, D61): the damage tensor H and its rate, tensor components in the
    ! crystal frame (zero for elements without irradiation).
    call sim_write_elts(sim, mesh, 'dmg', state%dmg_c(:, q, :), step)
    call sim_write_elts(sim, mesh, 'dmgrate', state%dmgrate_c(:, q, :), step)

    ! The step's end (HDF5: its row of /results/steps with the time, the XDMF file; step 9.1).
    call sim_end_step(sim, step, time)

  end subroutine write_step

  !> Run limits after an increment that did not complete the run (D44, v2-fixed C32):
  !> max_incr, max_total_time, and, when given, max_strain (largest |normal engineering
  !> strain| of the bounding box) and max_eqstrain (sqrt(2/3 sum e_ii^2), v2 driver_triaxcsr).
  !> A tripped limit writes the results of the increment as the current step (unless they were
  !> just written) and ends the run with an error.
  !>
  !> Arguments:
  !>   i_step   (in)  integer  current step
  !>   written  (in)  logical  the results of this increment were just written as the step's
  !> Locals:
  !>   strain   real(rk)(3)    normal engineering strains of the bounding box (periodic: of
  !>                           Fbar - I) [-]
  !>   fbar     real(rk)(3, 3) periodic: Fbar [-]
  !>   eqstrain real(rk)       equivalent strain [-]
  !>   reason   character(200) the tripped limit (empty: none)
  subroutine check_limits(i_step, written)

    integer, intent(in) :: i_step
    logical, intent(in) :: written

    real(rk) :: strain(3), eqstrain, fbar(3, 3)
    character(len=200) :: reason

    ! Increments and time.
    reason = ''
    if (incr >= config%max_incr) write (reason, '(a,i0,a)') 'max_incr ', config%max_incr, &
      & ' reached'
    if (time >= config%max_total_time) write (reason, '(a,es11.4,a)') 'max_total_time ', &
      & config%max_total_time, ' s reached'

    ! Macroscopic strains (only when a strain limit is given: collective).
    if (config%has_max_strain .or. config%has_max_eqstrain) then
      if (config%bc_type == 'periodic') then
        ! Periodic: the normal components of Fbar - I (step 8.3b).
        fbar = constraint_macro_rate(cons, reshape(coo, [3_ik*mesh%num_nodes_loc]))
        strain = [fbar(1, 1), fbar(2, 2), fbar(3, 3)] - 1.0_rk
      else
        strain = bc_lengths(mesh, coo)/lengths0 - 1.0_rk
      end if
      eqstrain = sqrt(2.0_rk/3.0_rk*sum(strain**2))
      if (config%has_max_strain) then
        if (maxval(abs(strain)) >= config%max_strain) write (reason, '(a,es11.4,a,es11.4)') &
          & 'max_strain ', config%max_strain, ' reached: strain ', maxval(abs(strain))
      end if
      if (config%has_max_eqstrain) then
        if (eqstrain >= config%max_eqstrain) write (reason, '(a,es11.4,a,es11.4)') &
          & 'max_eqstrain ', config%max_eqstrain, ' reached: equivalent strain ', eqstrain
      end if
    end if

    ! A tripped limit: the results of this increment as the current step, then stop.
    if (len_trim(reason) == 0) return
    if (.not. written) call write_step(i_step)
    call sim_close(sim)
    write (text, '(a,i0,a,i0,a,es11.4,a)') 'stopped in step ', step_offset + i_step, &
      & ' at increment ', incr, ' (t = ', time, ' s): '
    call fatal_error('fepx', trim(text) // ' ' // trim(reason))

  end subroutine check_limits

  !> The macroscopic row of periodic conditions (step 8.3b, Q-P11): Fbar from the drive pairs'
  !> current separations, Fbar-dot from their velocities (R-P4 Eq. 72) and the volume-averaged
  !> Cauchy stress (zero at step 0), appended to results/macro. Collective; nothing without
  !> periodic conditions.
  !>
  !> Arguments:
  !>   step   (in)  integer  step (0: the initial state)
  !>   incr_w (in)  integer  increment (cumulative)
  !> Locals:
  !>   fbar, fdot  real(rk)(3, 3)  Fbar [-] and Fbar-dot [1/s]
  !>   stress_avg  real(rk)(6)     volume-averaged Cauchy stress [MPa]
  subroutine write_macro(step, incr_w)

    integer, intent(in) :: step, incr_w

    real(rk) :: fbar(3, 3), fdot(3, 3), stress_avg(6)

    if (config%bc_type /= 'periodic') return
    fbar = constraint_macro_rate(cons, reshape(coo, [3_ik*mesh%num_nodes_loc]))
    fdot = constraint_macro_rate(cons, vel)
    stress_avg = 0.0_rk
    if (step > 0) stress_avg = periodic_stress_average(mesh, elt, coo, state%stress_s)
    call sim_write_macro(sim, step, incr_w, time, fbar, fdot, stress_avg)

  end subroutine write_macro

  !> The previous step of a step: step i_cur - 1, or, for step 1 of a restarted run, the
  !> restart point's step (step 8.4).
  !>
  !> Arguments:
  !>   i_cur  (in)  integer              the step (> 1, or 1 of a restarted run)
  !>   rec    (out, result) type(step_record_t)  its previous step
  function previous_step(i_cur) result(rec)

    integer, intent(in) :: i_cur
    type(step_record_t) :: rec

    if (i_cur == 1) then
      rec = sched%prev
    else
      rec = step_schedule_record(sched, config, i_cur - 1)
    end if

  end function previous_step

  !> Write the restart point at the end of step i_done (print restart, step 8.4,
  !> addendum_restart.md §4, §5): every rank its file restart.rank<r> in simulation.sim/restart,
  !> written as <file>.new and renamed over the previous one after every rank has written.
  !> Records, in order: the run (the step's number, increments, time, strains reached, the
  !> control-face force, successive subdivided increments), the step's schedule record, the
  !> output files' counts (broadcast from the root), the configuration's fingerprint, the
  !> triaxial and load-rate carried values, then the rank's mesh part (global element and node ids,
  !> original coordinates, phases, orientations: a check), the coordinates, the velocities, the
  !> whole material state, and the increments again (end marker). Collective.
  !>
  !> Arguments:
  !>   i_done  (in)  integer  the step just completed
  !> Locals:
  !>   rst      type(restart_file_t)       the file
  !>   name     character(:), allocatable  its name
  !>   counts   integer(ik)(7)             lines of the force, convergence and macroscopic
  !>                                       files, largest step written, rows of the HDF5 force,
  !>                                       convergence and macroscopic tables (root's)
  !>   conv_attempts, conv_bc  logical     the convergence file's columns (root's)
  !>   j        integer                    fingerprint line
  subroutine restart_write(i_done)

    integer, intent(in) :: i_done

    type(restart_file_t) :: rst
    character(len=:), allocatable :: name
    integer(ik) :: counts(7)
    logical :: conv_attempts, conv_bc
    integer :: j

    ! The output files' counts, known by the root.
    counts = [sim%force_lines, sim%conv_lines, sim%macro_lines, sim%num_steps, sim%force_rows, &
      & sim%conv_rows, sim%macro_rows]
    call par_bcast(counts)
    conv_attempts = sim%conv_attempts
    conv_bc = sim%conv_bc
    call par_bcast(conv_attempts)
    call par_bcast(conv_bc)

    ! This rank's file, as <name>.new.
    name = sim_restart_dir(sim) // '/restart.rank' // int_str(int(my_rank, ik))
    call restart_create(rst, name // '.new', num_ranks, my_rank)

    ! The run and the step just completed.
    call restart_put(rst, 'run%step', step_offset + i_done)
    call restart_put(rst, 'run%incr', incr)
    call restart_put(rst, 'run%time', time)
    call restart_put(rst, 'run%strain_now', strain_now)
    call restart_put(rst, 'run%force_prev', force_prev)
    call restart_put(rst, 'run%num_successive', num_successive)
    call restart_put_step(rst, step_schedule_record(sched, config, i_done))

    ! The output files.
    call restart_put(rst, 'sim%counts', counts)
    call restart_put(rst, 'sim%conv_attempts', conv_attempts)
    call restart_put(rst, 'sim%conv_bc', conv_bc)

    ! The configuration's fingerprint.
    call restart_put(rst, 'cfg%num_lines', size(config%fingerprint))
    do j = 1, size(config%fingerprint)
      call restart_put(rst, 'cfg%line', trim(config%fingerprint(j)))
    end do
    call restart_put(rst, 'cfg%strict', config%fingerprint_strict)

    ! Triaxial load control and load rates: the values carried from step to step.
    call restart_put(rst, 'bc%vel_face', vel_face)
    call restart_put(rst, 'bc%bc_delta', bc_delta)
    call restart_put(rst, 'bc%vel_dwell_load', vel_dwell_load)
    call restart_put(rst, 'bc%vel_dwell_unload', vel_dwell_unload)
    call restart_put(rst, 'bc%action', action)
    call restart_put(rst, 'bc%prev_action', prev_action)
    call restart_put(rst, 'bc%stress_max', stress_max)

    ! The rank's mesh part (checked at reading), then the fields and the state.
    call restart_put(rst, 'mesh%elt_id_glob', mesh%elt_id_glob)
    call restart_put(rst, 'mesh%node_id_glob', mesh%node_id_glob)
    call restart_put(rst, 'mesh%coo_s', mesh%coo_s)
    call restart_put(rst, 'mesh%elt_phase', mesh%elt_phase)
    call restart_put(rst, 'mesh%rot_c2s', mesh%rot_c2s)
    call restart_put(rst, 'coo', coo)
    call restart_put(rst, 'vel', vel)
    call restart_put_state(rst, state)
    call restart_put(rst, 'end%incr', incr)
    call restart_close(rst)

    ! Every rank has written: the new files replace the old ones.
    call par_barrier()
    call restart_rename(name // '.new', name)
    call par_barrier()
    write (text, '(a,i0,a,i0,a)') 'Restart point written (step ', step_offset + i_done, &
      & ', increment ', incr, ')'
    call log_info(trim(text), 2)

  end subroutine restart_write

  !> Read the run part of this rank's restart file (restart on, step 8.4): the restart point's
  !> counters and its step (prev_rec; into the rst_* variables), the output files' counts, and
  !> the configuration's fingerprint, compared with this run's
  !> (check_fingerprint). The file stays open for restart_read_fields. Collective (the restart
  !> point must be the same on every rank).
  !>
  !> Locals:
  !>   name      character(:), allocatable  the file
  !>   exists    logical                    it exists
  !>   counts    integer(ik)(7)             the output files' counts (restart_write)
  !>   num_lines, j  integer                fingerprint lines, index
  !>   lines     character(512)(:)          the restart point's fingerprint
  !>   strict    logical(:)                 its strict lines
  !>   line      character(:), allocatable  one line read
  !>   incr_min, incr_max  integer(ik)      the restart point's increments over the ranks
  subroutine restart_read_globals()

    character(len=:), allocatable :: name, line
    logical :: exists
    integer(ik) :: counts(7), incr_min, incr_max
    integer :: num_lines, j
    character(len=512), allocatable :: lines(:)
    logical, allocatable :: strict(:)

    ! This rank's file.
    name = 'simulation.sim/restart/restart.rank' // int_str(int(my_rank, ik))
    inquire (file=name, exist=exists)
    if (.not. exists) call fatal_error('fepx', 'restart on: no restart file ' // name // &
      & ' (the run to continue needs print restart)')
    call restart_open(rst_in, name, num_ranks, my_rank)
    restarted = .true.

    ! The run and the restart point's step.
    call restart_get(rst_in, 'run%step', step_offset)
    rst_incr = 0
    call restart_get(rst_in, 'run%incr', rst_incr)
    rst_time = 0.0_rk
    call restart_get(rst_in, 'run%time', rst_time)
    rst_strain_now = 0.0_rk
    call restart_get(rst_in, 'run%strain_now', rst_strain_now)
    rst_force_prev = 0.0_rk
    call restart_get(rst_in, 'run%force_prev', rst_force_prev)
    rst_num_successive = 0
    call restart_get(rst_in, 'run%num_successive', rst_num_successive)
    call restart_get_step(rst_in, prev_rec)

    ! The same restart point on every rank (an interrupted renaming would mix two).
    call par_min(rst_incr, incr_min)
    call par_max(rst_incr, incr_max)
    if (incr_min /= incr_max) call fatal_error('fepx', 'restart on: the restart files of ' // &
      & 'the ranks are of different restart points (increments ' // int_str(incr_min) // &
      & ' to ' // int_str(incr_max) // ')')

    ! The output files.
    counts = 0_ik
    call restart_get(rst_in, 'sim%counts', counts)
    rst_lines = counts(1:3)
    rst_num_steps = counts(4)
    rst_rows = counts(5:7)
    rst_conv = .false.
    call restart_get(rst_in, 'sim%conv_attempts', rst_conv(1))
    call restart_get(rst_in, 'sim%conv_bc', rst_conv(2))

    ! The configuration's fingerprint, compared with this run's.
    num_lines = 0
    call restart_get(rst_in, 'cfg%num_lines', num_lines)
    allocate (lines(num_lines), strict(num_lines))
    do j = 1, num_lines
      call restart_get(rst_in, 'cfg%line', line)
      lines(j) = line
    end do
    call restart_get(rst_in, 'cfg%strict', strict)
    call log_info('Restart: continuing from step ' // log_int_text(step_offset) // &
      & ' (increment ' // log_int_text(rst_incr) // ', t = ' // log_real_text(rst_time, 4) // &
      & ' s); the ' // log_count_text(config%num_steps, 'step') // ' of simulation.cfg ' // &
      & 'follow it', 1)
    call check_fingerprint(lines, strict)

  end subroutine restart_read_globals

  !> Compare the restart point's configuration fingerprint with this run's (step 8.4,
  !> addendum_restart.md §3 R1): a strict line (material, boundary conditions and constraints,
  !> control component, periodic kinematics, input files) of either that is not in the other is
  !> an error naming it; any other difference a warning (the continuation is then not that of
  !> an uninterrupted run).
  !>
  !> Arguments:
  !>   lines   (in)  character(*)(:)  the restart point's lines
  !>   strict  (in)  logical(:)       its strict lines
  !> Locals:
  !>   j       integer                line
  subroutine check_fingerprint(lines, strict)

    character(len=*), intent(in) :: lines(:)
    logical, intent(in) :: strict(:)

    integer :: j

    ! Lines of the restart point not in this configuration.
    do j = 1, size(lines)
      if (any(config%fingerprint == lines(j))) cycle
      if (strict(j)) call fatal_error('fepx', 'restart on: "' // trim(lines(j)) // &
        & '" of the restart point''s configuration is not in simulation.cfg (the material, ' // &
        & 'the boundary conditions, the input files and the output format must not change)')
      call log_warning('Restart on: "' // trim(lines(j)) // '" of the restart point''s ' // &
        & 'configuration is not in simulation.cfg.', 2)
    end do

    ! Lines of this configuration not in the restart point's.
    do j = 1, size(config%fingerprint)
      if (any(lines == config%fingerprint(j))) cycle
      if (config%fingerprint_strict(j)) call fatal_error('fepx', 'restart on: "' // &
        & trim(config%fingerprint(j)) // '" in simulation.cfg is not in the restart ' // &
        & 'point''s configuration (the material, the boundary conditions, the input files ' // &
        & 'and the output format must not change)')
      call log_warning('Restart on: "' // trim(config%fingerprint(j)) // '" in ' // &
        & 'simulation.cfg is not in the restart point''s configuration.', 2)
    end do

  end subroutine check_fingerprint

  !> Read the fields part of this rank's restart file (restart on, step 8.4) and close it: the
  !> triaxial and load-rate carried values, the rank's mesh part (it must equal this run's: the
  !> same mesh and partition), then the
  !> coordinates, the velocities and the material state; the carried values of
  !> restart_read_globals are set (time, increments, strains reached, the control-face force,
  !> successive subdivided increments, the perturbation floors: the larger of the restart
  !> point's and this configuration's).
  !>
  !> Locals:
  !>   ids      integer(ik)(:)       global ids read
  !>   reals2   real(rk)(:, :)       original coordinates read [mm]
  !>   reals3   real(rk)(:, :, :)    (3, 3, num_elts) mesh orientations read [-]
  !>   end_incr integer(ik)          the end marker
  !>   differs  logical              the rank's mesh part differs
  !>   any_differs  logical          some rank's differs
  subroutine restart_read_fields()

    integer(ik), allocatable :: ids(:)
    real(rk), allocatable :: reals2(:, :), reals3(:, :, :)
    integer(ik) :: end_incr
    logical :: differs, any_differs

    ! Triaxial load control and load rates: the carried values (after the program's set-up of
    ! triaxial load control, which starts them at zero).
    call restart_get(rst_in, 'bc%vel_face', vel_face)
    call restart_get(rst_in, 'bc%bc_delta', bc_delta)
    call restart_get(rst_in, 'bc%vel_dwell_load', vel_dwell_load)
    call restart_get(rst_in, 'bc%vel_dwell_unload', vel_dwell_unload)
    call restart_get(rst_in, 'bc%action', action)
    call restart_get(rst_in, 'bc%prev_action', prev_action)
    rst_stress_max = 0.0_rk
    call restart_get(rst_in, 'bc%stress_max', rst_stress_max)

    ! The rank's mesh part.
    differs = .false.
    allocate (ids(mesh%num_elts))
    call restart_get(rst_in, 'mesh%elt_id_glob', ids)
    differs = differs .or. any(ids /= mesh%elt_id_glob)
    deallocate (ids)
    allocate (ids(size(mesh%node_id_glob)))
    call restart_get(rst_in, 'mesh%node_id_glob', ids)
    differs = differs .or. any(ids /= mesh%node_id_glob)
    allocate (reals2(3, size(mesh%coo_s, 2)))
    call restart_get(rst_in, 'mesh%coo_s', reals2)
    differs = differs .or. any(abs(reals2 - mesh%coo_s) > 0.0_rk)
    deallocate (ids)
    allocate (ids(mesh%num_elts))
    call restart_get(rst_in, 'mesh%elt_phase', ids)
    differs = differs .or. any(ids /= mesh%elt_phase)
    allocate (reals3(3, 3, mesh%num_elts))
    call restart_get(rst_in, 'mesh%rot_c2s', reals3)
    differs = differs .or. any(abs(reals3 - mesh%rot_c2s) > 0.0_rk)
    call par_any(differs, any_differs)
    if (any_differs) call fatal_error('fepx', 'restart on: the mesh (nodes, ' // &
      & 'coordinates, phases or orientations) or its partition differs from the restart ' // &
      & 'point''s')

    ! Coordinates, velocities, state; the end marker.
    call restart_get(rst_in, 'coo', coo)
    call restart_get(rst_in, 'vel', vel)
    call restart_get_state(rst_in, state)
    end_incr = -1_ik
    call restart_get(rst_in, 'end%incr', end_incr)
    if (end_incr /= rst_incr) call fatal_error('fepx', 'restart on: incomplete restart file')
    call restart_close(rst_in)

    ! The carried values.
    time = rst_time
    incr = rst_incr
    strain_now = rst_strain_now
    force_prev = rst_force_prev
    num_successive = rst_num_successive
    stress_max = max(stress_max, rst_stress_max)

  end subroutine restart_read_fields

  !> The command line (step 9.7, ROADMAP D94, Q64; dev/plan/terminal_output.md §3.6), handled
  !> before MPI starts (so that it works without mpirun): "--help" writes the banner without the
  !> run lines and the usage, "--version" the version text (build_info_mod), each then ends the
  !> program with exit status 0; any other argument (or a second one) is a fatal error, which
  !> ends the program with exit status 1 (command_line_stop). Without arguments it returns.
  !>
  !> Locals:
  !>   num_args  integer                    number of arguments
  !>   length    integer                    length of an argument
  !>   arg       character(:), allocatable  the first argument, then the one in error
  subroutine command_line()

    integer :: num_args, length
    character(len=:), allocatable :: arg

    ! No argument: run the simulation.
    num_args = command_argument_count()
    if (num_args == 0) return
    call get_command_argument(1, length=length)
    allocate (character(len=length) :: arg)
    call get_command_argument(1, arg)

    ! --help and --version, alone.
    if (num_args == 1 .and. arg == '--help') then
      call write_banner(.false.)
      call log_text(repeat('=', 72))
      call log_text('Usage: mpirun -np <ranks> fepx    run the simulation in the current')
      call log_text('                                  directory (simulation.cfg,')
      call log_text('                                  simulation.msh, and the optional')
      call log_text('                                  simulation.ori, .phase, .opt)')
      call log_text('       fepx --help                print this help')
      call log_text('       fepx --version             print the version')
      call log_text('See also: man fepx, the manual <https://fepx.info>.')
      stop
    end if
    if (num_args == 1 .and. arg == '--version') then
      write (output_unit, '(a)') build_version_text()
      stop
    end if

    ! Any other argument is an error (the second one after --help or --version).
    if (arg == '--help' .or. arg == '--version') then
      deallocate (arg)
      call get_command_argument(2, length=length)
      allocate (character(len=length) :: arg)
      call get_command_argument(2, arg)
    end if
    call set_abort_handler(command_line_stop)
    call fatal_error('', 'unknown option "' // arg // '" (see fepx --help)')

  end subroutine command_line

  !> Abort handler of the command line (before MPI starts): ends the program with the given
  !> exit status, without the run-time library's message (fatal_error has written the error).
  !>
  !> Arguments:
  !>   exit_status  (in)  integer  exit status of the program
  subroutine command_line_stop(exit_status)

    integer, intent(in) :: exit_status

    ! End quietly with the status.
    stop exit_status, quiet=.true.

  end subroutine command_line_stop

  !> Write the banner (step 9.7, Q58, Q59; dev/plan/terminal_output.md §3.1), as Neper's line
  !> for line: the title rule, the description, the version with the git hash and the commit
  !> date, the build (compiler, MPI library, HDF5 when built with it, separated by "|"), the
  !> number of ranks (run only), the URL, the copyright, the closing line.
  !>
  !> Arguments:
  !>   run  (in)  logical  a run (with the line of the number of ranks), not --help
  !> Locals:
  !>   compiler  character(:), allocatable  the compiler as name and version
  !>   i         integer                    position of " version " in it
  subroutine write_banner(run)

    logical, intent(in) :: run

    character(len=:), allocatable :: compiler
    integer :: i

    ! The compiler as name and version ("GCC version 15.2.0" -> "GCC 15.2.0").
    compiler = compiler_version()
    i = index(compiler, ' version ')
    if (i > 0) compiler = compiler(1:i) // compiler(i + 9:)

    call log_text('========================    F   E   P   X    ===========================')
    call log_info('A finite element software package for polycrystal plasticity.')
    call log_info('Version ' // build_version_text())
    if (hdf5_available) then
      call log_info('Built with: ' // compiler // '|' // par_library_version() // '|' // &
        & h5_library_version())
    else
      call log_info('Built with: ' // compiler // '|' // par_library_version())
    end if
    if (run) call log_info('Running on ' // log_count_text(int(num_ranks), 'rank') // '.')
    call log_info('<https://fepx.info>')
    call log_info('Copyright (C) 1996-2026, and MIT licensed, by the FEPX authors.')
    call log_info('---------------------------------------------------------------')

  end subroutine write_banner

  !> Write the setup summary (step 9.7, ROADMAP D94, dev/plan/terminal_output.md §3.2
  !> revision 7): from the configuration, the mesh, its distribution and the step schedule,
  !> after the setup, with the log lines the setup kept for it (hold mode): the material
  !> parameters of each phase followed by its derived quantities; "- Loading" with the
  !> constraints, triaxial estimates and step-schedule lines, then the plan of each step known
  !> before the run; "- Restart" (restarted runs); the options other than the defaults and the
  !> output; "- Mesh parameters"; "- Partitioning", with D30's warning when the mesh's
  !> partitions are not the ranks (Q65), the elements, ghost nodes and neighbors per rank.
  !> Collective (the per-rank statistics).
  !>
  !> Locals:
  !>   j                  integer      phase or step
  !>   num_min, num_max, num_sum  integer(ik)  a statistic over the ranks
  subroutine log_setup_summary()

    integer :: j
    integer(ik) :: num_min, num_max, num_sum

    ! Material: each phase's parameters, then its derived quantities (kept at the setup).
    call config_log_material(config)
    do j = 1, config%num_phases
      call config_log_phase(config, j)
      call log_lines_write(phase_held(j)%text)
    end do

    ! Loading: the conditions, the kept lines (constraints, estimates), the planned steps.
    call log_info('Loading: boundary_conditions ' // trim(config%bc_type) // ', ' // &
      & log_count_text(config%num_steps, 'step'), 1)
    call log_lines_write(loading_held)
    do j = 1, sched%num_planned
      call log_step_plan(j, 2)
    end do

    ! Restart (its own lines, kept).
    call log_lines_write(restart_held)

    ! Options and output.
    call config_log_options(config)

    ! Mesh.
    call log_info('Mesh parameters:', 1)
    call log_info('Node number: ' // log_int_text(mesh%num_nodes_glob), 2)
    call log_info('Elt  number: ' // log_int_text(mesh%num_elts_glob) // ' (' // &
      & elt_type_name(mesh%elt_type) // ')', 2)
    call log_info('Elsets: ' // log_int_text(mesh%num_elsets) // '; partitions: ' // &
      & log_int_text(mesh%num_parts_file), 2)

    ! The ranks, the warning for a partitioned mesh on 1 rank (D108), then the statistics over
    ! the ranks.
    call log_info('Partitioning: ' // log_count_text(num_ranks, 'rank'), 1)
    if (num_ranks == 1_ik .and. mesh%num_parts_file > 1_ik) then
      call log_warning('The mesh has ' // log_count_text(mesh%num_parts_file, 'partition') // &
        & ', but 1 rank was started: the rank takes the whole mesh.', 2)
    end if
    call par_min(mesh%num_elts, num_min)
    call par_max(mesh%num_elts, num_max)
    call log_info('Elements per rank: ' // log_int_text(num_min) // ' to ' // &
      & log_int_text(num_max), 2)
    call par_min(mesh%num_nodes_ghost, num_min)
    call par_max(mesh%num_nodes_ghost, num_max)
    call par_sum(mesh%num_nodes_ghost, num_sum)
    call log_info('Ghost nodes per rank: ' // log_int_text(num_min) // ' to ' // &
      & log_int_text(num_max) // ' (' // log_int_text(num_sum) // ' in total)', 2)
    call par_min(mesh%comm%num_neighbors, num_min)
    call par_max(mesh%comm%num_neighbors, num_max)
    call log_info('Neighbors per rank: ' // log_int_text(num_min) // ' to ' // &
      & log_int_text(num_max), 2)

  end subroutine log_setup_summary

  !> Write the closing summary (step 9.7, Q63; dev/plan/terminal_output.md §3.4): the end
  !> line, the results file, the increments, global iterations (NR) and CG iterations of the run
  !> (after the initialization), the warnings issued (all ranks), the wall time per phase (rank
  !> 0's clock: input and setup, initialization, increments with their element and material
  !> work, CG and the rest, output) and the elapsed wall time, the times aligned on their "s",
  !> then the closing rule.
  !>
  !> Locals:
  !>   time_increments  real(rk)  wall time of the increments, output excluded [s]
  subroutine log_closing_summary()

    real(rk) :: time_increments

    call log_info('Final step terminated. Simulation completed successfully.')
    call log_file('o', "Wrote results to `simulation.sim'.")
    call log_info('Increments: ' // log_int_text(incr - incr_start), 1)
    call log_info('Global iterations: ' // log_int_text(num_iters_run) // ' (' // &
      & log_int_text(num_nr_run) // ' NR)', 1)
    call log_info('CG iterations: ' // log_int_text(num_cg_run), 1)
    call log_info('Warnings issued: ' // log_int_text(num_warnings), 1)
    call log_info('Wall time per phase:', 1)
    time_increments = time_loop_end - time_init_end - (time_output_run)
    call time_line(2, 'Input and setup:', time_setup_end - time_start - time_output_setup)
    call time_line(2, 'Initialization:', time_init_end - time_setup_end)
    call time_line(2, 'Increments:', time_increments)
    call time_line(3, 'Element and material work:', time_elements_run)
    call time_line(3, 'CG:', time_cg_run)
    call time_line(3, 'Other:', max(0.0_rk, time_increments - time_elements_run - time_cg_run))
    call time_line(2, 'Output:', time_output)
    call time_line(0, 'Elapsed time:', par_wtime() - time_start)
    call log_text(repeat('=', 72))

  end subroutine log_closing_summary

  !> One line of the wall-time table of the closing summary: the label, then the time in
  !> seconds (3 decimals) right-aligned so that its last digit is in column 58 and its "s" in
  !> column 60 (dev/plan/terminal_output.md, line rules).
  !>
  !> Arguments:
  !>   depth    (in)  integer    depth of the line (0, 2 or 3)
  !>   label    (in)  character  the label, with its colon
  !>   seconds  (in)  real(rk)   the time [s]
  !> Locals:
  !>   value    character(:), allocatable  the time as text
  !>   lead     integer                    width of the depth's indentation and marker
  subroutine time_line(depth, label, seconds)

    integer, intent(in) :: depth
    character(len=*), intent(in) :: label
    real(rk), intent(in) :: seconds

    character(len=:), allocatable :: value
    integer :: lead

    value = seconds_text(seconds)
    lead = 0
    if (depth > 0) lead = 2*depth + 2
    call log_info(label // repeat(' ', max(1, 58 - 9 - lead - len(label) - len(value))) // &
      & value // ' s', depth)

  end subroutine time_line

  !> A wall time in seconds as text, 3 decimals, with the leading zero ("0.812").
  !>
  !> Arguments:
  !>   seconds  (in)  real(rk)  the time [s]
  !> Result:
  !>   text     character(:), allocatable  the time as text
  !> Locals:
  !>   buffer   character(24)  the formatted time
  function seconds_text(seconds) result(text)

    real(rk), intent(in) :: seconds
    character(len=:), allocatable :: text

    character(len=24) :: buffer

    ! F with a width that keeps the leading zero, then without the blanks.
    write (buffer, '(f24.3)') seconds
    text = trim(adjustl(buffer))

  end function seconds_text

  !> An integer as text, without blanks (for messages).
  !>
  !> Arguments:
  !>   value  (in)  integer(ik)  the integer
  !>   text   (out, result) character(:), allocatable
  !> Locals:
  !>   buffer  character(24)  formatted value
  pure function int_str(value) result(text)

    integer(ik), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=24) :: buffer

    write (buffer, '(i0)') value
    text = trim(buffer)

  end function int_str

end program fepx
