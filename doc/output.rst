.. _simulation_results:

Simulation Results (:file:`simulation.cfg`)
===========================================

The results to write are selected with one ``print`` line per variable::

    print <variable>
    print <variable>
    ...

They are written in the simulation directory :file:`simulation.sim` (see
:ref:`simulation_directory`) at the initial state (step 0) and at the end of every step whose
``print_results`` is 1, one file per variable and step. Face forces are written at every
increment. With ``output_format hdf5`` they are written to one HDF5 file instead, or as well
(see :ref:`output_formats`).

.. confval:: print

   A result variable to write (one ``print`` line per variable; the variables are described
   below). ``print restart`` writes restart points instead (:ref:`restart`).

   .. include:: generated/keywords/print.inc


Elemental quantities are computed at every quadrature point of an element, and written once per
element, at the quadrature point at the element centroid.

Values are written unrounded, with 17 significant digits (an exact copy of the double-precision
value), in the format ``es24.16``.

Components are written as follows:

- vector :math:`v`: :math:`v_1, v_2, v_3`, in the sample frame;
- symmetric tensor :math:`t`: :math:`t_{11}, t_{22}, t_{33}, t_{23}, t_{13}, t_{12}` (tensor
  components, not engineering shears), in the sample frame;
- the rates of rotation: ``spinrate`` as the skew-tensor components
  :math:`w_{12}, w_{13}, w_{23}`, ``rotrate``, ``rotrate_spin`` and ``rotrate_slip`` as
  rotation vectors (see below);
- the velocity gradient :math:`\mathbf L` (``velgrad``): its nine components
  :math:`L_{11}, L_{12}, L_{13}, L_{21}, \ldots, L_{33}`, :math:`L_{ij} = \partial v_i /
  \partial x_j`, in the sample frame;
- equivalent values (``*_eq``): :math:`\sqrt{2/3\; \mathbf a' : \mathbf a'}` of the
  deviatoric part :math:`\mathbf a'` of a strain or rate, :math:`\sqrt{3/2\;
  \boldsymbol\sigma' : \boldsymbol\sigma'}` (von Mises) of the stress;
- slip-system values: one per slip system of the element's phase.

Nodal Results
-------------

``coo``
    Coordinates [mm].

``disp``
    Displacements from the initial coordinates [mm].

``vel``
    Velocities [mm/s].

Elemental Results
-----------------

``stress``
    Cauchy stress [MPa].

``strain_el``
    Elastic strain [-].

``ori``
    Crystal orientation, written with the descriptor and convention of the input (e.g.
    Rodrigues vector, Euler--Bunge angles in degrees, quaternion).

``stress_eq``
    Equivalent (von Mises) stress [MPa].

``strain``
    Total strain [-], symmetric tensor: :math:`\sum \mathbf d \, \Delta t`, the deformation
    rate of each increment times its time increment, summed in the sample frame.

``strain_eq``, ``strain_el_eq``
    Equivalent total and elastic strain [-].

``defrate``
    Deformation rate :math:`\mathbf d = \mathrm{sym}(\mathbf L)` [1/s], symmetric tensor
    (volumetric part included).

``defrate_eq``
    Equivalent deformation rate [1/s] (of the deviatoric part).

``velgrad``
    Velocity gradient :math:`\mathbf L` [1/s] (nine components, see above).

``workrate``
    Work rate :math:`\boldsymbol\sigma : \mathbf d` [MPa/s].

``work``
    Work per unit volume [MPa = mJ/mm³]: the time integral of the work rate
    :math:`\boldsymbol\sigma : \mathbf d` (Cauchy stress, deformation rate), by the trapezoidal
    rule over each increment, :math:`W_n = W_{n-1} + \tfrac{\Delta t}{2}(\dot W_{n-1} +
    \dot W_n)`.

``rotrate``
    Lattice rotation rate [1/s]: the rotation of the crystal lattice over the last increment,
    divided by its time increment, as the rotation vector :math:`\boldsymbol\omega` (axis times
    rate of the angle, :math:`\mathbf W \mathbf v = \boldsymbol\omega \times \mathbf v`),
    sample frame, components :math:`\omega_1, \omega_2, \omega_3`. With orientations
    :math:`\mathbf R_{n-1}` and :math:`\mathbf R_n` (crystal to sample) at the start and end of
    the increment, it is the rotation vector of :math:`\mathbf R_n \mathbf R_{n-1}^T` divided
    by :math:`\Delta t`. Zero at step 0.

Crystal Plasticity Results
~~~~~~~~~~~~~~~~~~~~~~~~~~

The following are written for elements of ``crystal_plasticity`` phases; for elements of
``elastic`` phases the slip-system results are empty lines and the plastic results zero.
Slip-system results have one value per slip system of the element's phase, in the slip-system
order of :ref:`crystal_plasticity`.

``crss``
    Slip-system strengths :math:`\tau_{\mathrm{crit}}^\alpha` [MPa] (the initial strengths at step
    0).

``rss``
    Resolved shear stresses :math:`\tau^\alpha = \mathbf s^\alpha \cdot \boldsymbol\tau
    \, \mathbf m^\alpha` [MPa] of the Kirchhoff stress :math:`\boldsymbol\tau` (the stress of
    the slip kinetics; :math:`\boldsymbol\tau = \beta \boldsymbol\sigma`, :math:`\beta` the
    elastic volume ratio).

``sliprate``
    Slip rates :math:`\dot\gamma^\alpha` [1/s], signed, from ``rss`` and ``crss`` by the slip
    law.

``slip``
    Accumulated slip :math:`\gamma^\alpha = \sum \dot\gamma^\alpha \Delta t` [-] (signed;
    the end-of-increment rates).

``defrate_pl``
    Plastic deformation rate :math:`\mathbf d^p = \sum_\alpha \dot\gamma^\alpha
    \,\mathrm{sym}(\mathbf s^\alpha \otimes \mathbf m^\alpha)` [1/s], symmetric tensor,
    sample frame.

``strain_pl``
    Plastic strain [-], symmetric tensor: :math:`\sum \mathbf d^p \Delta t` accumulated in the
    **crystal frame** and written in the sample frame of the current orientation:
    the plastic strain of the lattice, which rotates with it.

``strain_pl_eq``
    Equivalent plastic strain :math:`\sqrt{2/3\; \boldsymbol\varepsilon^p :
    \boldsymbol\varepsilon^p}` [-] of ``strain_pl`` (which is deviatoric).

``strain_pl_sample``
    Plastic strain [-], symmetric tensor: :math:`\sum \mathbf d^p \Delta t` accumulated in the
    **sample frame**: the plastic strain of the material in the specimen axes.
    It equals ``strain_pl`` as long as the lattice does not rotate, and differs from it once it
    does (a few percent of the value at a few percent strain).

``strain_pl_sample_eq``
    Equivalent plastic strain [-] of ``strain_pl_sample``.

``defrate_pl_eq``
    Equivalent plastic deformation rate :math:`\sqrt{2/3\; \mathbf d^p : \mathbf d^p}`
    [1/s].

``workrate_pl``
    Plastic work rate :math:`\boldsymbol\sigma : \mathbf d^p` [MPa/s].

``work_pl``
    Plastic work per unit volume [MPa]: the time integral of ``workrate_pl`` by the trapezoidal
    rule, as ``work``.

``spinrate``
    Plastic spin rate :math:`\mathbf w^p = \sum_\alpha \dot\gamma^\alpha
    \,\mathrm{skw}(\mathbf s^\alpha \otimes \mathbf m^\alpha)` [1/s], skew tensor, sample
    frame, written as its components :math:`w^p_{12}, w^p_{13}, w^p_{23}`. Note that
    ``rotrate`` is written as a rotation vector, whose components are
    :math:`(-w_{23}, w_{13}, -w_{12})` of its skew tensor.

``rotrate_spin``, ``rotrate_slip``
    The two parts of the lattice rotation rate [1/s], rotation vectors in the form of
    ``rotrate``, of the two spins that rotated the lattice over the increment:
    ``rotrate_spin`` of the spin :math:`\hat{\mathbf w}^p` of the stress solve (the continuum
    spin with the elastic-strain corrections of the recoupled kinematics), ``rotrate_slip`` of
    the slip spin :math:`\sum_\alpha \dot\gamma^\alpha \,\mathrm{skw}(\mathbf s^\alpha
    \otimes \mathbf m^\alpha)`, so that ``rotrate`` = ``rotrate_spin`` − ``rotrate_slip``.
    ``rotrate_slip`` is the tensor of ``spinrate``, to the tolerance of the material
    point solve (its slip rates are those of the orientation update, before the last strength
    update; ``sliprate`` and ``spinrate`` are those of the final strengths).

``km_rho``
    Network dislocation density :math:`\rho_n` [1/mm²] of ``hardening kocks_mecking``
    (:ref:`crystal_plasticity`), one value per element; 0 for elements of phases
    with another law. Requires a Kocks--Mecking phase.

``km_k2``
    The dynamic-recovery coefficient :math:`k_2 = k_{20} (\dot\gamma_{k0}/\dot\Gamma)^{1/n}`
    [-] of ``hardening kocks_mecking`` at the slip rates of the element's final stress (new in
    FEPX); with ``km_rho`` and ``km_k1`` it gives the balance of dislocation storage
    :math:`k_1 \sqrt{\rho_n}` and recovery :math:`k_2 \rho_n`. 0 where no slip occurs
    (:math:`\dot\Gamma = 0`, :math:`k_2` undefined) and for elements of phases with another
    law. Requires a Kocks--Mecking phase.

``km_rhorate``
    Rate of the network dislocation density, :math:`\dot\rho_n = \dot\Gamma (k_1 \sqrt{\rho_n}
    - k_2 \rho_n)` [1/(mm² s)], at the slip rates of the converged stress; one value per
    element. 0 where no slip occurs and for elements of phases with another hardening law.

``dmg``
    Damage tensor :math:`\mathbf H` [1/mm²] of irradiation softening (:ref:`crystal_plasticity`,
    new in FEPX), symmetric tensor, **in the crystal
    frame**, one per element; 0 for elements of phases without irradiation. Requires a phase
    with ``irradiation on``. It may become indefinite (a negative contact density
    :math:`\mathbf N^\alpha : \mathbf H`, see :ref:`crystal_plasticity`).

``dmgrate``
    Rate :math:`\dot{\mathbf H}` [1/(mm² s)] at the converged slip rates and :math:`\mathbf H`,
    same layout as ``dmg``; 0 for elements of phases without irradiation. Requires a phase
    with ``irradiation on``.

Face Results
------------

``forces``
    For each face set of the mesh (for a cubic domain ``x0``, ``x1``, ``y0``, ``y1``, ``z0``,
    ``z1``), a file :file:`results/forces/<faset>` with a header and one row per increment::

      <step> <incr> <force_x> <force_y> <force_z> <area> <time>

    where the forces [N] are integrated over the face, ``<area>`` [mm²] is the current face area,
    and ``<time>`` [s] the time at the end of the increment. Step 0 is the initial state.

With periodic conditions (:ref:`periodic_conditions`) the face forces file is not written.

Run Results
-----------

``macro`` (periodic conditions; written without a ``print`` line)
    The macroscopic quantities, in :file:`results/macro`: a header, then one row at step 0 and
    one per increment::

      <step> <incr> <time> <Fbar_11> <Fbar_12> ... <Fbar_33> <Fdot_11> ... <Fdot_33> <s_11> <s_22>
      <s_33> <s_23> <s_13> <s_12>

    with :math:`\bar{\mathbf F}` the macroscopic deformation gradient from the drive pairs'
    current separations (row by row, 11 12 13 21 ... 33), :math:`\dot{\bar{\mathbf F}}` its rate
    from their velocities [1/s], and the volume average of the Cauchy stress over the current
    configuration [MPa] (0 at step 0).

``convergence``
    The global (velocity) iterations of every increment, in the file
    :file:`results/convergence`, with a header and one row per iteration::

      <incr> <iter> <SA|NR> <r_norm> <rx_norm> <f_norm> <delu_norm> <delux_norm> <u_norm> <cg_iter>

    where ``<incr>`` is the increment (0: the initialization increment, see
    :ref:`solver_parameters`), ``<iter>`` the iteration within it, ``SA`` or ``NR`` its method
    (successive approximation or Newton--Raphson), ``<r_norm>`` and ``<rx_norm>`` the 2-norm
    and the largest absolute value of the residual [N], ``<f_norm>`` the 2-norm of the force
    vector [N], ``<delu_norm>`` and ``<delux_norm>`` the 2-norm and the largest absolute value
    of the velocity change divided by ``<u_norm>``, the 2-norm of the velocities before the
    change [mm/s], and ``<cg_iter>`` the conjugate-gradient iterations of the linear solve. The
    iteration has converged when ``<delu_norm>`` falls below ``global_tol`` (see
    :ref:`solver_parameters`). Rows are written as the run goes, so a run that stops on a failed
    increment shows its iterations.

    With ``subincrementation on`` (:ref:`solver_parameters`), the rows of every attempt of an
    increment, accepted or rejected, are written under the increment's number, with two more
    columns (named in the header)::

      ... <cg_iter> <attempt> <A|R>

    where ``<attempt>`` numbers the attempts of the increment in the order they were solved
    (1, 2, ...) and ``A`` or ``R`` says whether it was accepted or rejected. A subdivided
    increment thus shows its rejected attempt, then the attempts of its parts. A run with a
    load or stress target (:ref:`loading`) has the same two columns: an increment that passed
    its target by more than the tolerance and was solved again (landing) shows the discarded
    attempt (marked ``A``: its solve converged; the log names it), then the next attempts, the
    last being the increment written. A run with load-controlled components (triaxial load
    control, :ref:`loading`) writes every solve of the boundary-velocity iteration, with two
    more columns after these::

      ... <attempt> <A|R> <bc_iter> <solve>

    where ``<bc_iter>`` is the iteration (1, 2, ...) and ``<solve>`` is 0 for the solve at the
    iteration's face velocities and 11, 22 or 33 for the solve with that component's face
    velocity perturbed. At a load or stress rate (two loops) the iterations of the scaling loop
    come first and those of the Newton loop are numbered on; the solve at zero velocity at the
    start of a reload is iteration 0. Otherwise the file has only the ten columns of the first row
    form.


.. _output_formats:

Output Formats
--------------

The format of the results is chosen with::

    output_format <format> [ <format> ]

where ``<format>`` is ``ascii`` (the default) or ``hdf5``; ``output_format ascii hdf5`` writes
both.

``ascii`` writes the text files described above, which Neper reads.
``hdf5`` writes the same results to the HDF5 file :file:`simulation.sim/simulation.h5`, with the
XDMF file :file:`simulation.sim/simulation.xdmf` that lets ParaView open them; with ``hdf5``
alone no text result file is written (the directory holds
:file:`inputs/`, :file:`restart/` and these two files). The HDF5 output needs a build with it
(``-DFEPX_ENABLE_HDF5``, see :ref:`installation`); in a build without it, ``output_format
hdf5`` is an error. Neper does not read the HDF5 file: to visualize with ``neper -V``, use
``output_format ascii`` (or both). On a restart (:ref:`restart`) ``output_format`` must be
the same as in the run that is continued.

.. confval:: output_format

   The output formats: ``ascii`` (the text files of :file:`simulation.sim`), ``hdf5``
   (:file:`simulation.h5` and its XDMF file), or both.

   .. include:: generated/keywords/output_format.inc


The HDF5 file holds the values of the text files, unrounded (64-bit floats), with the same
names, component orders, frames and units. Every dataset is a table of rows, as the HDF5 tools
and Python show it: one row per node or element (in global id order) or per increment, one
column per value:

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Dataset
     - Shape
     - Content
   * - ``/results/nodes/<var>/<var>.step<N>``
     - (nodes, 3)
     - a node variable at step *N*
   * - ``/results/elts/<var>/<var>.step<N>``
     - (elements, *k*)
     - an element variable at step *N*: the *k* values of its text file's lines
   * - ``/results/elts/<var>/num_values``
     - (elements, 1)
     - slip-system variables: the number of values of each element (its phase's slip
       systems); the table has the largest number of columns, the values after an element's
       own number are NaN
   * - ``/results/forces/<faset>``
     - (rows, 7)
     - ``step incr Fx Fy Fz area time``, one row per increment
   * - ``/results/convergence``
     - (rows, 14)
     - ``incr iter nr r_norm rx_norm f_norm delu_norm delux_norm u_norm cg_iter attempt
       accepted bc_iter solve``: the columns of :file:`results/convergence`, with ``nr`` 1 for
       an NR iteration and 0 for SA, ``accepted`` 1 for ``A`` and 0 for ``R``, and every column
       always present (``attempt`` 1, ``accepted`` 1, ``bc_iter`` and ``solve`` 0 where the text
       file has no such column)
   * - ``/results/macro``
     - (rows, 27)
     - the columns of :file:`results/macro` (periodic conditions)
   * - ``/results/steps``
     - (rows, 2)
     - ``step time``: the steps written and their times [s]
   * - ``/mesh/coo``
     - (nodes, 3)
     - the initial coordinates [mm]
   * - ``/mesh/elts``
     - (elements, 10 or 20)
     - the nodes of each element (1-based ids, the mesh file's order)
   * - ``/mesh/elts_xdmf``
     - (elements, 10 or 20)
     - the same, 0-based, in the node order of XDMF and VTK (used by the XDMF file)
   * - ``/mesh/elset``, ``/mesh/phase``
     - (elements, 1)
     - the elset id (of the mesh file) and the phase of each element

Each table of rows has an attribute ``columns`` (the column names). The root group has the
attributes ``format`` (``fepx-hdf5 1``), ``num_nodes``, ``num_elts``, ``num_elsets``,
``num_ranks``, ``element_type`` (``tet10`` or ``hex20``), ``orientation`` (the orientation
descriptor and convention of ``ori``), ``inputs``, ``node_results``, ``elt_results`` and
``num_steps`` (the completed steps).

The XDMF file has one grid per step written, at the step's time: the elements, the step's
coordinates (``coo`` if printed, otherwise the initial ones), and as fields ``elset``,
``phase`` and every variable printed. Node variables are vectors; element variables are
scalars, vectors (``rotrate``, ``rotrate_spin``, ``rotrate_slip``) or tables of *k* components
in the order of the text files (for a symmetric tensor 11 22 33 23 13 12: ParaView numbers
them 0 to 5). The file is rewritten at every step, so the results of a run that stopped can be
opened.

Reading the file in Python, for example with h5py:

.. code-block:: python

  import h5py

  with h5py.File("simulation.sim/simulation.h5", "r") as h5:
      steps = h5["/results/steps"][:]                   # (step, time) rows
      stress = h5["/results/elts/stress/stress.step3"][:]  # (elements, 6)
      force = h5["/results/forces/z1"][:]               # step incr Fx Fy Fz area time
      print(h5.attrs["element_type"], stress[:, 2].mean(), force[-1, 4])

Restart Points
--------------

``print restart`` writes restart points (see :ref:`restart`), not a result.
