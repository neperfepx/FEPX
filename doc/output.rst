.. _simulation_results:

Simulation Results (:file:`simulation.cfg`)
==============================================

Simulation results include nodal results, elemental results, mesh face results, simulation restart information, and other, advanced output. The printing of the field variable files is defined as follows::

    print <result1> <result2> ...

where :data:`<result1>`, :data:`<result2>`, etc. are the particular result fields to be output.
Generally, and particularly for nodes and elements, outputs are provided per line (one node or element, etc. per line) and are space delimited within a line.  Note that, elemental quantities are internally calculated at each Gauss quadrature point within the element (15 total), but are only printed as one value per element, which is that associated with the quadrature point that falls at the element centroid.

The results are written in a :ref:`simulation_directory`, where each result is written in its own file, at each step (one file for each result and each step).  Results can have scalar values, vectorial values or tensorial values. In the result files, values for the different entities (nodes or elements) are written on successive lines, with components written on successive columns (space delimited). The components are written as follows:

    - scalar, :math:`a`: :math:`a`;
    - vector, :math:`v`: :math:`v_1`, :math:`v_2`, :math:`v_3`;
    - symmetrical tensor, :math:`t`: :math:`t_{11}`, :math:`t_{22}`, :math:`t_{33}`, :math:`t_{23}`, :math:`t_{31}`, :math:`t_{12}` (Voigt notation);
    - skew-symmetrical tensor, :math:`t`: :math:`t_{12}`, :math:`t_{13}`, :math:`t_{23}`;
    - non-symmetrical tensor, :math:`t`:  :math:`t_{11}`, :math:`t_{12}`, :math:`t_{13}`, :math:`t_{21}`, :math:`t_{22}`, :math:`t_{23}`, :math:`t_{31}`, :math:`t_{32}`, :math:`t_{33}`.

.. _nodal_results:

Nodal Results
-------------

Coordinates (:data:`coo`)
    Position vector, written in the sample basis.

Velocities (:data:`vel`)
    Velocity vector, written in the sample basis.

Displacements (:data:`disp`)
    Displacement vector, written in the sample basis.

.. _elemental_results:

Elemental Results
-----------------

Crystallographic Orientation (:data:`ori`)
  Orientation written using the same parameterization as that used in input, which can be

  - the Rodrigues vector, :math:`{\bf r} = {\bf t} \tan{(\omega / 2)}`;
  - the Euler-Bunge angles, :math:`(\varphi_1,\,\phi,\,\varphi_2)` expressed in degrees;
  - the Euler-Kocks angles, :math:`(\Psi,\,\Theta,\,\phi)` expressed in degrees;
  - the axis-angle pair, :math:`(t_1,\,t_2,\,t_3,\,\omega)`, where :math:`\bf{t}` is the (unit) axis of rotation and :math:`\omega` is the angle of rotation about said axis, expressed in degrees;
  - the quaternion, :math:`(q_0,\,q_1,\,q_2,\,q_3)`.

Slip System Resolved Shear Stress (:data:`rss`)
  :math:`\tau^\alpha` scalar values, printed in the order specified in :ref:`slip_systems`.

Slip System Critical Resolved Shear Stress(es) (:data:`crss`)
  :math:`g_0` scalar value(s), printed in the order specified in :ref:`slip_systems`.
  For isotropic hardening, only one value is printed.

.. _slip_system_shear_rate:

Slip System Shear Rate (:data:`sliprate`)
  :math:`\dot\gamma^\alpha` scalar values, printed in the order specified in :ref:`slip_systems`.

Slip System Slips (:data:`slip`)
  :math:`\gamma^\alpha` scalar values, printed in the order specified in :ref:`slip_systems`.

Stress Tensor (:data:`stress`)
  :math:`\bf \sigma` symmetrical tensor, written in the sample basis.

Equivalent Stress (:data:`stress_eq`)
  :math:`\sigma_{eq}` scalar value, calculated based on the deviatoric portion of the stress tensor, :math:`{\bf \sigma}^{\prime}`, via the tensor inner product: :math:`\sigma_{eq} = \sqrt{ {3 \over 2} {\bf \sigma}^{\prime} : {\bf \sigma}^{\prime}}`

Total Strain Tensor (:data:`strain`)
  :math:`\bf E` symmetrical tensor, written in the sample basis.

Equivalent Total Strain (:data:`strain_eq`)
  :math:`E_{eq}` scalar value, calculated based on the deviatoric portion of the total strain tensor, :math:`{\bf E}^{\prime}`, via the tensor inner product: :math:`E_{eq} = \sqrt{ {2 \over 3} {\bf E}^{\prime} : {\bf E}^{\prime}}`

Elastic Strain Tensor (:data:`strain_el`)
  :math:`\bf E^e` symmetrical tensor, written in the sample basis.

Equivalent Elastic Strain (:data:`strain_el_eq`)
  :math:`E^e_{eq}` scalar value, calculated based on the deviatoric portion of the elastic strain tensor, :math:`{\bf E}^{e \prime}`, via the tensor inner product: :math:`E^e_{eq} = \sqrt{ {2 \over 3} {\bf E}^{e \prime} : {\bf E}^{e \prime}}`

Plastic Strain Tensor (:data:`strain_pl`)
  :math:`\bf E^p` symmetrical tensor, written in the sample basis.

Equivalent Plastic Strain (:data:`strain_pl_eq`)
  :math:`E^p_{eq}` scalar value, calculated based on the deviatoric portion of the elastic strain tensor, :math:`{\bf E}^{p \prime}`, via the tensor inner product: :math:`E^p_{eq} = \sqrt{ {2 \over 3} {\bf E}^{p \prime} : {\bf E}^{p \prime}}`

The below variables are of secondary general interest:

Velocity Gradient Tensor (:data:`velgrad`)
  :math:`{\bf L}^p` general tensor, written in the sample basis.

Deformation Rate Tensor (:data:`defrate`)
  :math:`\bf D` symmetrical tensor, written in the sample basis.

Equivalent Deformation Rate (:data:`defrate_eq`)
  :math:`D_{eq}` scalar value, corresponding to the :math:`\bf D` tensor inner product: :math:`D = \sqrt{ {2 \over 3} {\bf D} : {\bf D} }`

Plastic Deformation Rate Tensor (:data:`defrate_pl`)
  :math:`{\bf D}^p` symmetrical tensor, written in the sample basis.

Equivalent Plastic Deformation Rate (:data:`defrate_pl_eq`)
  :math:`D^p_{eq}` scalar value, corresponding to the :math:`\bf D^p` tensor inner product: :math:`D = \sqrt{ {2 \over 3} {\bf D^p} : {\bf D^p} }`

Plastic Spin Rate Tensor (:data:`spinrate`)
  :math:`{\bf W}^p` skew-symmetric tensor, written in the sample basis.

Work (:data:`work`)
  :math:`W` Scalar value, calculated as the time integration of the tensor inner product of the deformation rate tensor and the Cauchy stress tensor: :math:`W = \int{  (\sigma : {\bf D}) }\Delta t`.

Plastic Work (:data:`work_pl`)
  :math:`W^p` Scalar value, calculated as the time integration of the tensor inner product of the plastic deformation rate tensor and the deviatoric portion of the Cauchy stress tensor: :math:`W^p = \int{ ( \sigma ^ \prime :  {\bf D }^ p ) \Delta t}`.

Work Rate (:data:`workrate`)
  :math:`\dot{W}` scalar value, calculated as the tensor inner product of the deformation rate tensor and the Cauchy stress tensor:
  :math:`\dot{W} = \sigma : {\bf D}`.

Plastic Work Rate (:data:`workrate_pl`)
  :math:`\dot{W}^p` scalar value, calculated as the tensor inner product of the plastic deformation rate tensor and the deviatoric portion of the Cauchy stress tensor:
  :math:`\dot{W}^p = \sigma ^ \prime : {\bf D }^ p`

Lattice Reorientation Rate Vector (:data:`rotrate`)
  (reorientation axis x reorientation velocity) vector, written in the sample basis.

Spin Part of the Lattice Reorientation Rate Vector (:data:`rotrate_spin`)
  (reorientation axis x reorientation velocity) vector, spin part, written in the sample basis.

Slip Part of the Lattice Reorientation Rate Vector (:data:`rotrate_slip`)
  (reorientation axis x reorientation velocity) vector, slip part, written in the sample basis.

.. _surface_results:

Surface Results
---------------

Integrated Forces (:data:`force`)
  Integrated forces (or "loads") for all surfaces, provided at each increment and written as :file:`post.force.*` files, which are constructed via the defined :data:`faset_label` strings in the :file:`simulation.msh` file (for a cubic domain, :data:`x0`, :data:`x1`, :data:`y0`, :data:`y1`, :data:`z0` and :data:`z1`, where :data:`0` and :data:`1` represent the minimum and maximum coordinates).  These files are generally tabulated with the given structure::

    <step> <incr> <force_x> <force_y> <force_z> <surf_area> <current_time>

  where :data:`<step>` is the prescribed load step, :data:`<incr>` is the total increment value being printed, :data:`<force_x>` is the surface-integrated force in the :data:`x` direction, :data:`<force_y>` is the surface-integrated force in the :data:`y` direction, :data:`<force_z>` is the surface-integrated force in the :data:`z` direction, :data:`<surf_area>` is the current surface area of the given face, and :data:`<current_time>` is the total simulated time at the time of printing.

  If multiaxial loading is utilized, an additional :data:`<length>` column will be appended to the right of :data:`<current_time`. The :data:`<length>` column contains the maximal coordinate values of the domain and these values are stored in their associated face files. For example, the maximal mesh coordinate value in the :data:`x` direction is stored in the :file:`post.force.x0` and :file:`post.force.x1` files accordingly.

.. _advanced_results:

Advanced Results
----------------

Convergence Information (:data:`convergence`)
  Convergence statistics, provided at each increment and written as a :file:`post.conv` file, which is tabulated with the given structure::

    <incr> <iter> <nr> <r_norm> <rx_norm> <f_norm> <delu_norm> <delux_norm> <u_norm> <cg_iter>

  where :data:`<incr>` is the total increment value being printed, :data:`<iter>` is a sub-increment iteration, :data:`<nr>` is a boolean that notifies if the given iteration was a Newton-Raphson iteration, :data:`<r_norm>` is a residual norm, :data:`<rx_norm>` is the maximum absolute value of :data:`<r_norm>`, :data:`<f_norm>` is a force norm, :data:`<delu_norm>` is the change in velocity norm, :data:`<delux_norm>` is the maximum absolute value of :data:`<delu_norm>`, :data:`<u_norm>` is the velocity norm, and :data:`<cg_iter>` is the number of iterations the conjugate gradient solver performed. All norms are computed as 2-norm or the square root of the inner product of a vector by itself.

Restart Information (:data:`restart`)
  Restart files, which allows to restart a simulation at the last completed step.  Restart files are only written at the end of a successful step and not at individual increments.

