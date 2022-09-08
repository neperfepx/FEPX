.. _simulation_output:

Simulation Output
=================

In this section the output files are described. In general, the output can be broken down into four types: nodal output (variables that are calculated and printed at the finite element nodes), elemental output (variables that are calculated and printed at the finite elements), restart output (which contains all variables necessary to restart a simulation at a specific load step), and miscellaneous output (macroscopic forces and various simulation statistics).

All desired output must be defined in the :file:`simulation.config` file (see: :ref:`Printing Results <printing_results>` for a description of the print options). A (small) :file:`post.report` file is systematically printed, which contains information necessary for post-processing with Neper. All output described in this section is *raw* simulation output and can be post-processed into a more human-readable format with Neper's :data:`-S` module. Refer to the `Neper reference manual <https://neper.info/docs/neper.pdf>`_, for a more complete description of the :data:`-S` module.

Nodal Output
------------

Nodal output prints a single variable per finite element node. Raw output is printed on a per-core basis, and the general file name structure is :file:`post.variable.core<#>`, where :data:`variable` is the variable being printed, and :data:`<#>` denotes the ID of the core on which the data is being printed. In general, the file structure is:

::

    <step_number> <core_start_dof> <core_end_dof>
    <node_N> <values>
    <node_N+1> <values>
    <node_N+2> <values>
    ...
    <step_number> <core_start_dof> <core_end_dof>
    <node_N> <values>
    <node_N+1> <values>
    <node_N+2> <values>
    ...

Here, a header line prints for each deformation step, which details the deformation step number (:data:`<step_number>`), the initial degree of freedom that prints from that core (:data:`<core_start_dof>`), and the final degree of freedom that prints from that core (:data:`<core_end_dof>`). Note that the number of nodes per-core is one third of the number of degrees of freedoms per-core. For each load step, multiple values associated with a variable are printed per line (:data:`<values>`), preceded by the corresponding node number (:data:`<node_*>`). When multiple values are printed per line, values are space delimited.

Below, the specific nodal values available for printing are explained in detail.

.. option:: Coordinates

    In the :file:`post.coo.core*` files, the coordinates are printed. Each nodal coordinate is described in the orthonormal Cartesian sample basis, and one coordinate is printed per line (3 values per line). The components of the coordinates are printed in the order:

    :math:`x, y, z`

.. option:: Displacements

    In the :file:`post.disp.core*` files, the displacements are printed. Each nodal displacement is described in the orthonormal Cartesian sample basis, and one displacement is printed per line (3 values per line). The components of the displacements are printed in the order:

    :math:`d_x, d_y, d_z`

.. option:: Velocities

    In the :file:`post.vel.core*` files, the velocities are printed. Each nodal velocity is described in the orthonormal Cartesian sample basis, and one velocity is printed per line (3 values per line). The components of the velocities are printed in the order:

    :math:`v_x, v_y, v_z`

Elemental Output
----------------

Elemental output prints a single variable per finite element. Raw output is printed on a per-core basis, and the general file name structure is :file:`post.variable.core<#>`, where :data:`variable` is the variable being printed, and :data:`<#>` denotes the ID of the core on which the data is being printed. In general, the file structure is:

::

    <step_number> <core_start_elt> <core_end_elt>
    <element_N> <value(s)>
    <element_N+1> <value(s)>
    <element_N+2> <value(s)>
    ...
    <step_number> <core_start_elt> <core_end_elt>
    <element_N> <value(s)>
    <element_N+1> <value(s)>
    <element_N+2> <value(s)>
    ...

Here, a header line prints for each deformation step, which details the deformation step number (:data:`<step_number>`), the initial element that prints from that core (:data:`<core_start_elt>`), and the final element that prints from that core (:data:`<core_end_elt>`). For each load step, either a single or multiple value (:data:`<value(s)>`) associated with a variable is printed per line (for variables that are printed singularly per element, such as scalars, or for variables that print multiple values per element, such as tensors, respectively), preceded by the corresponding element number (:data:`<element_*>`). When multiple values are printed per line, values are space delimited.

FEPX calculates elemental quantities at each Gauss quadrature point within the element (15 total). However, only one value is printed -- that associated with the quadrature point that falls at the element centroid.

Below, the specific nodal values available for printing are explained in detail.

.. option:: Critical Resolved Shear Stress

    In the :file:`post.crss.core*` files, the critical resolved shear stress is printed. For the isotropic hardening assumption (:ref:`Optional Input Parameters <optional_input_parameters>`), one value is printed per element. For anisotropic hardening assumptions, the critical resolved shear stress is printed for each slip system per element, one line of values per element.

    For body centered cubic crystal symmetry, values are printed in the order:

    :math:`(0 1 \bar 1)[1 1 1],\,(1 0 \bar 1)[1 1 1],\,(1 \bar 1 0)[1 1 1],\,(0 1 1)[1 1 \bar 1],\,(1 0 1)[1 1 \bar 1],\,(1 \bar 1 0)[1 1 \bar 1],`

    :math:`(0 1 1)[1 \bar 1 1],\,(1 0 \bar 1)[1 \bar 1 1],\,(1 1 0)[1 \bar 1 1],\,(0 1 \bar 1)[1 \bar 1 \bar 1],\,(1 0 1)[1 \bar 1 \bar 1],\,(1 1 0)[1 \bar 1 \bar 1]`.

    For face centered cubic crystal symmetry, values are printed in the order:

    :math:`(1 1 1)[0 1 \bar 1],\,(1 1 1)[1 0 \bar 1],\,(1 1 1)[1 \bar 1 0],\,(1 1 \bar 1)[0 1 1],\,(1 1 \bar 1)[1 0 1],\,(1 1 \bar 1)[1 \bar 1 0],\,`

    :math:`(1 \bar 1 1)[0 1 1],\,(1 \bar 1 1)[1 0 \bar 1],\,(1 \bar 1 1)[1 1 0],\,(1 \bar 1 \bar 1)[0 1 \bar 1],\,(1 \bar 1 \bar 1)[1 0 1],\,(1 \bar 1 \bar 1)[1 1 0]`.

    For hexagonal close packed crystal symmetry, values are printed in the order (corresponding to the 3 basal, 3 prismatic, and 12 pyramidal slip systems):

    :math:`(0 0 0 1)[2 \bar 1 \bar 1 0],\,(0 0 0 1)[\bar 1 2 \bar 1 0],\,(0 0 0 1)[\bar 1 \bar 1 2 0],\,(0 1 \bar 1 0)[2 \bar 1 \bar 1 0],\,(\bar 1 0 1 0)[\bar 1 2 \bar 1 0],\,(1 \bar 1 0 0)[\bar 1 \bar 1 2 0],\,`

    :math:`(1 0 \bar 1 1)[\bar 2 1 1 3],\,(1 0 \bar 1 1)[\bar 1 \bar 1 2 3],\,(0 1 \bar 1 1)[\bar 1 \bar 1 2 3],\,(0 1 \bar 1 1)[1 \bar 2 1 3],\,(\bar 1 1 0 1)[1 \bar 2 1 3],\,(\bar 1 1 0 1)[2 \bar 1 \bar 1 3],\,`

    :math:`(\bar 1 0 1 1)[2 \bar 1 \bar 1 3],\,(\bar 1 0 1 1)[1 1 \bar 2 3],\,(0 \bar 1 1 1)[1 1 \bar 2 3],\,(0 \bar 1 1 1)[\bar 1 2 \bar 1 3],\,(1 \bar 1 0 1)[\bar 1 2 \bar 1 3],\,(1 \bar 1 0 1)[\bar 2 1 1 3]`.

    For body centered tetragonal crystal symmetry, values are printed in the order:

    :math:`(1 0 0)[0 0 1],\,(0 1 0)[0 0 1],\,(1 1 0)[0 0 1],\,(1 \bar 1 0)[0 0 1],\,(1 0 0)[0 1 0],\,(0 1 0)[1 0 0],`

    :math:`(1 1 0)[1 \bar 1 1],\,(1 1 0)[\bar 1 1 1],\,(1 \bar 1 0)[1 1 1],\,(1 \bar 1 0)[\bar 1 \bar 1 1],\,(1 1 0)[\bar 1 1 0],\,(1 \bar 1 0)[1 1 0],`

    :math:`(1 0 0)[0 1 1],\,(1 0 0)[0 1 \bar 1],\,(0 1 0)[1 0 1],\,(0 1 0)[1 0 \bar 1],\,(0 0 1)[1 1 0],\,(0 0 1)[1 \bar 1 0],`

    :math:`(0 0 1)[1 1 0],\,(0 0 1)[1 \bar 1 0],\,(1 0 1)[1 0 \bar 1],\,(1 0 \bar 1)[1 0 1],\,(0 1 1)[0 1 \bar 1],\,(0 1 \bar 1)[0 1 1],`

    :math:`(1 2 1)[\bar 1 0 1],\,(\bar 1 2 1)[1 0 1],\,(\bar 1 \bar 2 1)[1 0 1],\,(1 \bar 2 1)[\bar 1 0 1],\,(2 1 1)[0 \bar 1 1],\,(\bar 2 1 1)[0 \bar 1 1],`

    :math:`(\bar 2 \bar 1 1)[0 1 1],\,(2 \bar 1 1)[0 1 1].`

.. option:: Deformation Rate Tensor

    In the :file:`post.defrate.core*` files, the deformation rate tensor is printed. Each tensor, :math:`\bf D`, is printed in the sample basis. The independent components are printed, one tensor per line (6 values per line). The components, :math:`D _{ij}`, are printed in the order:

    :math:`D_{11},\,D_{22},\,D_{33},\,D_{23},\,D_{13},\,D_{12}`

.. option:: Equivalent Deformation Rate

    In the :file:`post.defrate-eq.core*` files, the equivalent deformation rate is printed. One scalar value is printed per element. The equivalent deformation rate, :math:`D`, is calculated based on the deformation rate tensor, :math:`{\bf D}`, via the tensor inner product:

    :math:`D = \sqrt{ {2 \over 3} {\bf D} : {\bf D} }`

.. option:: Plastic Deformation Rate Tensor

    In the :file:`post.defrate-pl.core*` files, the deviatoric plastic deformation rate tensor is printed. Each tensor, :math:`{\bf D}^p`, is printed in the sample basis. The independent components are printed, one tensor per line (6 values per line). The components, :math:`D^p_{ij}`, are printed in the order:

    :math:`D^p_{11},\,D^p_{22},\,D^p_{33},\,D^p_{23},\,D^p_{13},\,D^p_{12}`

.. option:: Equivalent Plastic Deformation Rate

    In the :file:`post.defrate-pl-eq.core*` files, the equivalent plastic deformation rate is printed. One scalar value is printed per element. The equivalent plastic deformation rate, :math:`D^p`, is calculated based on the plastic deformation rate tensor, :math:`{\bf D}^p`, via the tensor inner product:

    :math:`D^p = \sqrt{ {2 \over 3} {\bf D}^p : {\bf D}^p }`

.. option:: Elemental Volume

    In the :file:`post.elt-vol.core*` files, the elemental volume is printed. One scalar value is printed per element. The elemental volume is calculated as the Gaussian integration of the determinant of the Jacobian matrix:

    :math:`V_{el} = \sum_{i=1}^{n_{qp}}{det(J_i) w_i }`

.. option:: Crystallographic Orientation

    In the :file:`post.ori.core*` files, the crystallographic orientation is printed. Depending on the orientation parameterization used as input, the orientation values may range from 3 values per element (when using Rodrigues vector, Euler-Bunge angles and Euler-Kocks angles parameterizations) or 4 values per element (when using axis-angle or quaternion parameterizations). One orientation is printed per line (3 or 4 values per line).

    For Rodrigues: :math:`r_1,\,r_2,\,r_3`, where the Rodrigues vector is :math:`{\bf r} = {\bf t} \tan{(\omega / 2)}`.

    For Euler-Bunge: :math:`\phi_1,\,\theta,\,\phi_2` (where :math:`\phi_1` is the rotation about the :math:`z` axis, :math:`\theta` is the rotation about the :math:`x^{\prime}` axis, and :math:`\phi_2` is the rotation about the :math:`z^{\prime \prime}` axis, all in degrees).

    For Euler-Kocks: :math:`\Psi,\,\Theta,\,\phi` (where :math:`\Psi` is the rotation about the :math:`z` axis, :math:`\Theta` is the rotation about the :math:`y^{\prime}` axis, and :math:`\phi` is the rotation about the :math:`z^{\prime \prime}` axis, all in degrees).

    For axis-angle: :math:`t_1,\,t_2,\,t_3,\,\omega` (where :math:`\bf{t}` is the normalized axis of rotation and :math:`\omega` is the angle of rotation about said axis, in degrees).

    For quaternion: :math:`q_0,\,q_1,\,q_2,\,q_3`, where :math:`q_0 = \cos{(\omega / 2)}` and :math:`q_i = t_i \sin{(\omega / 2)}` for :math:`i = 1,\,2,\,3`.

.. option:: Slip System Shear

    In the :file:`post.slip.core*` files, the accumulated slip system shear is printed. The slip system shear is printed for each slip system per element, one line of values per element.

    For body centered cubic crystal symmetry, values are printed in the order:

    :math:`(0 1 \bar 1)[1 1 1],\,(1 0 \bar 1)[1 1 1],\,(1 \bar 1 0)[1 1 1],\,(0 1 1)[1 1 \bar 1],\,(1 0 1)[1 1 \bar 1],\,(1 \bar 1 0)[1 1 \bar 1],`
    :math:`(0 1 1)[1 \bar 1 1],\,(1 0 \bar 1)[1 \bar 1 1],\,(1 1 0)[1 \bar 1 1],\,(0 1 \bar 1)[1 \bar 1 \bar 1],\,(1 0 1)[1 \bar 1 \bar 1],\,(1 1 0)[1 \bar 1 \bar 1]`.

    For face centered cubic crystal symmetry, values are printed in the order:

    :math:`(1 1 1)[0 1 \bar 1],\,(1 1 1)[1 0 \bar 1],\,(1 1 1)[1 \bar 1 0],\,(1 1 \bar 1)[0 1 1],\,(1 1 \bar 1)[1 0 1],\,(1 1 \bar 1)[1 \bar 1 0],\,`
    :math:`(1 \bar 1 1)[0 1 1],\,(1 \bar 1 1)[1 0 \bar 1],\,(1 \bar 1 1)[1 1 0],\,(1 \bar 1 \bar 1)[0 1 \bar 1],\,(1 \bar 1 \bar 1)[1 0 1],\,(1 \bar 1 \bar 1)[1 1 0]`.

    For hexagonal close packed crystal symmetry, values are printed in the order (corresponding to the 3 basal, 3 prismatic, and 12 pyramidal slip systems):

    :math:`(0 0 0 1)[2 \bar 1 \bar 1 0],\,(0 0 0 1)[\bar 1 2 \bar 1 0],\,(0 0 0 1)[\bar 1 \bar 1 2 0],\,(0 1 \bar 1 0)[2 \bar 1 \bar 1 0],\,(\bar 1 0 1 0)[\bar 1 2 \bar 1 0],\,(1 \bar 1 0 0)[\bar 1 \bar 1 2 0],\,`
    :math:`(1 0 \bar 1 1)[\bar 2 1 1 3],\,(1 0 \bar 1 1)[\bar 1 \bar 1 2 3],\,(0 1 \bar 1 1)[\bar 1 \bar 1 2 3],\,(0 1 \bar 1 1)[1 \bar 2 1 3],\,(\bar 1 1 0 1)[1 \bar 2 1 3],\,(\bar 1 1 0 1)[2 \bar 1 \bar 1 3],\,`
    :math:`(\bar 1 0 1 1)[2 \bar 1 \bar 1 3],\,(\bar 1 0 1 1)[1 1 \bar 2 3],\,(0 \bar 1 1 1)[1 1 \bar 2 3],\,(0 \bar 1 1 1)[\bar 1 2 \bar 1 3],\,(1 \bar 1 0 1)[\bar 1 2 \bar 1 3],\,(1 \bar 1 0 1)[\bar 2 1 1 3]`.

    For body centered tetragonal crystal symmetry, values are printed in the order:

    :math:`(1 0 0)[0 0 1],\,(0 1 0)[0 0 1],\,(1 1 0)[0 0 1],\,(1 \bar 1 0)[0 0 1],\,(1 0 0)[0 1 0],\,(0 1 0)[1 0 0],`
    :math:`(1 1 0)[1 \bar 1 1],\,(1 1 0)[\bar 1 1 1],\,(1 \bar 1 0)[1 1 1],\,(1 \bar 1 0)[\bar 1 \bar 1 1],\,(1 1 0)[\bar 1 1 0],\,(1 \bar 1 0)[1 1 0],`
    :math:`(1 0 0)[0 1 1],\,(1 0 0)[0 1 \bar 1],\,(0 1 0)[1 0 1],\,(0 1 0)[1 0 \bar 1],\,(0 0 1)[1 1 0],\,(0 0 1)[1 \bar 1 0],`
    :math:`(0 0 1)[1 1 0],\,(0 0 1)[1 \bar 1 0],\,(1 0 1)[1 0 \bar 1],\,(1 0 \bar 1)[1 0 1],\,(0 1 1)[0 1 \bar 1],\,(0 1 \bar 1)[0 1 1],`
    :math:`(1 2 1)[\bar 1 0 1],\,(\bar 1 2 1)[1 0 1],\,(\bar 1 \bar 2 1)[1 0 1],\,(1 \bar 2 1)[\bar 1 0 1],\,(2 1 1)[0 \bar 1 1],\,(\bar 2 1 1)[0 \bar 1 1],`
    :math:`(\bar 2 \bar 1 1)[0 1 1],\,(2 \bar 1 1)[0 1 1].`

.. option:: Slip System Shear Rate

    In the :file:`post.sliprate.core*` files, the slip system shear rate is printed. The slip system shear rate is printed for each slip system per element, one line of values per element.

    For body centered cubic crystal symmetry, values are printed in the order:

    :math:`(0 1 \bar 1)[1 1 1],\,(1 0 \bar 1)[1 1 1],\,(1 \bar 1 0)[1 1 1],\,(0 1 1)[1 1 \bar 1],\,(1 0 1)[1 1 \bar 1],\,(1 \bar 1 0)[1 1 \bar 1],`
    :math:`(0 1 1)[1 \bar 1 1],\,(1 0 \bar 1)[1 \bar 1 1],\,(1 1 0)[1 \bar 1 1],\,(0 1 \bar 1)[1 \bar 1 \bar 1],\,(1 0 1)[1 \bar 1 \bar 1],\,(1 1 0)[1 \bar 1 \bar 1]`.

    For face centered cubic crystal symmetry, values are printed in the order:

    :math:`(1 1 1)[0 1 \bar 1],\,(1 1 1)[1 0 \bar 1],\,(1 1 1)[1 \bar 1 0],\,(1 1 \bar 1)[0 1 1],\,(1 1 \bar 1)[1 0 1],\,(1 1 \bar 1)[1 \bar 1 0],\,`
    :math:`(1 \bar 1 1)[0 1 1],\,(1 \bar 1 1)[1 0 \bar 1],\,(1 \bar 1 1)[1 1 0],\,(1 \bar 1 \bar 1)[0 1 \bar 1],\,(1 \bar 1 \bar 1)[1 0 1],\,(1 \bar 1 \bar 1)[1 1 0]`.

    For hexagonal close packed crystal symmetry, values are printed in the order (corresponding to the 3 basal, 3 prismatic, and 12 pyramidal slip systems):

    :math:`(0 0 0 1)[2 \bar 1 \bar 1 0],\,(0 0 0 1)[\bar 1 2 \bar 1 0],\,(0 0 0 1)[\bar 1 \bar 1 2 0],\,(0 1 \bar 1 0)[2 \bar 1 \bar 1 0],\,(\bar 1 0 1 0)[\bar 1 2 \bar 1 0],\,(1 \bar 1 0 0)[\bar 1 \bar 1 2 0],\,`
    :math:`(1 0 \bar 1 1)[\bar 2 1 1 3],\,(1 0 \bar 1 1)[\bar 1 \bar 1 2 3],\,(0 1 \bar 1 1)[\bar 1 \bar 1 2 3],\,(0 1 \bar 1 1)[1 \bar 2 1 3],\,(\bar 1 1 0 1)[1 \bar 2 1 3],\,(\bar 1 1 0 1)[2 \bar 1 \bar 1 3],\,`
    :math:`(\bar 1 0 1 1)[2 \bar 1 \bar 1 3],\,(\bar 1 0 1 1)[1 1 \bar 2 3],\,(0 \bar 1 1 1)[1 1 \bar 2 3],\,(0 \bar 1 1 1)[\bar 1 2 \bar 1 3],\,(1 \bar 1 0 1)[\bar 1 2 \bar 1 3],\,(1 \bar 1 0 1)[\bar 2 1 1 3]`.

    For body centered tetragonal crystal symmetry, values are printed in the order:

    :math:`(1 0 0)[0 0 1],\,(0 1 0)[0 0 1],\,(1 1 0)[0 0 1],\,(1 \bar 1 0)[0 0 1],\,(1 0 0)[0 1 0],\,(0 1 0)[1 0 0],`
    :math:`(1 1 0)[1 \bar 1 1],\,(1 1 0)[\bar 1 1 1],\,(1 \bar 1 0)[1 1 1],\,(1 \bar 1 0)[\bar 1 \bar 1 1],\,(1 1 0)[\bar 1 1 0],\,(1 \bar 1 0)[1 1 0],`
    :math:`(1 0 0)[0 1 1],\,(1 0 0)[0 1 \bar 1],\,(0 1 0)[1 0 1],\,(0 1 0)[1 0 \bar 1],\,(0 0 1)[1 1 0],\,(0 0 1)[1 \bar 1 0],`
    :math:`(0 0 1)[1 1 0],\,(0 0 1)[1 \bar 1 0],\,(1 0 1)[1 0 \bar 1],\,(1 0 \bar 1)[1 0 1],\,(0 1 1)[0 1 \bar 1],\,(0 1 \bar 1)[0 1 1],`
    :math:`(1 2 1)[\bar 1 0 1],\,(\bar 1 2 1)[1 0 1],\,(\bar 1 \bar 2 1)[1 0 1],\,(1 \bar 2 1)[\bar 1 0 1],\,(2 1 1)[0 \bar 1 1],\,(\bar 2 1 1)[0 \bar 1 1],`
    :math:`(\bar 2 \bar 1 1)[0 1 1],\,(2 \bar 1 1)[0 1 1].`

.. option:: Plastic Spin Rate Tensor

    In the :file:`post.spinrate.core*` files, the skew-symmetric plastic spin rate tensor is printed. Each tensor, :math:`{\bf W}^p`, is printed in the sample basis. The independent components are printed, one tensor per line (3 values per line). The components, :math:`W^p_{ij}`, are printed in the order:

    :math:`W^p_{12}, W^p_{13}, W^p_{23}`

.. option:: Total Strain Tensor

    In the :file:`post.strain.core*` files, the total strain tensor is printed. Each tensor, :math:`\bf E`, is printed in the sample basis. The independent components are printed, one tensor per line (6 values per line). The components, :math:`E_{ij}`, are printed in the order:

    :math:`E_{11}, E_{22}, E_{33}, E_{23}, E_{13}, E_{12}`

.. option:: Equivalent Total Strain

    In the :file:`post.strain-eq.core*` files, the equivalent total strain is printed. One scalar value is printed per element. The equivalent total strain, :math:`E`, is calculated based on the deviatoric portion of the total strain tensor, :math:`{\bf E}^\prime`. via the tensor inner product:

    :math:`E = \sqrt{ {2 \over 3} {\bf E}^\prime : {\bf E}^\prime}`

.. option:: Elastic Strain Tensor

    In the :file:`post.strain-el.core*` files, the elastic strain tensor is printed. Each tensor, :math:`\bf {E}^e`, is printed in the sample basis. The independent components are printed, one tensor per line (6 values per line). The components, :math:`E^e_{ij}`, are printed in the order:

    :math:`E^e_{11}, E^e_{22}, E^e_{33}, E^e_{23}, E^e_{13}, E^e_{12}`

.. option:: Equivalent Elastic Strain

    In the :file:`post.strain-el-eq.core*` files, the equivalent elastic strain is printed. One scalar value is printed per element. The equivalent elastic strain, :math:`E^{e}`, is calculated based on the deviatoric portion of the elastic strain tensor, :math:`{\bf E}^{e \prime}`, via the tensor inner product:

    :math:`E^e = \sqrt{ {2 \over 3} {\bf E}^{e \prime} : {\bf E}^{e \prime}}`

.. option:: Plastic Strain Tensor

    In the :file:`post.strain-pl.core*` files, the plastic strain tensor is printed. Each tensor, :math:`\bf {E}^p`, is printed in the sample basis. The independent components are printed, one tensor per line (6 values per line). The components, :math:`E^p_{ij}`, are printed in the order:

    :math:`E^p_{11}, E^p_{22}, E^p_{33}, E^p_{23}, E^p_{13}, E^p_{12}`

.. option:: Equivalent Plastic Strain

    In the :file:`post.strain-pl-eq.core*` files, the equivalent plastic strain is printed. One scalar value is printed per element. The equivalent plastic strain, :math:`E^p`, is calculated based on the plastic strain tensor, :math:`{\bf E}^p`, via the tensor inner product:

    :math:`E^p = \sqrt{ {2 \over 3} {\bf E}^p : {\bf E}^p }`

.. option:: Stress Tensor

    In the :file:`post.stress.core*` files, the symmetric stress tensor is printed. Each tensor, :math:`\bf \sigma`, is printed in the sample basis. The independent components are printed, one tensor per line (6 values per line). The components, :math:`\sigma_{ij}`, are printed in the order:

    :math:`\sigma_{11}, \sigma_{22}, \sigma_{33}, \sigma_{23}, \sigma_{13}, \sigma_{12}`

.. option:: Equivalent Stress

    In the :file:`post.stress-eq.core*` files, the equivalent stress is printed. One scalar value is printed per element. The equivalent stress, :math:`\sigma`, is calculated based on the deviatoric portion of the stress tensor, :math:`{\bf \sigma}^{\prime}`, via the tensor inner product:

    :math:`\sigma = \sqrt{ {3 \over 2} {\bf \sigma}^{\prime} : {\bf \sigma}^{\prime}}`

.. option:: Velocity Gradient Tensor

    In the :file:`post.velgrad.core*` files, the velocity gradient tensor is printed. Each tensor, :math:`\bf L`, is printed in the sample basis. One tensor is printed per line (9 values per line). The components, :math:`L_{ij}`,  are printed in the order:

    :math:`L_{11}, L_{12}, L_{13}, L_{21}, L_{22}, L_{23}, L_{31}, L_{32}, L_{33}`

.. option:: Work

    In the :file:`post.work.core*` files, the work is printed. One scalar value is printed per element. The work is calculated as the time integration of the tensor inner product of the deformation rate tensor and the Cauchy stress tensor:

    :math:`W = \int{  (\sigma : {\bf D}) }\Delta t`

.. option:: Plastic Work

    In the :file:`post.work-pl.core*` files, the plastic work is printed. One scalar value is printed per element. The plastic work is calculated as the time integration of the tensor inner product of the plastic deformation rate tensor and the deviatoric portion of the Cauchy stress tensor:

    :math:`W^p = \int{ ( \sigma ^ \prime :  {\bf D }^ p ) \Delta t}`

.. option:: Work Rate

    In the :file:`post.workrate.core*` files, the work rate is printed. One scalar value is printed per element. The work rate is calculated as the tensor inner product of the deformation rate tensor and the Cauchy stress tensor:

    :math:`\dot{W} = \sigma : {\bf D}`

.. option:: Plastic Work Rate

    In the :file:`post.workrate-pl.core*` files, the plastic work rate is printed. One scalar value is printed per element. The plastic work rate is calculated as the tensor inner product of the plastic deformation rate tensor and the deviatoric portion of the Cauchy stress tensor:

    :math:`\dot{W}^p = \sigma ^ \prime : {\bf D }^ p`


.. _restart_output:

Restart Output
--------------

If the :data:`print restart` command is present in the :file:`simulation.config` file, a set of additional restart files will be generated from the simulation. These files are written at the end of each prescribed step and contain necessary information to restart a given simulation (:ref:`Restarting a Simulation <sim_restart>` for information on how to restart a simulation). Two types of restart files are generated, a control file, :file:`rst<N>.control`, and per-core field files, :file:`rst<N>.field.core*` (where :data:`<N>` indicates which simulation the files describe, 0 indexing). Both file types are unformatted (or binary) files and are generally unmodifiable. The structures of the data stored within both files for the various deformation modes follow.

Uniaxial Restart Control
~~~~~~~~~~~~~~~~~~~~~~~~

The :file:`rst<N>.control` file for uniaxial loading modes contains the following data in the given order:

::

    current_step <step>
    previous_load_array <load_x> <load_y> <load_z>
    step_complete_flag <logical>
    previous_timestep_value <time>
    current_incr <increment>
    current_time <time>
    surface_1_loads <load_x> <load_y> <load_z>
    ...
    surface_6_loads <load_x> <load_y> <load_z>
    previous_prescribed_load <load>
    current_surface_areas <area_surf_1> ... <area_surf_6>
    initial_surface_areas <area_surf_1> ... <area_surf_6>

Multiaxial CSR Restart Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :file:`rst<N>.control` file for multiaxial constant strain rate loading modes contains the following data in the given order:

::

    current_step <step>
    current_load_array <load_x> <load_y> <load_z>
    previous_load_array <load_x> <load_y> <load_z>
    step_complete_flag <logical>
    previous_timestep_value <time>
    current_incr <increment>
    current_time <time>
    surface_1_loads <load_x> <load_y> <load_z>
    ...
    surface_6_loads <load_x> <load_y> <load_z>
    current_surface_areas <area_surf_1> ... <area_surf_6>
    initial_surface_areas <area_surf_1> ... <area_surf_6>
    current_mesh_lengths <length_x> <length_y> <length_z>
    initial_mesh_lengths <length_x> <length_y> <length_z>
    current_control_velocity <vel_x> <vel_y> <vel_z>
    s_pert_mag <vel>
    t_pert_mag <vel>

Multiaxial CLR Restart Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :file:`rst<N>.control` file for multiaxial constant load rate loading modes contains the following data in the given order:

::

    current_step <step>
    current_load_array <load_x> <load_y> <load_z>
    previous_load_array <load_x> <load_y> <load_z>
    first_incr_in_step <logical>
    current_incr <increment>
    current_time <time>
    surface_1_loads <load_x> <load_y> <load_z>
    ...
    surface_6_loads <load_x> <load_y> <load_z>
    current_surface_areas <area_surf_1> ... <area_surf_6>
    initial_surface_areas <area_surf_1> ... <area_surf_6>
    current_mesh_lengths <length_x> <length_y> <length_z>
    initial_mesh_lengths <length_x> <length_y> <length_z>
    current_control_velocity <vel_x> <vel_y> <vel_z>
    previous_control_action <integer>
    current_control_action <integer>
    initial_load_dwell_velocity <vel_x> <vel_y> <vel_z>
    initial_unload_dwell_velocity <vel_x> <vel_y> <vel_z>

Restart Field Data
~~~~~~~~~~~~~~~~~~

All loading modes also write field data on a per-core basis to :file:`rst<N>.field.core*` files. These files contain the necessary field variable information in order to spatially define the total state of the virtual sample at the time of printing. The following field data arrays are written to the files in the given order:

::

    coords <coords>
    velocity <velocity>

    c0_angs <orientation>
    c_angs <orientation>
    rstar <rotation>
    rstar_n <rotation>
    wts <weight>
    crss <crss>
    crss_n <crss>

    gela_kk_bar <strain>
    gsig_vec_n <stress>
    pela_kk_bar <strain>
    psig_vec_n <stress>
    e_elas_kk_bar <strain>
    sig_vec_n <stress>

    eqstrain <strain>
    eqplstrain <strain>
    gamma <shear>

    el_work_n <work>
    el_workp_n <work>
    el_work_rate_n <work_rate>
    el_workp_rate_n <work_rate>

Miscellaneous Output
--------------------

In addition to nodal and elemental variable printing, miscellaneous output is available for printing and include simulation convergence data, surface-integrated forcing data, and a simulation report file. The optional input commands and output file formats are described in this section.

Convergence Statistics Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the :data:`print convergence` command is present in the :file:`simulation.config` file, various convergence statistics for the performed simulation will be output with convergence values provided at each increment. This :file:`post.conv` file is tabulated with the given structure:

::

    <incr> <iter> <nr> <r_norm> <rx_norm> <f_norm> <delu_norm> <delux_norm> <u_norm> <cg_iter>

where :data:`<incr>` is the total increment value being printed, :data:`<iter>` is a sub-increment iteration, :data:`<nr>` is a boolean that notifies if the given iteration was a Newton-Raphson iteration, :data:`<r_norm>` is a residual norm, :data:`<rx_norm>` is the maximum absolute value of :data:`<r_norm>`, :data:`<f_norm>` is a force norm, :data:`<delu_norm>` is the change in velocity norm, :data:`<delux_norm>` is the maximum absolute value of :data:`<delu_norm>`, :data:`<u_norm>` is the velocity norm, and :data:`<cg_iter>` is the number of iterations the conjugate gradient solver performed. All norms are computed as :math:`l^{2}`-norms or the square root of the inner product of a vector by itself.

Surface Forces Output
~~~~~~~~~~~~~~~~~~~~~

If the :data:`print forces` command is present in the :file:`simulation.config` file, loads for all surfaces in the performed simulation will be output with load values provided at each increment. The :file:`post.force.*` file names are constructed via the defined :data:`faset_label` strings in the :file:`simulation.msh` file. The default :data:`faset_label` order is :data:`x0, x1, y0, y1, z0, z1` which defines the six orthogonal and planar surfaces that bound a domain. For example, :data:`x0` refers to the face where the nodal coordinate component values in the :data:`x` direction are minimum and the file :file:`post.force.x0` would contain the surface-integrated forces on this face. These files are generally tabulated with the given structure:

::

    <step> <incr> <force_x> <force_y> <force_z> <surf_area> <current_time>


where :data:`<step>` is the prescribed load step, :data:`<incr>` is the total increment value being printed, :data:`<force_x>` is the surface-integrated force in the :data:`x` direction, :data:`<force_y>` is the surface-integrated force in the :data:`y` direction, :data:`<force_z>` is the surface-integrated force in the :data:`z` direction, :data:`<surf_area>` is the current surface area of the given face, and :data:`<current_time>` is the total simulated time at the time of printing.

If multiaxial loading is utilized, an additional :data:`<length>` column will be appended to the right of :data:`<current_time`. The :data:`<length>` column contains the maximal coordinate values of the domain and these values are stored in their associated face files. For example, the maximal mesh coordinate value in the :data:`x` direction is stored in the :file:`post.force.x0` and :file:`post.force.x1` files accordingly.

Simulation Report File
~~~~~~~~~~~~~~~~~~~~~~

The :file:`post.report` file is always printed for a simulation. The report file is for utilization with Neper and contains the following information:

::

    number_of_nodes <num_nodes>
    number_of_elements <num_elems>
    number_of_partitions <num_part>
    number_of_elements_byparition <part1_num_elems> ... <partN_num_elems>
    number_of_nodes_byparition <part1_num_nodes> ... <partN_num_nodes>
    number_of_slip_systems <num_slip_systems_for_crystal_type>
    orientation_definition <orientation_descriptor>:<orientation_convention>
    results_nodes <nodal_output_files>
    results_elements <elemental_output_files>
    number_of_steps <number_of_completed_steps>
