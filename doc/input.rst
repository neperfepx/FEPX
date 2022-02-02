.. _simulation_input:

Simulation Input
================

By default, a minimum of two files is necessary to completely define a simulation. These files are the configuration file, :file:`simulation.config`, and the mesh file, :file:`simulation.msh`. The configuration file defines the material parameters, control of the simulation (i.e., boundary conditions and loading history), printing of output files, and various optional input. The mesh file contains the polycrystal finite element mesh information (grain morphologies, phases and crystal orientations) as well as simulation-related information on the domain faces and the mesh partitions used for parallel simulations.

The Configuration File (:file:`simulation.config`)
--------------------------------------------------

This file contains the necessary definitions to run a simulation. It is structured in several, successive blocks that define different aspects of the simulation.  Each of these blocks is headed by a line starting by :data:`#`, which provides a short description of the block.  The structure for the :file:`simulation.config` file is as follows::

  # Optional Input
      <key_phrase> <value>
      ...

  # Material Parameters
      <key_phrase> <value>
      ...

  # Deformation History
      <key_phrase> <value>
      ...

  # Boundary Conditions
      <key_phrase> <value>
      ...

  # Printing Results
      <key_phrase> <value>
      ...

Configuration options within a block may generally be provided in any order; however, the material parameters and deformation history must follow specific orders. It is recommended that the overall structure of the :file:`simulation.config` file follow the example structure above. The :data:`Optional Input` block should always precede any others.

Any piece of text that is preceded by a :data:`#` is assumed to be a comment and is ignored (which also makes block headers optional), while :data:`<key_phrase> <value>` is the input structure of the file, where :data:`<key_phrase>` is the input command and :data:`<value>` are the associated parameters for the input command (as will be defined in the following sections). A single line in the file should only ever pertain to a single :data:`<key_phrase>`/:data:`<value>` pairing. All strings are interpreted literally and should be lowercase except where otherwise stated.

Material Description
~~~~~~~~~~~~~~~~~~~~

The model is described in :ref:`Model Description <model_description>`, and the corresponding input parameters are made clear and written in fixed font, most often as in: :math:`\gamma` (:data:`gamma`).  The specification of the input parameters in the :file:`simulation.config` file is detailed in :ref:`Input Parameters <input_parameters>`.

.. _model_description:

Model Description
^^^^^^^^^^^^^^^^^

The material is described via an elastic response (Hooke's law) and a plastic response (rate dependent plastic flow and hardening).

The stress, :math:`\sigma`, is related to the elastic strain, :math:`\epsilon`, via Hooke's law:

.. math::

    \sigma = \cal C \epsilon,

where :math:`\cal C` is the stiffness tensor. Written in Voigt notation, the above equation is expanded for cubic materials:

.. math::

    \begin{Bmatrix}
        \sigma_{11} \\
        \sigma_{22} \\
        \sigma_{33} \\
        \sigma_{23} \\
        \sigma_{13} \\
        \sigma_{12}
    \end{Bmatrix} =
    \begin{bmatrix}
        C_{11} & C_{12} & C_{12} & & & \\
        C_{12} & C_{11} & C_{12} & & & \\
        C_{12} & C_{12} & C_{11} & & & \\
        & & & C_{44} & & \\
        & & & & C_{44} & \\
        & & & & & C_{44}
    \end{bmatrix}
    \begin{Bmatrix}
        \epsilon_{11} \\
        \epsilon_{22} \\
        \epsilon_{33} \\
        2\epsilon_{23} \\
        2\epsilon_{13} \\
        2\epsilon_{12}
    \end{Bmatrix}.

Here, the strength of materials convention is utilized, where the shear factors of 2 are written in the strain vector. Special attention must be paid to ensure that the correct stiffness values are chosen, to align with the input convention used here.

For this convention, the Zener anisotropy ratio for cubic materials (which quantifies the level of elastic anisotropy, with 1 being perfectly isotropic) would be written as:

.. math::

    A = {2 C_{44} \over C_{11} - C_{12}}.

For example, Tungsten (W) is a nearly perfectly elastically isotropic cubic (BCC) material, with :math:`C_{11} = 522.4` GPa, :math:`C_{12} = 204.4` GPa, and :math:`C_{44} = 160.8` GPa. This would yield a Zener ratio of 1.01.

The elastic constitutive relation may also be expanded for hexagonal materials:

.. math::

    \begin{Bmatrix}
        \sigma_{11} \\
        \sigma_{22} \\
        \sigma_{33} \\
        \sigma_{23} \\
        \sigma_{13} \\
        \sigma_{12}
    \end{Bmatrix} =
    \begin{bmatrix}
        C_{11} & C_{12} & C_{13} & & & \\
        C_{12} & C_{11} & C_{13} & & & \\
        C_{13} & C_{13} & C_{33} & & & \\
        & & & C_{44} & & \\
        & & & & C_{44} & \\
        & & & & & \left( C_{11}-C_{12}\right)/2
    \end{bmatrix}
    \begin{Bmatrix}
        \epsilon_{11} \\
        \epsilon_{22} \\
        \epsilon_{33} \\
        2\epsilon_{23} \\
        2\epsilon_{13} \\
        2\epsilon_{12}
    \end{Bmatrix}.

Note that to allow for the decoupling of the hydrostatic and deviatoric portions of the elastic deformation, the following must be satisfied: :math:`C_{33} = C_{11} + C_{12} - C_{13}` (FEPX, thus, expects no input for :math:`C_{33}`).

For tetragonal materials, the elastic constitutive relation is expanded to:

.. math::

    \begin{Bmatrix}
        \sigma_{11} \\
        \sigma_{22} \\
        \sigma_{33} \\
        \sigma_{23} \\
        \sigma_{13} \\
        \sigma_{12}
    \end{Bmatrix} =
    \begin{bmatrix}
        C_{11} & C_{12} & C_{13} & & & \\
        C_{12} & C_{11} & C_{13} & & & \\
        C_{13} & C_{13} & C_{33} & & & \\
        & & & C_{44} & & \\
        & & & & C_{44} & \\
        & & & & & C_{66}
    \end{bmatrix}
    \begin{Bmatrix}
        \epsilon_{11} \\
        \epsilon_{22} \\
        \epsilon_{33} \\
        2\epsilon_{23} \\
        2\epsilon_{13} \\
        2\epsilon_{12}
    \end{Bmatrix}.

Again, to allow for the decoupling of the hydrostatic and deviatoric portions of the elastic deformation, :math:`C_{33} = C_{11} + C_{12} - C_{13}`.

Overall, the stiffness tensor is thus defined by the input parameters :data:`c11`, :data:`c12`, and :data:`c44` for cubic materials, :data:`c11`,  :data:`c12`,  :data:`c13`, and  :data:`c44` for hexagonal materials, and :data:`c11`, :data:`c12`, :data:`c13`, :data:`c44`, and :data:`c66` for tetragonal materials.

The kinematics of slip are described by a power law:

.. math::

    \dot{\gamma}^{\alpha} = \dot{\gamma}_{0} \left( \left| {\tau}^{\alpha} \right| \over g^{\alpha} \right)^{1/m} \rm sgn({\tau}^{\alpha}),

where :math:`\dot{\gamma}_0` (:data:`gammadot_0`) is the fixed-rate strain rate scaling coefficient (expressed in [1/s]), and :math:`m` (:data:`m`) is the rate sensitivity exponent.

.. note:: All variables presented are detailed by their dimensions (if applicable) instead of any specific unit. No unit system is inherently assumed by FEPX and the chosen unit system and value magnitudes should be consistent with the chosen length scale for the domain. For example, if it is assumed that the length scale is *mm* and SI units are to be used, then [force/area] will be understood to be *MPa*. The unit for time, however, is always assumed to be seconds (*s*).

For an isotropic hardening assumption, slip system strength evolution (hardening) is modeled by (note that for HCP and BCT materials, the implementation of an isotropic hardening assumption is  such that the shape of the single crystal yield surface is maintained. That is, the slip rates on the basal, prismatic, and pyramidal slip systems will harden such that the ratios of slip strengths remains constant. The consequence of this is that for HCP and BCT materials, this is not a "true" isotropic assumption, as the different slip families may harden at different rates, depending on the ratios of slip system strengths):

.. math ::

    \dot{g^{\alpha}} = h_{0} \left (g_{s0} - g^{\alpha} \over g_{s0} - g_{0} \right)^{n} \dot{\gamma},

where :math:`h_0` (:data:`h_0`) is the fixed-state hardening rate scaling coefficient, :math:`g_{s0}` (:data:`g_s0`) is the initial slip system saturation strength (expressed in [force/area]), :math:`g_0` (:data:`g_0`) is the initial slip system strength (expressed in [force/area]), and :math:`n` (:data:`n`) is the non-linear Voce hardening exponent. In the above equation, :math:`\dot{\gamma}` is calculated as:

.. math ::

    \dot{\gamma} = \sum_{\alpha} \left|\dot{\gamma}^{\alpha}\right|.

The slip system saturation strength may be evolved as a function of the slip activity. In this case, the hardening expression takes the form:

.. math ::

    \dot{g^{\alpha}} = h_{0} \left (g_{s}(\dot{\gamma}) - g^{\alpha} \over g_{s}(\dot{\gamma}) - g_{0} \right)^{n} \dot{\gamma},

where :math:`g_{s}(\dot{\gamma})` is the function for the saturation strength, which evolves via:

.. math ::

    g_{s}(\dot{\gamma}) = g_{s0} \left (\dot{\gamma} \over \dot{\gamma}_{s0} \right)^{m'},

where :math:`g_{s0}` (:data:`g_s0`) is the initial slip system saturation strength (expressed in [force/area]), :math:`m'` (:data:`m_prime`) is the saturation strength rate scaling exponent, and :math:`\dot{\gamma}_{s0}` (:data:`gammadot_s0`) is the initial saturation slip system shear rate. Again, in the above two equations, :math:`\dot{\gamma}` is calculated as the sum of the absolute value of the individual slip system shear rates, as defined above.

For a cyclic hardening assumption, the slip system strength evolution (hardening) is modeled by:

.. math ::

    \dot{g^{\alpha}} = h_{0} \left (g_{s}(\dot{\gamma}) - g^{\alpha} \over g_{s}(\dot{\gamma}) - g_{0} \right)^{n} f,

where :math:`f` is calculated as:

.. math ::

    f = \sum_{\beta = 0}^{n_{a}} \left|\dot{\gamma}^{\beta}\right| .

A slip system that contributes to hardening (:math:`n_{a}` total systems contributing to hardening) is that which has a change in shear greater than a critical value:

.. math ::

    \Delta\gamma_{crit} = a \left[\, g / g_{s}(\dot{\gamma}) \right]^{c},

where the material parameters here are :math:`a` (:data:`cyclic_a`) and :math:`c` (:data:`cyclic_c`). A more complete description can be found in Turkmen *et al.* [TURKMEN04]_. Note that minor differences exist between the implemented model described above and the formulation described in the paper.

For an anisotropic hardening assumption, slip system strength evolution (hardening) is modeled by:

.. math ::

    \dot{g^{\alpha}} = h_{0} \left (g_{s}(\dot{\gamma}) - g^{\alpha} \over g_{s}(\dot{\gamma}) - g_{0} \right)^{n} \dot{\gamma} h_{\alpha \beta},

where the model parameters are the same as the isotropic case described above, with the addition of :math:`h_{\alpha \beta}`, the slip interaction matrix. The slip interaction matrix only allows for interactions from direct and coplanar slip families. The slip interaction matrix is defined by the diagonal entry, :math:`d`, and the off-diagonal entries, :math:`h_{1},\dots, h_{n}`. These input parameters are defined by :data:`diag`, :data:`h1`, :data:`h2`, :data:`h3`, and :data:`h4` for FCC materials, :data:`diag`, :data:`h1`, :data:`h2`, :data:`h3`, :data:`h4`, :data:`h5`, and :data:`h6` for BCC materials, :data:`diag`, :data:`h1`, :data:`h2`, :data:`h3`, :data:`h4`, :data:`h5`, :data:`h6`, and :data:`h7` for HCP materials, and :data:`diag`, :data:`h1`, :data:`h2`, :data:`h3`, :data:`h4`, :data:`h5`, :data:`h6`, :data:`h7`, :data:`h8`, :data:`h9`, and :data:`h10` for BCT materials. A more complete description can be found in Carson *et al.* [CARSON17]_.

In any of the above hardening models, the base slip system strength may be modified to consider the effects of the presence of precipitates. This is performed via:

.. math ::

    g_{0} = g_{0} + a_{p} \left( f_{p} r_{p} \over b_{p}  \right)^{1 \over 2},

where :math:`a_{p}` (:data:`a_p`) is the precipitate hardening scaling coefficient, :math:`f_{p}` (:data:`f_p`) is the precipitate volume fraction, :math:`r_{p}` (:data:`r_p`) is the average precipitate diameter, and :math:`b_{p}` (:data:`b_p`) is the average Burgers' vector for the precipitate phase. Currently, the increase in strength due to the presence of precipitates is applied globally to all elements.

For a hexagonal material, the crystal c/a ratio (:data:`c_over_a`) must be defined. The initial slip strengths (:data:`g_0`) must be provided as a set of three values. The order of the values for :data:`g_0` correspond to the basal slip family strength, prismatic slip family strength, and pyramidal slip family strength.

For a tetragonal material, the crystal c/a ratio (:data:`c_over_a`) must be defined. The initial slip system strengths (:data:`g_0`) must be provided as a set of ten values. The order of the values for :data:`g_0` correspond to the slip families :math:`\left\{100\right)\left<001\right]`, :math:`\left\{110\right)\left<001\right]`, :math:`\left\{100\right)\left<010\right]`, :math:`\left\{110\right)\left<1 \bar 11\right]`, :math:`\left\{110\right)\left<1 \bar10\right]`, :math:`\left\{100\right)\left<011\right]`, :math:`\left\{001\right)\left<010\right]`, :math:`\left\{001\right)\left<110\right]`, :math:`\left\{011\right)\left<01 \bar 1\right]` and :math:`\left\{211\right)\left<01 \bar 1\right]`.

.. _input_parameters:

Input Parameters
^^^^^^^^^^^^^^^^

The material (as defined in the mesh file) can include one or several phases (to which grains are assigned), and the mechanical behavior of these phases must be defined accordingly. The number of phases must first be provided:

::

    number_of_phases <nphases>

The material parameters for a particular phase should be defined entirely for said phase before parameters for any subsequent phases are defined.

Each phase requires the specification of a consistent set of single-crystal material parameters, prefaced by

::

    phase <phase_id>

where :data:`<phase_id>`, the phase identification number, ranges from 1 to :data:`<nphases>`.

First, the crystal symmetry is defined by:

::

    crystal_type <ctype>

where the :data:`<ctype>` can be :data:`fcc`, :data:`bcc`, :data:`hcp`, and :data:`bct` for face-centered cubic, body-centered cubic, hexagonal close-packed, and body-centered tetragonal respectively.

The single-crystal elastic and plastic material parameters of the phase must be defined. Depending on the crystal symmetry, the total number of required parameters varies.

Anisotropic elastic constants are defined using the strength of materials convention, as described previously in :ref:`Model Description <model_description>`. The input is, for :data:`fcc` and :data:`bcc` crystal symmetry:

::

    c11 <modulus>
    c12 <modulus>
    c44 <modulus>

for :data:`hcp` crystal symmetry:

::

    c11 <modulus>
    c12 <modulus>
    c13 <modulus>
    c44 <modulus>

and for :data:`bct` crystal symmetry:

::

    c11 <modulus>
    c12 <modulus>
    c13 <modulus>
    c44 <modulus>
    c66 <modulus>

where :data:`<modulus>` are expressed in [force/area]. For :data:`hcp` and :data:`bct` materials, the :math:`C_{33}` (:data:`c33`) elastic constant is constrained by the other moduli and is not required as direct input.

For :data:`hcp` and :data:`bct` materials, an additional crystal parameter needs to be provided:

::

    c_over_a <ratio>

Crystallographic slip (plasticity) parameters are defined as:

::

    m <value(s)>
    gammadot_0 <value>
    h_0 <strength>
    g_0 <strength(s)>
    g_s0 <strength>
    n <value>

For hexagonal and tetragonal materials, multiple values may be provided for the rate sensitivity exponent, :data:`m`, and the initial slip system strengths, :data:`g_0`. If a single value is provided for :data:`m`, a constant rate sensitivity exponent is assumed across all slip families. Otherwise, three or ten values (for hexagonal or tetragonal, respectively) may be provided for :data:`m`, and the rate sensitivity exponents are applied on a per-family basis. Additionally, :data:`g_0` must be defined by three unique values for hexagonal materials and ten unique values for tetragonal materials. The order of the values for :data:`m` and :data:`g_0` for hexagonal materials correspond to the: basal slip family, prismatic slip family, and pyramidal slip family. For tetragonal materials, the ten values correspond to the slip families :math:`\left\{100\right)\left<001\right]`, :math:`\left\{110\right)\left<001\right]`, :math:`\left\{100\right)\left<010\right]`, :math:`\left\{110\right)\left<1 \bar 11\right]`, :math:`\left\{110\right)\left<1 \bar10\right]`, :math:`\left\{100\right)\left<011\right]`, :math:`\left\{001\right)\left<010\right]`, :math:`\left\{001\right)\left<110\right]`, :math:`\left\{011\right)\left<01 \bar 1\right]` and :math:`\left\{211\right)\left<01 \bar 1\right]`.

Saturation strength evolution (optional) is by default disabled, and may be enabled by defining both of the necessary parameters for saturation strength evolution. These parameters do not need to be defined if the saturation strength is not intended to evolve. To enable saturation strength evolution, define:

::

    m_prime <value>
    gammadot_s0 <value>

Cyclic hardening (optional) is by default disabled. It may be enabled via:

::

 hard_type cyclic_isotropic

If cyclic hardening is enabled, each phase requires the definition of two additional parameters by:

::

    cyclic_a <cyc_a>
    cyclic_c <cyc_c>

where both :data:`cyclic_a` and :data:`cyclic_c` values are model parameters for a critical value of accumulated shear strain used to modify the form of the Voce hardening law [TURKMEN04]_.

Anisotropic hardening (optional) is by default disabled. It may be enabled via:

::

    hard_type anisotropic

If anisotropic hardening is enabled, each phase requires the definition of slip interaction matrix values which vary based on crystal symmetry.

For :data:`fcc` crystal symmetry:

::

    latent_parameters <diag> <h1> <h2> <h3> <h4>

For :data:`bcc` crystal symmetry:

::

    latent_parameters <diag> <h1> <h2> <h3> <h4> <h5> <h6>

For :data:`hcp` crystal symmetry:

::

    latent_parameters <diag> <h1> <h2> <h3> <h4> <h5> <h6> <h7>

For :data:`bct` crystal symmetry:

::

    latent_parameters <diag> <h1> <h2> <h3> <h4> <h5> <h6> <h7> <h8> <h9> <h10>

where :data:`<diag>` is the diagonal coefficient and :data:`h{1-10}` are the in-plane interaction coefficients [CARSON17]_.

Strengthening due to the presence of precipitates is by default disabled, and may be enabled by defining the necessary parameters. These parameters do not need to be defined if precipitate strengthening is not intended to be considered. To enable precipitate strengthening, define:

::

    a_p <strength>
    f_p <volume_fraction>
    r_p <length>
    b_p <length>

.. _deformation_history:

Deformation History
~~~~~~~~~~~~~~~~~~~

A variety of deformation modes are available that are capable of reproducing various mechanical loading configurations. A deformation history is defined by both steps and increments, where steps are made of one or several increments. Steps define the strain or load targets that are to be reached during the simulation, and results can only be printed at the end of steps. One or several increments occur within each step to reach the prescribed step target while ensuring numerical stability. A step target can be expressed in terms of strain or load. The strain refers to the engineering strain as computed from the displacement of the loading surface and the initial sample length along the loading direction, and the load refers to the total force on the loading surface. Of course, the relative order of the steps defined for the deformation history matters and should be written in an ascending manner.

Deformation histories are divided into uniaxial loading and multiaxial loading. In general, the multiaxial loading definition is technically triaxial in nature; however, biaxial loading may be performed by zeroing one of the load columns accordingly. The available deformation history configuration options follow.

Uniaxial
^^^^^^^^

Uniaxial loading is always strain controlled (i.e., constant strain rate); however, either specific strain targets or specific load targets may be prescribed. For strain targeting, the number of increments for a given step must be provided as opposed to a time-step value. For load targeting, the bounds on the time-step value are provided in order to control both the accuracy and, indirectly, the number of increments taken per step. These time-step values should be defined relative to the :data:`strain_rate` value (:ref:`Boundary Conditions <boundary_conditions>`).

Strain targeting allows the definition of loading to specific uniaxial strain states. This deformation history is defined as follows:

::

    def_control_by uniaxial_strain_target
    number_of_strain_steps <nsteps>
    target_strain <target_val> <n_incr> <print_flag>
    ...

where :data:`<nsteps>` is the number of strain steps that are defined in the file after this line, :data:`<target_val>` is the desired strain value to be reached, :data:`<n_incr>` is the number of increments to be performed in order to complete the step, and :data:`<print_flag>` allows for the printing (or not) of specific steps. The options available for :data:`<print_flag>` are: :data:`print_data` or :data:`suppress_data`.

Load targeting allows the definition of loading to specific uniaxial load states. This deformation history is defined as follows:

::

    def_control_by uniaxial_load_target
    number_of_load_steps <nsteps>
    target_load <target_val> <dt_max> <dt_min> <print_flag>
    ...

where :data:`<nsteps>` is the number of load steps that are defined in the file after this line, :data:`<target_val>` is the desired load value to be reached, :data:`<dt_max>` is the maximum time-step value to be used for a given increment, :data:`<dt_min>` is the minimum time-step value to be used for a given increment, and :data:`<print_flag>` allows for the printing (or not) of specific steps. The options available for :data:`<print_flag>` are: :data:`print_data` or :data:`suppress_data`.

Strain rate jumps are also available for both uniaxial deformation modes and are defined by adding the following input to the block:

::

    number_of_strain_rate_jumps <njumps>
    strain_rate_jump <target_step> <new_strain_rate>
    ...

where :data:`<njumps>` is the number of strain rate jumps defined in the file after this line, :data:`<target_step>` defines which :data:`target_strain` step is assigned a new strain rate, and :data:`<new_strain_rate>` is the new strain rate to be assigned and has units of [1/s]. In general, and for numerical stability, the strain rate jumps should be of a similar magnitude to the :data:`strain_rate` defined previously.

Multiaxial
^^^^^^^^^^

Multiaxial loading is always strain controlled (internally) and operates at either a constant engineering strain rate or constant load rate; however, only specific load targets may be prescribed.
The principal loading directions must be aligned with the coordinate axes of the mesh and the surface face normals should likewise be coincident with the coordinate axis of the mesh. Symmetry boundary conditions (zero normal velocities) are enforced on the three faces of minimal coordinates (:data:`<*0>`), and, in the general case, non-zero normal velocities are applied to the faces of maximal coordinates (:data:`<*1>`). The velocity on the primary control surface is held constant through the simulation (except during a strain rate jump).

Multiaxial loading with a constant strain rate (CSR) is defined as follows:

::

    def_control_by triaxial_constant_strain_rate
    number_of_csr_load_steps <nsteps>
    target_csr_load <load_x> <load_y> <load_z> <dt_max> <dt_min> <print_flag>
    ...


where :data:`<nsteps>` is the number of CSR load steps that are defined in the file after this line, :data:`<load_x>` is the desired load value to be reached in the :data:`x` direction, :data:`<load_y>` is the desired load value to be reached in the :data:`y` direction, :data:`<load_z>` is the desired load value to be reached in the :data:`z` direction, :data:`<dt_max>` is the maximum time-step value to be used for a given increment, :data:`<dt_min>` is the minimum time-step value to be used for a given increment, and :data:`<print_flag>` allows for the printing (or not) of specific steps. The options available for :data:`<print_flag>` are: :data:`print_data` or :data:`suppress_data`.

Strain rate jumps are available for this deformation mode and are defined by adding the following input to the block:

::

    number_of_strain_rate_jumps <njumps>
    strain_rate_jump <target_step> <new_strain_rate>
    ...

where :data:`<njumps>` is the number of strain rate jumps defined in the file after this line, :data:`<target_step>` defines which :data:`target_csr_load` step is assigned a new strain rate, and :data:`<new_strain_rate>` is the new strain rate to be assigned and has units of [1/s].

Multiaxial loading with a constant load rate (CLR) is defined as follows:

::

    def_control_by triaxial_constant_load_rate
    number_of_clr_load_steps <nsteps>
    target_clr_load <load_x> <load_y> <load_z> <target_time_incr> <print_flag>
    ...

where :data:`<nsteps>` is the number of CLR load steps that are defined in the file after this line, :data:`<load_x>` is the desired load value to be reached in the :data:`x` direction, :data:`<load_y>` is the desired load value to be reached in the :data:`y` direction, :data:`<load_z>` is the desired load value to be reached in the :data:`z` direction, :data:`<target_time_incr>` is the physical time increment to be reached for the given :data:`target_clr_load` steps for a given load rate, and :data:`<print_flag>` allows for the printing (or not) of specific steps. The options available for :data:`<print_flag>` are: :data:`print_data` or :data:`suppress_data`.

Load rate jumps and dwell episodes are available for this deformation mode. A dwell episode maintains the macroscopic loads of the step in which it is defined, but holds the ramp rate at zero for the amount of time defined by :data:`<dwell_time>`. These options are defined as follows:

- For load rate jumps:

  ::

      number_of_load_rate_jumps <njumps>
      load_rate_jump <target_step> <new_ramp_rate>
      ...

  where :data:`<njumps>` is the number of load rate jumps defined in the file after this line, :data:`<target_step>` defines which :data:`target_clr_load` step is assigned a new load rate, and :data:`<new_load_rate>` is the new load rate to be assigned and has units of [force/s].


- For dwell episodes:

  ::

    number_of_dwell_episodes <nepisodes>
    dwell_episode <target_step> <dwell_time> <target_time_incr> <print_flag>
    ...

  where :data:`<nepisodes>` is the number of dwell episodes defined in the file after this line, :data:`<target_step>` defines which :data:`target_clr_load` step is assigned to dwell, :data:`<dwell_time>` is the physical amount of time in [s] for a given dwell episode, :data:`<target_time_incr>` is the physical time increment to be reached for the given dwell episode, and :data:`<print_flag>` allows for the printing (or not) of specific steps. The options available for :data:`<print_flag>` are: :data:`print_data` or :data:`suppress_data`.

.. _boundary_conditions:

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

Standard, simple boundary conditions are available for automatic definition with minimal input and are computed internally for each simulation based on the definitions in the :file:`simulation.config` file. This ensures that standard boundary conditions are consistently defined for all simulations and increases the portability of the :file:`simulation.config` file.
Alternatively, custom boundary conditions can be defined, as described separately, in :ref:`Externally Defined Boundary Conditions (Optional) <external_bcs>`.

Uniaxial
^^^^^^^^

Uniaxial definitions are available for three different constraint configurations. The available uniaxial constraint configuration options follow.

- *Grip boundary conditions* fully constrain two opposite faces in the spatial domain. The first face is fully fixed in all sample directions while the second face has a strain rate applied in the face normal direction while the other two sample directions are fully fixed. All other faces are unconstrained.

  .. figure:: images/gripbcs.png
     :width: 50%
     :align: center

     Simplified schematic of the applied velocities for grip boundary conditions. The loading face is :data:`z1` and the sample is being loading in the :data:`+Z` direction.

  Grip boundary conditions are defined as follows::

    boundary_conditions uniaxial_grip
    loading_direction <sample_dir>
    loading_face <face_label>
    strain_rate <strain_rate>


  where :data:`<sample_dir>` is the direction along the a positive sample axis in which the sample is loaded, :data:`<face_label>` is the face on which the loading is applied (the opposing face is fully fixed), and :data:`<strain_rate>` is the strain rate value in units of [1/s].


- *Symmetry boundary conditions* constrain four faces in the spatial domain. The three :data:`<*0>` faces are fixed in the face normal directions and unconstrained in the other two sample directions. The fourth :data:`<*1>` face has a strain rate applied in the face normal direction while the other two sample directions are fully fixed. The selection of the :data:`<*1>` face is based on the defined :data:`loading_direction`.

  .. figure:: images/symmbcs.png
     :width: 50%
     :align: center

     Simplified schematic of the applied velocities for symmetry boundary conditions. The sample is being loaded in the :data:`+Z` direction.

  Symmetry boundary conditions are defined as follows::

    boundary_conditions uniaxial_symmetry
    loading_direction <sample_dir>
    strain_rate <strain_rate>

  where :data:`<sample_dir>` is the direction along the a positive sample axis in which the sample is loaded, and :data:`<strain_rate>` is the strain rate value in units of [1/s].


- *Minimal boundary conditions* are a modification of grip boundary conditions that only constrain two opposite faces in the face normal directions and two corner nodes in the spatial domain. The selection of the constrained faces is based on the defined :data:`loading_direction`. The first node is always fully fixed where the :data:`<*0>` faces converge. The second node is defined relative to the defined :data:`loading_direction` and is constrained to prevent rigid body rotation about the :data:`loading_direction` axis.

  .. figure:: images/minimalbcs.png
     :width: 50%
     :align: center

     Simplified schematic of the applied velocities for minimal boundary conditions. The sample is being loaded in the :data:`+Z` direction. The two blue corner nodes are constrained to prevent rigid body translation and motion.

  Minimal boundary conditions are defined as follows::

    boundary_conditions uniaxial_minimal
    loading_direction <sample_dir>
    strain_rate <strain_rate>

  where :data:`<sample_dir>` is the direction along the a positive sample axis in which the sample is loaded, and :data:`<strain_rate>` is the strain rate value in units of [1/s].

Multiaxial
^^^^^^^^^^

Multiaxial boundary conditions are generally consistent across modes, however, the input rate type varies depending on the mode. For both modes, the :data:`loading_direction` defines the primary control direction in which the normal velocities are held constant throughout the simulation.

Multiaxial loading with a constant strain rate (CSR) is defined as follows:

::

    boundary_conditions triaxial
    loading_direction <sample_dir>
    strain_rate <strain_rate>

where :data:`<sample_dir>` is the direction along the a positive sample axis in which the sample is loaded, and :data:`<strain_rate>` is the strain rate value in units of [1/s].

Multiaxial loading with a constant load rate (CLR) is defined as follows:

::

    boundary_conditions triaxial
    loading_direction <sample_dir>
    load_rate <load_rate>

where :data:`<sample_dir>` is the direction along the a positive sample axis in which the sample is loaded, and :data:`<load_rate>` is the loading rate value in units of [force/s].

.. _external_bcs:

Externally Defined Boundary Conditions (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Boundary conditions different from the ones defined in :ref:`Boundary Conditions <boundary_conditions>` can be applied using an external :file:`simulation.bcs` containing per-node constraints. Of course, this file should be generated for a singular :file:`simulation.msh` file (if the finite element mesh in the :file:`simulation.msh` file changes, then the associated :file:`simulation.bcs` file will need to be generated anew).

To read in external boundary conditions, the following line must be added to the :file:`simulation.config` file:

::

    read_bcs_from_file

The per-node constraints are defined as follows, in the :file:`simulation.bcs` file:

::

    <node_id> <coord_index> <vel>
    ...

where :data:`<node_id>` is a unique 1-indexed identification number, :data:`<coord_index>` defines the sample axis the constraint is applied to, and :data:`<vel>` is the velocity being applied to the node in the constraint direction. The options for :data:`<coord_index>` are: :data:`x`, :data:`y`, or :data:`z`. A singular :data:`<coord_index>`/:data:`<vel>` pair should be defined per-line for a given :data:`<node_id>`.  The velocities should be prescribed relative to the mesh dimensions and time-step size in order to produce expected strain rates.

.. _external_ori_phase:

Externally Defined Orientations and Phases (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Crystallographic phase and orientations different from the ones defined in the :file:`simulation.msh` file can be defined in external files by adding appropriate commands to the configuration file.

To read in external orientations, the following line must be added to the :file:`simulation.config` file:

::

    read_ori_from_file

More detailed information on the structure of this external file can be found in :ref:`External Orientation Assignment (Optional) <external_oris>`.

To read in external grain/phase assignments, the following line must be added to the :file:`simulation.config` file:

::

    read_phase_from_file

More detailed information on the structure of this external file can be found in :ref:`External Crystallographic Phase Assignment (Optional) <external_phase>`.

.. _printing_results:

Printing Results
~~~~~~~~~~~~~~~~

Each field variable file to be output from a simulation must be individually defined. This includes nodal output, elemental output, simulation restart information, and other miscellaneous output (see :ref:`Simulation Output <simulation_output>` for a complete description of all output). The printing of a given field variable file is defined as follows:

::

    print <output_file_name>
    ...

where :data:`<output_file_name>` is the particular field variable file to be output. The available options for :data:`<output_file_name>` are:


- :data:`coo`: Nodal coordinates

- :data:`crss`: Critical resolved shear stress

- :data:`defrate`: Deformation rate tensor

- :data:`defrate-eq`: Equivalent deformation rate

- :data:`defrate-pl`: Plastic deformation rate tensor

- :data:`defrate-pl-eq`: Equivalent plastic deformation rate

- :data:`disp`: Nodal displacements

- :data:`elt-vol`: Elemental volume

- :data:`ori`: Crystallographic orientations

- :data:`slip`: Slip system shear

- :data:`sliprate`: Slip system shear rate

- :data:`spinrate`: Plastic spin rate tensor

- :data:`strain`: Total strain tensor

- :data:`strain-eq`: Equivalent total strain

- :data:`strain-el`: Elastic strain tensor

- :data:`strain-el-eq`: Equivalent elastic strain

- :data:`strain-pl`: Plastic strain tensor

- :data:`strain-pl-eq`: Equivalent plastic strain

- :data:`stress`: Stress tensor

- :data:`stress-eq`: Equivalent stress

- :data:`vel`: Nodal velocity

- :data:`velgrad`: Velocity gradient tensor

- :data:`work`: Work

- :data:`work-pl`: Plastic work

- :data:`workrate`: Work Rate

- :data:`workrate-pl`: Plastic work rate

- :data:`forces`: Surface forces

- :data:`convergence`: Simulation convergence statistics

- :data:`restart`: Simulation restart data

A full description of each output variable can be found in :ref:`Simulation Output <simulation_output>`.

.. _optional_input_parameters:

Optional Input Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

These options may pertain to specific deformation modes or control standard simulation behavior. All possible inputs presented in this section have default values already defined.

- :data:`max_incr \<increments\>` specifies the maximum number of increments (default: :data:`50000`).

- :data:`max_total_time \<time\>` specifies the maximum deformation time (default: :data:`12000.0`).

- :data:`check_necking {on,off}` specifies whether or not to terminate simulation when specimen begins to neck (default: :data:`off`).

- :data:`load_tol \<tolerance\>` is the the target load tolerance. A small positive load tolerance (e.g. 0.1 :math:`\times` control surface area) improves load control while reducing the number of small steps near target loads (default: :data:`0.0`).

- :data:`dtime_factor \<time\>` is a number greater than or equal to 1 which is used when calculating time increments near target loads (default: :data:`1.001`).

- :data:`hard_type {isotropic,anisotropic,cyclic_isotropic}` specifies the hardening model to use (default: :data:`isotropic`).

- :data:`max_bc_iter \<iterations\>` specifies the maximum number of boundary condition iterations (default: :data:`10`).

- :data:`min_pert_frac \<fraction\>` is the minimum fraction of the control velocity by which the secondary and tertiary surface velocities are perturbed during boundary condition iterations (default: :data:`0.001`).

- :data:`load_tol_abs \<tolerance\>` is the absolute tolerance on the secondary and tertiary loads. The absolute load criterion is that both loads are within the absolute load tolerance of the ideal load. Loads are considered to be within tolerance if either the absolute or relative criterion is satisfied (default: :data:`0.1`).

- :data:`load_tol_rel \<tolerance\>` is the relative load tolerance on the secondary and tertiary loads. It represents a fraction of the load in the control direction. The relative load criterion is that the difference between the load and ideal load, normalized by the load in the control direction, is less than the relative load tolerance. Loads are considered to be within tolerance if either the absolute or relative criterion is satisfied (default: :data:`0.001`).

- :data:`max_strain_incr \<increments\>` specifies the maximum strain increment for dwell episodes (default: :data:`0.001`).

- :data:`max_strain \<strain\>` specifies the maximum allowable macroscopic strain (default: :data:`0.2`).

- :data:`max_eqstrain \<strain\>` specifies the maximum allowable macroscopic equivalent strain (default: :data:`0.2`).

- :data:`max_iter_hard_limit \<iterations\>` specifies the maximum allowable iterations on the Backward Euler approximation used to update hardnesses (default: :data:`10`).

Additional input commands related to restart capabilities and external file read-in are presented separately in :ref:`Restarting a Simulation <sim_restart>`, :ref:`External Orientation Assignment (Optional) <external_oris>`, and :ref:`External Crystallographic Phase Assignment (Optional) <external_phase>`, respectively.

Optional Convergence Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These options modify the tolerances and general behavior of the solution algorithms and should only be modified by those who know what they are doing. All possible inputs presented in this section have default values already defined.

Velocity Convergence
^^^^^^^^^^^^^^^^^^^^

The velocity solver employs a hybrid successive-approximation/Newton-Raphson algorithm. Convergence of the velocity solution is based on a convergence parameter, which unless otherwise noted, is defined as the norm of the change in the velocity field, divided by the norm of the velocity field, :math:`||\Delta u||/||u||`. Other parameters are also used to assess the convergence of the velocity solution. The following parameters pertain to the convergence of the velocity solver:

- :data:`nl_max_iters \<iterations\>`: maximum allowable number of iterations of the nonlinear velocity solver (default: :data:`50`).

- :data:`nl_tol_strict \<tolerance\>`: desired tolerance on the elasto-viscoplastic velocity solution (default: :data:`5e-4`).

- :data:`nl_tol_loose \<tolerance\>`: acceptable level of convergence if the desired level of convergence cannot be reached via :data:`nl_tol_strict` (default: :data:`5e-4`).

- :data:`nl_tol_min \<tolerance\>`: tolerance on the norm of the change in velocity, divided by the number of degrees of freedom, :math:`(||u||/ \rm max(ndof))`. This parameter is useful for assessing convergence when the macroscopic velocity is near zero (default: :data:`1e-10`).

- :data:`nl_tol_switch_ref \<tolerance\>`: value of the convergence parameter at which the solution algorithm switches from successive-approximation to Newton-Raphson. To only use successive-approximations, set the value of :data:`nr_tol_switch_ref` equal to the value of :data:`nl_tol_strict` (default: :data:`1e-2`).

- :data:`nl_tol_conv \<tolerance\>`: parameter between 0 and 1 that is used to assess whether the Newton-Raphson algorithm is converging slowly (default: :data:`0.2`).

Conjugate Gradient Convergence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The solution of the linear system of equations :math:`[K]\{\Delta u\} = -\{R\}` is performed using a conjugate gradient solver. The following parameters pertain to the convergence of the conjugate gradient solver:

- :data:`cg_max_iters \<iterations\>`: maximum allowable number of iterations of the conjugate gradient solver (default: :data:`16000`).

- :data:`cg_tol \<tolerance\>`: desired tolerance on the conjugate gradient solver (default: :data:`1e-8`).

Material State Convergence
^^^^^^^^^^^^^^^^^^^^^^^^^^

The convergence of the material stress state for both the viscoplastic and elasto-viscoplastic solutions is assessed by the following parameters:

- :data:`sx_max_iters_state \<iterations\>`: maximum number of iterations on material state (default: :data:`100`).

- :data:`sx_max_iters_newton \<iterations\>`: maximum number of iterations of the Newton algorithm used to solve for crystal stress (default: :data:`100`).

- :data:`sx_tol \<tolerance\>`: tolerance on the stress solution (default: :data:`1e-4`).

The Mesh File (:file:`simulation.msh`)
--------------------------------------

This file contains the finite element mesh information along with phase assignments and crystal orientations. The mesh file is generally generated by Neper and not directly modified.  A brief description is provided below, which a more complete description can be found in the `Neper reference manual <https://neper.info>`_.  The file can be opened by Gmsh for interactive visualization.

The file is structured in several, successive fields that define different aspects of the mesh. Each of these fields is wrapped by :data:`$<Field>/$End<Field>` lines, where :data:`<Field>` is a short description of the information stored within the block. A typical :file:`simulation.msh` file will contain the following information:

- Mesh Format (:data:`$MeshFormat`),

- Mesh Version (:data:`$MeshVersion`),

- Nodes (:data:`$Nodes`),

- Elements (:data:`$Elements`),

- Surface Element Sets (:data:`$Fasets`),

- Crystal Orientations (:data:`$ElsetOrientations` or :data:`$ElementOrientations`),

- Grain/Phase Assignments (:data:`$Groups`).

Additionally, the :file:`simulation.msh` file may also include fields with partition information for both the nodes and elements if the domain is decomposed for parallel execution and surface node sets (:data:`$NSets`).

Embedded microstructural information (phases and orientations) with the :file:`simulation.msh` may be overridden by external files, :file:`simulation.ori` and :file:`simulation.phase`, if the appropriate commands are added to the :file:`simulation.config` file (:ref:`Externally Defined Orientations or Phase Assignments (Optional) <external_ori_phase>`).

.. _external_oris:

External Orientation Assignment (Optional, :file:`simulation.ori`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The embedded orientation assignments within the :file:`simulation.msh` may be overridden via an external :file:`simulation.ori` file. This file contains formatting identical to the associated fields in the mesh file and is defined as:

For per-grain (or :data:`Elset`) orientations:

::

    $ElsetOrientations
    <number_of_ori_entities> <orientation_descriptor>:<orientation_convention>
    <entity_id> <ori_des1> ...
    ...
    $EndElsetOrientations

For per-element orientations:

::

    $ElementOrientations
    <number_of_ori_entities> <orientation_descriptor>:<orientation_convention>
    <entity_id> <ori_des1> ...
    ...
    $EndElementOrientations

where :data:`<number_of_ori_entities>` is the number of unique orientations defined in the section, :data:`<orientation_descriptor>` is the parameterization for the orientations (see options below), :data:`<orientation_convention>` describes the basis transformation route for the orientations provided, :data:`<entity_id>` is a unique 1-indexed identification number, and :data:`<ori_des*>` are the components of the unique orientation. Available options for :data:`<orientation_convention>` are: :data:`active` or :data:`passive`. Following the usual terminology, an active orientation assumes that which describes a basis transformation from the sample basis to the crystal basis (sample-to-crystal), while a passive orientation convention assumes that which describes a basis transformation from the crystal basis to the sample basis (“crystal-to-sample”).

The following :data:`<orientation_descriptor>` types are available (associated per-line formats are also described):

- For :data:`rodrigues`, each orientation is described by :math:`r_1, r_2, r_3`, where :math:`{\bf r} = {\bf t} \tan{ (\omega / 2)}` (see: :data:`axis-angle` for definitions of :math:`{\bf t}` and :math:`\omega`). The per-line format is:

::

    <entity_id> <r_1> <r_2> <r_3>

- For :data:`euler-bunge`, each orientation is described by :math:`\phi_1, \theta, \phi_2`, where :math:`\phi_1` is the rotation about the :math:`z` axis, :math:`\theta` is the rotation about the :math:`z^{\prime}` axis, and :math:`\phi_2` is the rotation about the :math:`z^{\prime \prime}` axis, all in degrees). The per-line format is:

::

    <entity_id> <phi_1> <Phi> <phi_2>


- For :data:`euler-kocks`, each orientation is described by :math:`\Psi, \Theta, \phi`, where :math:`\Psi` is the rotation about the :math:`z` axis, :math:`\Theta` is the rotation about the :math:`y^{\prime}` axis, and :math:`\phi` is the rotation about the :math:`z^{\prime \prime}` axis, all in degrees). The per-line format is:

::

    <entity_id> <Psi> <Theta> <phi>


- For :data:`axis-angle`, each orientation is described by :math:`t_1, t_2, t_3, \omega`, where :math:`\bf{t}` is the normalized axis of rotation and :math:`\omega` is the angle of rotation about said axis, in degrees. The per-line format is:

::

    <entity_id> <t_1> <t_2> <t_3> <omega>

- For :data:`quaternion`, each orientation is described by :math:`q_0, q_1, q_2, q_3`, where :math:`q_0 = \cos{(\omega / 2)}` and :math:`q_i = t_i \sin{(\omega / 2)}` for :math:`i = 1, 2, 3` (see: :data:`axis-angle` for definitions of :math:`{\bf t}` and :math:`\omega`). The per-line format is:

::

    <entity_id> <q_0> <q_1> <q_2> <q_3>

.. _external_phase:

External Crystallographic Phase Assignment (Optional, :file:`simulation.phase`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The embedded grain/phase assignments within the :file:`simulation.msh` may be overridden via an external :file:`simulation.phase` file. This file contains formatting identical to the associated fields in the mesh file and is defined as::

    $Groups
    <group_entity>
    <number_of_group_entities>
    <entity_id> <group>
    ...
    $EndGroups

where :data:`<group_entity>` defines the phase assignment method and must always be defined as :data:`elset`, :data:`<number_of_group_entities>` is the number of grain/phase pairs defined in the field, :data:`<entity_id>` is a unique 1-indexed identification number, and :data:`<group>` is an 1-indexed value that defines the phase for a given grain.

References
----------

.. [CARSON17] R. Carson, M. Obstalecki, M. Miller, and P. R. Dawson. Characterizing heterogeneous intragranular deformations in polycrystalline solids using diffraction-based and mechanics-based metrics. *Modelling and Simulation in Materials Science and Engineering*, 25:055008, 2017.

.. [TURKMEN04] H. S. Turkmen, M. P. Miller, P. R. Dawson, and J. C. Moosbrugger. A slip-based model for strength evolution during cyclic loading. *Journal of Engineering Materials and Technology*, 126(4):329-338, 2004.
