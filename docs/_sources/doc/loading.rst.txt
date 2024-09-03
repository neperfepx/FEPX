.. _loading:

Loading (:file:`simulation.cfg`)
===================================

.. _boundary_conditions:

Boundary Conditions
-------------------

General boundary conditions can be applied (to general domains) by assiging velocities to nodes or node sets. Multi-point constraints are also available, to couple node velocities.  RVE-type boundary conditions can also by applied, which are based on RVE-type variables, such as a strain rate.

General Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

General boundary conditions are velocity conditions applied on specific nodes or node sets defined in the :ref:`mesh_file`.
Any number of conditions can be prescribed as follows::

    set_bc vel <nset> <direction> <vel> ...
    set_bc vel <nset> <direction> <vel>
    ...

where each :data:`set_bc` line defines a velocity boundary condition, :data:`<nset>` is the label of a node set (alternatively a node id), :data:`<direction>` is the direction of the constraint applied to the node set, and :data:`<vel>` is the velocity applied to the node set in the constraint direction.  The options for :data:`<direction>` are: :data:`<x>`, :data:`<y>` or :data:`<z>`. Up to three :data:`<direction>`/:data:`<vel>` pairs can be defined per-line for a given :data:`set_bc <nset>`.

.. _mpcs:

Multi-Point Constraints (MPCs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Multi-point constraints link the degrees of freedom of several nodes by a linear relation between nodal velocities of primary and secondary nodes [#mpc_note]_:

  .. math::
     v_{x,y,z}^s=v_{x,y,z}^p+D_{x,y,z}

where :math:`v_{x,y,z}^s` and :math:`v_{x,y,z}^p` are the velocity on the primary and secondary nodes, respectively and :math:`D_{x,y,z}` the velocity offset.

The following (and only) type of MPCs, :data:`mpc1`, is such that one or several velocity component(s) of all nodes belonging to a node set are equal::

    set_mpc1 vel <nset> <direction> ...

where :data:`<nset>` is the node set and :data:`<d>` is the direction  (among :data:`x`, :data:`y` and :data:`z`).  Up to three directions can be defined per line for a given :data:`set_mpc1 vel \<nset\>`.

Any number of conditions can be prescribed as follows::

    set_mpc1 vel <nset> <direction> ...
    set_mpc1 vel <nset> <direction> ...
    ...

MPCs can typically be used to impose that a face should remain planar during a simulation (but not necessarily at the same coordinate).  An example of use is provided in :ref:`cruciform_specimen_simulation`.

RVE-type Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

RVE-type boundary conditions apply to cubic (or right prisms) domains, for which the loading can typically be defined in terms of (components of the) velocity gradient or strain rate.  They may also be applied to to other domains, such as cylinders, or periodic tessellation.
For a cubic domain, the faces are named :data:`x0`, :data:`y0` and :data:`z0` at minimum (typically 0) coordinates, and :data:`x1`, :data:`y1`, :data:`z1` at maximum (typically 1) coordinates.  In this section, the minimum-coordinate face is called the *reference* face, while the maximum-coordinate face is called the *loading* face.

To apply strain rates along the sample axes, the input is::

    set_bc strainrate <component> <strainrate> [ <method> ]
    ...

To apply velocity gradient components (which allows for shear), the input is::

    set_bc velgrad <component> <velgrad> [ <method> ]
    ...

:data:`<component>` can be :data:`11`, :data:`22` or :data:`33` for :data:`strainrate`, and
:data:`11`, :data:`12`, :data:`13`, :data:`21`, :data:`22`, :data:`23`,
:data:`31`, :data:`32` or :data:`33` for :data:`velgrad`.
The (optional) method argument indicates how to translate the specified velocity gradient or strain rate into actual velocities, on the corresponding faces (other faces are unconstrained).  Options for :data:`<method>` are:

  - :data:`minimal` (the default): zero velocity along the loading direction is applied on the reference face, while prescribed velocity along the loading direction is applied on the loading face.  Other directions are unconstrained.

  - :data:`grip`: fully constrained mode, where zero velocity in all directions is applied on the reference face, while prescribed velocity in the loading direction is applied on the loading face, with 0 velocity in other directions.  Note that :data:`grip` boundary conditions make particular sense when the loading face is not expected to change shape (or surface area), like in pure shear; however, applying :data:`grip` boundary conditions in different directions will generally result in conflicting velocity conditions (and an error).

  - :data:`periodic`: to apply Periodic Boundary Conditions (PBCs) in the case of a periodic tessellation. To use with :data:`set_bc strainrate`. Available for triaxial traction-compression loading, i.e. non-zero components can be :data:`11`, :data:`22` or :data:`33`. The shear components :data:`12` (or :data:`21`), :data:`23` (or :data:`32`) and :data:`13` (or :data:`31`) must be specified and equal to :data:`0`. Otherwise explicitly indicated in the simulation, rigid-body motions do not need to be fixed by the user when at least one :data:`periodic` option is specified in :file:`simulation.cfg`.

.. note:: Several :data:`set_bc vel`, :data:`set_bc strainrate` and :data:`set_bc velgrad` inputs can be combined as long as they do not conflict.  One or several :data:`set_bc vel` inputs can for example be used to fix the rigid-body motions that may remain as a result of :data:`set_bc velgrad` or :data:`set_bc strainrate` boundary conditions.

.. _deformation_history:

Deformation History (Steps)
---------------------------

.. A variety of deformation modes are available that are capable of reproducing various mechanical loading configurations. A deformation history is defined by both steps and increments, where steps are made of one or several increments. Steps define the time, strain or load targets that are to be reached during the simulation, and results can only be printed at the end of steps. One or several increments occur within each step to reach the prescribed step target while ensuring numerical stability.  The time and strain targets are cumulative.  The strain refers to the engineering strain as computed from the displacement of the loading surface and the initial sample length along the loading direction, and the load refers to the total force on the loading surface. Of course, the relative order of the steps defined for the deformation history matters and should be written in an ascending manner.

A simulation is divided into a series of *steps*.  There can be any number of steps, and different types of *targets* can be specified for the steps; however, the same type must be used for all steps.  (If the number of steps is not defined, it is determined from the number of provided targets.)  The format for specifying steps is as follows:

  .. code::

     [ number_of_steps <n> ]

     [ target_time      <time_1>   ... <time_n>   ]
     [ target_strain<c> <strain_1> ... <strain_n> ]
     [ target_load<d>   <load_1>   ... <load_n>   ]

:data:`number_of_steps` defines the number of steps.
:data:`target_time`, :data:`target_strain<c>` and :data:`target_load<d>` (mutually exclusive) can be used to define the step targets,
where :data:`time` is the cumulative time,
:data:`strain<c>` is the :data:`<c>` component of the engineering strain (among :data:`11`, :data:`22` and :data:`33`),
and
:data:`load<c>` is the load alond direction :data:`<d>` (among :data:`1`, :data:`2` and :data:`3`).
If the number of steps is undefined as the targets are specified, it is set to the number of target values.

A step is divided into (time) *increments*.  There can be any number of increments within a step, and it can be defined from different variables; however, it must be defined using the same variable for all steps. (For conveninence, if the number of provided values for a parameter is smaller than the number of steps, then the last specified value is used for all next steps.)  The format for specifying increments is as follows.

  .. code::

     [ number_of_incrs <num_incrs_1> ... <num_incrs_n> ]
     [ dtime           <dtime_1>     ... <dtime_n>     ]
     [ dtime_min       <dtime_min_1> ... <dtime_min_n> ]
     [ dstrain         <dstrain_1>   ... <dstrain_n>   ]

All variables are mutually exclusive.
:data:`number_of_incrs` corresponds to the number of increments,
:data:`dtime` corresponds to a time increment,
:data:`dtime_min` corresponds to a minimum time increment (needed by :data:`target_load`),
and
:data:`dstrain` corresponds to a strain increment.

Finally, it is possible to specify whether the results should be printed (or not) at the end of the steps (optional, default :data:`1`):

  .. code::

     [ print_results <print_flag_1> ... <print_flag_n> ]

.. **Revise below**
.. 
.. Strain rate jumps are also available for both uniaxial deformation modes and are defined by adding the following input to the block:
.. 
.. ::
.. 
..     number_of_strain_rate_jumps <njumps>
..     strain_rate_jump <target_step> <new_strain_rate>
..     ...
.. 
.. where :data:`<njumps>` is the number of strain rate jumps defined in the file after this line, :data:`<target_step>` defines which :data:`target_strain` step is assigned a new strain rate, and :data:`<new_strain_rate>` is the new strain rate to be assigned and has units of [1/s]. In general, and for numerical stability, the strain rate jumps should be of a similar magnitude to the :data:`strain_rate` defined previously.

[LEGACY] Multiaxial
-------------------

.. note:: The features described in this section are legacy and will change in a future version.

Multiaxial loading operates on cubic domains, at either a constant engineering strain rate or constant load rate; however, only specific load targets may be prescribed.
The principal loading directions must be aligned with the coordinate axes of the mesh and the surface face normals should likewise be coincident with the coordinate axis of the mesh. Symmetry boundary conditions (zero normal velocities) are enforced on the three faces of minimal coordinates (:data:`<*0>`), and, in the general case, non-zero normal velocities are applied to the faces of maximal coordinates (:data:`<*1>`). The velocity on the primary control surface is held constant through the simulation (except during a strain rate jump).

Multiaxial loading with a constant strain rate (CSR) is defined as follows:

::

    def_control_by triaxial_constant_strain_rate
    number_of_steps <num_steps>
    target_csr_load <load_x> <load_y> <load_z> <dt_max> <dt_min> <print_flag>
    ...


where :data:`<num_steps>` is the number of CSR load steps that are defined in the file after this line, :data:`<load_x>` is the desired load value to be reached in the :data:`x` direction, :data:`<load_y>` is the desired load value to be reached in the :data:`y` direction, :data:`<load_z>` is the desired load value to be reached in the :data:`z` direction, :data:`<dt_max>` is the maximum time-step value to be used for a given increment, :data:`<dt_min>` is the minimum time-step value to be used for a given increment, and :data:`<print_flag>` allows for the printing (or not) of specific steps. The options available for :data:`<print_flag>` are: :data:`print_data` or :data:`suppress_data`.

Strain rate jumps are available for this deformation mode and are defined by adding the following input to the block:

::

    number_of_strain_rate_jumps <njumps>
    strain_rate_jump <target_step> <new_strain_rate>
    ...

where :data:`<njumps>` is the number of strain rate jumps defined in the file after this line, :data:`<target_step>` defines which :data:`target_csr_load` step is assigned a new strain rate, and :data:`<new_strain_rate>` is the new strain rate to be assigned and has units of [1/s].

Multiaxial loading with a constant load rate (CLR) is defined as follows:

::

    def_control_by triaxial_constant_load_rate
    number_of_steps <num_steps>
    target_clr_load <load_x> <load_y> <load_z> <target_time_incr> <print_flag>
    ...

where :data:`<num_steps>` is the number of CLR load steps that are defined in the file after this line, :data:`<load_x>` is the desired load value to be reached in the :data:`x` direction, :data:`<load_y>` is the desired load value to be reached in the :data:`y` direction, :data:`<load_z>` is the desired load value to be reached in the :data:`z` direction, :data:`<target_time_incr>` is the physical time increment to be reached for the given :data:`target_clr_load` steps for a given load rate, and :data:`<print_flag>` allows for the printing (or not) of specific steps. The options available for :data:`<print_flag>` are: :data:`print_data` or :data:`suppress_data`.

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

Multiaxial boundary conditions are generally consistent across modes, however, the input rate type varies depending on the mode. For both modes, the :data:`loading_direction` defines the primary control direction in which the normal velocities are held constant throughout the simulation.

Multiaxial loading with a constant strain rate (CSR) is defined as follows:

::

    boundary_conditions triaxial
    loading_direction <sample_dir>
    strain_rate <strain_rate>

where :data:`<sample_dir>` is the direction along the positive sample axis in which the sample is loaded, and :data:`<strain_rate>` is the strain rate value in units of [1/s].

Multiaxial loading with a constant load rate (CLR) is defined as follows:

::

    boundary_conditions triaxial
    loading_direction <sample_dir>
    load_rate <load_rate>

where :data:`<sample_dir>` is the direction along the positive sample axis in which the sample is loaded, and :data:`<load_rate>` is the loading rate value in units of [force/s].

.. [#mpc_note] MPCs are implemented by eliminating degrees of freedom between coupled nodes.
