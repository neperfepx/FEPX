.. _restart_tutorial:

[LEGACY] Restarting a Simulation with Appended Load Steps
---------------------------------------------------------

An example use case for the restart capabilities are to append additional loading steps to a completed simulation. This example covers the cyclic triaxial deformation of an FCC polycrystal at a constant strain rate. Each restart simulation adds an addition load-unload cycle. Loads are applied normal to the surface, maintaining proportional macroscopic load ratios of :math:`0`::math:`0`::math:`1` for the :math:`x`::math:`y`::math:`z` directions, respectively. Load tolerance options are prescribed and cyclic hardening is activated. The primary loading direction is set to be in the :math:`z` direction. Elemental equivalent plastic deformation rate, slip system shears, slip system shear rates, nodal coordinates, and restart files are output. Restart files are only printed on the first cycle. Material parameters are those for a AL6XN steel and are provided in the tables below. Illustrations of the results are provided in the figures below.

The included shell script will run the initial simulation in :file:`tutorials/restart/cycle1` which runs 2 load steps (a single load-unload cycle) and prints the files necessary to restart the simulation from the final state (using print option :data:`print restart`). After successful completion of the first cycle, the mesh file and restart files are copied into the secondary directory (:file:`tutorials/restart/cycle2`) and the simulation is performed again for another 2 load steps (a second load-unload cycle). The configuration file for the second cycle (:file:`simulation_cycle2.config`) contains the following input to allow for simulating additional load steps to those already completed:

::

    restart on

along with the additional load steps (as described in the :ref:`Deformation History <deformation_history>` section). The restarted simulation will continue with the load steps as defined in :file:`simulation_cycle2.config`. Restart control information will print to the console upon the execution of the second cycle to briefly assess the state of the sample when the simulation is restarted. Material parameters are defined as:

.. list-table:: Single crystal elastic constants.
    :widths: 25 25 25 25 25
    :align: center
    :header-rows: 1

    * - Phase
      - Type
      - :math:`C_{11}` [MPa]
      - :math:`C_{12}` [MPa]
      - :math:`C_{44}` [MPa]
    * - :math:`\gamma`
      - FCC
      - :math:`204.6 \times 10^3`
      - :math:`137.7 \times 10^3`
      - :math:`126.2 \times 10^3`

.. list-table:: Initial slip system strengths and other plasticity parameters.
    :widths: 10 10 10 10 10 10 10
    :align: center
    :header-rows: 1

    * - Phase
      - :math:`m` [-]
      - :math:`\dot{\gamma_{0}}` [1/s]
      - :math:`h_{0}` [MPa]
      - :math:`g_{0}` [MPa]
      - :math:`g_{s0}` [MPa]
      - :math:`n` [-]
    * - :math:`\alpha`
      - 0.02
      - 1.0
      - 375.0
      - 160.0
      - 1000.0
      - 1.0

.. list-table:: Cyclic hardening parameters.
    :widths: 25 25 25
    :align: center
    :header-rows: 1

    * - Phase
      - :math:`a` [-]
      - :math:`c` [-]
    * - :math:`\alpha`
      - 0.05
      - 3.50

.. figure:: restart/4_all.png
    :width: 75%
    :align: center

    Sample after the second cycle is completed, (left) colored by accumulated slip shear on the :math:`(1 \bar 1 1)[0 1 1]` system and (right) colored by accumulated slip shear on the :math:`(1 \bar 1 1)[1 0 \bar 1]` system.
