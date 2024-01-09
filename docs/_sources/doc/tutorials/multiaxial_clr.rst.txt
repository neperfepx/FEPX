.. _multiaxial_clr:

[LEGACY] Multiaxial Control with Constant Load Rate
===================================================

This example (:file:`tutorials/triaxCLR`) covers the triaxial deformation of a dual phase FCC/BCC polycrystal (phase map shown below) at a constant load rate followed by a dwell episode and subsequent unloading. Loads are applied normal to the surface, maintaining proportional macroscopic load ratios of :math:`-0.375`::math:`-0.625`::math:`1` for the :math:`x`::math:`y`::math:`z` directions, respectively. Load tolerance options are prescribed. The primary loading direction is set to be in the :math:`z` direction, and a dwell episode is initiated on the second step. Elemental critical resolved shear stresses and equivalent strains are output, along with the nodal coordinates and surface-integrated forces. Material parameters are those for the austenitic (:math:`\gamma`) and ferritic (:math:`\alpha`) phases of an LDX-2101 steel, and are provided in the tables below. Illustrations of the results are provided in the figures below.

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
    * - :math:`\alpha`
      - BCC
      - :math:`236.9 \times 10^3`
      - :math:`140.6 \times 10^3`
      - :math:`116.0 \times 10^3`

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
    * - :math:`\gamma,\alpha`
      - 0.02
      - 1.0
      - 391.9
      - 200.0
      - 335.0
      - 1.0

.. figure:: multiaxial_clr/3_all.png
    :width: 75%
    :align: center

    (left) Grain and phase assignment distribution in the virtual sample. Red-colored grains are :math:`\gamma`-phase and green-colored grains are :math:`\alpha`-phase. (right) Elastically unloaded sample colored by critical resolved shear stress.

.. figure:: multiaxial_clr/3_all2.png
    :width: 100%
    :align: center

    (left) Evolution of the macroscopic normal strains, and (right) evolution of the macroscopic stress and strain on the :math:`z` surface. Note that the load rate in the :math:`z` direction is always held constant during the simulation (except during the dwell episode) while the other two are automatically modified to maintain load proportionality throughout the simulation.
