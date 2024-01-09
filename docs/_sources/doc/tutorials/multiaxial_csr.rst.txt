.. _multiaxial_csr:

[LEGACY] Multiaxial Control with Constant Strain Rate
=====================================================

This example (:file:`tutorials/triaxSLR`) covers the biaxial deformation of an HCP polycrystal at a constant strain rate. Loads are applied normal to the surface, maintaining proportional macroscopic load ratios of :math:`-1`::math:`0`::math:`1` for the :math:`x`::math:`y`::math:`z` directions, respectively. Load tolerance options are prescribed and latent hardening is enabled, as well as saturation strength evolution. The primary loading direction is set to be in the :math:`x` direction, and the strain rate is doubled on the second step. Elemental values for the equivalent plastic strain and plastic work are output, along with the nodal coordinates and surface-integrated forces. Material parameters are those for the :math:`\alpha` phase of Ti-6Al-4V and are provided in the below tables. The latent parameters are input values to the hardening interaction matrix [CARSON17]_. Illustrations of the results are provided in the below figures.

.. list-table:: Single crystal elastic constants.
    :widths: 25 25 25 25 25 25
    :align: center
    :header-rows: 1

    * - Phase
      - Type
      - :math:`C_{11}` [MPa]
      - :math:`C_{12}` [MPa]
      - :math:`C_{13}` [MPa]
      - :math:`C_{44}` [MPa]
    * - :math:`\alpha`
      - HCP
      - :math:`169.66 \times 10^3`
      - :math:`88.66 \times 10^3`
      - :math:`61.66 \times 10^3`
      - :math:`42.50 \times 10^3`

.. list-table:: Plasticity parameters.
    :widths: 10 10 10 10 10 10 10 10 10
    :align: center
    :header-rows: 1

    * - Phase
      - :math:`m` [-]
      - :math:`\dot{\gamma_{0}}` [1/s]
      - :math:`h_{0}` [MPa]
      - :math:`g_{s0}` [MPa]
      - :math:`m^\prime` [-]
      - :math:`\dot{\gamma_{s}}` [1/s]
      - :math:`n` [-]
      - :math:`c/a` [-]
    * - :math:`\alpha`
      - 0.01
      - 1.0
      - 190.0
      - 530.0
      - 1.1
      - 1.0
      - 1.0
      - 1.587

.. list-table:: Initial slip system strengths and hardening parameters.
    :widths: 10 10 10 10 10 10
    :align: center
    :header-rows: 1

    * - Phase
      - :math:`g_0` (basal) [MPa]
      - :math:`g_0` (prismatic) [MPa]
      - :math:`g_0` (pyramidal) [MPa]
      - :math:`h_{diag}` [-]
      - :math:`h_{1}-h_{7}` [-]
    * - :math:`\alpha`
      - 390.0
      - 468.0
      - 663.0
      - 1.0
      - 1.4

.. figure :: multiaxial_csr/2_all.png
    :width: 75%
    :align: center

    Deformed sample at the end of the second load step (deformation field is exaggerated 10x for illustrative purposes), colored by (left) plastic work (:math:`W^{p}`) and (right) equivalent plastic strain (:math:`\bar\epsilon^{P}`). Note that, unlike the deformed sample in :ref:`uniaxial_ex`, a multiaxial simulation will maintain the orthogonal, planar surfaces throughout the simulation.

.. figure :: multiaxial_csr/2_normalstraintime.png
    :width: 70%
    :align: center

    Evolution of the macroscopic normal strains. Note the strain rate increase corresponding to the strain-rate jump defined for step 2.
