.. _material_behavior:

Crystal Behavior (:file:`simulation.cfg`)
============================================

.. and the corresponding input parameters are made clear and written in fixed font, most often as in: :math:`\gamma` (:data:`gamma`).

The mechanical behavior of the crystal implies the definition of the crystal structure and the definition of the mechanical behavior *per se*.  If the polycrystal (as defined in the mesh file) contains several phases, then the properties of each of the phases must be defined, successively.
The format is as follows::

  # Material Description

    number_of_phases <nphases>

    phase 1

        crystal_type <ctype>
        [ c_over_a <ratio> ]

        { parameter definition }

    ...

    phase <nphases>

        crystal_type <ctype>
        [ c_over_a <ratio> ]

        { parameter definition }


where the crystal type (:data:`<ctype>`) can be :data:`fcc`, :data:`bcc`, :data:`hcp` and :data:`bct`, for face-centered cubic, body-centered cubic, hexagonal close-packed and body-centered tetragonal, respectively. For :data:`hcp` and :data:`bct`, the :math:`c/a` ratio (:data:`c_over_a`) also needs to be provided.
:data:`{ parameter definition }` represents a block of text that defines a mechanical behavior.  The mechanical behavior of a crystal includes :ref:`aniso_elasticity` and :ref:`crystal_plasticity`.

.. note:: All variables presented are detailed by their dimensions (if applicable) instead of any specific unit. No unit system is inherently assumed by FEPX, and the chosen unit system and value magnitudes should be consistent with the chosen length scale for the domain. For example, if it is assumed that the length scale is *mm* and SI units are to be used, then [force/area] will be understood to be *MPa*. The unit for time, however, is always assumed to be seconds (*s*).

.. _aniso_elasticity:

Anisotropic Elasticity
-----------------------

The stress, :math:`\sigma`, is related to the elastic strain, :math:`\epsilon`, via Hooke's law:

.. math::
    \sigma = \cal C \, \epsilon,

where :math:`\cal C` is the stiffness tensor. Using Voigt notation and the strength of materials convention, :math:`\sigma` and :math:`\epsilon` are written as vectors of components :math:`\sigma_{11}`, :math:`\sigma_{22}`, :math:`\sigma_{33}`, :math:`\sigma_{23}`, :math:`\sigma_{13}`, :math:`\sigma_{12}`, and  :math:`\epsilon_{11}`, :math:`\epsilon_{22}`, :math:`\epsilon_{33}`, :math:`2\,\epsilon_{23}`, :math:`2\,\epsilon_{13}`, :math:`2\,\epsilon_{12}`.  :math:`\cal C` is symmetrical and has 21 independent components; however, its structure is greatly simplified by the crystal symmetry,  and the number of independent components (hence necessary parameters) is greatly reduced.  The moduli are expressed in [force/area].

Cubic symmetry
^^^^^^^^^^^^^^

For *cubic* symmetry (:data:`bcc` and :data:`fcc` crystals), Hooke's law expands as:

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

The corresponding input is::

    c11 <modulus>
    c12 <modulus>
    c44 <modulus>

Hexagonal symmetry
^^^^^^^^^^^^^^^^^^

For *hexagonal* symmetry (:data:`hcp` crystals), Hooke's law expands as:

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
    \end{Bmatrix},

and the following must additionally be satisfied: :math:`C_{33} = C_{11} + C_{12} - C_{13}` (which is why no input is expected for :math:`C_{33}`) [#decoupling]_.

The corresponding input is::

    c11 <modulus>
    c12 <modulus>
    c13 <modulus>
    c44 <modulus>

Tetragonal symmetry
^^^^^^^^^^^^^^^^^^^

For *tetragonal* symmetry (:data:`bct` crystals), Hooke's law expands as:

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
    \end{Bmatrix},

and the following must additionally be satisfied: :math:`C_{33} = C_{11} + C_{12} - C_{13}` (which is why no input is expected for :math:`C_{33}`) [#decoupling]_.

The corresponding input is::

    c11 <modulus>
    c12 <modulus>
    c13 <modulus>
    c44 <modulus>
    c66 <modulus>

.. note:: Special attention must be paid to ensure that the correct stiffness values are chosen, to align with the input convention used here.  For this convention, the Zener anisotropy ratio for cubic materials (which quantifies the level of elastic anisotropy, with 1 being perfectly isotropic) would be written as:

  .. math::

      A = {2 C_{44} \over C_{11} - C_{12}}.

  For example, Tungsten (W) is a nearly perfectly elastically isotropic cubic-symmetry material, with :math:`C_{11} = 522.4` GPa, :math:`C_{12} = 204.4` GPa, and :math:`C_{44} = 160.8` GPa, which corresponds to a Zener ratio of :math:`\simeq 1.01`.

.. _crystal_plasticity:

Crystal Plasticity
------------------

Slip kinetics
^^^^^^^^^^^^^

The kinematics of slip are described by a power law:

.. math::

    \dot{\gamma}^{\alpha} = \dot{\gamma}_{0} \left( \left| {\tau}^{\alpha} \right| \over g^{\alpha} \right)^{1/m} \rm sgn({\tau}^{\alpha}),

where :math:`\dot{\gamma}_0` (:data:`gammadot_0`) is the fixed-rate strain rate scaling coefficient (typically 1 by convention, and expressed in [1/s]), and :math:`m` (:data:`m`) is the rate sensitivity exponent.

The corresponding input is::

    m <rate_sensitivity(ies)>
    gammadot_0 <slip rate(s)>

One :data:`m` value may be provided for all slip systems (or families), or several :data:`m` values may be provided, each applying to a specific family of slip systems (see :ref:`slip_systems`).

Hardening
^^^^^^^^^

Hardening models of increasing complexity are available, but are all evolutions of a *base model*, which corresponds to a Voce-type hardening law with isotropic hardening (same behavior on all systems).  Each evolution of this model enrich one aspect of the behavior and are introduced compared to this model, but the evolutions can be combined with each other (leading to more general hardening laws). An example will be provided.

.. _base_model:

Base Model
""""""""""

The slip system strength evolution is modeled by:

.. math ::

    \dot{g^{\alpha}} = h_{0} \left (g_s - g^{\alpha} \over g_s - g_0 \right)^{n} \dot{\gamma},
    \quad \hbox{with} \quad
    \dot{\gamma} = \sum_{\alpha} \left|\dot{\gamma}^{\alpha}\right|,

where
:math:`g_0` (:data:`g_0`) is the slip system initial strength (expressed in [force/area]),
:math:`g_s` (:data:`g_s`) is the slip system saturation strength (expressed in [force/area]),
:math:`h_0` (:data:`h_0`) is the fixed-state hardening rate scaling coefficient, and
:math:`n` (:data:`n`) is the non-linear Voce hardening exponent.
:math:`g_0` may be defined by one value for all systems, one value per slip family, or one value per slip system (see :ref:`slip_systems`).  Note that this law implied *isotropic* hardening, while :ref:`anisotropic_hardening` is described in the following.

The corresponding input is:


.. code::

    [ hardening saturation[,isotropic] ]

    g_0 <strength(s)>
    g_s <strength>
    h_0 <strength>
    n   <hardening exponent>

.. _saturation_strength:

Saturation Strength Evolution
"""""""""""""""""""""""""""""

Compared to the base model, the slip system saturation strength evolves as a function of the slip activity (constant :math:`g_s` becomes :math:`g_s(\dot\gamma)`):

.. math ::

    \dot{g^{\alpha}} = h_{0} \left (g_{s}(\dot{\gamma}) - g^{\alpha} \over g_{s}(\dot{\gamma}) - g_{0} \right)^{n} \dot{\gamma},
    \quad \hbox{with} \quad
    \dot{\gamma} = \sum_{\alpha} \left|\dot{\gamma}^{\alpha}\right|,

where :math:`g_{s}(\dot{\gamma})` evolves via:

.. math ::

    g_{s}(\dot{\gamma}) = g_{s0} \left (\dot{\gamma} \over \dot{\gamma}_{s0} \right)^{m'},

which introduces three new variables:
:math:`g_s` (:data:`g_s`) is the initial slip system saturation strength (expressed in [force/area]),
:math:`m'` (:data:`m_prime`) is the saturation strength rate scaling exponent, and
:math:`\dot{\gamma}_{s0}` (:data:`gammadot_s0`) is the initial saturation slip system shear rate (typically 1 by convention, and expressed in [1/s]).

Compared to the base model, replace :data:`g_s` by the following input:

.. code::

    hardening saturation_evolution

    g_s0        <strength>
    gammadot_s0 <slip rate(s)>
    m_prime     <strength>

.. _cyclic_hardening:

Cyclic Hardening
""""""""""""""""

Compared to the base model, :math:`\dot\gamma` becomes :math:`f`, defined as follows [TURKMEN04]_:

.. math ::

    \dot{g^{\alpha}} = h_{0} \left (g_s - g^{\alpha} \over g_s - g_{0} \right)^{n} f,
    \quad \hbox{with} \quad
    f = \sum_{\beta = 0}^{n_{a}} \left|\dot{\gamma}^{\beta}\right|.

So, a slip system that contributes to hardening (:math:`n_{a}` systems in total) is that which has a change in shear greater than a critical value:

.. math ::

    \Delta\gamma_{crit} = a \, \left( \frac{g}{g_s} \right)^c.

:math:`a` (:data:`cyclic_a`) and :math:`c` (:data:`cyclic_c`) are model parameters for a critical value of accumulated shear strain used to modify the form of the Voce hardening law.

Compared to the base model, modify the hardening type and provide the :math:`a` and :math:`c` values::

    hardening cyclic

    cyclic_a    <strength(s)>
    cyclic_c    <strength>

.. _precipitation_hardening:

Precipitation Hardening
"""""""""""""""""""""""

The base slip system strength is modified to consider the effects of the presence of precipitates. This is performed by replacing :math:`g_{0}` by :math:`g_{0} + g_{p}`, where :math:`g_p` is the additional strength contribution due to a precipitate phase:

.. math ::

    \dot{g^{\alpha}} = h_{0} \left (g_s - g^{\alpha} \over g_s - (g_0 + g_p) \right)^{n} \dot{\gamma},
    \quad \hbox{with} \quad
    \dot{\gamma} = \sum_{\alpha} \left|\dot{\gamma}^{\alpha}\right|,

There are two options for precipitate strengthening [COURTNEY90]_. Below a critical average precipitate radius (the "cutting" regime), the contribution to strength is calculated via:

.. math ::

	g_p = a_p \left(\frac{f_p \, r_p}{b_p}\right)^{\frac{1}{2}},

where :math:`f_{p}` (:data:`f_p`) is the precipitate volume fraction, :math:`r_p` (:data:`r_p`) is the average precipitate radius, :math:`b_p` :data:`b_p` is the precipitate's Burgers' vector, and :math:`a_p` (:data:`a_p`) is the cutting scaling coefficient.

Above a critical average precipitate radius (the "bowing" regime), the contribution to strength is calculated via:

.. math ::

	g_p = \frac{c_p \, G_m \, b_m}{L - 2 \, r_p},

where :math:`G_m` is the shear modulus of the matrix, :math:`c_p` (:data:`c_p`) the bowing scaling coefficient, and  :math:`L` is the average center to center spacing of the precipitates, calculated as:

.. math ::

	L = \frac{r_p}{\sqrt{f_p}}.

The increase in strength due to the presence of precipitates is applied globally equally to all elements.

.. _input_parameters:

Compared to the base model, provide the following additional input::

    hardening precipitation

    a_p <strength>
    f_p <volume_fraction>
    r_p <length>
    b_p <length>
    c_p <coefficient>

To disable precipitate cutting, omit :data:`a_p`. To disable precipitate bowing, omit :data:`c_p`.

.. _anisotropic_hardening:

Anisotropic Hardening
"""""""""""""""""""""

Contrary to the case of isotropic hardening, for while that slip on given system generates the same hardening on all systems, and so the single crystal yield surface retains the same shape and grows isotropically [#isotropic_hardening]_, anisotropic hardening is such that slip on a given system generates different hardenings on the systems.

Compared to the base model, the hardening coefficient (:math:`h_0`) is multiplied by the slip interaction matrix, :math:`h_{\alpha\beta}`:

.. math ::

    \dot{g^{\alpha}} = h_{0} \, h_{\alpha \beta} \, \left (g_s - g^{\alpha} \over g_s - g_{0} \right)^{n} \dot{\gamma},

The diagonal entries of the matrix correspond to *self-hardening*, while the off-diagonal entries correspond to *latent hardening*.  The special case of isotropic hardening corresponds to equal self-hardening and latent hardening, i.e. :math:`h_{\alpha \beta} = 1 \, \forall \, \alpha,\beta`.
The coefficients generally depend on the type of interactions between systems, but two specific cases are implemented:

- 2-parameter interaction matrix, where the parameters are the diagonal (self-hardening) value, :math:`d`, and the off-diagonal (latent hardening) value, :math:`h`.

  The corresponding input is::

    hardening anisotropic

    interaction_matrix_parameters <diag> <h>

- Self-hardening + co-planar latent hardening, for which latent hardening is limited to coplanar slip systems (other interaction parameters are 0) [CARSON17]_.  In this case, the slip interaction matrix is defined by the diagonal entry, :math:`d`, and the off-diagonal entries, :math:`h_{1},\dots, h_{n}` (depending on the crystal symmetry / number of slip planes).

  The corresponding input is::

    hardening anisotropic

    interaction_matrix_parameters <diag> <h1> <h2> <h3> <h4>                                            (for fcc)
    interaction_matrix_parameters <diag> <h1> <h2> <h3> <h4> <h5> <h6>                                  (for bcc)
    interaction_matrix_parameters <diag> <h1> <h2> <h3> <h4> <h5> <h6> <h7> <h8> <h9> <h10> <h11> <h12> (for bcc with 112 slip systems considered)
    interaction_matrix_parameters <diag> <h1> <h2> <h3> <h4> <h5> <h6> <h7>                             (for hcp)
    interaction_matrix_parameters <diag> <h1> <h2> <h3> <h4> <h5> <h6> <h7> <h8> <h9> <h10>             (for bct)

.. _coupling_model_evolutions:

Coupling Model Evolutions
"""""""""""""""""""""""""

An example of using several evolutions of the :ref:`base_model` is as follows.  We introduce both :ref:`saturation_strength` and :ref:`anisotropic_hardening`.  The slip law becomes:

.. math::

    \dot{g^{\alpha}} = h_{0} h_{\alpha\beta} \left (g_{s}(\dot{\gamma}) - g^{\alpha} \over g_{s}(\dot{\gamma}) - g_{0} \right)^{n} \dot{\gamma},
    \quad \hbox{with} \quad
    \dot{\gamma} = \sum_{\alpha} \left|\dot{\gamma}^{\alpha}\right|.

where :math:`g_{s}(\dot{\gamma})` is the function for the saturation strength, which evolves via:

.. math ::

    g_{s}(\dot{\gamma}) = g_{s0} \left (\dot{\gamma} \over \dot{\gamma}_{s0} \right)^{m'}.

The corresponding input is::

    hardening saturation_evolution,anisotropic

    g_0                          <strength(s)>
    g_s0                         <strength>
    h_0                          <strength>
    n                            <hardening exponent>
    gammadot_s0                  <slip rate(s)>
    m_prime                      <strength>
    interaction_matrix_parameters <diag> <h>

.. [CARSON17] R. Carson, M. Obstalecki, M. Miller, and P.R. Dawson. Characterizing heterogeneous intragranular deformations in polycrystalline solids using diffraction-based and mechanics-based metrics. *Modelling and Simulation in Materials Science and Engineering*, 25:055008, 2017.

.. [TURKMEN04] H.S. Turkmen, M.P. Miller, P.R. Dawson, and J.C. Moosbrugger. A slip-based model for strength evolution during cyclic loading. *Journal of Engineering Materials and Technology*, 126(4):329-338, 2004. (Note that minor differences exist between the implemented model described above and the formulation described in the paper).

.. [COURTNEY90] Thomas H. Courtney. *Mechanical behavior of materials*. eng. 2. ed. Long Grove, Illinois: Waveland Press, 2005. isbn: 9781577664253.

.. [#decoupling] This allows for the decoupling of the hydrostatic and deviatoric portions of the elastic deformation, in the problem resolution.

.. [#isotropic_hardening] Note that, for crystals with different families of slip systems (HCP and BCT materials), this implies faster (absolute) hardening on the slip systems of larger slip strengths (specifically, for an HCP material, the slip rates on the basal, prismatic, and pyramidal slip systems will harden such that the ratios of slip strengths remain constant).
