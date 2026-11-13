.. _material_behavior:

Crystal Behavior (:file:`simulation.cfg`)
=========================================

The mechanical behavior of a crystal is defined by its crystal structure and its constitutive
model. If the polycrystal (as defined in the mesh file) contains several phases, each phase is
defined in its own block, opened by a ``phase`` line::

  number_of_phases <nphases>

  phase 1
    crystal_type <ctype>
    model <model>
    [ c_over_a <ratio> ]
    { parameters }

  ...

  phase <nphases>
    ...

The crystal type ``<ctype>`` is ``fcc``, ``bcc``, ``hcp`` or ``bct`` (face-centered cubic,
body-centered cubic, hexagonal close-packed, body-centered tetragonal). The model ``<model>`` is
``elastic`` (anisotropic elasticity only) or ``crystal_plasticity`` (anisotropic elasticity
and crystal plasticity). For ``hcp`` and ``bct``, the :math:`c/a` ratio ``c_over_a`` defines the
slip systems; it is required with ``crystal_plasticity`` and not used by the elastic model.

Units are MPa for moduli and strengths and s for time.

.. confval:: number_of_phases

   The number of phases of the polycrystal. Each phase from 1 to this number is defined by a
   ``phase`` block, and the mesh assigns a phase to every element (or :file:`simulation.phase`
   does, with ``read_phase_from_file``).

   .. include:: generated/keywords/number_of_phases.inc

.. confval:: phase

   Opens the block of phase *i*: the keywords that follow, up to the next ``phase`` line, define
   that phase. Every phase from 1 to ``number_of_phases`` has one block.

   .. include:: generated/keywords/phase.inc

.. confval:: crystal_type

   The crystal structure of the phase: ``fcc`` (face-centered cubic), ``bcc`` (body-centered
   cubic), ``hcp`` (hexagonal close-packed) or ``bct`` (body-centered tetragonal). It sets the
   elastic constants to give and the slip systems (:ref:`crystal`).

   .. include:: generated/keywords/crystal_type.inc

.. confval:: model

   The constitutive model of the phase: ``elastic`` (anisotropic elasticity only) or
   ``crystal_plasticity`` (anisotropic elasticity and crystal plasticity,
   :ref:`crystal_plasticity`).

   .. include:: generated/keywords/model.inc

.. confval:: c_over_a

   The lattice ratio :math:`c/a` of an ``hcp`` or ``bct`` phase, which defines the directions
   and planes of its slip systems. Required with ``model crystal_plasticity``; an ``elastic``
   phase does not use it.

   .. include:: generated/keywords/c_over_a.inc


.. _aniso_elasticity:

Anisotropic Elasticity
----------------------

The Kirchhoff stress :math:`\boldsymbol{\tau}` is related to the elastic strain
:math:`\boldsymbol{\epsilon}^e` by Hooke's law,

.. math::
    \boldsymbol{\tau} = \mathcal{C} \, \boldsymbol{\epsilon}^e,

where :math:`\mathcal{C}` is the stiffness tensor, written in the crystal frame and rotated to
the sample frame with the crystal orientation. The elastic constants are given in the usual
Voigt form: stress components :math:`\sigma_{11}, \sigma_{22}, \sigma_{33}, \sigma_{23},
\sigma_{13}, \sigma_{12}` and strain components :math:`\epsilon_{11}, \epsilon_{22},
\epsilon_{33}, 2\epsilon_{23}, 2\epsilon_{13}, 2\epsilon_{12}`. FEPX uses the full
anisotropic stiffness: the hydrostatic and deviatoric parts are not decoupled. The constants
must form a positive definite stiffness (checked).

Cubic Symmetry
~~~~~~~~~~~~~~

For ``fcc`` and ``bcc`` crystals:

.. math::

    \begin{Bmatrix}
        \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12}
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
        \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\
        2\epsilon_{23} \\ 2\epsilon_{13} \\ 2\epsilon_{12}
    \end{Bmatrix}.

The corresponding input is::

    c11 <modulus>
    c12 <modulus>
    c44 <modulus>

Hexagonal Symmetry
~~~~~~~~~~~~~~~~~~

For ``hcp`` crystals (c-axis along the crystal axis 3):

.. math::

    \begin{Bmatrix}
        \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12}
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
        \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\
        2\epsilon_{23} \\ 2\epsilon_{13} \\ 2\epsilon_{12}
    \end{Bmatrix}.

The corresponding input is::

    c11 <modulus>
    c12 <modulus>
    c13 <modulus>
    c33 <modulus>
    c44 <modulus>

.. note::

   :math:`C_{33}` is an independent constant: the full anisotropic stiffness is used.

Tetragonal Symmetry
~~~~~~~~~~~~~~~~~~~

For ``bct`` crystals:

.. math::

    \begin{Bmatrix}
        \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12}
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
        \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\
        2\epsilon_{23} \\ 2\epsilon_{13} \\ 2\epsilon_{12}
    \end{Bmatrix}.

The corresponding input is::

    c11 <modulus>
    c12 <modulus>
    c13 <modulus>
    c33 <modulus>
    c44 <modulus>
    c66 <modulus>

.. note::

   With this input convention, the Zener anisotropy ratio of a cubic crystal is
   :math:`A = 2 C_{44} / (C_{11} - C_{12})`. For example, tungsten, a nearly isotropic cubic
   crystal, has :math:`C_{11} = 522.4` GPa, :math:`C_{12} = 204.4` GPa and :math:`C_{44} =
   160.8` GPa, hence :math:`A \simeq 1.01`.

The elastic constants:

.. confval:: c11

   Elastic constant :math:`C_{11}` [MPa], in the Voigt form above. Every crystal type.

   .. include:: generated/keywords/c11.inc

.. confval:: c12

   Elastic constant :math:`C_{12}` [MPa]. Every crystal type.

   .. include:: generated/keywords/c12.inc

.. confval:: c13

   Elastic constant :math:`C_{13}` [MPa]. ``hcp`` and ``bct`` only (an error for cubic
   crystals).

   .. include:: generated/keywords/c13.inc

.. confval:: c33

   Elastic constant :math:`C_{33}` [MPa], independent of the others. ``hcp`` and ``bct`` only.

   .. include:: generated/keywords/c33.inc

.. confval:: c44

   Elastic constant :math:`C_{44}` [MPa]. Every crystal type.

   .. include:: generated/keywords/c44.inc

.. confval:: c66

   Elastic constant :math:`C_{66}` [MPa]. ``bct`` only (for ``hcp``, :math:`C_{66} = (C_{11} -
   C_{12})/2`).

   .. include:: generated/keywords/c66.inc


.. _crystal_plasticity:

Crystal Plasticity
------------------

.. note::

   The keywords ``slip_families``, ``m``, ``gammadot_0``, ``tau_0`` and the hardening keywords
   are an error with ``model elastic``.

Slip Systems and Families
~~~~~~~~~~~~~~~~~~~~~~~~~

The slip systems of a crystal are grouped in *families* (see :ref:`crystal`). The active
families are chosen with::

    [ slip_families <family1> <family2> ... ]

where each family is written in Miller notation (Miller--Bravais for ``hcp``) as
``{plane}<direction>``, exactly as listed in :ref:`crystal` (no blanks, same signs and digit
order), or ``all`` for every family of the crystal type. A family may be given once; ``all``
cannot be combined with other families. Without ``slip_families``, the defaults are: ``fcc`` and
``hcp`` all families, ``bcc`` the :math:`\{110\}\langle111\rangle` family only, ``bct`` all
families. Only the systems of the active families are modelled.

**Order.** The slip systems are numbered family after family, in the order in which the
families are listed in ``slip_families``; within a family, in the order of :ref:`crystal`. The
per-family parameters below and the per-system output (e.g. ``crss``, ``sliprate``) follow
this order. With the default families or ``all``, the order is that of :ref:`crystal`.

**Per-family parameters.** The parameters that may differ between slip systems, ``m``,
``gammadot_0`` and ``tau_0`` (and ``h_0``, ``tau_s`` and ``n`` with ``hardening
saturation_anisotropic``), are given either

- with **one value per active family**, in the order of ``slip_families`` (or of the default
  families when it is not given): the first value for the first family listed, and so on; or
- with **one value per slip system**, in the slip-system order.

The rows of ``interaction_matrix`` follow the slip-system order too.

A single value is accepted only when exactly one family is active (it is then the one value per
family). With several families, a value is given for each family, even if they are equal. Any
other number of values is an error that lists the accepted counts:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Active families
     - One per family
     - One per system
   * - ``fcc`` (one family)
     - 1
     - 12
   * - ``bcc`` :math:`\{110\}\langle111\rangle` only
     - 1
     - 12
   * - ``bcc`` both families
     - 2
     - 24
   * - ``hcp`` all families
     - 3
     - 18
   * - ``bct`` all families
     - 10
     - 32

(With fewer active families, the counts are those of the active families.)

For example, a ``bcc`` phase with both families:

.. literalinclude:: ../tests/regression/cases/sx_bcc_112/simulation.cfg
   :start-at: slip_families
   :end-before: hardening
   :dedent: 2

numbers the :math:`\{110\}` systems 1--12 and the :math:`\{112\}` systems 13--24, with initial
strengths 200 MPa and 240 MPa; with ``slip_families {112}<111> {110}<111>``, the
:math:`\{112\}` systems would be 1--12.

.. confval:: slip_families

   The active slip families of the phase: ``all``, or one or more of the families of its
   crystal type below, each written exactly as shown (``{plane}<direction>``) and given once.
   Only the slip systems of these families are modelled. Their order matters: it numbers the
   slip systems (those of the first family first), and the per-family values of :confval:`m`,
   :confval:`gammadot_0` and :confval:`tau_0` (and of :confval:`h_0`, :confval:`tau_s` and
   :confval:`n` with ``saturation_anisotropic``) are given in the same order, one value per
   family listed here. For example, ``slip_families {112}<111> {110}<111>`` with ``tau_0 240
   200`` gives the :math:`\{112\}` systems 240 MPa and the :math:`\{110\}` systems 200 MPa.

   .. include:: generated/keywords/slip_families.inc


Slip Kinetics
~~~~~~~~~~~~~

The slip rate on system :math:`\alpha` follows a power law,

.. math::

    \dot{\gamma}^{\alpha} = \dot{\gamma}_{0}^{\alpha}
    \left( \frac{\left| \tau^{\alpha} \right|}{\tau_{\mathrm{crit}}^{\alpha}} \right)^{1/m^{\alpha}}
    \operatorname{sgn}(\tau^{\alpha}),

where :math:`\tau^{\alpha}` is the resolved shear stress, :math:`\tau_{\mathrm{crit}}^{\alpha}` the
slip system strength (``crss`` in the output), :math:`\dot{\gamma}_{0}^{\alpha}` (``gammadot_0``)
the reference slip rate [1/s] and :math:`m^{\alpha}` (``m``) the rate sensitivity. Both are
required, strictly positive, and given per family or per system as above. A slip system does not
slip while :math:`|\tau^{\alpha}|/\tau_{\mathrm{crit}}^{\alpha}` is at or below a small cutoff that
depends on :math:`m^{\alpha}` (for :math:`m = 0.01` it is about 0.057). The input is::

    m          <rate_sensitivity(ies)>
    gammadot_0 <slip_rate(s)>

.. confval:: m

   The rate sensitivity :math:`m^\alpha` of the slip law: one value per active slip family,
   in the order of :confval:`slip_families`, or one per slip system.

   .. include:: generated/keywords/m.inc

.. confval:: gammadot_0

   The reference slip rate :math:`\dot\gamma_0^\alpha` [1/s] of the slip law: one value per
   active slip family, in the order of :confval:`slip_families`, or one per slip system.

   .. include:: generated/keywords/gammadot_0.inc


Slip System Strengths
~~~~~~~~~~~~~~~~~~~~~

The slip system strengths with no hardening contribution, :math:`\tau_0^{\alpha}` [MPa], are
given per family or per system as above (required, strictly positive)::

    tau_0 <strength(s)>

``tau_0`` is one keyword for every hardening law. With the saturation laws and ``none`` it is the
initial strength. With ``kocks_mecking`` it is the lattice (virgin) strength, the strength with
no dislocations, and the initial strength is :math:`\tau_0^\alpha + b \mu_0 \sqrt{h_n \rho_0}`
(`Kocks--Mecking hardening`_).

.. confval:: tau_0

   The slip system strengths with no hardening contribution, :math:`\tau_0^\alpha` [MPa]: one
   value per active slip family, in the order of :confval:`slip_families`, or one per slip
   system. They are the initial strengths with the saturation laws and ``none``, the lattice
   strengths with ``kocks_mecking``.

   .. include:: generated/keywords/tau_0.inc


Hardening
~~~~~~~~~

The evolution of the strengths is chosen per phase with::

    [ hardening saturation_isotropic|saturation_anisotropic|kocks_mecking|none ]

- ``saturation_isotropic`` (the default) and ``saturation_anisotropic`` are modified Voce laws
  (`Saturation laws`_), which the `Saturation-law modifiers`_ can change;
- ``kocks_mecking`` is a dislocation-density law (`Kocks--Mecking hardening`_), to which
  `Irradiation softening`_ can be added;
- ``none`` keeps the strengths constant (perfect plasticity); ``h_0``, ``tau_s``, ``n`` and
  ``interaction_matrix`` are then an error.

Notation: :math:`\tau^\alpha` is the resolved shear stress of slip system :math:`\alpha` and
:math:`\tau_{\mathrm{crit}}^\alpha` its current strength (critical resolved shear stress, the
output ``crss``); :math:`\tau_0`, :math:`\tau_s`, :math:`\tau_{s0}` and :math:`\tau_p` are
the strengths below.

.. confval:: hardening

   The hardening law of a ``crystal_plasticity`` phase: ``saturation_isotropic`` (the default),
   ``saturation_anisotropic``, ``kocks_mecking`` or ``none`` (constant strengths), described
   below.

   .. include:: generated/keywords/hardening.inc

Saturation Laws
^^^^^^^^^^^^^^^

Both saturation laws are modified Voce laws with the parameters::

    h_0   <rate_coefficient(s)>     # MPa
    tau_s <saturation_strength(s)>  # MPa
    n     <exponent(s)>

.. confval:: h_0

   The hardening rate coefficient :math:`h_0` [MPa] of the saturation laws: one value with
   ``saturation_isotropic``; with ``saturation_anisotropic``, one value (for all systems), one
   per active slip family in the order of :confval:`slip_families`, or one per slip system.

   .. include:: generated/keywords/h_0.inc

.. confval:: tau_s

   The saturation strength :math:`\tau_s` [MPa] of the saturation laws (:math:`\tau_{s0}` with
   ``saturation_evolution on``), given like :confval:`h_0`; it must be greater than ``tau_0`` (plus
   :math:`\tau_p` with precipitates) on every system.

   .. include:: generated/keywords/tau_s.inc

.. confval:: n

   The hardening exponent :math:`n` of the saturation laws, given like :confval:`h_0`.

   .. include:: generated/keywords/n.inc

Isotropic Saturation
""""""""""""""""""""

``saturation_isotropic`` (the default) gives every slip system the same
strength increment, computed from slip system 1:

.. math::

    \dot{\tau}_{\mathrm{crit}}^{\alpha} = h_0 \left( \frac{\tau_s - \tau_{\mathrm{crit}}^{1}}{\tau_s
    - \tau^{1}_0} \right)^{n} \sum_{\beta} \left| \dot{\gamma}^{\beta} \right| \quad \text{for every
      } \alpha,

which stops once :math:`\tau_{\mathrm{crit}}^1` reaches :math:`\tau_s`. ``h_0``, ``tau_s`` and ``n``
take one value each, and ``tau_s`` must be greater than ``tau_0`` of slip system 1. The initial
strengths ``tau_0`` may still differ between families or systems; the systems keep their
differences.

Anisotropic Saturation
""""""""""""""""""""""

``saturation_anisotropic`` evolves each system with its own
strength and parameters, driven by the slip on all systems through the interaction matrix
:math:`H`:

.. math::

    \dot{\tau}_{\mathrm{crit}}^{\alpha} = h_0^{\alpha}
    \left( \frac{\tau_s^{\alpha} - \tau_{\mathrm{crit}}^{\alpha}}{\tau_s^{\alpha} - \tau^{\alpha}_0}
    \right)^{n^{\alpha}} \sum_{\beta} H_{\alpha\beta} \left| \dot{\gamma}^{\beta} \right|,

each system stopping once it reaches its own :math:`\tau_s^{\alpha}`. ``h_0``, ``tau_s`` and ``n``
take one value (for all systems), one per active family, or one per slip system, like
``tau_0``; ``tau_s`` must be greater than ``tau_0`` on every system. :math:`H` is
:math:`N \times N` for the :math:`N` active slip systems: by default 1 on the diagonal (self
hardening) and 1.4 off it (latent hardening), or given in full, one line per row in the
slip-system order (see `Slip systems and families`_), with :math:`N` values :math:`\ge 0` per
line::

    interaction_matrix <H_11> <H_12> ... <H_1N>
    interaction_matrix <H_21> <H_22> ... <H_2N>
    ...

.. confval:: interaction_matrix

   One row of the interaction matrix :math:`H` of ``saturation_anisotropic``: :math:`N` values
   :math:`\ge 0`, on :math:`N` lines for the :math:`N` active slip systems, rows and columns in
   the slip-system order (that of :confval:`slip_families`). By default 1 on the diagonal and
   1.4 off it.

   .. include:: generated/keywords/interaction_matrix.inc


With :math:`H` all ones and uniform parameters and initial strengths, the
anisotropic law gives the isotropic one. ``interaction_matrix`` is an error with
``saturation_isotropic``.

Saturation-Law Modifiers
^^^^^^^^^^^^^^^^^^^^^^^^

Three modifiers change the saturation laws (each ``off`` by default, per phase; an error with
``hardening none``, ``hardening kocks_mecking`` or ``model elastic``). Their parameters are required
when the modifier is ``on`` and are an error when it is ``off``. They can be combined.

Saturation-Strength Evolution
"""""""""""""""""""""""""""""

With ``saturation_evolution on``, the saturation strength evolves with the total slip rate :math:`\dot\Gamma = \sum_\beta |\dot\gamma^\beta|`, one value
for all systems:

.. math::

    \tau_s = \tau_{s0} \left( \frac{\dot\Gamma}{\dot\gamma_{s0}} \right)^{m'}
    \quad (\tau_s = \tau_{s0} \text{ when } \dot\Gamma = 0),

::

    saturation_evolution on
    gammadot_s0 <rate>       # 1/s, > 0
    m_prime     <exponent>   # > 0

.. confval:: saturation_evolution

   ``on``: the saturation strength evolves with the total slip rate (above).

   .. include:: generated/keywords/saturation_evolution.inc

.. confval:: gammadot_s0

   The reference slip rate :math:`\dot\gamma_{s0}` [1/s] of the saturation-strength evolution.

   .. include:: generated/keywords/gammadot_s0.inc

.. confval:: m_prime

   The exponent :math:`m'` of the saturation-strength evolution.

   .. include:: generated/keywords/m_prime.inc


``tau_s`` (:math:`\tau_{s0}`) then takes one value. A system stops hardening while its strength is
at or above the current :math:`\tau_s` (it never softens), so at low slip rates (:math:`\tau_s`
below :math:`\tau_0`) nothing hardens.

Cyclic Hardening
""""""""""""""""

With ``cyclic_hardening on`` (after Turkmen et al. 2004), only slip systems whose slip accumulated since the last load reversal, :math:`A^\beta = \sum |\dot\gamma^\beta|
\Delta t`, has reached a critical value drive hardening, and their drive is kept between
:math:`10^{-6}` and :math:`10` times the equivalent deformation rate :math:`d_{eq}` of the point::

    cyclic_hardening on
    cyclic_a <a>   # 0 < a < 1
    cyclic_c <c>   # > 0

.. confval:: cyclic_hardening

   ``on``: only the slip systems whose slip since the last load reversal has reached a critical
   value drive hardening (below).

   .. include:: generated/keywords/cyclic_hardening.inc

.. confval:: cyclic_a

   The coefficient :math:`a` of the critical accumulated slip :math:`a
   (\tau_{\mathrm{crit}}/\tau_s)^c` (:math:`0 < a < 1`).

   .. include:: generated/keywords/cyclic_a.inc

.. confval:: cyclic_c

   The exponent :math:`c` of the critical accumulated slip.

   .. include:: generated/keywords/cyclic_c.inc


With ``saturation_isotropic``, :math:`\dot\Gamma` is replaced by

.. math::

    f = \operatorname{clamp}\Big( \sum_{\beta:\, A^\beta \ge a (\tau_{\mathrm{crit}}^1/\tau_s)^c}
    |\dot\gamma^\beta| ,\ 10^{-6} d_{eq},\ 10\, d_{eq} \Big);

with ``saturation_anisotropic``, each system is judged
by its own strength and contributes to the interaction sum only past its threshold:

.. math::

    f^\alpha = \operatorname{clamp}\Big( \sum_{\beta:\, A^\beta \ge a
    (\tau_{\mathrm{crit}}^\beta/\tau_s^\beta)^c} H_{\alpha\beta} |\dot\gamma^\beta| ,\ 10^{-6}
    d_{eq},\ 10\, d_{eq} \Big).

:math:`A^\beta` is reset to zero at the start of every step that reverses the strain rate of the
control component.

Precipitate Strengthening
"""""""""""""""""""""""""

With ``precipitate_strengthening on``, a strength :math:`\tau_p` is added to the initial strength :math:`\tau_0` of every system (also in the hardening laws), from
precipitate cutting, bowing, or the smaller of the two when both parameter sets are given::

    precipitate_strengthening on
    f_p <volume_fraction>   # 0 < f_p <= 1, required
    r_p <radius>            # mm, required
    a_p <coefficient>       # MPa   } cutting
    b_p <burgers_vector>    # mm    }
    c_p <coefficient>       # -     } bowing
    b_m <burgers_vector>    # mm    }

.. confval:: precipitate_strengthening

   ``on``: a precipitate strength :math:`\tau_p` (cutting, bowing, or the smaller of the two) is
   added to the initial strength of every system (below).

   .. include:: generated/keywords/precipitate_strengthening.inc

.. confval:: f_p

   The volume fraction :math:`f_p` of the precipitates.

   .. include:: generated/keywords/f_p.inc

.. confval:: r_p

   The precipitate radius :math:`r_p`.

   .. include:: generated/keywords/r_p.inc

.. confval:: a_p

   The coefficient :math:`a_p` [MPa] of precipitate cutting (given with ``b_p``).

   .. include:: generated/keywords/a_p.inc

.. confval:: b_p

   The Burgers vector :math:`b_p` of precipitate cutting (given with ``a_p``).

   .. include:: generated/keywords/b_p.inc

.. confval:: c_p

   The coefficient :math:`c_p` of precipitate bowing (given with ``b_m``).

   .. include:: generated/keywords/c_p.inc

.. confval:: b_m

   The Burgers vector :math:`b_m` of the matrix, for precipitate bowing (given with ``c_p``).

   .. include:: generated/keywords/b_m.inc


.. math::

    \tau_p^{cut} = a_p \left( \frac{f_p\, r_p}{b_p} \right)^{1/2}, \qquad
    \tau_p^{bow} = \frac{c_p\, G_m\, b_m}{L - 2 r_p}, \quad L = r_p \sqrt{\pi / f_p},

with :math:`G_m` the Voigt--Reuss--Hill average shear modulus of the phase's elastic constants
(Hill 1952). Only ratios of lengths enter, so any length unit may be used consistently.
``tau_s`` must exceed ``tau_0`` plus :math:`\tau_p`. Strengths read from
:file:`simulation.opt` are used as given, without :math:`\tau_p` (`Per-element values`_).

Kocks--Mecking Hardening
^^^^^^^^^^^^^^^^^^^^^^^^

``hardening kocks_mecking`` evolves one network dislocation density
:math:`\rho_n` per material point, in the Kocks--Mecking form of Barton et al. (2013, J. Mech.
Phys. Solids 61, 341, Eqs. 3 and 6, without the irradiation term) written in the physical density
(Beyerlein and Tomé 2008, Int. J. Plast. 24, 867, Eq. 3.6), and computes the strengths from it:

.. math::

    \tau_{\mathrm{crit}}^{\alpha} = \tau_0^{\alpha} + b\, \mu_0 \sqrt{h_n\, \rho_n}, \qquad
    \dot\rho_n = \dot\Gamma \left( k_1 \sqrt{\rho_n} - k_2\, \rho_n \right), \qquad
    k_2 = k_{20} \left( \frac{\dot\gamma_{k0}}{\dot\Gamma} \right)^{1/n},

with :math:`\dot\Gamma = \sum_\beta |\dot\gamma^\beta|` and :math:`\rho_n = \rho_0` at the
start. The first term of :math:`\dot\rho_n` is dislocation storage, the second dynamic recovery;
the density tends to :math:`(k_1/k_2)^2`, from below or above. The same strength increment is
added to every system; ``tau_0`` (:math:`\tau_0^\alpha`, here the lattice or virgin strength:
the strength with no dislocations) may differ by family or system. The
parameters, one value each per phase, all required::

    hardening kocks_mecking
    km_b           <burgers_vector>      # b, mm, > 0
    km_h_n         <coefficient>         # h_n, -, >= 0
    km_rho_0       <initial_density>     # rho_0, 1/mm^2, > 0
    km_k1          <storage>             # k_1, 1/mm, >= 0
    km_k20         <recovery>            # k_20, -, >= 0
    km_gammadot_k0 <reference_rate>      # 1/s, > 0
    km_n           <rate_exponent>       # n (k_2 uses 1/n), > 0

.. confval:: km_b

   The Burgers vector :math:`b` [mm] of the Kocks--Mecking strength.

   .. include:: generated/keywords/km_b.inc

.. confval:: km_h_n

   The coefficient :math:`h_n` of the network-density term of the strength.

   .. include:: generated/keywords/km_h_n.inc

.. confval:: km_rho_0

   The initial network dislocation density :math:`\rho_0` [1/mm²].

   .. include:: generated/keywords/km_rho_0.inc

.. confval:: km_k1

   The storage coefficient :math:`k_1` [1/mm].

   .. include:: generated/keywords/km_k1.inc

.. confval:: km_k20

   The recovery coefficient :math:`k_{20}` of :math:`k_2`.

   .. include:: generated/keywords/km_k20.inc

.. confval:: km_gammadot_k0

   The reference slip rate :math:`\dot\gamma_{k0}` [1/s] of :math:`k_2`.

   .. include:: generated/keywords/km_gammadot_k0.inc

.. confval:: km_n

   The rate exponent :math:`n` of :math:`k_2` (which uses :math:`1/n`).

   .. include:: generated/keywords/km_n.inc


``km_rho_0`` is the dislocation density the run starts from, so the initial strengths are
:math:`\tau_0^\alpha + b \mu_0 \sqrt{h_n \rho_0}`, not ``tau_0``. :math:`\mu_0` is not an input:
it is the Voigt--Reuss--Hill average shear modulus of the phase's elastic constants (Hill 1952),
written to the log with the initial strength increment. :math:`k_1` carries a length unit:
in mm\ :sup:`-1` it is 1/1000 of its value in m\ :sup:`-1`.

Barton et al. write the law in :math:`h = \rho_n/\rho_0^B` with a reference density
:math:`\rho_0^B` ("an arbitrary dislocation density in specified units") and
:math:`\dot h = \dot\Gamma (k_1^B \sqrt{h} - k_2 h)`. Their parameters convert as
:math:`k_1 = k_1^B \sqrt{\rho_0^B}` and :math:`\rho_0 = \rho_0^B h(0)`; :math:`k_2` is
unchanged. Their Table 1 and §3.2 (iron) give, in FEPX's units (:math:`b` = 2.48 Å;
:math:`k_1 = k_1^B \sqrt{\rho_0^B} = 450 \sqrt{10^6\ \mathrm{mm}^{-2}}`; :math:`\rho_0` = 2e13 m\ :sup:`-2`, Barton's start,
20 times his :math:`\rho_0^B` of 1e12 m\ :sup:`-2`):

.. literalinclude:: ../tests/regression/cases/km_sx_bcc/simulation.cfg
   :start-at: km_b
   :end-at: km_rho_0
   :dedent: 2

(the exponent :math:`1/n` of :math:`k_2` is missing from the printed Eq. 6; Table 1 lists
:math:`n`). No modifier may be used with this law, nor ``h_0``, ``tau_s``, ``n`` or
``interaction_matrix``. :math:`\rho_0` may be given per element set or per element
(`Per-element values`_). The outputs ``km_rho``, ``km_k2`` and ``km_rhorate``
(:ref:`simulation_results`) give the density and the recovery coefficient.

Irradiation Softening
^^^^^^^^^^^^^^^^^^^^^

``irradiation on`` (only with ``hardening kocks_mecking``) adds the
defect term of Barton et al. (2013, J. Mech. Phys. Solids 61, 341, Eqs. 3--5) to the
Kocks--Mecking strengths. Irradiation produces dislocation loops; they are described by a
symmetric damage tensor :math:`\mathbf H` [1/mm²] per material point, a density of loop
perimeter, which raises the strengths of the slip systems whose planes the loops cross. Slip
sweeps the loops away (annihilation), so the strengths drop as the deformation proceeds: the
stress--strain curve shows a peak, a drop and, once the defects are cleared, renewed hardening
by the network density.

.. math::

    \tau_{\mathrm{crit}}^{\alpha} = \tau_0^{\alpha} + b\, \mu_0 \left( \sqrt{h_n\, \rho_n}
      + \sqrt{h_d \left\langle \mathbf N^{\alpha} : \mathbf H \right\rangle_+} \right), \qquad
    \dot{\mathbf H} = -\eta \sum_\alpha \left( \mathbf N^{\alpha} : \mathbf H \right)
      \mathbf N^{\alpha} \left| \dot\gamma^{\alpha} \right|,

with :math:`\mathbf N^\alpha = \mathbf n^\alpha \otimes \mathbf n^\alpha` (the unit normal of
the slip plane of system :math:`\alpha`) and :math:`\rho_n` evolving as above. The loops move
with the lattice: :math:`\mathbf H` is kept in the crystal frame. At the start,

.. math::

    \mathbf H_0 = 3\, d_l\, \rho_l\; \frac{1}{N_\text{dirs}}
      \sum_{i=1}^{N_\text{dirs}} \left( \mathbf I - \mathbf n_i \otimes \mathbf n_i \right)

(Barton Eq. 4: hexagonal loops of diameter :math:`d_l`, perimeter
:math:`3 d_l`, number density :math:`\rho_l`, equally distributed over the
:math:`N_\text{dirs}` loop normals :math:`\mathbf n_i` of the listed families). The parameters,
one value each per phase, all required with ``irradiation on`` and rejected without it::

    irradiation        on
    irr_h_d            <coefficient>         # h_d, defect hardening coefficient, -, >= 0
    irr_eta            <efficiency>          # eta, loop annihilation efficiency, -, >= 0
    irr_loop_density   <number_density>      # rho_l, 1/mm^3, >= 0
    irr_loop_diameter  <diameter>            # d_l, mm, > 0
    irr_loop_normals   <family> ...          # <111>, <100>, <110> (one or more)

.. confval:: irradiation

   ``on``: the irradiation defect term (dislocation loops) is added to the Kocks--Mecking
   strengths (above).

   .. include:: generated/keywords/irradiation.inc

.. confval:: irr_h_d

   The defect hardening coefficient :math:`h_d`.

   .. include:: generated/keywords/irr_h_d.inc

.. confval:: irr_eta

   The loop annihilation efficiency :math:`\eta`.

   .. include:: generated/keywords/irr_eta.inc

.. confval:: irr_loop_density

   The number density :math:`\rho_l` [1/mm³] of the loops.

   .. include:: generated/keywords/irr_loop_density.inc

.. confval:: irr_loop_diameter

   The loop diameter :math:`d_l` [mm].

   .. include:: generated/keywords/irr_loop_diameter.inc

.. confval:: irr_loop_normals

   The families of the loop normals: one or more of ``<111>``, ``<100>``, ``<110>``.

   .. include:: generated/keywords/irr_loop_normals.inc


Barton's iron (Table 1, §3.2) in FEPX's units: ``irr_h_d 0.675``, ``irr_eta 100``,
``irr_loop_diameter 2.48e-5`` (100 b), ``irr_loop_normals <111>``, and loop densities of
8.15e11 to 1.63e13 mm\ :sup:`-3` (8.15e20 to 1.63e22 m\ :sup:`-3`). The loop-normal families
are written in Miller notation, each normal counted once for a :math:`\pm` pair (4, 3 and 6
normals); every listed normal is equally populated. For cubic crystals every family gives an
isotropic :math:`\mathbf H_0 = 2 \rho_l d_l \mathbf I`, so the initial defect
term is the same on every system. Irradiation is available for fcc and bcc phases; for hcp and
bct it is an error (their loop normals are not defined yet).

Barton's evolution law does not keep every contact density :math:`\mathbf N^\alpha : \mathbf H`
positive: slip on one plane also lowers the contact density of every plane not orthogonal to
it, and under multiple slip that of a latent plane can pass zero. FEPX keeps the law as
published and uses the positive part :math:`\langle x \rangle_+ = \max(x, 0)` in the strength:
a system whose contact density has reached zero has no defect term (FEPX's reading of Barton's
Eq. 3, which is undefined there). The output ``dmg`` shows :math:`\mathbf H` as computed.

:math:`\rho_l` and :math:`d_l` may be given per element set or per element
(`Per-element values`_). A zero loop density gives exactly Kocks--Mecking hardening. The
outputs ``dmg`` and ``dmgrate`` (:ref:`simulation_results`) give :math:`\mathbf H` and its rate;
the defect term of a system follows from ``crss`` and ``km_rho`` as
:math:`\tau_{\mathrm{crit}}^\alpha - \tau_0^\alpha - b \mu_0 \sqrt{h_n \rho_n}`.

Time Integration
^^^^^^^^^^^^^^^^

Every hardening law is integrated the same way. A forward-Euler estimate is followed, in each
iteration of the point solve, by a backward-Euler update at the current slip rates, solved with
Newton's method. The interaction between strengths, stress and slip rates is resolved by the
point solve's iterations (see :ref:`solver_parameters`, ``point_max_iters_hardening``). The
scheme is first-order accurate in the time increment, for every law (verification test
``hardening_integration``).

**Saturation laws**: the update is in the strengths :math:`\tau_{\mathrm{crit}}^\alpha`.

**Kocks--Mecking**: the update is in the variable :math:`s = \sqrt{\rho_n}`, in which the law is
linear: :math:`\dot s = (\dot\Gamma/2)(k_1 - k_2 s)`. The update is then exact in one Newton step
and the density stays positive for any time increment (in :math:`\rho_n` itself, Newton's method
fails for low densities and slip increments of 0.5--2 %). The stored and written quantity is
:math:`\rho_n`. The update's tolerance is ``point_tol_hardening`` times :math:`\sqrt{\rho_n}` at
the start of the increment.

**Irradiation**: :math:`\mathbf H` is updated with :math:`\rho_n`. The update is linear in
:math:`\mathbf H`, so it is exact in one Newton step for any time increment; its tolerance is
``point_tol_hardening`` times the largest component of :math:`\mathbf H` at the start of the
increment. The forward-Euler estimate that starts the point solve is replaced by the start
value where :math:`\Delta t\, \eta\, \dot\Gamma \ge 1` (forward Euler would overshoot the
contact density of an active system past zero and start the iteration far from the solution).
The state loop of a phase with irradiation is damped (:ref:`solver_parameters`).

.. _per_element_values:

Per-Element Values
~~~~~~~~~~~~~~~~~~

With :confval:`read_from_opt_file`, some parameters of the crystal-plasticity phases are read per
element set or per element from :file:`simulation.opt` (format: :ref:`opt_file`). They replace,
element by element, the phase's values given in the configuration, which stay required. The
variables read depend on the phase's hardening law:

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Section of :file:`simulation.opt`
     - Replaces
     - Read for phases with
   * - ``Crss``
     - :confval:`tau_0`
     - ``saturation_isotropic``, ``saturation_anisotropic``
   * - ``CrssSat``
     - :confval:`tau_s` (:math:`\tau_{s0}` with ``saturation_evolution on``)
     - ``saturation_isotropic``, ``saturation_anisotropic``
   * - ``HardCoeff``
     - :confval:`h_0`
     - ``saturation_isotropic``, ``saturation_anisotropic``
   * - ``Rho0``
     - :confval:`km_rho_0`
     - ``kocks_mecking``
   * - ``IrrLoopDensity``
     - :confval:`irr_loop_density`
     - ``kocks_mecking`` with ``irradiation on``
   * - ``IrrLoopDiameter``
     - :confval:`irr_loop_diameter`
     - ``kocks_mecking`` with ``irradiation on``

Each section comes per element set (``$Elset<variable>``) or per element
(``$Element<variable>``), not both. Values given for phases that do not read the variable
(for example ``Crss`` for a ``kocks_mecking`` phase, or any value for an elastic phase) are
ignored, with a warning. The irradiation variables are optional: a variable not in the file takes
the configuration value for every element.

Strengths read from the file are used as given, without the precipitate strength
:math:`\tau_p` (a warning is written). With ``initialization isotropic_viscoplastic``, the
phase's yield-surface vertices are computed from the configuration's ``tau_0``, an approximate
start for elements with other strengths.

.. _crystal:

Crystal Structures
------------------

The slip families and slip systems of each crystal type (``crystal_type``), as used by
``slip_families`` and numbered in the per-system input and output.

.. _slip_systems:

Slip Families
~~~~~~~~~~~~~

The slip systems are grouped in families. The tokens below are those accepted by
``slip_families`` (see :ref:`crystal_plasticity`), written exactly as shown; each names its
family in Miller (``hcp``: Miller--Bravais) notation. The default families are those used
without ``slip_families``.

.. list-table::
   :header-rows: 1
   :widths: 12 52 18 18

   * - Crystal
     - Families (tokens, in this order)
     - Systems
     - Default
   * - ``fcc``
     - ``{111}<110>``
     - 12
     - all
   * - ``bcc``
     - ``{110}<111>``, ``{112}<111>``, ``{123}<111>``
     - 12 + 12 + 24
     - ``{110}<111>``
   * - ``hcp``
     - ``{0001}<11-20>`` (basal), ``{10-10}<11-20>`` (prismatic), ``{10-11}<11-23>``
       (pyramidal :math:`\langle c+a \rangle`)
     - 3 + 3 + 12
     - all
   * - ``bct``
     - ``{100}<001>``, ``{110}<001>``, ``{100}<010>``, ``{110}<1-11>``, ``{110}<1-10>``,
       ``{100}<011>``, ``{001}<100>``, ``{001}<110>``, ``{101}<10-1>``, ``{121}<-101>``
     - 2 + 2 + 2 + 4 + 2 + 4 + 2 + 2 + 4 + 8
     - all

The ``bct`` families are those of :math:`\beta`-Sn. The ``bcc`` family ``{123}<111>`` has the
24 :math:`\{123\}` planes, each with the one :math:`\langle 111 \rangle` direction it contains
(pencil glide on all three families, as in Barton et al., 2013, J. Mech. Phys.
Solids 61, 341).

Slip Systems
~~~~~~~~~~~~

The slip systems of every crystal type, family after family in the order of the table above,
with their plane and direction in Miller (``hcp``: Miller--Bravais) indices. These are the
systems and their order (``bcc`` systems 25--48: the ``{123}<111>`` family). With ``slip_families``,
only the systems of the listed families are used, family after family in the order of
``slip_families``, each family in the order below; the numbering of the per-system input and output
follows (see :ref:`crystal_plasticity`).

.. include:: slip_systems_table.rst
