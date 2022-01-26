---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**

Provide a clear and concise description of what the bug is.

**To Reproduce**

Provide a *minimal working example* that can be used to reproduce the bug, that is, *complete* but *excluding unnecessary options*, *large inputs* and *excessive computation times*.  The example should include all necessary input files, and the terminal output.

*Formatting is important*.  As example should be written as *plain text* (no screenshot, etc.) and provided between triple backticks (Markdown formatting). Input or output text files that are too large to be provided in-text can be attached (GitHub requires a .txt extension).  Images can be used if necessary.

To preview your message, use "Preview".  To edit your message after it is published, use "Edit".

Here is an example of properly formatted example:

```
==========================     F   E   P   X   ==========================
Info   : A finite element software package for  polycrystal plasticity.
Info   : Version 1.2.1
Info   : Running on 4 cores.
Info   : <https://fepx.info>
Info   : Copyright (C) 1996-2021, DPLab, ACME Lab.
Info   : ---------------------------------------------------------------
Info   : Start time: 2022-1-7 at 16:16
Info   : Loading simulation...
Info   :   [i] Parsing file `simulation.config'...
Info   :   - Material parameters:
Info   :     > Number of phases: 1
Info   :     > phase 1 - crystal type:   FCC
Info   :     > m:               0.500000E-01
Info   :     > gammadot_0:      0.100000E+01
Info   :     > h_0:             0.200000E+03
Info   :     > g_0:             0.210000E+03
Info   :     > g_s0:            0.330000E+03
Info   :     > n:               0.100000E+01
Info   :     > c11:             0.245000E+06
Info   :     > c12:             0.155000E+06
Info   :     > c44:             0.625000E+05
Info   :   - Boundary conditions:
Info   :     > Uniaxial grip
Info   :     > Strain targeting, constant strain rate
Info   :     > Strain rate:   0.100000E-01
Info   :     > Loading direction: z
Info   :     > Loading face: z1
Info   :   - Deformation history:
Info   :     > Number of strain steps: 2
Info   :   [i] Parsed file `simulation.config'.
Info   :   [i] Parsing file `simulation.msh'...
Info   :   - Mesh parameters:
Info   :     > Node number: 8737
Info   :     > Elt  number: 5641
Info   :   [i] Parsed file `simulation.msh'.
Info   : Initializing simulation...
Info   :   - Initializing parallel gather/scatter routines...
Info   :   - Initializing fields from isotropic viscoplastic solution
Info   :     > ITMETHOD_VP: Iteration 1
Info   :       . Solving NL iteration... R = 0.5282E+00 (254 iters)
Info   :     > ITMETHOD_VP: Iteration 2
Info   :       . Solving NL iteration... R = 0.8981E-01 (288 iters)
Info   :     > Converged in 2 iterations
Info   :     > Ready for anisotropic elasto-viscoplastic simulation.
Info   : Running step 1...
Info   :   - Increment 1: t = 0.2500 secs, dt = 0.2500 secs
Info   :     > ITMETHOD_EVPS: Iteration 1
Info   :       . Solving SA iteration... R = 0.2249E+00 (384 iters)
Info   :     > ITMETHOD_EVPS: Iteration 2
Info   :       . Solving SA iteration... R = 0.2161E-01 (361 iters)
Info   :     > ITMETHOD_EVPS: Iteration 3
[...]
Info   : Elapsed time:   148.898 secs.
Info   : Final step terminated. Simulation completed successfully.
========================================================================
```

You may cut parts of the terminal output as long as all information needed for debugging is not removed.
