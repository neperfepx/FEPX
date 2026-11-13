.. _tutorials:

Tutorials
=========

The tutorials below are complete simulations. Each page shows the mesh, the material, the input
files (to download), the run and its terminal output, the simulation directory, and the
results, with the commands that make the figures: field images by Neper (``neper -V``, with
POV-Ray) and curves by Python (NumPy and Matplotlib).

The inputs are in the source distribution, :file:`doc/tutorials/<name>/`, and the meshes
in :file:`tests/regression/cases/meshes/` (Neper meshes, msh file version 2.4, see
:ref:`file_formats`). Every tutorial but the cruciform specimen (about 1.5 h on 4 processes) is
run by the test suite, so the inputs shown here are tested; the terminal outputs, directory
trees and figures of the pages are made from the same inputs by the script
:file:`tutorials/tools/make_tutorial_pages.sh` of the source distribution.

Each tutorial is run in its own directory holding :file:`simulation.cfg` and
:file:`simulation.msh` (save the downloaded mesh under that name), and the other input files
named. Length units are mm, stresses MPa, times s.

.. note:: The crystal types of the tutorials' inputs are written in upper case (``FCC``,
   ``BCC``, ``HCP``): FEPX reads them in any case, but the post-processing of Neper 5.0.1
   (``neper -S``) stops on a lower-case crystal type (and on ``bct``); later versions read
   both.

.. toctree::
   :maxdepth: 1

   tutorials/simple_simulation.rst
   tutorials/simple_rve.rst
   tutorials/periodic_rve.rst
   tutorials/cruciform.rst
   tutorials/overwriting.rst
   tutorials/multiaxial_csr.rst
   tutorials/multiaxial_clr.rst
   tutorials/restart.rst
   tutorials/mixed_targets.rst
   tutorials/hcp_single_crystal.rst
   tutorials/postprocessing.rst

Running on a cluster with a job scheduler is described in :ref:`running_cluster`.
