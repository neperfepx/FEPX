.. _tutorials:

Tutorials
=========

.. Several example simulations come pre-packaged to get you started with running simulations. These examples are reference cases to show how a simulation should be built and how FEPX can interface with Neper in order to prepare mesh files as well as post-process a simulation directory. All examples contain the necessary configuration (:file:`simulation.cfg`) and mesh (:file:`simulation.msh`) files, along with a shell script to generate the mesh file directly from Neper. In order to run the provided shell scripts, you must have a configured installation of Neper present on your system. In the following, visualizations of the undeformed and deformed mesh are generated with Neper while graphs are generated with Gnuplot.

.. All examples can be run either in serial or parallel (:ref:`Running a Simulation <running_a_simulation>`), but the included scripts are pre-set to run in parallel on 4 cores, and parallel execution with OpenMPI via :data:`mpirun` is assumed.

.. A polycrystal containing 50 grains is generated via Voronoi tessellation for all examples. Each cell in the tessellation represents a discrete grain in the domain and all grains are volumetrically discretized into finite elements. Visualizations of the tessellated domain (morphology) and the associated finite element mesh are shown below. Length units are assumed to be [mm], thus, all pressure units assumed to be [MPa] for the simulation (including input parameters).

.. toctree::
   :maxdepth: 1

   tutorials/simple_simulation.rst
   tutorials/simple_rve_simulation.rst
   tutorials/periodic_rve_simulation.rst
   tutorials/cruciform_specimen_simulation.rst
   tutorials/overwriting.rst
   tutorials/hpc.rst

Legacy capabilities / tutorials:

.. toctree::
   :maxdepth: 1

   tutorials/restart.rst
   tutorials/multiaxial_csr.rst
   tutorials/multiaxial_clr.rst
