.. _running_a_simulation:

Running a Simulation
====================

First, a simulation may be run serially -- that is, on a single core -- by executing the binary from the terminal:

.. code-block:: console

  $ fepx

However, simulations will generally require more computational resources to run in reasonable time, which can be done using a parallel computer architecture.

A simulation may be run in parallel utilizing MPI:

.. code-block:: console

  $ mpirun -np <N> fepx

where :data:`<N>` refers to the number of MPI processes desired. Note that your local installation of MPI may not utilize :program:`mpirun` and instead an alternative MPI command may be required.

Submitting FEPX to a Job Scheduling Program
-------------------------------------------

Performing simulations on high performance computing clusters typically requires interfacing with a job scheduling program. These programs have a number of directives that are too numerous to define here. For sake of illustration, however, generic submission scripts for the Slurm and Torque job scheduling programs are provided. The generation of these scripts is highly dependent on the local configuration, and you are encouraged to work with the system administrator of your cluster if you are unsure on how to properly build a submission script. Both example scripts below are designed to submit a parallel job of FEPX to 4 nodes with 16 cores per node (i.e., 64 total tasks) to the :data:`main` queue.

A generic submission script, :file:`runslurm.sh`, for submitting to a Slurm scheduler:

.. code-block:: bash

    #!/bin/bash
    #SBATCH -J fepxjob
    #SBATCH -q main
    #SBATCH --ntasks 64
    #SBATCH --ntasks-per-node 16
    #SBATCH -o output.%A
    #SBATCH -e errors.%A

    srun --mpi=pmi2 fepx

would be run by entering the following into the terminal from within the simulation directory:

.. code-block:: console

  $ sbatch runslurm.sh

A generic submission script, :file:`runtorque.sh`, for submitting to a Torque scheduler:

.. code-block:: bash

    #!/bin/bash
    #PBS -N fepxjob
    #PBS -q main
    #PBS -l nodes=4:ppn=16
    #PBS -k oe
    #PBS -j oe

    # Change to the current working directory
    cd $PBS_O_WORKDIR

    # Calculate the total number of cores requested
    NP=`cat $PBS_NODEFILE | wc -l`

    mpirun -np $@{NP@} fepx

would be run by entering the following into the terminal from within the simulation directory:

.. code-block:: console

  $ qsub runtorque.sh

.. _sim_restart:

Restarting a Simulation
-----------------------

A simulation may be restarted only if the restart files were printed as simulation output on the previous run (:ref:`Restart Output <>`). Printing restart files outputs a single :file:`rst<N>.control` file and a :file:`rst<N>.field.core@var<#>` file for each individual core, where :data:`<N>` refers to the restart ID (0 indexing), and :data:`<#>` denotes the ID of the core on which the data is being printed. These restart files must be included in the simulation directory along with the configuration file, the mesh file, and any external files included with the simulation.

A simulation may be restarted by adding the following line to the :file:`simulation.config` file::

  restart on

When a simulation restarts, it will attempt to find the simulation restart files with the highest index, :data:`<N>`. It will write output variable data to a new set of files, :file:`post.<var>.rst<N+1>.core*`, where :data:`<var>` is the requested output variable name and :data:`rst<N+1>` is the restart label applied to all new output variable files. Likewise, if restart files are again printed, their index will increase to :data:`<N+2>` (the previous restart's files will not be overwritten).

A simulation restart must be performed with the same number of cores that were used to run the original simulation. Restart files are always written at the end of a successful step and not at individual increments.

Each restarted simulation is considered a new simulation, albeit with initialized field variables as written in the restart files. Step and increment indices, as well as the simulation time, are all reset to 0.

When restarting a simulation, the prescribed deformation history (:ref:`Deformation History <deformation_history>`) should include only additional steps. The restarted simulation will attempt to follow the entire deformation history as prescribed in the :file:`simulation.config` file, and will not consider steps that were completed in the initial simulation.
