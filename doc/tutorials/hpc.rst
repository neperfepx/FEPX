.. _hpc:

Submitting FEPX to a Job Scheduling Program
-------------------------------------------

Performing simulations on high performance computing clusters typically requires interfacing with a job scheduling program. These programs have a number of directives that are too numerous to define here. For sake of illustration, however, generic submission scripts for the Slurm and Torque job scheduling programs are provided.  Both example scripts below are designed to submit a parallel job of FEPX to 4 nodes with 16 cores per node (i.e., 64 total tasks) to the :data:`main` queue.

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
