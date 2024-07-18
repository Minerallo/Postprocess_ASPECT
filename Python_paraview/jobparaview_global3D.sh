#!/bin/bash

## SBATCH -A bbkponsm
## #SBATCH -A bbp00039
#SBATCH -A bbk00014

 
module load gcc/9.2.0 openmpi/gcc.9 anaconda3 llvm/9.0.0 paraview

# # srun --partition=large96 --nodes=1 --ntasks-per-node=96 --mem=747000mb --pty bash


# Ensure the cpus-per-task option is propagated to srun commands
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# srun --distribution=block:block --hint=nomultithread pvbatch test_simple_sphere.py
# mpirun --map-by socket:pe=$OMP_NUM_THREADS pvbatch /scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/test_simple_sphere.py
mpirun --map-by socket:pe=$OMP_NUM_THREADS pvbatch /scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/Model_extract_global_3D_auto.py



