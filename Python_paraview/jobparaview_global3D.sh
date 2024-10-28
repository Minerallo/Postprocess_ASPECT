#!/bin/bash

#SBATCH -A bbk00014
module load HRZIBenv sw.clx.el9 slurm anaconda3/2023.09

module unload openmpi/gcc/5.0.3

ml impi

/sw/comm/impi/mpi/2021.13/bin/mpirun -np 96 /sw/viz/paraview/x86_64.el9/ParaView-headless-5.13.1-egl-MPI-Linux-Python3.10-x86_64/bin/pvbatch /scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/auto_V07d/Model_extract_global_3D_auto.py