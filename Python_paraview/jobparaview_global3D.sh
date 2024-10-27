#!/bin/bash

#SBATCH -A bbk00014

module unload openmpi/gcc/5.0.3

ml impi
mpirun -np --map-by socket:pe=$OMP_NUM_THREADS /sw/viz/paraview/x86_64.el9/ParaView-headless-5.13.1-egl-MPI-Linux-Python3.10-x86_64/bin/pvbatch /scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/Model_extract_global_3D_auto.py