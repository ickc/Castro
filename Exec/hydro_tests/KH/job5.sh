#!/bin/bash
#SBATCH -N 32
#SBATCH -C haswell
#SBATCH -p debug
#SBATCH -J Castro-KH-problem-5
#SBATCH -t 00:30:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
srun -n 1024 -c 2 --cpu_bind=cores ./Castro2d.gnu.MPI.ex inputs.2d.p5
