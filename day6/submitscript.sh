#!/bin/bash
#SBATCH -N 1 
#SBATCH --ntasks-per-node=8
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
module load cuda/11.2
module load fftw/3.3.5
##module load nvhpc/22.9
module load spack/0.17
module load nvhpc/22.1-gcc-11.2.0-axka
###module --ignore-cache load "nvhpc/22.9"
###module load nvhpc/22.9
cd $SLURM_SBUMIT_DIR 
burger/burger.exe
