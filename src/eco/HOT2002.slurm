#!/bin/bash
#SBATCH --partition=compute 
#SBATCH --job-name=Palmer 
#SBATCH --mail-type=END # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1 
#SBATCH --mem=1gb 
#SBATCH --time=00:20:00 
#SBATCH --output=Palmer_%j.out

module purge
module load pgi
module load openmpi
module load netcdf
#make clean
#make
./test_adjoint <input
./adjoint_driver <input
./hessian_driver <input_hess
