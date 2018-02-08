#!/bin/bash -l
#SBATCH -n 1
#SBATCH --job-name="diffusion"
#SBATCH --output out_diffusion_%j
#SBATCH --partition=slurm_shortgpu,slurm_courtesy
hostname
cd $SLURM_SUBMIT_DIR
./diffusion
