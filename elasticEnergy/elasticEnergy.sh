#!/bin/bash -l
#SBATCH -n 1
#SBATCH --job-name="elasticEnergy"
#SBATCH -t 0-0:20:0
#SBATCH --output elasticEnergy_out%j
#SBATCH --partition=slurm_shortgpu,slurm_courtesy
cd $SLURM_SUBMIT_DIR
module load cuda/8.0 
./elasticEnergy
rm elasticEnergy_out*
