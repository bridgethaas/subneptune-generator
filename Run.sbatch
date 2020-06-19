#!/bin/bash
#SBATCH --account=pi-larogers
#SBATCH --job-name=HaasTest
#SBATCH --time=36:00:00
##SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bhaas@uchicago.edu
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module load Anaconda3/5.1.0

export OMP_NUM_THREADS=16

#./clean
./clear_everything
#./mk
python run_orion.py

