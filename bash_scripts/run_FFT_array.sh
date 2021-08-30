#!/bin/bash

#SBATCH --job-name=FFT_array
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jsunseri@ufl.edu
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --output=FFT_%A_%a.log    # Standard output and error log
#SBATCH --array=0-3
#SBATCH --qos=narayanan

module load python
python FFT_4pcf_script.py ${SLURM_ARRAY_TASK_ID} 5 256 2 1 128 
