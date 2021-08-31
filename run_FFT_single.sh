#!/bin/bash

#SBATCH --job-name=FFT_single
#SBATCH --output=logs/FFT_single.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jsunseri@ufl.edu
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=40gb
#SBATCH --qos=narayanan

module load python
python FFT_4pcf_script_odd_parity.py 0 10 256 2 1 128
echo Finished!
