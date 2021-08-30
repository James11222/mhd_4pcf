#!/bin/bash

#SBATCH --job-name=FFT_all
#SBATCH --output=logs/FFT_all.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jsunseri@ufl.edu
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=40gb
#SBATCH --qos=narayanan

module load python
python FFT_4pcf_script.py 0 10 256 2 1 128
echo Finished the first density field
python FFT_4pcf_script.py 1 10 256 2 1 128
echo Finished the second density field
python FFT_4pcf_script.py 2 10 256 2 1 128
echo Finished the third density field
python FFT_4pcf_script.py 3 10 256 2 1 128
echo Finished!
