#!/bin/bash

#SBATCH --job-name=FFT_all
#SBATCH --output=logs/FFT_all.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jsunseri@ufl.edu
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=40gb
#SBATCH --qos=narayanan

module load python
echo running b1p1 simulation
python FFT_4pcf_script.py 0 0 5 256 0 1 128
echo Finished the 1st density field
python FFT_4pcf_script.py 0 1 5 256 0 1 128
echo Finished the 2nd density field
python FFT_4pcf_script.py 0 2 5 256 0 1 128
echo Finished the 3rd density field
python FFT_4pcf_script.py 0 3 5 256 0 1 128
echo Finished the 4th density field
python FFT_4pcf_script.py 0 4 5 256 0 1 128
echo Finished the 5th density field
python FFT_4pcf_script.py 0 5 5 256 0 1 128
echo Finished the 6th density field
python FFT_4pcf_script.py 0 6 5 256 0 1 128
echo Finished the 7th density field
python FFT_4pcf_script.py 0 7 5 256 0 1 128
echo Finished the 8th density field
python FFT_4pcf_script.py 0 8 5 256 0 1 128
echo Finished 9th and final density field!
