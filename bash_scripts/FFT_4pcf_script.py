import numpy as np
from subprocess import call
import astropy.io.fits as pyf
import sys
import os
import _4PCF_Algorithm as fourpcf

print("Importing necessary libraries and packages")

############################################################################

# an example of what you should type in the terminal is
#python3 FFT_4pcf_script.py run_name_index nbins resolution ell_max bin_min bin_max

############################################################################

# run_name = 'dens_t800' #take in from command line arg
# file_name = run_name +'.fits.gz'

file_list = ['dens_t800_b.1p.1', 'dens_t800_b.1p1', 'dens_t800_b1p.1', 'dens_t800_b1p1']
# file_list = ['dens_t800_0', 'dens_t800_1', 'dens_t800_2', 'dens_t800_3']
# data_dir = 'MHD_Sims_Data/'
run_name = file_list[int(sys.argv[1])]
print(sys.argv[1] + '--' + run_name)
file_name =os.getcwd() + '/data/' + run_name +'.fits.gz'

hdulist = pyf.open(file_name)
data = hdulist[0].data.astype(np.float64)

#we only care about excess and deficit
data = data - np.mean(data)

############################################################################


nbins, resolution, ell_max, bin_min, bin_max = int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])

# make sure to create an FFT_Files directory inside your current directory of running this
save_dir = os.getcwd() + '/FFT_Files/'

save_name = 'FFT_' + run_name + '_' 

#create _4pcf object using my fourpcf code
_4pcf_object = fourpcf.measure_4pcf(density_field_data=data,
                                save_dir=save_dir, save_name = save_name, 
                                nbins=nbins, ld_one_d=resolution, ell_max=ell_max, bin_min=bin_min, bin_max=bin_max)


############################################################################

#run the algorithm
_4pcf_object.run_all(verbose_flag=True)

#save the zeta in a 6d array .npy file
_4pcf_object.save_zeta()


