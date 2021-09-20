import numpy as np
from subprocess import call
import astropy.io.fits as pyf
import sys
import os
import _4PCF_Algorithm as fourpcf

print("Importing necessary libraries and packages")

############################################################################

# an example of what you should type in the terminal is
#python3 FFT_4pcf_script.py mhd_index time_index  nbins resolution ell_max bin_min bin_max

############################################################################

mhd_folders = ['b1p1/', 'b1p.1/', 'b.1,p1/', 'b.1p.1/'] #mhd_index 0-3
time_folders = ['t_500/','t_550/','t_600/','t_650/','t_700/','t_750/','t_800/','t_850/','t_900/'] #time_index 0-8
time_files = ['t500','t550','t600','t650','t700','t750','t800','t850','t900'] #time_index 0-8
file_names = ['dens_' + time_files[i]+'.fits.gz' for i in range(len(time_folders))]

mhd_index = int(sys.argv[1])
time_index = int(sys.argv[2])

file_name =os.getcwd() + '../data/users.flatironinstitute.org/~bburkhart/data/CATS/MHD/256/' + mhd_folders[mhd_index] + time_folders[time_index] + file_names[time_index]

print(file_name)
hdulist = pyf.open(file_name)
data = hdulist[0].data.astype(np.float64)

#we only care about excess and deficit
data = data - np.mean(data)

############################################################################


nbins, resolution, ell_max, bin_min, bin_max = int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])

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


