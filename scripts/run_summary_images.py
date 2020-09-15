import os
import re
from astropy.io import fits 
import scipy
import scipy.signal
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from  scipy.signal import medfilt2d

#include_path='/Users/simon/common/python/include/'
#sys.path.append(include_path)

import KineSummaryCompact

sourcedir='/Users/simon/common/ppdisks/HD163296/kine/data_1/dgaussmoments_tclean_HD_163296briggs0.5_12CO_auto_wide/'
filename_source=sourcedir+'im_g_v0.fits'
#workdir='work_dgaussmoments_tclean_robust0.5_pix4_Merid_Filtered_fixPAinc/'
workdir='work_dgaussmoments_LOSNoise_tclean_briggs0.5_pix2_Merid_fixPAinc/'

#workdir='workd_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_pix2_Merid/'
#workdir='work_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_pix2_Merid_Filtered/'
#workdir='workd_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_Merid_Filtered/'


delta_planet = 2.3
PA_planet=-3. * np.pi / 180.
x_planet= delta_planet * np.sin(PA_planet)
y_planet= delta_planet * np.cos(PA_planet)
RegionOfInterest=[x_planet,y_planet]

        
file_continuum='/Users/simon/common/ppdisks/HD163296/DSHARP_continuum/guvmem_runs/mem_lS0.0_lL0.0_nogrid/mod_out.fits'
KineSummaryCompact.exec_summary_allrads(workdir,filename_source,file_continuum=file_continuum,vsyst=5.768,RegionOfInterest=RegionOfInterest)



delta_planet = 2.3
PA_planet=-3. * np.pi / 180.
x_planet= delta_planet * np.sin(PA_planet)
y_planet= delta_planet * np.cos(PA_planet)


PA=312.379492
inc=2.330638

PArad=PA*np.pi/180.
x_planetp=x_planet*np.cos(PArad)-y_planet*np.sin(PArad)
y_planetp=x_planet*np.sin(PArad)+y_planet*np.cos(PArad)

x_planetpp = x_planetp / np.fabs(np.cos(inc))
y_planetpp = y_planetp

RegionOfInterest=[x_planetpp,y_planetpp]

print("RegionOfInterest",RegionOfInterest)
file_continuum='/Users/simon/common/ppdisks/HD163296/DSHARP_continuum/polarmaps/polarmaps_modout_default/mod_out_z_stretched.fits'
KineSummaryCompact.exec_summary_faceon(workdir,filename_source,file_continuum=file_continuum,vsyst=5.768,RegionOfInterest=RegionOfInterest)

