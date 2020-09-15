import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#from pylab import *
import matplotlib.colors as colors
import re
from astropy import constants as const
import os    
import matplotlib.gridspec as gridspec


import PlotRotorient

workdir='prev_work/work_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_pix2_Merid_Filtered/'

#workdir='work_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_pix2_Merid_Filtered/'

workdir='work_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_pix4_Merid_Filtered/'

workdir='work_dgaussmoments_LOSNoise_model_cube_lS0.0005_lL1e-05_xwide_pix4_Merid_Filtered/'

#workdir='work_dgaussmoments_tclean_robust0.5_pix2_Merid/'

PlotRotorient.execfig(workdir,a_min=0.2,a_max=3.0)


