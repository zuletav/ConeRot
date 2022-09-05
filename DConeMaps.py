import sys
import numpy as np
import os
import os.path
from scipy import ndimage
from astropy.io import fits as pf
import re
from copy import deepcopy
from astropy.wcs import WCS
from scipy import optimize
import time
from time import gmtime,strftime

#import matplotlib as plt
#import matplotlib.pyplot as plt
#from pylab import *
#import matplotlib.colors as colors


include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
import ConeRot.funcs_DConeMaps 


class Model():
    def __init__(self,
                 filename_source='',
                 filename_errormap='',
                 workdir='',
                 PA=95.,  # pass degrees
                 inc=40.*np.pi/180., # pass radians
                 tanpsi=0.05,
                 RA=False,
                 DEC=False,
                 x_center=0.,
                 y_center=0.,
                 # ErrorMap=False,
                 DoErrorMap=False,
                 XCheckInv=False,
                 ComputeSkyImages=True,
                 DoAzimuthalProfile=False,
                 PlotRadialProfile=True,
                 vsyst=0.,
                 sigma_vsyst=0.,
                 ComputeSystVelo=False,
                 DoAccr=False,
                 DoMerid=False,
                 a_min=0.2,
                 a_max=0.6,
                 LegBkgd=-1,
                 fieldscale=1.0,
                 pixscale_factor=1.,
                 unitscale=1.,
                 typicalerror=1E-2,
                 InjectNoise=True,
                 PrintOptimStatus=True,
                 DumpAllFitsFiles=False,
                 DoDCone=False,
                 domain=(),
                 RunMCMC=False,
                 BlindMCMC=False,
                 Hdu=False,
                 Hduw=False,
                 Hducentered=False,
                 Hduwcentered=False,
                 Hduregion=False,
                 Hdudiff=False,
                 Hdumoddrot=False,
                 Hduregion_faceon=False,
                 Verbose=False,
                 VerboseInit=False,
                 fout=False,
                 filelog='log_output.txt',
                 StoreRegions=False,
                 DoConjGrad=False,
                 DoMinuit=False,
                 ClearWorkDir=False,
                 InheritGlobalInit=False, # to force same initial conditions for all regions
                 PA0=0.,
                 inc0=10.,
                 tanpsi0=0.1,
                 a_min_regions=0.2,
                 a_max_regions=0.3,
                 iregion = 0, # region index
                 n_abins=12,
                 InheritMumap=False,  # pass mumap from a previous orientation - used as weights in KepAmps
                 mumap=None,
                 Hdumumap=None,
                 HduDConemoddrot=None,
                 HdudiffDConemoddrot=None,
                 Hdudiff_faceon=None,
                 Hduresamp_faceon=None, 
                 Hdurrs=None,
                 Hduphis=None,
                 Hdurrs_faceon=None,
                 Hduphis_faceon=None,
                 chi2_prev=0.,
                 polarchi2=0.,
                 skychi2=0.,
                 velodev_med=0.,
                 velodev_std=0.,
                 velodev_std2=0.,
                 DConeDeltaChi2=False,
                 RadialProfile=None,
                 bmaj=0.07, # arcsec
                 bmin=0.07, # arcsec
                 Ncorr=1.,  # dimensionless
                 DoFarSideOnly=False,
                 ExtendRegions=False, # True: extends inner and outer regions to whole radial domain - for plots only
                 RestrictAvToRadialDomain=False,
                 TriangleFile='triangle.png',
                 Diskgeometry={}):

        


        initlocals=locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            if VerboseInit:
                print( "DConeMaps setting ",a_attribute," to ",initlocals[a_attribute])
            setattr(self,a_attribute,initlocals[a_attribute])
        


    def conicpolar_expansions(self):
        return ConeRot.funcs_DConeMaps.exec_conicpolar_expansions(self)
            
    def prep_files(self):
        ConeRot.funcs_DConeMaps.exec_prep_files(self)

    def grid_4center(self):
        ConeRot.funcs_DConeMaps.exec_grid_4center(self)
