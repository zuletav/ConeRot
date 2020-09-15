import sys
import numpy as np
import re

import os
HOME=os.environ.get('HOME')
include_path='/Users/simon/common/python/include/'
#include_path=HOME+'/common/python/conemaps-git/'
sys.path.append(include_path)

import ConeRot.MasterDConeMaps as MasterDConeMaps

######################################################################

#sourcedir='/Users/simon/common/ppdisks/HD163296/kine/data_1/dgaussmoments_LOSNoise_restored_cube_lS0.0005_lL1e-05_xwide_robust0.5/'

sourcedir='/Users/simon/common/ppdisks/HD163296/kine/data_1/dgaussmoments_restored_cube_lS0.0005_lL1e-05_xwide_robust0.5/'
sourcedir='/Users/simon/common/ppdisks/HD163296/kine/data_1/dgaussmoments_tclean_HD_163296briggs0.5_12CO_auto_wide/'
sourcedir='/Users/simon/common/ppdisks/HD163296/kine/data_1/dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide/'


RunMCMCmaster=False

S=MasterDConeMaps.Setup(
    filename_source=sourcedir+'im_g_v0.fits',
    filename_errormap=sourcedir+'im_g_v0_errormap.fits',
    #filename_errormap=sourcedir+'im_g_v0_e.fits',
    #workdir='work_dgaussmoments_restored_cube_lS0.0005_lL1e-05_xwide_robust0.5_pix6_Accr/',  # with trailing back slash
    #workdir='work_dgaussmoments_restored_cube_lS0.0005_lL1e-05_xwide_robust0.5_Merid/',  # with trailing back slash
    #workdir='work_dgaussmoments_restored_cube_lS0.0005_lL1e-05_xwide_robust0.5_pix2_Merid_Filtered/',  # with trailing back slash
    #workdir='work_dgaussmoments_LOSNoise_restored_cube_lS0.0005_lL1e-05_xwide_robust0.5_pix6_Merid/',  # with trailing back slash
    #workdir='work_dgaussmoments_tclean_robust0.5_pix6_Merid/',  # with trailing back slash
    #workdir='work_dgaussmoments_tclean_robust0.5_pix2_Merid/',  # with trailing back slash
    #workdir='work_dgaussmoments_tclean_robust0.5_pix4_Merid_Filtered/',  # with trailing back slash
    #workdir='work_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_pix2_Merid/',  # with trailing back slash
    workdir='work_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_pix2_Merid_Filtered/',  # with trailing back slash
    #workdir='work_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_Merid_Filtered/',  # with trailing back slash
    #workdir='work_dgaussmoments_model_cube_lS0.0005_lL1e-05_xwide_Merid_Filtered/',  # with trailing back slash
    #workdir='work_dgaussmoments_restored_cube_lS0.0005_lL1e-05_xwide_robust0.5_pix3/',  # with trailing back slash
    DoErrorMap=True,
    typicalerror=0.1, #  km/s

    ComputeSystVelo=False,  # best run this only once, then pass value in vsyst
    #vsyst= 5.76636 ,  # 5.76860412669,
    #vsyst=5.748031675193103,
    #vsyst= 5.74850163679975,
    #vsyst=5.749784298010025,
    vsyst=5.7454612183344285,
    
    fieldscale=1.,
    #pixscale_factor=6.0,   #6.0
    pixscale_factor=2.0,   #6.0
    #pixscale_factor=3.0,   #6.0
    unitscale=1.,

    PA=312.379492,
    inc=2.330638,
    tanpsi=-0.255051,

    # using global PA: 312.732843
    # using global inc: 2.330638
    # using global tanpsi: -0.255051

    rangePA=40.,
    rangeinc=15.*np.pi/180.,
    rangetanpsi=0.8,

    a_min=0.4,
    a_max=3.0,

    DoRegions=True,
    #a_min_regions=0.35,
    #a_min_regions=0.4,
    #a_max_regions=4.5,
    a_min_regions=0.05,
    a_max_regions=3.5,
    #n_abins=14, #  #minimum 3 for overlap
    #n_abins=15, #  #minimum 3 for overlap
    n_abins=15, #  #minimum 3 for overlap

    DoAccr=False,
    DoAccr_fixPAinc=False,
    DoMerid_fixPAinc=True,
    ClearWorkDir=True,
    DoExec=True,  # Execute the full optimization
    DoFixOrient=True, # Execute the restricted optimization, with fixed orientation

    DumpAllFitsFiles=False,

    x_center=0.,
    y_center=0.,

    #strelka10:32:55~/common/ppdisks/HD163296/kine/data_1$beam.py dgaussmoments_restored_cube_lS0.0005_lL1e-05_xwide_robust0.5/im_g_v0.fits 
    #0.083 x 0.083 / 0.1 
    #bmaj=0.083, # arcsec
    #bmin=0.083, # arcsec
    bmaj=0.091, # arcsec
    bmin=0.091, # arcsec

    DoConjGrad=True,
    DoMinuit=False,  # BROKEN 

    RunMCMC=RunMCMCmaster,
    RecoverMCMC=RunMCMCmaster, # RunMCMC
    n_cores_MCMC=30, #30
    Nit=120,
    burn_in=50,
    nwalkers= 20)

S.domain=( ('PA',(S.PA-S.rangePA/2.,S.PA+S.rangePA/2.)), ('inc',(S.inc-S.rangeinc/2.,S.inc+S.rangeinc/2.)),('tanpsi',(S.tanpsi-S.rangetanpsi/2.,S.tanpsi+S.rangetanpsi/2.)))


if S.DoExec:
    S.Run()

if S.DoFixOrient:
    S.RunFixOrient()

import PlotRotorient

PlotRotorient.execfig(S.workdir)

