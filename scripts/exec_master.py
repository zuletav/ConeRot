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

sourcedir='/Users/simon/common/ppdisks/HD135344B/DMoments/12CO_sgaussmoments_constub_z/' # im_g_v0.fits 

#sourcedir='/Users/simon/common/ppdisks/HD135344B/DMoments/12CO_sgaussminuit_contsub_guvmem_lS0.003_lL0.0_smooth/' # im_g_v0.fits 

#sourcedir='/Users/simon/common/ppdisks/HD135344B/linefit/output_iminuit_multiso_dev/' # im_g_v0.fits 


S=MasterDConeMaps.Setup(
    filename_source=sourcedir+'im_g_v0.fits',
    filename_errormap=sourcedir+'im_g_v0_errormap.fits',

    # workdir='work_gmom_8_werrmap/',  # with trailing back slash
    # workdir='work_guvmem_werrmap/',  # with trailing back slash
    # workdir='work_guvmem_werrmap_MCMC/',  # with trailing back slash
    #workdir='work_guvmem_werrmap_MCMC_dev/',  # with trailing back slash
    #workdir='work_guvmem_werrmap_MCMC/',  # with trailing back slash
    #workdir='work_linefit_dev_werrmap_MCMC/',  # with trailing back slash
    #workdir='work_sgauss_contsub_MCMC/',  # with trailing back slash
    workdir='work_sgauss_contsub_dev/',  # with trailing back slash
    #workdir='work_test_vdisp/',  # with trailing back slash
    #workdir='work_test_sK/',  # with trailing back slash
    #workdir='work_test2/',  # with trailing back slash

    DoErrorMap=True,
    typicalerror=0.1, #  km/s

    ComputeSystVelo=False,  # best run this only once, then pass value in vsyst
    vsyst= 7.108,

    fieldscale=1.,
    pixscale_factor=1.0,   #6.0
    unitscale=1.,

    PA=243.,  #Stolker + PA=62
    inc=20.*np.pi/180. ,  #Stolker + inc=11
    tanpsi=0.3, 

    rangePA=40.,
    rangeinc=40.*np.pi/180.,
    rangetanpsi=0.6,

    a_min=0.15,
    a_max=0.35,

    DoRegions=True,
    a_min_regions=0.15,
    a_max_regions=0.35,
    n_abins=4, #minimum 3 for overlap

    DoAccr=False,
    DoAccr_fixPAinc=True,

    ClearWorkDir=True,
    DoExec=True,  # Execute the full optimization
    DoFixOrient=True, # Execute the restricted optimization, with fixed orientation

    DumpAllFitsFiles=False,

    x_center=0.,
    y_center=0.,

    #python ~/common/python/simon_examplescripts/beam.py    /Users/simon/common/ppdisks/HD163296/kine/data/HD163296_CO.fits
    #0.104 x 0.095 / -80.2 
    bmaj=0.054, # arcsec
    bmin=0.054, # arcsec

    DoConjGrad=True,
    DoMinuit=False,  # BROKEN 

    RunMCMC=True,
    RecoverMCMC=True, # RunMCMC
    n_cores_MCMC=30, #30
    Nit=120,
    burn_in=50,
    nwalkers= 20)

S.domain=( ('PA',(S.PA-S.rangePA/2.,S.PA+S.rangePA/2.)), ('inc',(S.inc-S.rangeinc/2.,S.inc+S.rangeinc/2.)),('tanpsi',(S.tanpsi-S.rangetanpsi/2.,S.tanpsi+S.rangetanpsi/2.)))


if S.DoExec:
    S.Run()

if S.DoFixOrient:
    S.RunFixOrient()

