import sys
import numpy as np
import re
from copy import copy,deepcopy
import os


from optparse import OptionParser

HOME=os.environ.get('HOME')
include_path='/Users/simon/common/python/include/'
#include_path=HOME+'/common/python/conemaps-git/'
sys.path.append(include_path)

import ConeRot.MasterDConeMaps as MasterDConeMaps


parser = OptionParser()
parser.add_option("-r", "--retrograde", action="store_true", dest="RetroGrade", default=False, help="toggle retrograde orientation")
parser.add_option("-f", "--forceorient", action="store_true", dest="ForceOrient", default=False, help="toggle force input orientation in FixPAinc run")

#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

print("options.RetroGrade:", options.RetroGrade)
print("options.ForceOrient:", options.ForceOrient)


######################################################################

exec_master_script=sys.argv[0]


#sourcedir='/Users/simon/common/ppdisks/HD100546/kine_2020/Dmoments_2020/dgaussmoments_LOSNoise_calibrated_12co_concatall_avg_selfcal_contsub_briggs1.0.image/'
sourcedir='./dgaussmoments_LOSNoise_image_12CO_2-1_contsub_wnoise_smooth/'


RunMCMCmaster=False
RunMaster=True

######################################################################
#RunMCMCmaster=False

if RunMaster:
    ClearWorkDir=True
else:
    ClearWorkDir=False


workdir='work_Merid_varchi2'

if options.ForceOrient:
    workdir += '_ForceOrient'
workdir += '/'

if options.RetroGrade:
    PA=320.-180.
    inc=(180.-40.)*np.pi/180.
else:
    PA=320.
    inc=40.*np.pi/180.
    

S=MasterDConeMaps.Setup(
    filename_source=sourcedir+'im_g_v0.fits',
    filename_errormap=sourcedir+'im_velo_errormap.fits',
    #filename_errormap=sourcedir+'im_g_v0_e.fits',
    #workdir='dgaussmoments_LOSNoise_12co_pix2_wide_Merid_Filtered/'
    #workdir='dgaussmoments_LOSNoise_12co_pix2_wide_Merid_Filtered/',
    #workdir='dgaussmoments_LOSNoise_calibrated_12co_pix2_Merid/',
    #workdir='dgaussmoments_LOSNoise_12co_pix3_wide_Merid_Filtered_MCMC/',
    #workdir='dgaussmoments_LOSNoise_12co_pix2_wide_Merid_Filtered_MCMC/',
    #workdir='work_dgaussmoments_LOSNoise_tclean_briggs0.5_pix3_FloorNoise_Merid_Filtered_MCMCb/',
    #workdir='work_dgaussmoments_LOSNoise_tclean_briggs0.5_pix3_FloorNoise_Merid_Filtered_MCMC_noNcorr/',
    #workdir='work_dgaussmoments_LOSNoise_tclean_briggs0.5_pix3_FloorNoise_Merid_Filtered_MCMC_invNcorr/',
    workdir=workdir, 
    #workdir='dgaussmoments_LOSNoise_12co_pix2_wide_Merid_Filtered_MCMC/',
    DoErrorMap=True,
    typicalerror=0.001, #  km/s

    ComputeSystVelo=False,  # best run this only once, then pass value in vsyst
    #vsyst=5.6763827,
    vsyst=0., 
    
    fieldscale=1., #1.
    #pixscale_factor=6.0,   #6.0
    pixscale_factor=1.0,   #6.0
    #pixscale_factor=3.0,   #6.0
    unitscale=1.,

    PA=PA,
    inc=inc, 
    tanpsi=0.,

    # using global PA: 312.732843
    # using global inc: 2.330638
    # using global tanpsi: -0.255051

    rangePA=40.,
    rangeinc=30.*np.pi/180.,
    #rangeinc=0.1*np.pi/180.,
    rangetanpsi=0.8,

    #a_min=0.15,
    #a_max=0.75,
    a_min=0.25,
    a_max=0.5,

    DoRegions=True,
    #a_min_regions=0.35,
    a_min_regions=0.15,
    a_max_regions=0.9,
    #n_abins=17, #minimum 3 for overlap
    #n_abins=17, #minimum 3 for overlap
    n_abins=11, #minimum 3 for overlap

    DoAccr=False,
    DoAccr_fixPAinc=False,
    DoMerid_fixPAinc=True,
    ClearWorkDir=ClearWorkDir,
    DoExec=RunMaster,  # Execute the full optimization
    DoFixOrient=RunMaster, # Execute the restricted optimization, with fixed orientation

    DumpAllFitsFiles=False,

    x_center=0.,
    y_center=0.,

    #strelka10:32:55~/common/ppdisks/HD163296/kine/data_1$beam.py dgaussmoments_restored_cube_lS0.0005_lL1e-05_xwide_robust0.5/im_g_v0.fits 
    #python ~/common/python/simon_examplescripts/beam.py    /Users/simon/common/ppdisks/HD163296/kine/data/HD163296_CO.fits
    #0.104 x 0.095 / -80.2 
    bmaj=59E-3, # arcsec
    bmin=42E-3, # arcsec

    DoConjGrad=True,
    DoMinuit=False,  # BROKEN 

    DoFarSideOnly=False,
    RunMCMC=RunMCMCmaster,
    RecoverMCMC=RunMCMCmaster, # RunMCMC
    n_cores_MCMC=30, #30
    Nit=200,
    burn_in=100,
    nwalkers= 20,
    exec_master_script=exec_master_script)

S.domain=( ('PA',(S.PA-S.rangePA/2.,S.PA+S.rangePA/2.)), ('inc',(S.inc-S.rangeinc/2.,S.inc+S.rangeinc/2.)),('tanpsi',(S.tanpsi-S.rangetanpsi/2.,S.tanpsi+S.rangetanpsi/2.)))


if S.DoExec:
    S.Run()

SFixOrient=copy(S)
if S.DoFixOrient:

    SFixOrient.RunFixOrient(ForceGlobalOrient=options.ForceOrient, Force_allradsPA=S.PA, Force_allradsinc=S.inc)
else:
    SFixOrient.workdir=re.sub('/$','_fixPAinc/',SFixOrient.workdir)

import ConeRot.RotOrient.PlotRotorient

#ConeRot.RotOrient.PlotRotorient.execfig(S.workdir)
ConeRot.RotOrient.PlotRotorient.execfig(S.workdir,SFixOrient.filename_source, ForceGlobalOrient=options.ForceOrient, Force_allradsPA=S.PA, Force_allradsinc=S.inc,WithComparData=False,WithComparRadTWind=True)


import ConeRot.KineSummaryCompact

#KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,SFixOrient.filename_source,file_continuum='',vsyst=S.vsyst)



file_continuum='image_out_wl880_smooth.fits'


ConeRot.KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,SFixOrient.filename_source,file_continuum=file_continuum,vsyst=S.vsyst)



ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,SFixOrient.filename_source,file_continuum=file_continuum,vsyst=S.vsyst)

#ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,SFixOrient.filename_source,file_continuum=file_continuum,vsyst=S.vsyst,Zoom=True,side=1.2)
#KineSummaryCompact.exec_summary_faceon(workdir,filename_source,file_continuum=file_continuum,vsyst=5.768)


