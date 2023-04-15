import sys
import numpy as np
import re
from copy import copy,deepcopy
import os
from optparse import OptionParser

import ConeRot.MasterDConeMaps as MasterDConeMaps

sourcedir='/home/simon/rsynccommon/ppdisks/HD163296/kine/data_1/dgaussmoments_LOSNoise_tclean_HD_163296briggs0.5_12CO_auto_wide/'
workdir='work_tclean_pix3_var_chi2'

sourcedir='/home/simon/rsynccommon/ppdisks/HD163296/kine_2021/DMoments/DSHARP_12CO_dgaussmoments_z/'
workdir='work_DSHARP'
workdir='work_DSHARP_pix2'

distance=101.5 #pc

inc_RichTeague=180.-46.8

#a_min=0.4
a_min=0.6
a_max=2.5

a_min_regions=0.3
a_max_regions=4.

a_min_plot=a_min_regions
a_max_plot=a_max_regions

######################################################################

parser = OptionParser()
parser.add_option("-r", "--retrograde", action="store_true", dest="RetroGrade", default=False, help="toggle retrograde orientation (RT trials only)")
parser.add_option("-f", "--forceorient", action="store_true", dest="ForceOrient", default=False, help="toggle force input orientation in FixPAinc run")
parser.add_option("-F", "--farside", action="store_true", dest="DoFarSideOnly", default=False, help="toggle far side only")
parser.add_option("-M", "--MCMC", action="store_true", dest="RunMCMCmaster", default=False, help="toggle MCMC optim")
parser.add_option("-d", "--dry-run", action="store_false", dest="RunMaster", default=True, help="toggle dry run")
parser.add_option("-o", "--NoVarOrient", action="store_false", dest="DoVarOrient", default=True, help="no variable PA, inc profile, use with --forceorient")
parser.add_option("-R", "--Regions", action="store_true", dest="Regions", default=False, help="use regions")

#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

print("options.RetroGrade:", options.RetroGrade)
print("options.ForceOrient:", options.ForceOrient)
print("options.DoFarSideOnly:", options.DoFarSideOnly)
print("options.RunMCMCmaster:", options.RunMCMCmaster)
print("options.RunMaster:", options.RunMaster)
print("options.DoVarOrient:", options.DoVarOrient)


######################################################################

exec_master_script=sys.argv[0]

RunMCMCmaster=options.RunMCMCmaster
RunMaster=options.RunMaster
Regions=options.Regions

######################################################################
#RunMCMCmaster=False

if RunMaster:
    ClearWorkDir=True
else:
    ClearWorkDir=False

DoExec=False
PlotVarPAinc=False
if options.DoVarOrient:
    DoExec=RunMaster
    PlotVarPAinc=True
    
    

if options.ForceOrient:
    workdir += '_ForceOrient'
if options.DoFarSideOnly:
    workdir += '_FarSide'
if options.RunMCMCmaster:
    workdir += '_MCMC'
    
workdir += '/'

if options.RetroGrade:
    PA=320.-180.
    inc=(180.-40.)*np.pi/180.
else:
    PA=312.8
    inc=inc_RichTeague*np.pi/180.


    
S=MasterDConeMaps.Setup(
    filename_source=sourcedir+'im_g_v0.fits',
    filename_errormap=sourcedir+'im_velo_errormap.fits',
    workdir=workdir, 
    DoErrorMap=True,
    #DoErrorMap=False,
    typicalerror=0.1, #  km/s

    ComputeSystVelo=True,  # best run this only once, then pass value in vsyst
    vsyst=5.761962,  # 5.7454612183344285,
    #USED VSYST= 5.76188229931048

    fieldscale=1., #1.
    pixscale_factor=2.0,   #6.0
    unitscale=1.,

    PA=PA,
    inc=inc, 
    tanpsi=0.,

    rangePA=40.,
    rangeinc=15.*np.pi/180.,
    #rangeinc=0.1*np.pi/180.,
    rangetanpsi=1.,

    a_min=a_min, 
    a_max=a_max, 

    DoRegions=Regions,
    a_min_regions=a_min_regions,
    a_max_regions=a_max_regions,
    n_abins=13,# 

    DoAccr=False,
    DoAccr_fixPAinc=False,
    DoMerid_fixPAinc=True,
    
    ClearWorkDir=ClearWorkDir,
    DoExec=DoExec,  # Execute the full optimization
    DoFixOrient=RunMaster, # Execute the restricted optimization, with fixed orientation

    DumpAllFitsFiles=False,

    x_center=0., # from the continuum
    y_center=0.,

    bmaj=0.091, # arcsec
    bmin=0.078, # arcsec
    DoConjGrad=True,
    DoMinuit=False,  # BROKEN 

    DoFarSideOnly=options.DoFarSideOnly,
    RunMCMC=RunMCMCmaster,
    RecoverMCMC=RunMCMCmaster, # RunMCMC
    n_cores_MCMC=30, #30
    #Nit=300,
    Nit=200,
    #burn_in=150,
    burn_in=100,
    #nwalkers= 20,
    nwalkers= 10,
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

#ConeRot.RotOrient.PlotRotorient.execfig(S.workdir,SFixOrient.filename_source, ForceGlobalOrient=False, Force_allradsPA=S.PA, Force_allradsinc=S.inc, WithComparData='rteague_compardata/velocity_profiles_a.npy', WithComparRadTWind=False, rgaps=[0.4433,  0.8575, 1.3923, 2.3])

ConeRot.RotOrient.PlotRotorient.execfig(S.workdir,SFixOrient.filename_source, distance=distance, ForceGlobalOrient=options.ForceOrient, Force_allradsPA=S.PA, Force_allradsinc=S.inc, WithComparData= 'rteague_compardata/velocity_profiles_a.npy', WithComparRadTWind=False, PlotVarPAinc=PlotVarPAinc, rgaps=[0.4433,  0.8575, 1.3923, 2.3])


import ConeRot.KineSummaryCompact

#KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,SFixOrient.filename_source,file_continuum='',vsyst=S.vsyst)






delta_planet = 2.3
PA_planet=-3. * np.pi / 180.
x_planet= delta_planet * np.sin(PA_planet)
y_planet= delta_planet * np.cos(PA_planet)


PA=SFixOrient.PA
inc=SFixOrient.inc

PArad=PA*np.pi/180.
x_planetp=x_planet*np.cos(PArad)-y_planet*np.sin(PArad)
y_planetp=x_planet*np.sin(PArad)+y_planet*np.cos(PArad)

x_planetpp = x_planetp / np.fabs(np.cos(inc))
y_planetpp = y_planetp

RegionOfInterest=[x_planet,y_planet]

RegionOfInterest_faceon=[x_planetpp,y_planetpp]



file_continuum='/home/simon/rsynccommon/ppdisks/HD163296/DSHARP_continuum/guvmem_runs/mem_lS0.0_lL0.0_nogrid/mod_out.fits'

ConeRot.KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,SFixOrient.filename_source,file_continuum=file_continuum,vsyst=S.vsyst,RegionOfInterest=RegionOfInterest)

file_continuum='/home/simon/rsynccommon/ppdisks/HD163296/DSHARP_continuum/polarmaps/polarmaps_modout_default/mod_out_z_stretched.fits'

ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,SFixOrient.filename_source,file_continuum=file_continuum,vsyst=S.vsyst,RegionOfInterest=RegionOfInterest_faceon)

ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,SFixOrient.filename_source,file_continuum=file_continuum,vsyst=S.vsyst,Zoom=True)

#ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,SFixOrient.filename_source,file_continuum=file_continuum,vsyst=S.vsyst,Zoom=True,side=1.2)
#KineSummaryCompact.exec_summary_faceon(workdir,filename_source,file_continuum=file_continuum,vsyst=5.768)


