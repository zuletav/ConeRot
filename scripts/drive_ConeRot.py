import sys
import numpy as np
import re
from copy import copy,deepcopy
import os
from optparse import OptionParser

HOME=os.environ.get('HOME')
include_path=HOME+'/common/python/include/'
sys.path.append(include_path)
import ConeRot.MasterDConeMaps as MasterDConeMaps


sourcedir='hydro/12CO_dgauss_commonsigma/'

workdir='work_pix2_Smom8'
workdir='work_pix2_Smom8_noMerid'

sourcedir='/strelka_ssd/simon/C8_simus/AKIRA/hydro/13CO_sgauss_12co32_planetvortex_p05J/'
workdir='work_13co32_p05J_sgauss'


sourcedir='/strelka_ssd/simon/C8_simus/AKIRA/hydro/12CO_sgauss_12co32_planetvortex_p10J/'
workdir='work_12co32_p10J_sgauss'


sourcedir='/strelka_ssd/simon/C8_simus/AKIRA/hydro/13CO_sgauss_12co32_planetvortex_p10J/'
workdir='work_13co32_p10J_sgauss'

distance=140. #pc


#a_min=0.4
a_min=0.3
a_max=1.3

a_min_regions=0.15
a_max_regions=1.5

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
parser.add_option("-m", "--Merid", action="store_true", dest="DoMerid", default=False, help="use meridional flows")

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
print("options.DoMerid:", options.DoMerid)


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
    
    
if not options.DoMerid:
    workdir += '_nomerid'   
if Regions:
    workdir += '_Regions'
if options.ForceOrient:
    workdir += '_ForceOrient'
if options.DoFarSideOnly:
    workdir += '_FarSide'
if options.RunMCMCmaster:
    workdir += '_MCMC'
    
workdir += '/'

if options.RetroGrade:
    PA=320.-180.
    inc=(180.-30.)*np.pi/180.
else:
    PA=45.+180.
    inc=(33.)*np.pi/180.

    
S=MasterDConeMaps.Setup(
    filename_source=sourcedir+'im_g_v0.fits',
    #filename_source=sourcedir+'Smom_8.fits',
    filename_errormap=sourcedir+'im_g_v0_e.fits',
    workdir=workdir, 
    #DoErrorMap=False,
    DoErrorMap=False,
    typicalerror=0.02, #  km/s

    ComputeSystVelo=True,  # best run this only once, then pass value in vsyst
    vsyst=0.,  # 5.7454612183344285,
    #USED VSYST= 5.76188229931048
    #USED VSYST= 5.767168500867334
    #Rich Teague: 5.763

    fieldscale=1., #1.
    pixscale_factor=1.0,   #6.0
    unitscale=1.,

    PA=PA,
    inc=inc, 
    tanpsi=0.,

    rangePA=50.,
    rangeinc=15.*np.pi/180.,
    #rangeinc=0.1*np.pi/180.,
    rangetanpsi=1.,

    a_min=a_min, 
    a_max=a_max, 

    DoRegions=Regions,
    a_min_regions=a_min_regions,
    a_max_regions=a_max_regions,
    n_abins=10,# 

    DoAccr=False,
    DoAccr_fixPAinc=False,
    DoMerid_fixPAinc=options.DoMerid, 
    
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
    #nwalkers= 10, #leave as default
    exec_master_script=exec_master_script)



S.domain=( ('PA',(S.PA-S.rangePA/2.,S.PA+S.rangePA/2.)), ('inc',(S.inc-S.rangeinc/2.,S.inc+S.rangeinc/2.)),('tanpsi',(S.tanpsi-S.rangetanpsi/2.,S.tanpsi+S.rangetanpsi/2.)))


if S.DoExec:
    S.Run()

SFixOrient=copy(S)
if S.DoFixOrient:
    SFixOrient.Nit=100
    SFixOrient.burn_in=50
    SFixOrient.RunFixOrient(ForceGlobalOrient=options.ForceOrient, Force_allradsPA=S.PA, Force_allradsinc=S.inc)
else:
    SFixOrient.workdir=re.sub('/$','_fixPAinc/',SFixOrient.workdir)

import ConeRot.RotOrient.PlotRotorient


vsys=ConeRot.RotOrient.PlotRotorient.execfig(S.workdir,SFixOrient.filename_source, distance=distance, ForceGlobalOrient=options.ForceOrient, Force_allradsPA=S.PA, Force_allradsinc=S.inc, PlotVarPAinc=PlotVarPAinc, title='',alabel='')

vsys=ConeRot.RotOrient.PlotRotorient.execfig(S.workdir,SFixOrient.filename_source, distance=distance, ForceGlobalOrient=options.ForceOrient, Force_allradsPA=S.PA, Force_allradsinc=S.inc, PlotVarPAinc=PlotVarPAinc, title='',alabel='',PlotVarOrient=False)

SFixOrient.vsyst=vsys

import ConeRot.KineSummaryCompact


ConeRot.KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,SFixOrient.filename_source,vsyst=SFixOrient.vsyst,UseScatter=False)
ConeRot.KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,SFixOrient.filename_source,vsyst=SFixOrient.vsyst,UseScatter=True)

#file_continuum='/home/simon/rsynccommon/ppdisks/HD163296/DSHARP_continuum/polarmaps/polarmaps_modout_default/mod_out_z_stretched.fits'

ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,SFixOrient.filename_source,vsyst=SFixOrient.vsyst,UseScatter=False)

#ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,SFixOrient.filename_source,vsyst=SFixOrient.vsyst,Zoom=True)

#ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,SFixOrient.filename_source,file_continuum=file_continuum,vsyst=S.vsyst,Zoom=True,side=1.2)
#KineSummaryCompact.exec_summary_faceon(workdir,filename_source,file_continuum=file_continuum,vsyst=5.768)


