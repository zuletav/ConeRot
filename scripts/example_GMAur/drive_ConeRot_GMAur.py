import sys
import numpy as np
import re
from copy import copy, deepcopy
import os
from pprint import pprint

from ConeRot import MasterDConeMaps

distance = 121.246

sourcedir = '/home/simon/MAPS_kine/GMAur/data/GMAur_13CO21_modout_lS0.00032_lL0.0_sgauss-mask_zoom-restored_new/'
workdir = 'conerot_output'

a_min = 1.0
a_max = 1.5
a_min_regions = 0.4
a_max_regions = 1.9

a_min_plot = a_min_regions
a_max_plot = a_max_regions

PA = 57.2  
inc = 53.2 * np.pi / 180.  #

######################################################################


parser = MasterDConeMaps.gen_scipt_options_argparse()
#(options, args) = parser.parse_args()
options  = parser.parse_args()
pprint(options)

######################################################################

exec_master_script = sys.argv[0]

RunMCMCmaster = options.RunMCMCmaster
RunMaster = options.RunMaster
Regions = options.Regions

######################################################################

if RunMaster:
    DoExec = True
    ClearWorkDir = True
else:
    DoExec = False
    ClearWorkDir = False

PlotVarPAinc = False
if options.DoVarOrient:
    PlotVarPAinc = True

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

print("workdir>>>> ", workdir)

S = MasterDConeMaps.Setup(
    filename_source=sourcedir + 'im_g_v0.fits',
    filename_errormap=sourcedir + 'im_g_v0_e.fits',
    workdir=workdir,
    DoErrorMap=True,
    typicalerror=0.1,  #  km/s
    ComputeSystVelo=True,  # best run this only once, then pass value in vsyst
    vsyst=5.6,
    #fieldscale=1., #1.
    fieldscale=1.0,  #1.
    pixscale_factor=2.0,  #3.
    unitscale=1.,
    PA=PA,
    inc=inc,
    tanpsi=0.1,
    rangePA=20.,
    rangeinc=30. * np.pi / 180.,
    rangetanpsi=0.2,
    a_min=a_min,  #0.17
    a_max=a_max,  #0.27
    DoRegions=Regions,
    RestrictAvToRadialDomain=False,
    a_min_regions=a_min_regions,
    a_max_regions=a_max_regions,
    n_abins=7,  # 6   #minimum 3 for overlap
    DoAccr=False,
    DoAccr_fixPAinc=False,
    DoMerid_fixPAinc=options.DoMerid,
    ClearWorkDir=ClearWorkDir,
    DoExec=DoExec,  # Execute the full optimization
    DumpAllFitsFiles=False,
    #x_center=0.002,  # from the continuum
    #y_center=0.012,
    x_center=0.0,
    y_center=0.0,
    bmaj=204E-3,  # arcsec
    bmin=155E-3,  # arcsec
    DoConjGrad=True,
    DoMinuit=False,  # BROKEN 
    DoFarSideOnly=options.DoFarSideOnly,
    RunMCMC=RunMCMCmaster,
    RecoverMCMC=RunMCMCmaster,  # RunMCMC
    n_cores_MCMC=30,  #30
    Nit=400,
    burn_in=250,
    VerboseInit=False,
    exec_master_script=exec_master_script)

S.domain = (('PA', (S.PA - S.rangePA / 2., S.PA + S.rangePA / 2.)),
            ('inc', (S.inc - S.rangeinc / 2., S.inc + S.rangeinc / 2.)),
            ('tanpsi', (S.tanpsi - S.rangetanpsi / 2.,
                        S.tanpsi + S.rangetanpsi / 2.)))

if S.DoExec:
    S.Run()

SFixOrient = copy(S)
if Regions:
    if S.DoFixOrient and S.DoExec:
        SFixOrient.RunFixOrient(ForceGlobalOrient=options.ForceOrient,
                                Force_allradsPA=S.PA,
                                Force_allradsinc=S.inc)
    else:
        SFixOrient.workdir = re.sub('/$', '_fixPAinc/', SFixOrient.workdir)


        
"""
plotting results
"""

import ConeRot.RotOrient.PlotRotorient

rgaps = [[0.49, 0.52], [0.79, 0.81], [1.69, 1.71]]

if Regions:
    vsys = ConeRot.RotOrient.PlotRotorient.execfig(
        S.workdir,
        SFixOrient.filename_source,
        distance=distance,
        ForceGlobalOrient=options.ForceOrient,
        Force_allradsPA=S.PA,
        Force_allradsinc=S.inc,
        WithComparData=False,
        WithComparRadTWind=False,
        PlotVarPAinc=PlotVarPAinc,
        rgaps=rgaps,
        title='GMAur',
        DoAUBar=False,
        #alabel='',
        PlotVarOrient=True)

    #vsys = ConeRot.RotOrient.PlotRotorient.execfig(
    #    S.workdir,
    #    SFixOrient.filename_source,
    #    distance=distance,
    #    ForceGlobalOrient=options.ForceOrient,
    #    Force_allradsPA=S.PA,
    #    Force_allradsinc=S.inc,
    #    WithComparData=False,
    #    WithComparRadTWind=False,
    #    PlotVarPAinc=PlotVarPAinc,
    #    rgaps=rgaps,
    #    title='AS209',
    #    DoAUBar=False,
    #    alabel='',
    #    PlotVarOrient=True)
    print("returned from execfig vsys", vsys)

    #vsys=ConeRot.RotOrient.PlotRotorient.execfig(S.workdir,SFixOrient.filename_source, distance=distance, ForceGlobalOrient=options.ForceOrient, Force_allradsPA=S.PA, Force_allradsinc=S.inc, WithComparData= False, WithComparRadTWind=False, PlotVarPAinc=PlotVarPAinc, rgaps=rgaps,title='HD100546',DoAUBar=False,alabel='',PlotVarOrient=False)
    print("returned from execfig vsys", vsys)

else:
    vsys = ConeRot.RotOrient.PlotRotorient.execfig(
        S.workdir,
        SFixOrient.filename_source,
        distance=distance,
        ForceGlobalOrient=options.ForceOrient,
        Force_allradsPA=S.PA,
        Force_allradsinc=S.inc,
        WithComparData=False,
        WithComparRadTWind=False,
        PlotVarPAinc=PlotVarPAinc,
        VarOrient=False,
        a_min=a_min_plot,
        a_max=a_max_plot,
        Plot_vRot_Global=True,
        Plot_vRot_VarOrient=False,
        Plot_vRot_VarOrient_FixIncPA=False,
        rgaps=rgaps)

SFixOrient.vsyst = vsys
S.vsyst = vsys

from ConeRot import KineSummaryCompact

file_continuum = False  # './continuum/median_restored_finetav_fullim.fits'

KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,
                                        SFixOrient.filename_source,
                                        file_continuum=file_continuum,
                                        vsyst=S.vsyst,
                                        AllRads=False,
                                        a_min=a_min_plot,
                                        a_max=a_max_plot)

KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,
                                        SFixOrient.filename_source,
                                        file_continuum=file_continuum,
                                        vsyst=S.vsyst,
                                        AllRads=True,
                                        a_min=a_min_plot,
                                        a_max=a_max_plot)

file_continuum = False  # './continuum/median_restored_finetav_z_stretched.fits'

KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,
                                       SFixOrient.filename_source,
                                       file_continuum=file_continuum,
                                       vsyst=S.vsyst,
                                       AllRads=False,
                                       a_min=a_min_plot,
                                       a_max=a_max_plot,
                                       Zoom=True,
                                       side=1.5)

KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,
                                       SFixOrient.filename_source,
                                       file_continuum=file_continuum,
                                       vsyst=S.vsyst,
                                       AllRads=True,
                                       a_min=a_min_plot,
                                       a_max=a_max_plot,
                                       Zoom=False,
                                       side=1.5)

KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,
                                       SFixOrient.filename_source,
                                       file_continuum=file_continuum,
                                       vsyst=S.vsyst,
                                       AllRads=True,
                                       a_min=a_min_plot,
                                       a_max=a_max_plot,
                                       Zoom=True,
                                       side=3.,
                                       UseScatter=False)
