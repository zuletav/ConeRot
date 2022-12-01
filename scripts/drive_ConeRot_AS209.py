import sys
import numpy as np
import re
from copy import copy, deepcopy
import os
from optparse import OptionParser

HOME = os.environ.get('HOME')
include_path = '/home/simon/common/python/include/'
#include_path=HOME+'/common/python/conemaps-git/'
sys.path.append(include_path)
import ConeRot.MasterDConeMaps as MasterDConeMaps

distance = 121.246

sourcedir = '/home/simon/AS209/guvmem_runs/12CO21/momentmaps/AS209_CO21_modout_lS0.00032_lL0.0_dgauss/'
workdir = 'work_modout_lS0.00032_lL0.0_dgauss_numba_c'

a_min = 0.7
a_max = 1.0
a_min_regions = 0.3
a_max_regions = 2.3

PA = 86.7  # continuum
inc = (180. - 35.3) * np.pi / 180.  # teague
tanpsi = 0.

#PA (85.74220998465579, 0.13529773191532968, 0.1327572377106918) 85.76418586465633
#inc (145.09924886452643, 0.07437000793240145, 0.08022007127306097) 145.10917890648307
#dra_off (0.0003049257218234983, 0.0002611812938863823, 0.00023822235134958746) 0.00028508059219256446
#ddec_off (0.0006638206724156879, 0.0002150304853814315, 0.0001963689464940382) 0.0006907247091780622

# #####################################################################
# #####################################################################

a_min_plot = a_min_regions
a_max_plot = a_max_regions

parser = OptionParser()
parser.add_option("-r",
                  "--retrograde",
                  action="store_true",
                  dest="RetroGrade",
                  default=False,
                  help="toggle retrograde orientation (RT trials only)")
parser.add_option("-f",
                  "--forceorient",
                  action="store_true",
                  dest="ForceOrient",
                  default=False,
                  help="toggle force input orientation in FixPAinc run")
parser.add_option("-F",
                  "--farside",
                  action="store_true",
                  dest="DoFarSideOnly",
                  default=False,
                  help="toggle far side only")
parser.add_option("-M",
                  "--MCMC",
                  action="store_true",
                  dest="RunMCMCmaster",
                  default=False,
                  help="toggle MCMC optim")
parser.add_option("-d",
                  "--dry-run",
                  action="store_false",
                  dest="RunMaster",
                  default=True,
                  help="toggle dry run")
parser.add_option("-o",
                  "--NoVarOrient",
                  action="store_false",
                  dest="DoVarOrient",
                  default=True,
                  help="no variable PA, inc profile, use with --forceorient")
parser.add_option("-R",
                  "--Regions",
                  action="store_true",
                  dest="Regions",
                  default=False,
                  help="use regions")
parser.add_option("-m",
                  "--Merid",
                  action="store_true",
                  dest="DoMerid",
                  default=False,
                  help="use meridional flows")

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

exec_master_script = sys.argv[0]

RunMCMCmaster = options.RunMCMCmaster
RunMaster = options.RunMaster
Regions = options.Regions

######################################################################

if re.match('^(.*)\/$', workdir):
    workdir = m.group(1)

if RunMaster:
    ClearWorkDir = True
else:
    ClearWorkDir = False

DoExec = False
PlotVarPAinc = False
if options.DoVarOrient:
    DoExec = RunMaster
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
    filename_errormap=sourcedir + 'im_g_v0_errormap.fits',
    workdir=workdir,
    DoErrorMap=True,
    typicalerror=0.1,  #  km/s
    ComputeSystVelo=True,  # best run this only once, then pass value in vsyst
    vsyst=4.67,
    #fieldscale=1., #1.
    fieldscale=1.5,  #1.
    pixscale_factor=3.0,  #3.
    unitscale=1.,
    PA=PA,
    inc=inc,
    tanpsi=-0.3,
    rangePA=20.,
    rangeinc=30. * np.pi / 180.,
    rangetanpsi=0.6,
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
    bmaj=134E-3,  # arcsec
    bmin=100E-3,  # arcsec
    DoConjGrad=True,
    DoMinuit=False,  # BROKEN 
    DoFarSideOnly=options.DoFarSideOnly,
    RunMCMC=RunMCMCmaster,
    RecoverMCMC=RunMCMCmaster,  # RunMCMC
    n_cores_MCMC=90,  #30
    Nit=400,
    burn_in=250,
    exec_master_script=exec_master_script)

S.domain = (('PA', (S.PA - S.rangePA / 2., S.PA + S.rangePA / 2.)),
            ('inc', (S.inc - S.rangeinc / 2., S.inc + S.rangeinc / 2.)),
            ('tanpsi', (S.tanpsi - S.rangetanpsi / 2.,
                        S.tanpsi + S.rangetanpsi / 2.)))

if S.DoExec:
    S.Run()

SFixOrient = copy(S)
if Regions:
    if S.DoFixOrient:
        SFixOrient.RunFixOrient(ForceGlobalOrient=options.ForceOrient,
                                Force_allradsPA=S.PA,
                                Force_allradsinc=S.inc)
    else:
        SFixOrient.workdir = re.sub('/$', '_fixPAinc/', SFixOrient.workdir)

import ConeRot.RotOrient.PlotRotorient

rgaps = [[0.18, 0.4], [1.6, 1.7]]
#pdl> p ( (190.-75/2.)/$d)
#1.48058252427184
#pdl> p ( (190.+75/2.)/$d)
#2.20873786407767
# Walsk+ 2014

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
        title='AS209',
        DoAUBar=True,
        alabel='',
        PlotVarOrient=True)

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
        title='AS209',
        DoAUBar=False,
        alabel='',
        PlotVarOrient=True)
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

import ConeRot.KineSummaryCompact

file_continuum = False  # './continuum/median_restored_finetav_fullim.fits'

ConeRot.KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,
                                                SFixOrient.filename_source,
                                                file_continuum=file_continuum,
                                                vsyst=S.vsyst,
                                                AllRads=False,
                                                a_min=a_min_plot,
                                                a_max=a_max_plot)

ConeRot.KineSummaryCompact.exec_summary_allrads(SFixOrient.workdir,
                                                SFixOrient.filename_source,
                                                file_continuum=file_continuum,
                                                vsyst=S.vsyst,
                                                AllRads=True,
                                                a_min=a_min_plot,
                                                a_max=a_max_plot)

file_continuum = False  # './continuum/median_restored_finetav_z_stretched.fits'

ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,
                                               SFixOrient.filename_source,
                                               file_continuum=file_continuum,
                                               vsyst=S.vsyst,
                                               AllRads=False,
                                               a_min=a_min_plot,
                                               a_max=a_max_plot,
                                               Zoom=True,
                                               side=1.5)

ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,
                                               SFixOrient.filename_source,
                                               file_continuum=file_continuum,
                                               vsyst=S.vsyst,
                                               AllRads=True,
                                               a_min=a_min_plot,
                                               a_max=a_max_plot,
                                               Zoom=False,
                                               side=1.5)

ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,
                                               SFixOrient.filename_source,
                                               file_continuum=file_continuum,
                                               vsyst=S.vsyst,
                                               AllRads=True,
                                               a_min=a_min_plot,
                                               a_max=a_max_plot,
                                               Zoom=True,
                                               side=3.,
                                               UseScatter=False)
