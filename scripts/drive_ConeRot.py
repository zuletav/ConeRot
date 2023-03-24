import sys
import numpy as np
from copy import copy
import re
from optparse import OptionParser

from ConeRot import MasterDConeMaps

distance = 157.3 # pc GAIA 2018
incl = 180.-28
PAs  = 340.-180.

sourcedir='./model_sgauss/'
workdir = f'work_outer_inc{incl}_PA{PAs}'

title = r'HD142527 $^{12}$CO 6-5' # Title for the plot

a_min = 0.1
a_max = 1.
a_min_regions = 0.15
a_max_regions = 0.9

PA = PAs  # continuum
inc = (incl) * np.pi / 180.
tanpsi = -0.2

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
    filename_source=sourcedir + 'im_g_v0.fits', # Check the correct cube
    filename_errormap=sourcedir + 'im_g_v0_e.fits',
    workdir=workdir,
    DoErrorMap=True,
    typicalerror=0.1,  #  km/s
    ComputeSystVelo=False,  # best run this only once, then pass value in vsyst
    vsyst=3.72,
    #fieldscale=1., #1.
    fieldscale=1.0,  #1.
    pixscale_factor=1.0,  #3.
    unitscale=1.,
    PA=PA,
    inc=inc,
    tanpsi=-0.2,
    rangePA=30.,
    rangeinc=140. * np.pi / 180.,
    rangetanpsi=0.4,
    a_min=a_min,
    a_max=a_max,
    DoRegions=Regions,
    RestrictAvToRadialDomain=False,
    a_min_regions=a_min_regions,
    a_max_regions=a_max_regions,
    n_abins=2,
    DoAccr=False,
    DoAccr_fixPAinc=False,
    DoMerid_fixPAinc=options.DoMerid,
    ClearWorkDir=ClearWorkDir,
    DoExec=DoExec,  # Execute the full optimization
    DumpAllFitsFiles=False,
    #x_center=0.002,
    #y_center=0.012,
    x_center=0.0,
    y_center=0.0,
    bmaj=7.2E-5,  # arcsec
    bmin=5.7E-5,  # arcsec
    DoConjGrad=True,
    DoMinuit=False,  # BROKEN 
    DoFarSideOnly=options.DoFarSideOnly,
    RunMCMC=RunMCMCmaster,
    RecoverMCMC=RunMCMCmaster,  # RunMCMC
    n_cores_MCMC=4,
    Nit=10,
    burn_in=10,
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
        

from ConeRot.RotOrient import PlotRotorient

rgaps = False

if Regions:
    vsys = PlotRotorient.execfig(
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
        title=title,
        DoAUBar=True,
        alabel='',
        PlotVarOrient=True)

    vsys = PlotRotorient.execfig(
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
        title=title,
        DoAUBar=False,
        alabel='',
        PlotVarOrient=True)
    print("returned from execfig vsys", vsys)

else:
    vsys = PlotRotorient.execfig(
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

file_continuum = False

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

file_continuum = False

ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,
                                               SFixOrient.filename_source,
                                               file_continuum=file_continuum,
                                               vsyst=S.vsyst,
                                               AllRads=False,
                                               a_min=a_min_plot,
                                               a_max=a_max_plot,
                                               Zoom=True,
                                               side=1.0)

ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,
                                               SFixOrient.filename_source,
                                               file_continuum=file_continuum,
                                               vsyst=S.vsyst,
                                               AllRads=True,
                                               a_min=a_min_plot,
                                               a_max=a_max_plot,
                                               Zoom=False,
                                               side=1.0)

ConeRot.KineSummaryCompact.exec_summary_faceon(SFixOrient.workdir,
                                               SFixOrient.filename_source,
                                               file_continuum=file_continuum,
                                               vsyst=S.vsyst,
                                               AllRads=True,
                                               a_min=a_min_plot,
                                               a_max=a_max_plot,
                                               Zoom=True,
                                               side=1.0,
                                               UseScatter=False)
