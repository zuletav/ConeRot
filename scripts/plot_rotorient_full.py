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

HOME=os.environ.get('HOME')
print("HOME: ",HOME)
include_path=HOME+'/common/python/include/'
sys.path.append(include_path)

import ConeRot.RotOrient.LogResults as  LogResults
import ConeRot.RotOrient.StellarMass as StellarMass 
import ConeRot.RotOrient.RotCurve  as RotCurve
import ConeRot.RotOrient.Orient as Orient

#workdir='work_gmom_8_werrmap/'
#workdir='work_guvmem_werrmap_MCMC/'
#workdir='work_linefit_dev_werrmap_MCMC/'#
#workdir='work_sgauss_contsub_MCMC/'
#workdir='work_sgauss_contsub_dev/'
workdir='work_test_sK/'
#workdir='work_test/'

bmaj=0.054 # arcsec
distance=135.77 # GAIA  2018
nametag='im_g_v0'

Plot_vRot_Global=False  # if true will plot rotation curve for global (fixed) orientation  (so not averaged over Regions)
Plot_vRot_VarOrient=False  # if true will plot rotation curve for variable orientation too (so averaged over Regions)
VarOrient=True # will plot the variable orientation obtained from the Regions
#
Plot_vRot_VarOrient_FixIncPA = True
XCheckFixIncPA=False # cross check that the rot curve is the same for global optim  and global init optim for  fix PA and inc (should be the same PA, inc, psi)

rgaps=False #[0.4433,  0.8575, 1.3923, 2.3]

#######################################################################
# no further user input
######################################################################
 
#matplotlib.rc('text', usetex=True) 
#matplotlib.rc('font', family='arial') 
#matplotlib.rc('font', family='helvetica') 
matplotlib.rc('font', family='arial') 
matplotlib.rcParams.update({'font.size': 12})



workdir_fixincPA=re.sub('/$','_fixPAinc/',workdir)
DoFixIncPA=False
if os.path.isdir(workdir_fixincPA):
    DoFixIncPA=True

file_profile=workdir+nametag+'_radial_profile.dat'
file_profile_allrads=workdir+nametag+'_allrads_radial_profile.dat'
file_profile_fixincPA=workdir_fixincPA+nametag+'_radial_profile.dat'
file_profile_fixincPA_allrads=workdir_fixincPA+nametag+'_allrads_radial_profile.dat'


fileout_fig=workdir_fixincPA+'fig_rotorient_surfandmid_linear_dev_full.pdf'
print("loading file_profiles",file_profile," ",file_profile_fixincPA)

allprofiles=np.loadtxt(file_profile,unpack=True)
print("length allprofiles",len(allprofiles))
if (len(allprofiles) > 3):
    (rrs, v_Phi_prof, sv_Phi_prof, v_R_prof, sv_R_prof) = allprofiles # np.loadtxt(file_profile,unpack=True)
    DoAccr=True
    print("WARNING: found DoAccr optimization for global optim (no regions)  but will plot v_R averaged over regions (allrads)")
else:
    (rrs, v_Phi_prof, sv_Phi_prof) = np.loadtxt(file_profile,unpack=True)
    DoAccr=False
    
if VarOrient:
    allprofiles_allrads=np.loadtxt(file_profile_allrads,unpack=True)
    print("length allprofiles",len(allprofiles_allrads))
    if (len(allprofiles_allrads) > 3):
        (rrs_allrads, v_Phi_prof_allrads, sv_Phi_prof_allrads, v_R_prof_allrads, sv_R_prof_allrads) = allprofiles_allrads
        DoAccrAllRads=True
        print("WARNING: found DoAccr optimization with variable PA and inc ")
    else:
        (rrs_allrads, v_Phi_prof_allrads, sv_Phi_prof_allrads) = np.loadtxt(file_profile_allrads,unpack=True)
        DoAccrAllRads=False


DoAccr_fixIncPA=False
if DoFixIncPA:
    allprofiles_fixincPA=np.loadtxt(file_profile_fixincPA,unpack=True)
    print("length allprofiles fixincPA",len(allprofiles_fixincPA))
    if (len(allprofiles_fixincPA) > 3):
        (rrs_fixincPA, v_Phi_prof_fixincPA, sv_Phi_prof_fixincPA, v_R_prof_fixincPA, sv_R_prof_fixincPA) = allprofiles_fixincPA 
        DoAccr_fixIncPA=True
    else:
        (rrs_fixincPA, v_Phi_prof_fixincPA, sv_Phi_prof_fixincPA) = allprofiles_fixincPA 
        DoAccr_fixIncPA=False

    if VarOrient:
        allprofiles_fixincPA_allrads=np.loadtxt(file_profile_fixincPA_allrads,unpack=True)
        print("length allprofiles fixincPA",len(allprofiles_fixincPA_allrads))
        if (len(allprofiles_fixincPA_allrads) > 3):
            (rrs_fixincPA_allrads, v_Phi_prof_fixincPA_allrads, sv_Phi_prof_fixincPA_allrads, v_R_prof_fixincPA_allrads, sv_R_prof_fixincPA_allrads) = allprofiles_fixincPA_allrads
            DoAccrAllRads_fixIncPA=True
        else:
            (rrs_fixincPA_allrads, v_Phi_prof_fixincPA_allrads, sv_Phi_prof_fixincPA_allrads) = allprofiles_fixincPA_allrads 
            DoAccrAllRads_fixIncPA=False



nplotsy=0
if Plot_vRot_Global:
    nplotsy+=1    
    
if VarOrient:
    nplotsy+=1
    
if Plot_vRot_VarOrient:
    nplotsy+=1

if XCheckFixIncPA:
    nplotsy+=1
        
if DoAccr:
    nplotsy+=1
    
if DoAccr_fixIncPA:
    nplotsy+=1

if DoAccrAllRads and Plot_vRot_VarOrient:
    nplotsy+=1

if DoAccrAllRads and Plot_vRot_VarOrient_FixIncPA:
    nplotsy+=1

if Plot_vRot_VarOrient_FixIncPA:
    nplotsy+=1
    
        
print("nplotsy = ",nplotsy)    
figysize=nplotsy*4.
figxsize=9.

fig = plt.figure(constrained_layout=True,figsize=(figxsize,figysize))
gs = fig.add_gridspec(nplotsy, 1) #, width_ratios=[2., 1., 1.], height_ratios=[1.,1.])

    
(RunMCMC, proflist) = LogResults.load(workdir,FixPAInc=False)                

if RunMCMC:
    [r1s, r2s, rregions, incs, psis, PAs, allradsPA, allradsinc, allradspsi, PAs_MCMC, tanpsis_MCMC, incs_MCMC, allradsPAMCMC, allradsincMCMC,allradstanpsiMCMC] = proflist
else:
    [r1s, r2s, rregions, incs, psis, PAs, allradsPA, allradsinc, allradspsi] = proflist


a_min=np.min(r1s)
a_max=np.max(r2s)
print( "a_min ",a_min," a_max", a_max)



print( "psis",psis)
print( "allradspsi",allradspsi)
print( "incs",incs)

psierrors=np.zeros((2,len(psis)))
if RunMCMC:
    for ir,atanpsi in enumerate(tanpsis_MCMC):
        psierrors[1,ir]=180.*(np.arctan(atanpsi[1]+atanpsi[0])-np.arctan(atanpsi[0]))/np.pi
        psierrors[0,ir]=180.*(np.arctan(atanpsi[0])-np.arctan(atanpsi[0]-atanpsi[2]))/np.pi


incerrors=np.zeros((2,len(incs)))
if RunMCMC:
    for ir,ainc in enumerate(incs_MCMC):
        print( "ainc",ainc)
        incerrors[1,ir]=1.*ainc[1]*180./np.pi
        incerrors[0,ir]=1.*ainc[2]*180./np.pi

PAerrors=np.zeros((2,len(PAs)))
if RunMCMC:
    for ir,aPA in enumerate(PAs_MCMC):
        PAerrors[1,ir]=1.*aPA[1]
        PAerrors[0,ir]=1.*aPA[2]





if VarOrient:
    axprofile = fig.add_subplot(gs[(nplotsy-1),0])
    (ymin, ymax) = Orient.PlotOrientProfile(axprofile,rregions, PAs, allradsPA, PAerrors, incs, allradsinc,incerrors, psis, allradspsi, psierrors)
    

sigma_PA=np.std(PAs)
sigma_inc=np.std(incs)
sigma_psi=np.std(psis)
print( "Psi= %.1f +- %.1f deg" % (allradspsi,sigma_psi))
print( "Inc= %.2f +- %.2f deg" % (allradsinc,sigma_inc))
print( "PA= %.1f +- %.1f deg" % (allradsPA,sigma_PA))


print( "USING inclination ", allradsinc," for mass estimates")
cosi=np.fabs(np.cos(allradsinc*np.pi/180.))


#######################################################################
######################################################################

if DoFixIncPA:

    (RunMCMC, proflist) = LogResults.load(workdir_fixincPA,RunMCMC=RunMCMC,FixPAInc=True) 


    if RunMCMC:
        [r1s,r2s, rregions_fixincPA, psis_fixincPA, allradspsi_fixincPA, tanpsis_fixincPA_MCMC]=proflist
    else:
        [r1s,r2s, rregions_fixincPA, psis_fixincPA, allradspsi_fixincPA] = proflist


    a_min=np.min(r1s)
    a_max=np.max(r2s)
    print( "a_min ",a_min," a_max", a_max)

    print( "psis_fixincPA", psis_fixincPA)
    print( "allradspsi_fixincPA",allradspsi_fixincPA)

    psierrors_fixincPA=np.zeros((2,len(psis_fixincPA)))
    if RunMCMC:
        print( "psis_fixincPA MCMC", tanpsis_fixincPA_MCMC)
        for ir,atanpsi in enumerate(tanpsis_fixincPA_MCMC):
            psierrors_fixincPA[1,ir]=180.*(np.arctan(atanpsi[1]+atanpsi[0])-np.arctan(atanpsi[0]))/np.pi
            psierrors_fixincPA[0,ir]=180.*(np.arctan(atanpsi[0])-np.arctan(atanpsi[0]-atanpsi[2]))/np.pi

    if VarOrient:
        (ymin_fixincPA,ymax_fixincPA)=Orient.PlotOrientProfile_fixincPA(axprofile,rregions_fixincPA, psis_fixincPA, allradspsi_fixincPA, psierrors_fixincPA)

    sigma_psi_fixincPA=np.std(psis_fixincPA)
    print( "Psi fixincPA= %.1f +- %.1f deg" % (allradspsi_fixincPA,sigma_psi_fixincPA))

    
if VarOrient:

    ymin=min(ymin,ymin_fixincPA)
    ymax=max(ymax,ymax_fixincPA)

    axprofile.set_xlim(a_min,a_max)
    axprofile.set_ylim(ymin, ymax)
    axprofile.set_ylabel(r'deg')
    axprofile.set_xlabel(r'$r$ / arcsec')

    #axprofile.legend(fontsize=16)
    axprofile.legend()

    #axprofile.legend()
    axprofile.tick_params(axis='both', length = 8,  direction='in', pad=10)
    axprofile.tick_params(top='on',right='on', direction='in')
    axprofile.tick_params(which='minor', top='on', length = 4, right='on', direction='in')
    axprofile.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axprofile.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())


print( "USING inclination ", allradsinc," for mass estimates")
#cosi=np.cos(40.06*np.pi/180.)
cosi=np.fabs(np.cos(allradsinc*np.pi/180.))



######################################################################
## Correct v_phi for height over midplane 


if VarOrient:
    hhs=np.interp(rrs,rregions,np.tan(psis*np.pi/180.))
    hhs_fixincPA=np.interp(rrs_fixincPA,rregions_fixincPA,np.tan(psis_fixincPA*np.pi/180.))
else:
    hhs=np.tan(allradspsi*np.pi/180.)*np.ones(rrs.shape)

correct4midplane = ( (1. + hhs**2)**(3./4.) )
v_Phi_prof_mid = v_Phi_prof * correct4midplane
if DoFixIncPA:
    correct4midplane_fixincPA = ( (1. + hhs_fixincPA**2)**(3./4.))
    v_Phi_prof_mid_fixincPA = v_Phi_prof_fixincPA * correct4midplane_fixincPA 
    
jpos=0
if Plot_vRot_Global:
    axprofile = fig.add_subplot(gs[jpos,0])
    RotCurve.PlotV_phi(axprofile,rrs,a_min,a_max,v_Phi_prof,sv_Phi_prof,v_Phi_prof_mid,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=rgaps,label=r'global')
    jpos+=1

    
if DoFixIncPA and XCheckFixIncPA:
    axprofile = fig.add_subplot(gs[jpos,0])
    RotCurve.PlotV_phi(axprofile,rrs_fixincPA,a_min,a_max,v_Phi_prof_fixincPA,sv_Phi_prof_fixincPA,v_Phi_prof_mid_fixincPA,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=rgaps,label=r'fix $i$, PA global')
    jpos+=1


if Plot_vRot_VarOrient:
    axprofile = fig.add_subplot(gs[jpos,0])
    v_Phi_prof_mid_allrads = v_Phi_prof_allrads * correct4midplane   
    RotCurve.PlotV_phi(axprofile,rrs_allrads,a_min,a_max,v_Phi_prof_allrads,sv_Phi_prof_allrads,v_Phi_prof_mid_allrads,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=rgaps,label=r'region av.')
    jpos += 1

if DoFixIncPA and Plot_vRot_VarOrient_FixIncPA:
    axprofile = fig.add_subplot(gs[jpos,0])
    v_Phi_prof_mid_fixincPA_allrads = v_Phi_prof_fixincPA_allrads * correct4midplane_fixincPA 
    RotCurve.PlotV_phi(axprofile,rrs_allrads,a_min,a_max,v_Phi_prof_fixincPA_allrads,sv_Phi_prof_fixincPA_allrads,v_Phi_prof_mid_fixincPA_allrads,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=rgaps,label=r'fix $i$, PA region av.')
    jpos+=1


    
if (DoAccr and Plot_vRot_Global):
    axprofile = fig.add_subplot(gs[jpos,0])
    RotCurve.PlotV_R(axprofile,rrs_fixincPA,a_min,a_max,v_R_prof_fixincPA,sv_R_prof_fixincPA,ContinuumGaps=rgaps,label=r'global')
    jpos+=1

if (DoAccr_fixIncPA and Plot_vRot_Global):
    axprofile = fig.add_subplot(gs[jpos,0])
    RotCurve.PlotV_R(axprofile,rrs_fixincPA,a_min,a_max,v_R_prof_fixincPA,sv_R_prof_fixincPA,ContinuumGaps=rgaps,label=r'fix $i$, PA global')
    jpos+=1


if DoAccrAllRads and Plot_vRot_VarOrient:
    axprofile = fig.add_subplot(gs[jpos,0])
    RotCurve.PlotV_R(axprofile,rrs_allrads,a_min,a_max,v_R_prof_allrads,sv_R_prof_allrads,ContinuumGaps=rgaps,label=r'region av.')
    jpos += 1

if DoAccrAllRads_fixIncPA and Plot_vRot_VarOrient_FixIncPA:
    axprofile = fig.add_subplot(gs[jpos,0])
    RotCurve.PlotV_R(axprofile,rrs_fixincPA_allrads,a_min,a_max,v_R_prof_fixincPA_allrads,sv_R_prof_fixincPA_allrads,ContinuumGaps=rgaps,label=r'fix $i$, PA region av.')
    jpos += 1


plt.subplots_adjust(hspace=0.)

print( fileout_fig)
fig.savefig(fileout_fig, bbox_inches='tight')


