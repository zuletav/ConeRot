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


rgaps=[0.4433,  0.8575, 1.3923, 2.3]


#Plot_vRot_Global_FixIncPA = False,

def execfig(workdir, filename_source, bmaj=0.083, distance=101.50, a_min=-1,a_max=-1,WithComparData=False,WithComparRadTWind=False,rgaps=False,fileout_fig='default',Plot_vRot_Global=False, Plot_vRot_VarOrient=False, VarOrient=True, Plot_vRot_VarOrient_FixIncPA = True, PlotVarPAinc=True, ForceGlobalOrient=False, Force_allradsPA=0., Force_allradsinc=0.):
    
    XCheckFixIncPA=False # cross check that the rot curve is the same for global optim  and global init optim for  fix PA and inc (should be the same PA, inc, psi)

    Regions=True
    
    #nametag='im_g_v0'
    
    inbasename=os.path.basename(filename_source)
    inbasename=re.sub('\.fits','',inbasename)
    print("inbasename",inbasename)
    nametag=inbasename

    #workdir='work_gmom_8_werrmap_MCMC_pix2b/'
    #bmaj=0.104 # arcsec

    #Plot_vRot_Global=False  # if true will plot rotation curve for global (fixed) orientation  (so not averaged over Regions)
    #Plot_vRot_VarOrient=False  # if true will plot rotation curve for variable orientation too (so averaged over Regions)
    #VarOrient=True # will plot the variable orientation obtained from the Regions
    ##
    #Plot_vRot_Global_FixIncPA = False
    #Plot_vRot_VarOrient_FixIncPA = True
    #PlotVarPAinc=True


    XCheckFixIncPA=False # cross check that the rot curve is the same for global optim  and global init optim for  fix PA and inc (should be the same PA, inc, psi)

    if WithComparData:
        ######################################################################
        # Rich Teague data
        vcube_RT=np.load('data_0/rteague/velocity_profiles_a.npy')
        print("vcube_RT",vcube_RT.shape)      
        #vcube_RT (4, 3, 159)

        inc_RichTeague=(180.-46.7)*np.pi/180.

        
        rrs_RT=vcube_RT[0,0,:]
        v_phi_RT=1E-3*vcube_RT[1,0,:]/np.sin(inc_RichTeague)
        v_rad_RT=-1E-3*vcube_RT[2,0,:]/np.sin(inc_RichTeague)
        v_z_RT=(1E-3*vcube_RT[3,0,:]-5.763)/np.cos(inc_RichTeague)
        
        s_v_phi_RT_lo=1E-3*vcube_RT[1,1,:]/np.sin(inc_RichTeague)
        s_v_rad_RT_lo=1E-3*vcube_RT[2,1,:]/np.sin(inc_RichTeague)
        s_v_z_RT_lo=1E-3*vcube_RT[3,1,:]/np.fabs(np.cos(inc_RichTeague))
        
        s_v_phi_RT_up=1E-3*vcube_RT[1,2,:]/np.sin(inc_RichTeague)
        s_v_rad_RT_up=1E-3*vcube_RT[2,2,:]/np.sin(inc_RichTeague)
        s_v_z_RT_up=1E-3*vcube_RT[3,2,:]/np.fabs(np.cos(inc_RichTeague))
        ######################################################################

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


    if (fileout_fig == 'default'):
        fileout_fig=workdir_fixincPA+'fig_rotorient_surfandmid_linear_full.pdf'
        
    file_profile=workdir+nametag+'_radial_profile.dat'
    file_profile_allrads=workdir+nametag+'_allrads_radial_profile.dat'
    file_profile_fixincPA=workdir_fixincPA+nametag+'_radial_profile.dat'
    file_profile_fixincPA_allrads=workdir_fixincPA+nametag+'_allrads_radial_profile.dat'

    print("loading file_profiles",file_profile," ",file_profile_fixincPA)

    DoMerid=False
    DoAccr=False
    if os.path.isfile(file_profile):
        
        allprofiles=np.loadtxt(file_profile,unpack=True)
        print("length allprofiles",len(allprofiles))
        
        if (len(allprofiles) > 5):
            (rrs, v_Phi_prof, sv_Phi_prof, v_R_prof, sv_R_prof, v_z_prof, sv_z_prof) = allprofiles # np.loadtxt(file_profile,unpack=True)
            DoMerid=True
            print("WARNING: found DoMerid optimization for global optim (no regions)  but will plot v_R averaged over regions (allrads)")
        elif (len(allprofiles) > 3):
            (rrs, v_Phi_prof, sv_Phi_prof, v_R_prof, sv_R_prof) = allprofiles # np.loadtxt(file_profile,unpack=True)
            DoAccr=True
            print("WARNING: found DoAccr optimization for global optim (no regions)  but will plot v_R averaged over regions (allrads)")
        else:
            (rrs, v_Phi_prof, sv_Phi_prof) = np.loadtxt(file_profile,unpack=True)

    DoMeridAllRads=False
    DoAccrAllRads=False
    if os.path.isfile(file_profile_allrads):

        if VarOrient:
            allprofiles_allrads=np.loadtxt(file_profile_allrads,unpack=True)
            print("length allprofiles",len(allprofiles_allrads))
            if (len(allprofiles) > 5):
                (rrs_allrads, v_Phi_prof_allrads, sv_Phi_prof_allrads, v_R_prof_allrads, sv_R_prof_allrads, v_z_prof_allrads, sv_z_prof_allrads) = allprofiles_allrads # np.loadtxt(file_profile,unpack=True)
                DoMeridAllRads=True
                print("WARNING: found DoMerid optimization with variable PA and inc, could be degenerate")
            elif (len(allprofiles_allrads) > 3):
                (rrs_allrads, v_Phi_prof_allrads, sv_Phi_prof_allrads, v_R_prof_allrads, sv_R_prof_allrads) = allprofiles_allrads
                DoAccrAllRads=True
                print("WARNING: found DoAccr optimization with variable PA and inc, could be degenerate")
            else:
                (rrs_allrads, v_Phi_prof_allrads, sv_Phi_prof_allrads) = np.loadtxt(file_profile_allrads,unpack=True)
                DoAccrAllRads=False


    DoAccr_fixIncPA=False
    DoMerid_fixIncPA=False
    DoAccrAllRads_fixIncPA=False
    DoMeridAllRads_fixIncPA=False
    if os.path.isfile(file_profile_fixincPA):
        if DoFixIncPA:
            allprofiles_fixincPA=np.loadtxt(file_profile_fixincPA,unpack=True)
            print("length allprofiles fixincPA",len(allprofiles_fixincPA))
            if (len(allprofiles_fixincPA) > 5):
                (rrs_fixincPA, v_Phi_prof_fixincPA, sv_Phi_prof_fixincPA, v_R_prof_fixincPA, sv_R_prof_fixincPA, v_z_prof_fixincPA, sv_z_prof_fixincPA) = allprofiles_fixincPA # np.loadtxt(file_profile,unpack=True)
                DoMerid_fixIncPA=True
            elif (len(allprofiles_fixincPA) > 3):
                (rrs_fixincPA, v_Phi_prof_fixincPA, sv_Phi_prof_fixincPA, v_R_prof_fixincPA, sv_R_prof_fixincPA) = allprofiles_fixincPA 
                DoAccr_fixIncPA=True
            else:
                (rrs_fixincPA, v_Phi_prof_fixincPA, sv_Phi_prof_fixincPA) = allprofiles_fixincPA 
                DoAccr_fixIncPA=False

            if VarOrient:
                allprofiles_fixincPA_allrads=np.loadtxt(file_profile_fixincPA_allrads,unpack=True)
                print("length allprofiles fixincPA",len(allprofiles_fixincPA_allrads))
                if (len(allprofiles_fixincPA_allrads) > 5):
                    (rrs_fixincPA_allrads, v_Phi_prof_fixincPA_allrads, sv_Phi_prof_fixincPA_allrads, v_R_prof_fixincPA_allrads, sv_R_prof_fixincPA_allrads, v_z_prof_fixincPA_allrads, sv_z_prof_fixincPA_allrads) = allprofiles_fixincPA_allrads
                    DoMeridAllRads_fixIncPA=True
                elif (len(allprofiles_fixincPA_allrads) > 3):
                    (rrs_fixincPA_allrads, v_Phi_prof_fixincPA_allrads, sv_Phi_prof_fixincPA_allrads, v_R_prof_fixincPA_allrads, sv_R_prof_fixincPA_allrads) = allprofiles_fixincPA_allrads
                    DoAccrAllRads_fixIncPA=True
                else:
                    (rrs_fixincPA_allrads, v_Phi_prof_fixincPA_allrads, sv_Phi_prof_fixincPA_allrads) = allprofiles_fixincPA_allrads 


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

    if DoAccr_fixIncPA and Plot_vRot_Global:
        nplotsy+=1

    if DoMerid and Plot_vRot_Global:
        nplotsy+=2
    elif DoAccr and Plot_vRot_Global:
        nplotsy+=1

    if DoMerid_fixIncPA and Plot_vRot_Global:
        nplotsy+=2
    elif DoAccr_fixIncPA and Plot_vRot_Global:
        nplotsy+=1


        
    if DoMeridAllRads and Plot_vRot_VarOrient:
        nplotsy+=2
    elif DoAccrAllRads and Plot_vRot_VarOrient:
        nplotsy+=1

    if DoMeridAllRads_fixIncPA and Plot_vRot_VarOrient_FixIncPA:
        nplotsy+=2
    elif DoAccrAllRads_fixIncPA and Plot_vRot_VarOrient_FixIncPA:
        nplotsy+=1

    if Plot_vRot_VarOrient_FixIncPA:
        nplotsy+=1


    print("nplotsy = ",nplotsy)    
    figysize=nplotsy*4*7/9.
    figxsize=7.

    fig = plt.figure(constrained_layout=True,figsize=(figxsize,figysize))
    gs = fig.add_gridspec(nplotsy, 1) #, width_ratios=[2., 1., 1.], height_ratios=[1.,1.])

    if VarOrient:
        axprofile = fig.add_subplot(gs[(nplotsy-1),0])

    RunMCMC=False
    GlobalOrientProf=False
    if os.path.isdir(workdir):

        GlobalOrientProf=True
        
        (RunMCMC, proflist) = LogResults.load(workdir,FixPAInc=False)                
        
        if RunMCMC:
            print("picked up MCMC run")
            [r1s, r2s, rregions, incs, psis, PAs, allradsPA, allradsinc, allradspsi, PAs_MCMC, tanpsis_MCMC, incs_MCMC, allradsPAMCMC, allradsincMCMC,allradspsiMCMC] = proflist

            allradsinc=allradsincMCMC
            allradsPA=allradsPAMCMC
            allradspsi=allradspsiMCMC

        else:
            [r1s, r2s, rregions, incs, psis, PAs, allradsPA, allradsinc, allradspsi] = proflist

    
        if ForceGlobalOrient:
            print(">>>>>> Plotting  Fix Orient, with forced global orientation, at PA=",Force_allradsPA," inc=",Force_allradsinc)
            allradsPA=Force_allradsPA
            allradsinc=Force_allradsinc*180./np.pi

        if (a_min < 0):
            a_min=np.min(r1s)
        if (a_max <0):
            a_max=np.max(r2s)

        print( "a_min ",a_min," a_max", a_max)



        #print( "psis",psis)
        #print( "allradspsi",allradspsi)
        #print( "incs",incs)

        psierrors=np.zeros((2,len(psis)))
        if RunMCMC:
            for ir,atanpsi in enumerate(tanpsis_MCMC):
                psis[ir]=180.*np.arctan(atanpsi[0])/np.pi

                psierrors[1,ir]=180.*(np.arctan(atanpsi[1]+atanpsi[0])-np.arctan(atanpsi[0]))/np.pi
                psierrors[0,ir]=180.*(np.arctan(atanpsi[0])-np.arctan(atanpsi[0]-atanpsi[2]))/np.pi


        incerrors=np.zeros((2,len(incs)))
        if RunMCMC:
            for ir,ainc in enumerate(incs_MCMC):
                #print( "ainc",ainc)
                incs[ir]=ainc[0]*180./np.pi
                incerrors[1,ir]=1.*ainc[1]*180./np.pi
                incerrors[0,ir]=1.*ainc[2]*180./np.pi

        PAerrors=np.zeros((2,len(PAs)))
        if RunMCMC:
            for ir,aPA in enumerate(PAs_MCMC):
                PAs[ir]=aPA[0]
                PAerrors[1,ir]=1.*aPA[1]
                PAerrors[0,ir]=1.*aPA[2]




        if PlotVarPAinc:
            print("rregions",rregions)
            print("PAs",PAs)

            (ymin, ymax) = Orient.PlotOrientProfile(axprofile,rregions, PAs, allradsPA, PAerrors, incs, allradsinc,incerrors, psis, allradspsi, psierrors)


        sigma_PA=np.std(PAs)
        sigma_inc=np.std(incs)
        sigma_psi=np.std(psis)
        print( "Psi= %.1f +- %.1f deg" % (allradspsi,sigma_psi))
        print( "Inc= %.2f +- %.2f deg" % (allradsinc,sigma_inc))
        print( "PA= %.1f +- %.1f deg" % (allradsPA,sigma_PA))


        print( "USING inclination ", allradsinc," for mass estimates")
        cosi=np.fabs(np.cos(allradsinc*np.pi/180.))

        BackSide=False
        if (allradspsi < 0.):
            BackSide=True

    #######################################################################
    ######################################################################

    if DoFixIncPA:

        (RunMCMC, proflist) = LogResults.load(workdir_fixincPA,RunMCMC=RunMCMC,FixPAInc=True) 

        print("proflist",proflist)

        if RunMCMC:
            [r1s,r2s, rregions_fixincPA, psis_fixincPA, allradspsi_fixincPA, tanpsis_fixincPA_MCMC]=proflist
        else:
            [r1s,r2s, rregions_fixincPA, psis_fixincPA, allradspsi_fixincPA] = proflist

        if (len(rregions_fixincPA)==0):
            Regions=False
            VarOrient=False
            if ( (a_min < 0) or (a_max <0)):
                sys.exit("need to specify a_min and a_max in call to execfig")
            
        if (a_min < 0):
            a_min=np.min(r1s)
        if (a_max <0):
            a_max=np.max(r2s)
         
        print( "a_min ",a_min," a_max", a_max)

        #a_min=np.min(r1s)
        #a_max=np.max(r2s)
        print( "a_min ",a_min," a_max", a_max)

        #print( "psis_fixincPA", psis_fixincPA)
        #print( "allradspsi_fixincPA",allradspsi_fixincPA)

        psierrors_fixincPA=np.zeros((2,len(psis_fixincPA)))
        if RunMCMC:
            #print( "psis_fixincPA MCMC", tanpsis_fixincPA_MCMC)
            for ir,atanpsi in enumerate(tanpsis_fixincPA_MCMC):
                psis_fixincPA[ir]=180.*np.arctan(atanpsi[0])/np.pi
                psierrors_fixincPA[1,ir]=180.*(np.arctan(atanpsi[1]+atanpsi[0])-np.arctan(atanpsi[0]))/np.pi
                psierrors_fixincPA[0,ir]=180.*(np.arctan(atanpsi[0])-np.arctan(atanpsi[0]-atanpsi[2]))/np.pi

        if VarOrient:
            (ymin_fixincPA,ymax_fixincPA)=Orient.PlotOrientProfile_fixincPA(axprofile,rregions_fixincPA, psis_fixincPA, allradspsi_fixincPA, psierrors_fixincPA)

        sigma_psi_fixincPA=np.std(psis_fixincPA)
        print( "Psi fixincPA= %.1f +- %.1f deg" % (allradspsi_fixincPA,sigma_psi_fixincPA))


    if VarOrient:

        if not PlotVarPAinc:
            #ymin=min(ymin,ymin_fixincPA)
            #ymax=max(ymax,ymax_fixincPA)
            #else:
            ymin=ymin_fixincPA
            ymax=ymax_fixincPA

        else:

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


    if ForceGlobalOrient:
        print(">>>>>> Plotting  Fix Orient, with forced global orientation, at PA=",Force_allradsPA," inc=",Force_allradsinc)
        allradsPA=Force_allradsPA
        allradsinc=Force_allradsinc*180./np.pi

    print( "USING inclination ", allradsinc," for mass estimates")
    #cosi=np.cos(40.06*np.pi/180.)
    cosi=np.fabs(np.cos(allradsinc*np.pi/180.))

    BackSide=False
    if (allradspsi_fixincPA < 0.):
        BackSide=True

    ######################################################################
    ## Correct v_phi for height over midplane 


    if VarOrient:
        if GlobalOrientProf:
            hhs=np.interp(rrs,rregions,np.tan(psis*np.pi/180.))
        if DoFixIncPA:
            hhs_fixincPA=np.interp(rrs_fixincPA,rregions_fixincPA,np.tan(psis_fixincPA*np.pi/180.))
    else:
        hhs_fixincPA=np.tan(allradspsi_fixincPA*np.pi/180.)*np.ones(rrs_fixincPA.shape)

    if GlobalOrientProf:
        correct4midplane = ( (1. + hhs**2)**(3./4.) )
        v_Phi_prof_mid = v_Phi_prof * correct4midplane
        
    if DoFixIncPA:
        correct4midplane_fixincPA = ( (1. + hhs_fixincPA**2)**(3./4.))
        v_Phi_prof_mid_fixincPA = v_Phi_prof_fixincPA * correct4midplane_fixincPA 

    jpos=0
    if Plot_vRot_Global:
        axprofile = fig.add_subplot(gs[jpos,0])
        #RotCurve.PlotV_phi(axprofile,rrs_fixincPA,a_min,a_max,v_Phi_prof,sv_Phi_prof,v_Phi_prof_mid,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=rgaps,label=r'global')
        RotCurve.PlotV_phi(axprofile,rrs_fixincPA,a_min,a_max,v_Phi_prof_fixincPA,sv_Phi_prof_fixincPA,v_Phi_prof_mid_fixincPA,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=rgaps,label=r'')
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
        alabel=r'fix $i$, PA region av.'
        #alabel=''
        v_Phi_prof_mid_fixincPA_allrads = v_Phi_prof_fixincPA_allrads * correct4midplane_fixincPA 
        (vymin,vymax)=RotCurve.PlotV_phi(axprofile,rrs_fixincPA,a_min,a_max,v_Phi_prof_fixincPA_allrads,sv_Phi_prof_fixincPA_allrads,v_Phi_prof_mid_fixincPA_allrads,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=rgaps,label=alabel)

        #print("v_Phi_prof_fixincPA_allrads*np.sqrt(rrs_allrads)/np.sqrt(a_max)",v_Phi_prof_fixincPA_allrads*np.sqrt(rrs_allrads)/np.sqrt(a_max))
                
        if WithComparData:
            # Rich Teague data


            #print("v_phi_RT*np.sqrt(rrs_RT)/np.sqrt(a_max)",v_phi_RT*np.sqrt(rrs_RT)/np.sqrt(a_max))
            axprofile.plot(rrs_RT,v_phi_RT*np.sqrt(rrs_RT)/np.sqrt(a_max),color='C1',linewidth=1.5,linestyle='solid',label=r'$v_\phi$ eddy')
            axprofile.fill_between(rrs_RT,(v_phi_RT+s_v_phi_RT_up)*np.sqrt(rrs_RT)/np.sqrt(a_max),(v_phi_RT-s_v_phi_RT_lo)*np.sqrt(rrs_RT)/np.sqrt(a_max),  lw=0.1,color='C1', alpha=0.2, interpolate=True) #, step='mid'

            maskrange=( (rrs_RT > a_min) & (rrs_RT < a_max))
            vcomparmin=np.min((v_phi_RT[maskrange]-s_v_phi_RT_lo[maskrange])*np.sqrt(rrs_RT[maskrange])/np.sqrt(a_max))
            vcomparmax=np.max((v_phi_RT[maskrange]+s_v_phi_RT_up[maskrange])*np.sqrt(rrs_RT[maskrange])/np.sqrt(a_max))
            vymin=min(vymin,vcomparmin)
            vymax=max(vymax,vcomparmax)
            axprofile.set_ylim(vymin,vymax)

            
        axprofile.legend(loc='upper right')
        
        jpos+=1

    if (DoMerid and Plot_vRot_Global):
        axprofile = fig.add_subplot(gs[jpos,0])
        RotCurve.PlotV_z(axprofile,rrs_fixincPA,a_min,a_max,v_z_prof_fixincPA,sv_z_prof_fixincPA,BackSide=BackSide,ContinuumGaps=rgaps,label=r'global')
        jpos+=1
        axprofile = fig.add_subplot(gs[jpos,0])
        RotCurve.PlotV_R(axprofile,rrs_fixincPA,a_min,a_max,v_R_prof_fixincPA,sv_R_prof_fixincPA,ContinuumGaps=rgaps,label=r'global')
        jpos+=1
    elif (DoAccr and Plot_vRot_Global):
        axprofile = fig.add_subplot(gs[jpos,0])
        RotCurve.PlotV_R(axprofile,rrs_fixincPA,a_min,a_max,v_R_prof_fixincPA,sv_R_prof_fixincPA,ContinuumGaps=rgaps,label=r'global')
        jpos+=1


    if (DoMerid_fixIncPA and Plot_vRot_Global):
        axprofile = fig.add_subplot(gs[jpos,0])
        RotCurve.PlotV_z(axprofile,rrs_fixincPA,a_min,a_max,v_z_prof_fixincPA,sv_z_prof_fixincPA,BackSide=BackSide,ContinuumGaps=rgaps,label=r'')
        jpos+=1
        axprofile = fig.add_subplot(gs[jpos,0])
        RotCurve.PlotV_R(axprofile,rrs_fixincPA,a_min,a_max,v_R_prof_fixincPA,sv_R_prof_fixincPA,ContinuumGaps=rgaps,label=r'')
        jpos+=1
    elif (DoAccr_fixIncPA and Plot_vRot_Global):
        axprofile = fig.add_subplot(gs[jpos,0])
        RotCurve.PlotV_R(axprofile,rrs_fixincPA,a_min,a_max,v_R_prof_fixincPA,sv_R_prof_fixincPA,ContinuumGaps=rgaps,label=r'')
        jpos+=1


    if DoMeridAllRads and Plot_vRot_VarOrient:
        axprofile = fig.add_subplot(gs[jpos,0])
        RotCurve.PlotV_z(axprofile,rrs_allrads,a_min,a_max,v_z_prof_allrads,sv_z_prof_allrads,BackSide=BackSide,ContinuumGaps=rgaps,label=r'region av.')
        jpos += 1
        axprofile = fig.add_subplot(gs[jpos,0])
        RotCurve.PlotV_R(axprofile,rrs_allrads,a_min,a_max,v_R_prof_allrads,sv_R_prof_allrads,ContinuumGaps=rgaps,label=r'region av.')
        jpos += 1
    elif DoAccrAllRads and Plot_vRot_VarOrient:
        axprofile = fig.add_subplot(gs[jpos,0])
        RotCurve.PlotV_R(axprofile,rrs_allrads,a_min,a_max,v_R_prof_allrads,sv_R_prof_allrads,ContinuumGaps=rgaps,label=r'region av.')
        jpos += 1


    if DoMeridAllRads_fixIncPA and Plot_vRot_VarOrient_FixIncPA:
        axprofile = fig.add_subplot(gs[jpos,0])
        alabel=r'fix $i$, PA region av.'
        #alabel=''
        (vymin,vymax)=RotCurve.PlotV_z(axprofile,rrs_fixincPA_allrads,a_min,a_max,v_z_prof_fixincPA_allrads,sv_z_prof_fixincPA_allrads,BackSide=BackSide,ContinuumGaps=rgaps,label=alabel)
        if WithComparData:
            # Rich Teague data
            axprofile.plot(rrs_RT,v_z_RT*np.sqrt(rrs_RT)/np.sqrt(a_max),color='C1',linewidth=1.5,linestyle='solid',label=r'$v_z$ eddy')
            axprofile.fill_between(rrs_RT,(v_z_RT+s_v_z_RT_up)*np.sqrt(rrs_RT)/np.sqrt(a_max),(v_z_RT-s_v_z_RT_lo)*np.sqrt(rrs_RT)/np.sqrt(a_max),  lw=0.1,color='C1', alpha=0.2, interpolate=True) #, step='mid'


            maskrange=( (rrs_RT > a_min) & (rrs_RT < a_max))
            vcomparmin=np.min((v_z_RT[maskrange]-s_v_z_RT_lo[maskrange])*np.sqrt(rrs_RT[maskrange])/np.sqrt(a_max))
            vcomparmax=np.max((v_z_RT[maskrange]+s_v_z_RT_up[maskrange])*np.sqrt(rrs_RT[maskrange])/np.sqrt(a_max))
            vymin=min(vymin,vcomparmin)
            vymax=max(vymax,vcomparmax)

            
            print(">>>>> compar v_z ::",vymin,vymax)
            

            
            axprofile.set_ylim(vymin,vymax)


        if WithComparRadTWind:
            # Comparing with wind model

            zzs=hhs*rrs
            v_z_RadTWind=np.zeros(rrs.shape)
            v_R_RadTWind=np.zeros(rrs.shape)
            distance=100.
            Rhopp=rrs*distance
            Zpp=zzs*distance
            h_c_outer=0.06
            Rcavdust=20.
            v_term=1E4 # cm/s
            flaring_outer=0.15
            h_outer = h_c_outer * (rrs/Rcavdust)**(flaring_outer)
            H_outer=h_outer*Rhopp
            for iRhopp in list(range(len(Rhopp))):
                aRhopp=Rhopp[iRhopp]
                aZpp=Zpp[iRhopp]
                v_z_RadTWind[iRhopp]=0.
                v_R_RadTWind[iRhopp]=0.
                aH_outer=H_outer[iRhopp]
                print("aRhopp ",aRhopp,"aZpp",aZpp,"aH_outer",aH_outer)
                if ((aRhopp > Rcavdust) and (np.fabs(aZpp) > aH_outer)):
                    vwind=v_term * (1. - (Rcavdust/aRhopp)**2) * (1.-(aH_outer/aZpp)**2)
                    avrpp = vwind*np.cos(np.pi/4.) # *np.cos(np.pi/3.) #cylindrical
                    if (aZpp > 0.):
                        avzpp = vwind*np.sin(np.pi/4.)
                    else:
                        avzpp = -vwind*np.sin(np.pi/4.)

                    v_z_RadTWind[iRhopp]=avzpp/1.0E5 # km/s
                    v_R_RadTWind[iRhopp]=avrpp/1.0E5 # km/s

            if BackSide:
                v_z_RadTWind *= -1
                
            axprofile.plot(rrs,v_z_RadTWind*np.sqrt(rrs)/np.sqrt(a_max),color='C3',linewidth=1.5,linestyle='solid',label=r'$v_z$ RT')


        axprofile.legend(loc='lower right')
                
            
        jpos += 1

        axprofile = fig.add_subplot(gs[jpos,0])
        alabel=r'fix $i$, PA region av.'
        # alabel=''
        (vymin,vymax)=RotCurve.PlotV_R(axprofile,rrs_fixincPA_allrads,a_min,a_max,v_R_prof_fixincPA_allrads,sv_R_prof_fixincPA_allrads,ContinuumGaps=rgaps,label=alabel)
        if WithComparData:
            # Rich Teague data
            axprofile.plot(rrs_RT,v_rad_RT*np.sqrt(rrs_RT)/np.sqrt(a_max),color='C1',linewidth=1.5,linestyle='solid',label=r'$v_R$ eddy')
            axprofile.fill_between(rrs_RT,(v_rad_RT+s_v_rad_RT_up)*np.sqrt(rrs_RT)/np.sqrt(a_max),(v_rad_RT-s_v_rad_RT_lo)*np.sqrt(rrs_RT)/np.sqrt(a_max),  lw=0.1,color='C1', alpha=0.2, interpolate=True) #, step='mid'
            
            maskrange=( (rrs_RT > a_min) & (rrs_RT < a_max))
            vcomparmin=np.min((v_rad_RT[maskrange]-s_v_rad_RT_lo[maskrange])*np.sqrt(rrs_RT[maskrange])/np.sqrt(a_max))
            vcomparmax=np.max((v_rad_RT[maskrange]+s_v_rad_RT_up[maskrange])*np.sqrt(rrs_RT[maskrange])/np.sqrt(a_max))
            vymin=min(vymin,vcomparmin)
            vymax=max(vymax,vcomparmax)

            axprofile.set_ylim(vymin,vymax)

        if WithComparRadTWind:
            axprofile.plot(rrs,v_R_RadTWind*np.sqrt(rrs)/np.sqrt(a_max),color='C4',linewidth=1.5,linestyle='solid',label=r'$v_z$ RT')


            
        axprofile.legend(loc='lower right')
        jpos += 1
    elif DoAccrAllRads_fixIncPA and Plot_vRot_VarOrient_FixIncPA:
        axprofile = fig.add_subplot(gs[jpos,0])
        alabel=r'fix $i$, PA region av.'
        #alabel=''
        RotCurve.PlotV_R(axprofile,rrs_fixincPA_allrads,a_min,a_max,v_R_prof_fixincPA_allrads,sv_R_prof_fixincPA_allrads,ContinuumGaps=rgaps,label=alabel)
        if WithComparData:
            # Rich Teague data
            axprofile.plot(rrs_RT,v_rad_RT*np.sqrt(rrs_RT)/np.sqrt(a_max),color='C1',linewidth=1.5,linestyle='solid',label=r'$v_R$ eddy')
            axprofile.fill_between(rrs_RT,(v_rad_RT+s_v_rad_RT_up)*np.sqrt(rrs_RT)/np.sqrt(a_max),(v_rad_RT-s_v_rad_RT_lo)*np.sqrt(rrs_RT)/np.sqrt(a_max),  lw=0.1,color='C1', alpha=0.2, interpolate=True) #, step='mid'

        axprofile.legend()
        jpos += 1


    plt.subplots_adjust(hspace=0.)

    print( fileout_fig)
    fig.savefig(fileout_fig, bbox_inches='tight')

