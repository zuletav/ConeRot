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

HOME = os.environ.get('HOME')
print("HOME: ", HOME)
include_path = HOME + '/common/python/include/'
sys.path.append(include_path)

import ConeRot.RotOrient.LogResults as LogResults
import ConeRot.RotOrient.StellarMass as StellarMass
import ConeRot.RotOrient.RotCurve as RotCurve
import ConeRot.RotOrient.Orient as Orient

rgaps = [0.4433, 0.8575, 1.3923, 2.3]

#Plot_vRot_Global_FixIncPA = False,


def execfig(workdir,
            filename_source,
            bmaj=0.083,
            distance=101.50,
            a_min=-1,
            a_max=-1,
            WithComparData=False,
            WithComparRadTWind=False,
            rgaps=False,
            fileout_fig='default',
            Plot_vRot_Global=False,
            Plot_vRot_VarOrient=False,
            VarOrient=True,
            PlotVarOrient=True,
            Plot_vRot_VarOrient_FixIncPA=True,
            PlotVarPAinc=True,
            ForceGlobalOrient=False,
            Force_allradsPA=0.,
            Force_allradsinc=0.,
            alabel='default',
            RadialScaling=True,
            title='',
            DoAUBar=False,
            MidPlaneExtrapol=True):

    XCheckFixIncPA = False  # cross check that the rot curve is the same for global optim  and global init optim for  fix PA and inc (should be the same PA, inc, psi)

    Regions = True

    #nametag='im_g_v0'

    inbasename = os.path.basename(filename_source)
    inbasename = re.sub('\.fits', '', inbasename)
    print("inbasename", inbasename)
    nametag = inbasename

    #workdir='work_gmom_8_werrmap_MCMC_pix2b/'
    #bmaj=0.104 # arcsec

    #Plot_vRot_Global=False  # if true will plot rotation curve for global (fixed) orientation  (so not averaged over Regions)
    #Plot_vRot_Global_FixIncPA  # same as Plot_vRot_Global but for a run with fixed inc. and PA
    #Plot_vRot_VarOrient=False  # if true will plot rotation curve for variable orientation too (so averaged over Regions)
    #Plot_vRot_VarOrient_FixIncPA=False same as Plot_vRot_VarOrient but for fixed PA inc

    #VarOrient=True # will load the rotation curve data for  the variable orientation PA, inc, tanpsi obtained from the Regions - allrads labels
    #PlotVarOrient=True # will plot the variable orientation PA, inc, tanpsi obtained from the Regions
    #PlotVarPAinc=True # will plot the variable PA and inc obtained from the Regions

    if not PlotVarOrient:
        PlotVarPAinc = False

    XCheckFixIncPA = False  # cross check that the rot curve is the same for global optim  and global init optim for  fix PA and inc (should be the same PA, inc, psi)

    if WithComparData:
        ######################################################################
        # Rich Teague data 'data_0/rteague/velocity_profiles_a.npy'
        vcube_ComparData = np.load(WithComparData)
        print("vcube_ComparData", vcube_ComparData.shape)
        #vcube_ComparData (4, 3, 159)

        inc_RichTeague = (180. - 46.7) * np.pi / 180.

        rrs_ComparData = vcube_ComparData[0, 0, :]
        v_phi_ComparData = 1E-3 * vcube_ComparData[1, 0, :] / np.sin(
            inc_RichTeague)
        v_rad_ComparData = -1E-3 * vcube_ComparData[2, 0, :] / np.sin(
            inc_RichTeague)
        v_z_ComparData = (1E-3 * vcube_ComparData[3, 0, :] -
                          5.763) / np.cos(inc_RichTeague)

        s_v_phi_ComparData_lo = 1E-3 * vcube_ComparData[1, 1, :] / np.sin(
            inc_RichTeague)
        s_v_rad_ComparData_lo = 1E-3 * vcube_ComparData[2, 1, :] / np.sin(
            inc_RichTeague)
        s_v_z_ComparData_lo = 1E-3 * vcube_ComparData[3, 1, :] / np.fabs(
            np.cos(inc_RichTeague))

        s_v_phi_ComparData_up = 1E-3 * vcube_ComparData[1, 2, :] / np.sin(
            inc_RichTeague)
        s_v_rad_ComparData_up = 1E-3 * vcube_ComparData[2, 2, :] / np.sin(
            inc_RichTeague)
        s_v_z_ComparData_up = 1E-3 * vcube_ComparData[3, 2, :] / np.fabs(
            np.cos(inc_RichTeague))
        ######################################################################

    #######################################################################
    # no further user input
    ######################################################################

    #matplotlib.rc('text', usetex=True)
    #matplotlib.rc('font', family='arial')
    #matplotlib.rc('font', family='helvetica')
    matplotlib.rc('font', family='arial')
    matplotlib.rcParams.update({'font.size': 12})

    workdir_fixincPA = re.sub('/$', '_fixPAinc/', workdir)
    DoFixIncPA = False

    print("testing workdir_fixincPA", workdir_fixincPA)

    if os.path.isdir(workdir_fixincPA):
        DoFixIncPA = True
    else:
        print("not found, no workdir_fixincPA", workdir_fixincPA)
        DoFixIncPA = False

    VisibleXaxis_V_z = False
    VisibleXaxis_V_R = True
    if PlotVarOrient:
        VisibleXaxis_V_z = False
        VisibleXaxis_V_R = False

    if (fileout_fig == 'default'):
        if PlotVarOrient:
            if DoFixIncPA:
                fileouttag = workdir_fixincPA + 'fig_rotorient_surfandmid_linear_full'
            else:
                fileouttag = workdir + 'fig_rotorient_surfandmid_linear_full'
        else:
            if DoFixIncPA:
                fileouttag = workdir_fixincPA + 'fig_rot_surfandmid_linear_full'
            else:
                fileouttag = workdir + 'fig_rot_surfandmid_linear_full'

        if DoAUBar:
            fileouttag += '_wAUbar'
        if RadialScaling:
            fileout_fig = fileouttag + '.pdf'
        else:
            fileout_fig = fileouttag + '_noscale.pdf'

    file_profile = workdir + nametag + '_radial_profile.dat'
    file_profile_allrads = workdir + nametag + '_allrads_radial_profile.dat'
    file_profile_fixincPA = workdir_fixincPA + nametag + '_radial_profile.dat'
    file_profile_fixincPA_allrads = workdir_fixincPA + nametag + '_allrads_radial_profile.dat'  ### _allrads_ : averaged over regions

    DoMerid = False
    DoAccr = False
    if os.path.isfile(file_profile):

        allprofiles = np.loadtxt(file_profile, unpack=True)
        print("length allprofiles", len(allprofiles))

        if (len(allprofiles) > 5):
            (rrs, v_Phi_prof, sv_Phi_prof, v_R_prof, sv_R_prof, v_z_prof,
             sv_z_prof) = allprofiles  # np.loadtxt(file_profile,unpack=True)
            DoMerid = True
            print(
                "WARNING: found DoMerid optimization for global optim (no regions)  but will plot v_R averaged over regions (allrads)"
            )
        elif (len(allprofiles) > 3):
            (rrs, v_Phi_prof, sv_Phi_prof, v_R_prof,
             sv_R_prof) = allprofiles  # np.loadtxt(file_profile,unpack=True)
            DoAccr = True
            print(
                "WARNING: found DoAccr optimization for global optim (no regions)  but will plot v_R averaged over regions (allrads)"
            )
        else:
            (rrs, v_Phi_prof, sv_Phi_prof) = np.loadtxt(file_profile,
                                                        unpack=True)

    DoMeridAllRads = False  # averaged over regions
    DoAccrAllRads = False
    if os.path.isfile(file_profile_allrads):

        if VarOrient:
            allprofiles_allrads = np.loadtxt(file_profile_allrads, unpack=True)
            print("length allprofiles", len(allprofiles_allrads))
            if (len(allprofiles) > 5):
                (
                    rrs_allrads, v_Phi_prof_allrads, sv_Phi_prof_allrads,
                    v_R_prof_allrads, sv_R_prof_allrads, v_z_prof_allrads,
                    sv_z_prof_allrads
                ) = allprofiles_allrads  # np.loadtxt(file_profile,unpack=True)
                DoMeridAllRads = True
                print(
                    "WARNING: found DoMerid optimization with variable PA and inc, could be degenerate - AllRads"
                )
            elif (len(allprofiles_allrads) > 3):
                (rrs_allrads, v_Phi_prof_allrads, sv_Phi_prof_allrads,
                 v_R_prof_allrads, sv_R_prof_allrads) = allprofiles_allrads
                DoAccrAllRads = True
                print(
                    "WARNING: found DoAccr optimization with variable PA and inc, could be degenerate"
                )
            else:
                (rrs_allrads, v_Phi_prof_allrads,
                 sv_Phi_prof_allrads) = np.loadtxt(file_profile_allrads,
                                                   unpack=True)
                DoAccrAllRads = False

    DoAccr_fixIncPA = False
    DoMerid_fixIncPA = False
    DoAccrAllRads_fixIncPA = False
    DoMeridAllRads_fixIncPA = False
    if os.path.isfile(file_profile_fixincPA):
        if DoFixIncPA:
            print("loading file_profile", file_profile_fixincPA)

            allprofiles_fixincPA = np.loadtxt(file_profile_fixincPA,
                                              unpack=True)
            print("length allprofiles fixincPA", len(allprofiles_fixincPA))
            if (len(allprofiles_fixincPA) > 5):
                (
                    rrs_fixincPA, v_Phi_prof_fixincPA, sv_Phi_prof_fixincPA,
                    v_R_prof_fixincPA, sv_R_prof_fixincPA, v_z_prof_fixincPA,
                    sv_z_prof_fixincPA
                ) = allprofiles_fixincPA  # np.loadtxt(file_profile,unpack=True)
                DoMerid_fixIncPA = True
            elif (len(allprofiles_fixincPA) > 3):
                (rrs_fixincPA, v_Phi_prof_fixincPA, sv_Phi_prof_fixincPA,
                 v_R_prof_fixincPA, sv_R_prof_fixincPA) = allprofiles_fixincPA
                DoAccr_fixIncPA = True
            else:
                (rrs_fixincPA, v_Phi_prof_fixincPA,
                 sv_Phi_prof_fixincPA) = allprofiles_fixincPA
                DoAccr_fixIncPA = False

            if VarOrient:
                print("loading file_profile", file_profile_fixincPA_allrads)

                allprofiles_fixincPA_allrads = np.loadtxt(
                    file_profile_fixincPA_allrads, unpack=True)
                print("length allprofiles fixincPA",
                      len(allprofiles_fixincPA_allrads))
                if (len(allprofiles_fixincPA_allrads) > 5):
                    (rrs_fixincPA_allrads, v_Phi_prof_fixincPA_allrads,
                     sv_Phi_prof_fixincPA_allrads, v_R_prof_fixincPA_allrads,
                     sv_R_prof_fixincPA_allrads, v_z_prof_fixincPA_allrads,
                     sv_z_prof_fixincPA_allrads) = allprofiles_fixincPA_allrads
                    DoMeridAllRads_fixIncPA = True
                elif (len(allprofiles_fixincPA_allrads) > 3):
                    (rrs_fixincPA_allrads, v_Phi_prof_fixincPA_allrads,
                     sv_Phi_prof_fixincPA_allrads, v_R_prof_fixincPA_allrads,
                     sv_R_prof_fixincPA_allrads) = allprofiles_fixincPA_allrads
                    DoAccrAllRads_fixIncPA = True
                else:
                    (rrs_fixincPA_allrads, v_Phi_prof_fixincPA_allrads,
                     sv_Phi_prof_fixincPA_allrads
                     ) = allprofiles_fixincPA_allrads

    nplotsy = 0
    if Plot_vRot_Global:
        nplotsy += 1

    if PlotVarOrient:
        nplotsy += 1

    if Plot_vRot_VarOrient:
        nplotsy += 1

    if XCheckFixIncPA:
        nplotsy += 1

    if DoAccr:
        nplotsy += 1

    if DoAccr_fixIncPA and Plot_vRot_Global:
        nplotsy += 1

    if DoMerid and Plot_vRot_Global:
        nplotsy += 2
    elif DoAccr and Plot_vRot_Global:
        nplotsy += 1

    if DoMerid_fixIncPA and Plot_vRot_Global:
        nplotsy += 2
    elif DoAccr_fixIncPA and Plot_vRot_Global:
        nplotsy += 1

    if DoMeridAllRads and Plot_vRot_VarOrient:
        nplotsy += 2
    elif DoAccrAllRads and Plot_vRot_VarOrient:
        nplotsy += 1

    if DoMeridAllRads_fixIncPA and Plot_vRot_VarOrient_FixIncPA:
        nplotsy += 2
    elif DoAccrAllRads_fixIncPA and Plot_vRot_VarOrient_FixIncPA:
        nplotsy += 1

    if Plot_vRot_VarOrient_FixIncPA:
        nplotsy += 1

    print("Plot_vRot_VarOrient", Plot_vRot_VarOrient)
    print("nplotsy = ", nplotsy)
    figysize = nplotsy * 4 * 7 / 9.
    figxsize = 7.
    figysize = nplotsy * 2.7
    figxsize = 6.
    #figysize=nplotsy*2.25
    figxsize = 5.

    fig = plt.figure(constrained_layout=True, figsize=(figxsize, figysize))
    gs = fig.add_gridspec(
        nplotsy, 1)  #, width_ratios=[2., 1., 1.], height_ratios=[1.,1.])

    if PlotVarOrient:
        axprofile = fig.add_subplot(gs[(nplotsy - 1), 0])

    RunMCMC = False
    GlobalOrientProf = False

    print("testing workdir exists?:", workdir)

    if os.path.isdir(workdir):

        GlobalOrientProf = True

        (RunMCMC, proflist) = LogResults.load(workdir, FixPAInc=False)

        if RunMCMC:
            print("picked up MCMC run")
            [
                r1s, r2s, rregions, incs, psis, PAs, allradsPA, allradsinc,
                allradspsi, PAs_MCMC, tanpsis_MCMC, incs_MCMC, allradsPAMCMC,
                allradsincMCMC, allradspsiMCMC, vsys
            ] = proflist

            allradsinc = allradsincMCMC
            allradsPA = allradsPAMCMC
            allradspsi = allradspsiMCMC

        else:
            [
                r1s, r2s, rregions, incs, psis, PAs, allradsPA, allradsinc,
                allradspsi, vsys
            ] = proflist

        if ForceGlobalOrient:
            print(
                ">>>>>> Plotting  Fix Orient, with forced global orientation, at PA=",
                Force_allradsPA, " inc=", Force_allradsinc)
            allradsPA = Force_allradsPA
            allradsinc = Force_allradsinc * 180. / np.pi

        if (a_min < 0):
            a_min = np.min(r1s)
        if (a_max < 0):
            a_max = np.max(r2s)

        print("a_min ", a_min, " a_max", a_max)

        #print( "psis",psis)
        #print( "allradspsi",allradspsi)
        #print( "incs",incs)

        psierrors = np.zeros((2, len(psis)))
        if RunMCMC:
            for ir, atanpsi in enumerate(tanpsis_MCMC):
                psis[ir] = 180. * np.arctan(atanpsi[0]) / np.pi

                psierrors[1, ir] = 180. * (np.arctan(atanpsi[1] + atanpsi[0]) -
                                           np.arctan(atanpsi[0])) / np.pi
                psierrors[0, ir] = 180. * (np.arctan(
                    atanpsi[0]) - np.arctan(atanpsi[0] - atanpsi[2])) / np.pi

        incerrors = np.zeros((2, len(incs)))
        if RunMCMC:
            for ir, ainc in enumerate(incs_MCMC):
                #print( "ainc",ainc)
                incs[ir] = ainc[0] * 180. / np.pi
                incerrors[1, ir] = 1. * ainc[1] * 180. / np.pi
                incerrors[0, ir] = 1. * ainc[2] * 180. / np.pi

        PAerrors = np.zeros((2, len(PAs)))
        if RunMCMC:
            for ir, aPA in enumerate(PAs_MCMC):
                PAs[ir] = aPA[0]
                PAerrors[1, ir] = 1. * aPA[1]
                PAerrors[0, ir] = 1. * aPA[2]

        if PlotVarPAinc:
            print("rregions", rregions)
            print("PAs", PAs)

            (ymin, ymax) = Orient.PlotOrientProfile(axprofile, rregions, PAs,
                                                    allradsPA, PAerrors, incs,
                                                    allradsinc, incerrors,
                                                    psis, allradspsi,
                                                    psierrors)

            DoAUBar = False
            if DoAUBar:
                barlength = 10. / distance
                xxs = [
                    a_max - (a_max - a_min) * 0.05,
                    a_max - (a_max - a_min) * 0.05 - barlength
                ]
                yys = [
                    ymin + (ymax - ymin) * 0.08, ymin + (ymax - ymin) * 0.08
                ]
                axprofile.text(xxs[0], yys[0] + (ymax - ymin) * 0.03, '10 au')
                axprofile.plot(xxs, yys, color='C5')

            save_prof = np.zeros((len(rregions), 10))
            save_prof[:, 0] = rregions
            save_prof[:, 1] = PAs
            save_prof[:, 2] = PAerrors[0, :]
            save_prof[:, 3] = PAerrors[1, :]
            save_prof[:, 4] = incs
            save_prof[:, 5] = incerrors[0, :]
            save_prof[:, 6] = incerrors[1, :]
            save_prof[:, 7] = psis
            save_prof[:, 8] = psierrors[0, :]
            save_prof[:, 9] = psierrors[1, :]
            fileout_orientprofile = workdir + 'orient_profile.dat'
            np.savetxt(fileout_orientprofile,
                       save_prof)  # x,y,z equal sized 1D arrays

        sigma_PA = np.std(PAs)
        sigma_inc = np.std(incs)
        sigma_psi = np.std(psis)
        print("Psi= %.1f +- %.1f deg" % (allradspsi, sigma_psi))
        print("Inc= %.2f +- %.2f deg" % (allradsinc, sigma_inc))
        print("PA= %.1f +- %.1f deg" % (allradsPA, sigma_PA))

        print("USING inclination ", allradsinc, " for mass estimates")
        cosi = np.fabs(np.cos(allradsinc * np.pi / 180.))

        BackSide = False
        if (allradspsi < 0.):
            BackSide = True

    #######################################################################
    ######################################################################

    if DoFixIncPA:

        (RunMCMC, proflist) = LogResults.load(workdir_fixincPA,
                                              RunMCMC=RunMCMC,
                                              FixPAInc=True)

        print("proflist", proflist)
        if RunMCMC:
            [
                r1s, r2s, rregions_fixincPA, psis_fixincPA,
                allradspsi_fixincPA, tanpsis_fixincPA_MCMC, initPA, initinc,
                vsys
            ] = proflist
        else:
            [
                r1s, r2s, rregions_fixincPA, psis_fixincPA,
                allradspsi_fixincPA, initPA, initinc, vsys
            ] = proflist

        if (len(rregions_fixincPA) == 0):
            Regions = False
            VarOrient = False
            if ((a_min < 0) or (a_max < 0)):
                sys.exit("need to specify a_min and a_max in call to execfig")

        if (a_min < 0):
            a_min = np.min(r1s)
        if (a_max < 0):
            a_max = np.max(r2s)

        print("a_min ", a_min, " a_max", a_max)

        #a_min=np.min(r1s)
        #a_max=np.max(r2s)
        print("a_min ", a_min, " a_max", a_max)

        #print( "psis_fixincPA", psis_fixincPA)
        #print( "allradspsi_fixincPA",allradspsi_fixincPA)

        psierrors_fixincPA = np.zeros((2, len(psis_fixincPA)))
        if RunMCMC:
            #print( "psis_fixincPA MCMC", tanpsis_fixincPA_MCMC)
            for ir, atanpsi in enumerate(tanpsis_fixincPA_MCMC):
                psis_fixincPA[ir] = 180. * np.arctan(atanpsi[0]) / np.pi
                psierrors_fixincPA[
                    1, ir] = 180. * (np.arctan(atanpsi[1] + atanpsi[0]) -
                                     np.arctan(atanpsi[0])) / np.pi
                psierrors_fixincPA[0, ir] = 180. * (np.arctan(
                    atanpsi[0]) - np.arctan(atanpsi[0] - atanpsi[2])) / np.pi

        if PlotVarOrient:
            (ymin_fixincPA, ymax_fixincPA) = Orient.PlotOrientProfile_fixincPA(
                axprofile,
                rregions_fixincPA,
                psis_fixincPA,
                allradspsi_fixincPA,
                psierrors_fixincPA,
                allradsinc=initinc,
                allradsPA=initPA)
            if DoAUBar:
                barlength = 10. / distance
                xxs = [
                    a_max - (a_max - a_min) * 0.05,
                    a_max - (a_max - a_min) * 0.05 - barlength
                ]
                yys = [
                    ymin_fixincPA + (ymax_fixincPA - ymin_fixincPA) * 0.08,
                    ymin_fixincPA + (ymax_fixincPA - ymin_fixincPA) * 0.08
                ]
                axprofile.text(xxs[1],
                               yys[0] + (ymax_fixincPA - ymin_fixincPA) * 0.05,
                               '10 au')

                axprofile.plot(xxs, yys, color='C5')

            save_prof = np.zeros((len(rregions_fixincPA), 4))
            save_prof[:, 0] = rregions_fixincPA
            save_prof[:, 1] = psis_fixincPA
            save_prof[:, 2] = psierrors_fixincPA[0, :]
            save_prof[:, 3] = psierrors_fixincPA[1, :]

            fileout_orientprofile = workdir_fixincPA + 'psi_profile.dat'
            np.savetxt(fileout_orientprofile,
                       save_prof)  # x,y,z equal sized 1D arrays

        sigma_psi_fixincPA = np.std(psis_fixincPA)
        print("Psi fixincPA= %.1f +- %.1f deg" %
              (allradspsi_fixincPA, sigma_psi_fixincPA))

    if PlotVarOrient:

        if not PlotVarPAinc:
            #ymin=min(ymin,ymin_fixincPA)
            #ymax=max(ymax,ymax_fixincPA)
            #else:
            ymin = ymin_fixincPA
            ymax = ymax_fixincPA

        elif DoFixIncPA:

            ymin = min(ymin, ymin_fixincPA)
            ymax = max(ymax, ymax_fixincPA)

        axprofile.set_xlim(a_min, a_max)
        axprofile.set_ylim(ymin, ymax)

        axprofile.set_ylabel(r'deg')
        #axprofile.legend(loc='upper left')
        axprofile.legend()

        axprofile.set_xlabel(r'$R$ / arcsec')

        axprofile.tick_params(axis='both', length=8, direction='in', pad=10)
        axprofile.tick_params(top='on', right='on', direction='in')
        axprofile.tick_params(which='minor',
                              top='on',
                              length=4,
                              right='on',
                              direction='in')
        axprofile.xaxis.set_major_formatter(
            matplotlib.ticker.ScalarFormatter())
        axprofile.xaxis.set_minor_formatter(
            matplotlib.ticker.ScalarFormatter())

    if ForceGlobalOrient:
        print(
            ">>>>>> Plotting  Fix Orient, with forced global orientation, at PA=",
            Force_allradsPA, " inc=", Force_allradsinc)
        allradsPA = Force_allradsPA
        allradsinc = Force_allradsinc * 180. / np.pi

    print("USING inclination ", allradsinc, " for mass estimates")
    #cosi=np.cos(40.06*np.pi/180.)
    cosi = np.fabs(np.cos(allradsinc * np.pi / 180.))

    BackSide = False
    if DoFixIncPA:
        if (allradspsi_fixincPA < 0.):
            BackSide = True
    else:
        if (allradspsi < 0.):
            BackSide = True

    ######################################################################
    ## Correct v_phi for height over midplane

    if VarOrient:
        if GlobalOrientProf:
            hhs = np.interp(rrs, rregions, np.tan(psis * np.pi / 180.))
        if DoFixIncPA:
            hhs_fixincPA = np.interp(rrs_fixincPA, rregions_fixincPA,
                                     np.tan(psis_fixincPA * np.pi / 180.))
    else:
        if DoFixIncPA:
            hhs_fixincPA = np.tan(allradspsi_fixincPA * np.pi /
                                  180.) * np.ones(rrs_fixincPA.shape)
        else:
            hhs = np.tan(allradspsi * np.pi / 180.) * np.ones(rrs.shape)

    if GlobalOrientProf:
        correct4midplane = ((1. + hhs**2)**(3. / 4.))
        v_Phi_prof_mid = v_Phi_prof * correct4midplane

    if DoFixIncPA:
        correct4midplane_fixincPA = ((1. + hhs_fixincPA**2)**(3. / 4.))
        v_Phi_prof_mid_fixincPA = v_Phi_prof_fixincPA * correct4midplane_fixincPA

    jpos = 0
    if Plot_vRot_Global:
        axprofile = fig.add_subplot(gs[jpos, 0])
        #RotCurve.PlotV_phi(axprofile,rrs_fixincPA,a_min,a_max,v_Phi_prof,sv_Phi_prof,v_Phi_prof_mid,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=rgaps,label=r'global')
        if DoFixIncPA:
            (ymin, ymax, voff, scale_radprofile) = RotCurve.PlotV_phi(
                axprofile,
                rrs_fixincPA,
                a_min,
                a_max,
                v_Phi_prof_fixincPA,
                sv_Phi_prof_fixincPA,
                v_Phi_prof_mid_fixincPA,
                distance,
                cosi,
                bmaj,
                DoStellarMass=True,
                ContinuumGaps=rgaps,
                label=r'',
                RadialScaling=RadialScaling,
                title=title,
                MidPlaneExtrapol=MidPlaneExtrapol)
        else:
            (ymin, ymax, voff, scale_radprofile) = RotCurve.PlotV_phi(
                axprofile,
                rrs,
                a_min,
                a_max,
                v_Phi_prof,
                sv_Phi_prof,
                v_Phi_prof_mid,
                distance,
                cosi,
                bmaj,
                DoStellarMass=True,
                ContinuumGaps=rgaps,
                label=r'',
                RadialScaling=RadialScaling,
                title=title,
                MidPlaneExtrapol=MidPlaneExtrapol)

        jpos += 1

    if DoFixIncPA and XCheckFixIncPA:
        axprofile = fig.add_subplot(gs[jpos, 0])
        (ymin, ymax, voff, scale_radprofile) = RotCurve.PlotV_phi(
            axprofile,
            rrs_fixincPA,
            a_min,
            a_max,
            v_Phi_prof_fixincPA,
            sv_Phi_prof_fixincPA,
            v_Phi_prof_mid_fixincPA,
            distance,
            cosi,
            bmaj,
            DoStellarMass=True,
            ContinuumGaps=rgaps,
            label=r'fix $i$, PA global',
            RadialScaling=RadialScaling,
            title=title,
            MidPlaneExtrapol=MidPlaneExtrapol)
        jpos += 1

    if Plot_vRot_VarOrient:
        axprofile = fig.add_subplot(gs[jpos, 0])
        v_Phi_prof_mid_allrads = v_Phi_prof_allrads * correct4midplane
        (ymin, ymax, voff, scale_radprofile) = RotCurve.PlotV_phi(
            axprofile,
            rrs_allrads,
            a_min,
            a_max,
            v_Phi_prof_allrads,
            sv_Phi_prof_allrads,
            v_Phi_prof_mid_allrads,
            distance,
            cosi,
            bmaj,
            DoStellarMass=True,
            ContinuumGaps=rgaps,
            label=r'region av.',
            RadialScaling=RadialScaling,
            title=title,
            MidPlaneExtrapol=MidPlaneExtrapol)
        jpos += 1

    if DoFixIncPA and Plot_vRot_VarOrient_FixIncPA:
        print("plotting v_phi  with variable psi but fixed PA & inc")
        axprofile = fig.add_subplot(gs[jpos, 0])
        if (alabel == 'default'):
            alabel = r'fix $i$, PA region av.'
            alabel = r'$i$=%.1f PA=%d $\psi(R)$' % (allradsinc, allradsPA)

        v_Phi_prof_mid_fixincPA_allrads = v_Phi_prof_fixincPA_allrads * correct4midplane_fixincPA
        (vymin, vymax, voff, scale_radprofile) = RotCurve.PlotV_phi(
            axprofile,
            rrs_fixincPA,
            a_min,
            a_max,
            v_Phi_prof_fixincPA_allrads,
            sv_Phi_prof_fixincPA_allrads,
            v_Phi_prof_mid_fixincPA_allrads,
            distance,
            cosi,
            bmaj,
            DoStellarMass=True,
            ContinuumGaps=rgaps,
            label=alabel,
            RadialScaling=RadialScaling,
            title=title,
            MidPlaneExtrapol=MidPlaneExtrapol)

        #print("v_Phi_prof_fixincPA_allrads*np.sqrt(rrs_allrads)/np.sqrt(a_max)",v_Phi_prof_fixincPA_allrads*np.sqrt(rrs_allrads)/np.sqrt(a_max))

        if WithComparData:
            # Rich Teague data

            #print("v_phi_ComparData*np.sqrt(rrs_ComparData)/np.sqrt(a_max)",v_phi_ComparData*np.sqrt(rrs_ComparData)/np.sqrt(a_max))

            VKepNorm = True
            if (VKepNorm & RadialScaling):

                v_phi_ComparData_resamp = np.interp(rrs_fixincPA,
                                                    rrs_ComparData,
                                                    v_phi_ComparData)
                s_v_phi_ComparData_up_resamp = np.interp(
                    rrs_fixincPA, rrs_ComparData, s_v_phi_ComparData_up)
                s_v_phi_ComparData_lo_resamp = np.interp(
                    rrs_fixincPA, rrs_ComparData, s_v_phi_ComparData_lo)

                maskrange = ((rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))

                axprofile.plot(rrs_fixincPA[maskrange],
                               (v_phi_ComparData_resamp[maskrange] - voff) /
                               scale_radprofile,
                               color='C1',
                               linewidth=1.5,
                               linestyle='solid',
                               alpha=0.5,
                               label=r'eddy')
                axprofile.fill_between(
                    rrs_fixincPA[maskrange],
                    (v_phi_ComparData_resamp[maskrange] +
                     s_v_phi_ComparData_up_resamp[maskrange] - voff) /
                    scale_radprofile,
                    (v_phi_ComparData_resamp[maskrange] -
                     s_v_phi_ComparData_lo_resamp[maskrange] - voff) /
                    scale_radprofile,
                    lw=0.1,
                    color='C1',
                    alpha=0.2,
                    interpolate=True)  #, step='mid'

                vcomparmin = np.min(
                    (v_phi_ComparData_resamp[maskrange] -
                     s_v_phi_ComparData_lo_resamp[maskrange] - voff) /
                    scale_radprofile)
                vcomparmax = np.max(
                    (v_phi_ComparData_resamp[maskrange] +
                     s_v_phi_ComparData_up_resamp[maskrange] - voff) /
                    scale_radprofile)
                vymin = min(vymin, vcomparmin)
                vymax = max(vymax, vcomparmax)
                axprofile.set_ylim(vymin, vymax)

            elif RadialScaling:
                axprofile.plot(rrs_ComparData,
                               v_phi_ComparData * np.sqrt(rrs_ComparData) /
                               np.sqrt(a_max),
                               color='C1',
                               linewidth=1.5,
                               linestyle='solid',
                               alpha=0.5,
                               label=r'$v_\phi$ eddy')
                axprofile.fill_between(
                    rrs_ComparData,
                    (v_phi_ComparData + s_v_phi_ComparData_up) *
                    np.sqrt(rrs_ComparData) / np.sqrt(a_max),
                    (v_phi_ComparData - s_v_phi_ComparData_lo) *
                    np.sqrt(rrs_ComparData) / np.sqrt(a_max),
                    lw=0.1,
                    color='C1',
                    alpha=0.2,
                    interpolate=True)  #, step='mid'

                maskrange = ((rrs_ComparData > a_min) &
                             (rrs_ComparData < a_max))
                vcomparmin = np.min(
                    (v_phi_ComparData[maskrange] -
                     s_v_phi_ComparData_lo[maskrange]) *
                    np.sqrt(rrs_ComparData[maskrange]) / np.sqrt(a_max))
                vcomparmax = np.max(
                    (v_phi_ComparData[maskrange] +
                     s_v_phi_ComparData_up[maskrange]) *
                    np.sqrt(rrs_ComparData[maskrange]) / np.sqrt(a_max))
                vymin = min(vymin, vcomparmin)
                vymax = max(vymax, vcomparmax)
                axprofile.set_ylim(vymin, vymax)

        #axprofile.legend(loc='upper right')
        axprofile.legend()

        jpos += 1

    if (DoMerid and Plot_vRot_Global):
        axprofile = fig.add_subplot(gs[jpos, 0])
        if DoFixIncPA:
            RotCurve.PlotV_z(axprofile,
                             rrs_fixincPA,
                             a_min,
                             a_max,
                             v_z_prof_fixincPA,
                             sv_z_prof_fixincPA,
                             BackSide=BackSide,
                             ContinuumGaps=rgaps,
                             label=r'global',
                             RadialScaling=False,
                             VisibleXaxis=VisibleXaxis_V_z)
        else:
            RotCurve.PlotV_z(axprofile,
                             rrs,
                             a_min,
                             a_max,
                             v_z_prof,
                             sv_z_prof,
                             BackSide=BackSide,
                             ContinuumGaps=rgaps,
                             label=r'global',
                             RadialScaling=False,
                             VisibleXaxis=VisibleXaxis_V_z)

        jpos += 1
        axprofile = fig.add_subplot(gs[jpos, 0])
        RadialScalingVR = scale_radprofile
        if DoFixIncPA:
            RotCurve.PlotV_R(axprofile,
                             rrs,
                             a_min,
                             a_max,
                             v_R_prof,
                             sv_R_prof,
                             ContinuumGaps=rgaps,
                             label=r'global',
                             VisibleXaxis=VisibleXaxis_V_R,
                             RadialScaling=RadialScalingVR)
        else:
            RotCurve.PlotV_R(axprofile,
                             rrs,
                             a_min,
                             a_max,
                             v_R_prof,
                             sv_R_prof,
                             ContinuumGaps=rgaps,
                             label=r'global',
                             VisibleXaxis=VisibleXaxis_V_R,
                             RadialScaling=RadialScalingVR)

        jpos += 1
    elif (DoAccr and Plot_vRot_Global):
        axprofile = fig.add_subplot(gs[jpos, 0])
        RadialScalingVR = scale_radprofile
        RotCurve.PlotV_R(axprofile,
                         rrs_fixincPA,
                         a_min,
                         a_max,
                         v_R_prof_fixincPA,
                         sv_R_prof_fixincPA,
                         ContinuumGaps=rgaps,
                         label=r'global',
                         VisibleXaxis=VisibleXaxis_V_R,
                         RadialScaling=RadialScalingVR)
        jpos += 1

    if (DoMerid_fixIncPA and Plot_vRot_Global):
        axprofile = fig.add_subplot(gs[jpos, 0])
        RotCurve.PlotV_z(axprofile,
                         rrs_fixincPA,
                         a_min,
                         a_max,
                         v_z_prof_fixincPA,
                         sv_z_prof_fixincPA,
                         BackSide=BackSide,
                         ContinuumGaps=rgaps,
                         label=r'',
                         RadialScaling=False,
                         VisibleXaxis=VisibleXaxis_V_z)
        jpos += 1
        axprofile = fig.add_subplot(gs[jpos, 0])
        RadialScalingVR = scale_radprofile
        RotCurve.PlotV_R(axprofile,
                         rrs_fixincPA,
                         a_min,
                         a_max,
                         v_R_prof_fixincPA,
                         sv_R_prof_fixincPA,
                         ContinuumGaps=rgaps,
                         label=r'',
                         VisibleXaxis=VisibleXaxis_V_R,
                         RadialScaling=RadialScalingVR)
        jpos += 1
    elif (DoAccr_fixIncPA and Plot_vRot_Global):
        axprofile = fig.add_subplot(gs[jpos, 0])
        RadialScalingVR = scale_radprofile
        RotCurve.PlotV_R(axprofile,
                         rrs_fixincPA,
                         a_min,
                         a_max,
                         v_R_prof_fixincPA,
                         sv_R_prof_fixincPA,
                         ContinuumGaps=rgaps,
                         label=r'',
                         VisibleXaxis=VisibleXaxis_V_R,
                         RadialScaling=RadialScalingVR)
        jpos += 1

    if DoMeridAllRads and Plot_vRot_VarOrient:
        axprofile = fig.add_subplot(gs[jpos, 0])
        RotCurve.PlotV_z(axprofile,
                         rrs_allrads,
                         a_min,
                         a_max,
                         v_z_prof_allrads,
                         sv_z_prof_allrads,
                         BackSide=BackSide,
                         ContinuumGaps=rgaps,
                         label=r'region av.',
                         RadialScaling=False,
                         VisibleXaxis=VisibleXaxis_V_z)
        jpos += 1
        axprofile = fig.add_subplot(gs[jpos, 0])
        RadialScalingVR = scale_radprofile
        RotCurve.PlotV_R(axprofile,
                         rrs_allrads,
                         a_min,
                         a_max,
                         v_R_prof_allrads,
                         sv_R_prof_allrads,
                         ContinuumGaps=rgaps,
                         label=r'region av.',
                         VisibleXaxis=VisibleXaxis_V_R,
                         RadialScaling=RadialScalingVR)
        jpos += 1
    elif DoAccrAllRads and Plot_vRot_VarOrient:
        axprofile = fig.add_subplot(gs[jpos, 0])
        RadialScalingVR = scale_radprofile
        RotCurve.PlotV_R(axprofile,
                         rrs_allrads,
                         a_min,
                         a_max,
                         v_R_prof_allrads,
                         sv_R_prof_allrads,
                         ContinuumGaps=rgaps,
                         label=r'region av.',
                         VisibleXaxis=VisibleXaxis_V_R,
                         RadialScaling=RadialScalingVR)
        jpos += 1

    if DoMeridAllRads_fixIncPA and Plot_vRot_VarOrient_FixIncPA:

        print("plotting v_R v_Z  with variable psi but fixed PA & inc")
        axprofile = fig.add_subplot(gs[jpos, 0])
        if (alabel == 'default'):
            alabel = r'fix $i$, PA region av.'
            alabel = r'$i$=%.1f PA=%d' % (allradsinc, allradsPA)
            alabel = r'$i$=%.1f PA=%d $\psi(R)$' % (allradsinc, allradsPA)

        #alabel=''
        (vymin, vymax) = RotCurve.PlotV_z(axprofile,
                                          rrs_fixincPA_allrads,
                                          a_min,
                                          a_max,
                                          v_z_prof_fixincPA_allrads,
                                          sv_z_prof_fixincPA_allrads,
                                          BackSide=BackSide,
                                          ContinuumGaps=rgaps,
                                          label=alabel,
                                          RadialScaling=False,
                                          VisibleXaxis=VisibleXaxis_V_z)
        if WithComparData:
            # Rich Teague data
            axprofile.plot(rrs_ComparData,
                           v_z_ComparData * np.sqrt(rrs_ComparData) /
                           np.sqrt(a_max),
                           color='C1',
                           linewidth=1.5,
                           alpha=0.5,
                           linestyle='solid',
                           label=r'eddy')
            axprofile.fill_between(rrs_ComparData,
                                   (v_z_ComparData + s_v_z_ComparData_up) *
                                   np.sqrt(rrs_ComparData) / np.sqrt(a_max),
                                   (v_z_ComparData - s_v_z_ComparData_lo) *
                                   np.sqrt(rrs_ComparData) / np.sqrt(a_max),
                                   lw=0.1,
                                   color='C1',
                                   alpha=0.2,
                                   interpolate=True)  #, step='mid'

            maskrange = ((rrs_ComparData > a_min) & (rrs_ComparData < a_max))
            vcomparmin = np.min(
                (v_z_ComparData[maskrange] - s_v_z_ComparData_lo[maskrange]) *
                np.sqrt(rrs_ComparData[maskrange]) / np.sqrt(a_max))
            vcomparmax = np.max(
                (v_z_ComparData[maskrange] + s_v_z_ComparData_up[maskrange]) *
                np.sqrt(rrs_ComparData[maskrange]) / np.sqrt(a_max))
            vymin = min(vymin, vcomparmin)
            vymax = max(vymax, vcomparmax)

            print(">>>>> compar v_z ::", vymin, vymax)

            axprofile.set_ylim(vymin, vymax)

        if WithComparRadTWind:
            # Comparing with wind model

            zzs = hhs * rrs
            v_z_RadTWind = np.zeros(rrs.shape)
            v_R_RadTWind = np.zeros(rrs.shape)
            distance = 100.
            Rhopp = rrs * distance
            Zpp = zzs * distance
            h_c_outer = 0.06
            Rcavdust = 20.
            v_term = 1E4  # cm/s
            flaring_outer = 0.15
            h_outer = h_c_outer * (rrs / Rcavdust)**(flaring_outer)
            H_outer = h_outer * Rhopp
            for iRhopp in list(range(len(Rhopp))):
                aRhopp = Rhopp[iRhopp]
                aZpp = Zpp[iRhopp]
                v_z_RadTWind[iRhopp] = 0.
                v_R_RadTWind[iRhopp] = 0.
                aH_outer = H_outer[iRhopp]
                print("aRhopp ", aRhopp, "aZpp", aZpp, "aH_outer", aH_outer)
                if ((aRhopp > Rcavdust) and (np.fabs(aZpp) > aH_outer)):
                    vwind = v_term * (1. - (Rcavdust / aRhopp)**2) * (
                        1. - (aH_outer / aZpp)**2)
                    avrpp = vwind * np.cos(
                        np.pi / 4.)  # *np.cos(np.pi/3.) #cylindrical
                    if (aZpp > 0.):
                        avzpp = vwind * np.sin(np.pi / 4.)
                    else:
                        avzpp = -vwind * np.sin(np.pi / 4.)

                    v_z_RadTWind[iRhopp] = avzpp / 1.0E5  # km/s
                    v_R_RadTWind[iRhopp] = avrpp / 1.0E5  # km/s

            if BackSide:
                v_z_RadTWind *= -1

            axprofile.plot(rrs,
                           v_z_RadTWind * np.sqrt(rrs) / np.sqrt(a_max),
                           color='C3',
                           linewidth=1.5,
                           linestyle='solid',
                           label=r'$v_z$ RT')

        #axprofile.legend(loc='lower right')
        axprofile.legend()

        jpos += 1

        axprofile = fig.add_subplot(gs[jpos, 0])
        if (alabel == 'default'):
            alabel = r'fix $i$, PA region av.'
            alabel = r'$i$=%.1f PA=%d' % (allradsinc, allradsPA)
            alabel = r'$i$=%.1f PA=%d $\psi(R)$' % (allradsinc, allradsPA)

        # alabel=''
        RadialScalingVR = scale_radprofile
        (vymin, vymax) = RotCurve.PlotV_R(axprofile,
                                          rrs_fixincPA_allrads,
                                          a_min,
                                          a_max,
                                          v_R_prof_fixincPA_allrads,
                                          sv_R_prof_fixincPA_allrads,
                                          ContinuumGaps=rgaps,
                                          label=alabel,
                                          RadialScaling=RadialScalingVR,
                                          VisibleXaxis=VisibleXaxis_V_R)
        if WithComparData:
            # Rich Teague data

            VKepNorm = True
            if (VKepNorm & RadialScaling):

                v_rad_ComparData_resamp = np.interp(rrs_fixincPA,
                                                    rrs_ComparData,
                                                    v_rad_ComparData)
                s_v_rad_ComparData_up_resamp = np.interp(
                    rrs_fixincPA, rrs_ComparData, s_v_rad_ComparData_up)
                s_v_rad_ComparData_lo_resamp = np.interp(
                    rrs_fixincPA, rrs_ComparData, s_v_rad_ComparData_lo)

                maskrange = ((rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))

                axprofile.plot(rrs_fixincPA[maskrange],
                               v_rad_ComparData_resamp[maskrange] /
                               scale_radprofile,
                               color='C1',
                               linewidth=1.5,
                               linestyle='solid',
                               alpha=0.5,
                               label=r'eddy')
                axprofile.fill_between(
                    rrs_fixincPA[maskrange],
                    (v_rad_ComparData_resamp[maskrange] +
                     s_v_rad_ComparData_up_resamp[maskrange]) /
                    scale_radprofile,
                    (v_rad_ComparData_resamp[maskrange] -
                     s_v_rad_ComparData_lo_resamp[maskrange]) /
                    scale_radprofile,
                    lw=0.1,
                    color='C1',
                    alpha=0.2,
                    interpolate=True)  #, step='mid'

                vcomparmin = np.min((v_rad_ComparData_resamp[maskrange] -
                                     s_v_rad_ComparData_lo_resamp[maskrange]) /
                                    scale_radprofile)
                vcomparmax = np.max((v_rad_ComparData_resamp[maskrange] +
                                     s_v_rad_ComparData_up_resamp[maskrange]) /
                                    scale_radprofile)
                vymin = min(vymin, vcomparmin)
                vymax = max(vymax, vcomparmax)

                axprofile.set_ylim(vymin, vymax)

            else:
                axprofile.plot(rrs_ComparData,
                               v_rad_ComparData * np.sqrt(rrs_ComparData) /
                               np.sqrt(a_max),
                               color='C1',
                               linewidth=1.5,
                               linestyle='solid',
                               label=r'$v_R$ eddy')
                axprofile.fill_between(
                    rrs_ComparData,
                    (v_rad_ComparData + s_v_rad_ComparData_up) *
                    np.sqrt(rrs_ComparData) / np.sqrt(a_max),
                    (v_rad_ComparData - s_v_rad_ComparData_lo) *
                    np.sqrt(rrs_ComparData) / np.sqrt(a_max),
                    lw=0.1,
                    color='C1',
                    alpha=0.2,
                    interpolate=True)  #, step='mid'

                maskrange = ((rrs_ComparData > a_min) &
                             (rrs_ComparData < a_max))
                vcomparmin = np.min(
                    (v_rad_ComparData[maskrange] -
                     s_v_rad_ComparData_lo[maskrange]) *
                    np.sqrt(rrs_ComparData[maskrange]) / np.sqrt(a_max))
                vcomparmax = np.max(
                    (v_rad_ComparData[maskrange] +
                     s_v_rad_ComparData_up[maskrange]) *
                    np.sqrt(rrs_ComparData[maskrange]) / np.sqrt(a_max))
                vymin = min(vymin, vcomparmin)
                vymax = max(vymax, vcomparmax)

                axprofile.set_ylim(vymin, vymax)

        if WithComparRadTWind:
            axprofile.plot(rrs,
                           v_R_RadTWind * np.sqrt(rrs) / np.sqrt(a_max),
                           color='C4',
                           linewidth=1.5,
                           linestyle='solid',
                           label=r'$v_z$ RT')

        #axprofile.legend(loc='lower right')
        axprofile.legend()
        jpos += 1
    elif DoAccrAllRads_fixIncPA and Plot_vRot_VarOrient_FixIncPA:
        axprofile = fig.add_subplot(gs[jpos, 0])
        if (alabel == 'default'):
            alabel = r'fix $i$, PA region av.'
            alabel = r'$i$=%.1f PA=%d' % (allradsinc, allradsPA)
            alabel = r'$i$=%.1f PA=%d $\psi(R)$' % (allradsinc, allradsPA)

        #alabel=''
        RadialScalingVR = scale_radprofile

        RotCurve.PlotV_R(
            axprofile,
            rrs_fixincPA_allrads,
            a_min,
            a_max,
            v_R_prof_fixincPA_allrads,
            sv_R_prof_fixincPA_allrads,
            ContinuumGaps=rgaps,
            label=alabel,
            RadialScaling=RadialScalingVR,
            VisibleXaxis=VisibleXaxis_V_R,
        )
        if WithComparData:
            # Rich Teague data
            axprofile.plot(rrs_ComparData,
                           v_rad_ComparData * np.sqrt(rrs_ComparData) /
                           np.sqrt(a_max),
                           color='C1',
                           linewidth=1.5,
                           linestyle='solid',
                           label=r'$v_R$ eddy')
            axprofile.fill_between(rrs_ComparData,
                                   (v_rad_ComparData + s_v_rad_ComparData_up) *
                                   np.sqrt(rrs_ComparData) / np.sqrt(a_max),
                                   (v_rad_ComparData - s_v_rad_ComparData_lo) *
                                   np.sqrt(rrs_ComparData) / np.sqrt(a_max),
                                   lw=0.1,
                                   color='C1',
                                   alpha=0.2,
                                   interpolate=True)  #, step='mid'

        axprofile.legend()
        jpos += 1

    if (not PlotVarOrient):
        axprofile = fig.add_subplot(gs[-1, 0])

        #plt.setp(axprofile.get_xticklabels(),visible=True) #, fontsize=6)
        axprofile.set_xlim(a_min, a_max)
        print("LAST a_min ", a_min, " a_max", a_max)
        print("VisibileXaxis_V_z", VisibleXaxis_V_z)
        print("axprofile.get_xticklabels()", axprofile.get_xticklabels())
        axprofile.set_xlabel(r'$R$ / arcsec')

        axprofile.tick_params(axis='both', length=8, direction='in', pad=10)
        axprofile.tick_params(top='on', right='on', direction='in')
        axprofile.tick_params(which='minor',
                              top='on',
                              length=4,
                              right='on',
                              direction='in')
        axprofile.xaxis.set_major_formatter(
            matplotlib.ticker.ScalarFormatter())
        axprofile.xaxis.set_minor_formatter(
            matplotlib.ticker.ScalarFormatter())

    plt.subplots_adjust(hspace=0.)

    print(fileout_fig)
    fig.savefig(fileout_fig, bbox_inches='tight')

    return vsys
