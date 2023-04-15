import numpy as np
import matplotlib.pyplot as plt
import sys

import ConeRot.RotOrient.StellarMass as StellarMass


def drawgaps(axprofile, ContinuumGaps, ymin, ymax):
    for argap in ContinuumGaps:
        if (isinstance(argap, list)):
            argap1 = argap[0]
            argap2 = argap[1]

            #axprofile.fill_between(argap, [ymin,ymin], [ymax,ymax], lw=0.1, alpha=0.08, fc='c',hatch='//',zorder=2) #, step='mid'
            axprofile.fill_between(argap, [ymin, ymin], [ymax, ymax],
                                   lw=0.1,
                                   alpha=0.1,
                                   fc='c',
                                   zorder=2)  #, step='mid'
        else:
            axprofile.plot([argap, argap], [ymin, ymax],
                           color='C4',
                           linewidth=0.5,
                           linestyle='dotted')


def PlotV_phi(axprofile,
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
              ContinuumGaps=False,
              label='',
              RadialScaling=True,
              title='',
              MidPlaneExtrapol=True):

    ######################################################################
    # Plot rotation curves

    if (DoStellarMass):
        [vnorm, Mstar] = StellarMass.KepMass(axprofile,
                                             rrs_fixincPA,
                                             cosi,
                                             bmaj,
                                             v_Phi_prof_fixincPA,
                                             sv_Phi_prof_fixincPA,
                                             distance,
                                             a_min,
                                             a_max,
                                             RadialScaling=RadialScaling)

    rmax = np.max(rrs_fixincPA)

    axprofile.set_xlim(a_min, a_max)
    #maskrange=np.where( (rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))
    plotmask = np.where((rrs_fixincPA >= a_min) & (rrs_fixincPA <= a_max))
    maskrange = plotmask

    scale_radprofile = 1.
    if RadialScaling:
        if DoStellarMass:
            scale_radprofile = vnorm[maskrange]
        else:
            scale_radprofile = (np.sqrt(rrs_fixincPA[maskrange]) /
                                np.sqrt(a_max))
        voff = scale_radprofile
    else:
        voff = 0.

    ymin = 0.99 * np.min(
        ((v_Phi_prof_fixincPA[maskrange] - sv_Phi_prof_fixincPA[maskrange]) -
         voff) / scale_radprofile)

    dup = (
        (v_Phi_prof_fixincPA[maskrange] + sv_Phi_prof_fixincPA[maskrange]) -
        voff
    ) / scale_radprofile  # (np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max))

    dup_corr = (
        v_Phi_prof_mid_fixincPA[maskrange] - voff
    ) / scale_radprofile  # (np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max))

    if MidPlaneExtrapol:
        ymax = 1.005 * np.max(np.concatenate((dup, dup_corr)))
    else:
        ymax = 1.005 * np.max(dup)

    axprofile.set_ylim(ymin, ymax)

    axprofile.fill_between(
        rrs_fixincPA[plotmask],
        (((v_Phi_prof_fixincPA[plotmask] + sv_Phi_prof_fixincPA[plotmask]) -
          voff) / scale_radprofile),
        ((v_Phi_prof_fixincPA[plotmask] - sv_Phi_prof_fixincPA[plotmask] -
          voff) / scale_radprofile),
        lw=0.1,
        color='grey',
        alpha=0.2,
        interpolate=True)  #, step='mid'

    VKepNorm = True
    if RadialScaling:
        if VKepNorm:
            axprofile.set_ylabel('')
            linelabel = r'$(\overline{v}_{\phi} - v_K)/v_K$'
        else:
            axprofile.set_ylabel(
                r'$\sqrt{R/' + str(a_max) +
                '} \\times \\overline{v}_{\phi}(R)$ / km s$^{-1}$')
            linelabel = r'$\overline{v}_\phi$'
    else:
        axprofile.set_ylabel(r'$\overline{v}_{\phi}$ / km s$^{-1}$')
        linelabel = r'$\overline{v}_\phi$'

    axprofile.plot(rrs_fixincPA[plotmask],
                   ((v_Phi_prof_fixincPA[plotmask] - voff) / scale_radprofile),
                   color='grey',
                   linewidth=1.5,
                   linestyle='solid',
                   label=linelabel)

    if (ContinuumGaps):
        drawgaps(axprofile, ContinuumGaps, ymin, ymax)

    #print( "ymax v_Phi_prof_mid_fixincPA", np.max(v_Phi_prof_mid_fixincPA[plotmask]*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max)))
    #print( "ymax v_Phi_prof_mid_fixincPA", np.max(v_Phi_prof_mid_fixincPA[plotmask]*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max)))

    if MidPlaneExtrapol:
        linelabel += ' mid.'
        axprofile.plot(rrs_fixincPA[plotmask],
                       (v_Phi_prof_mid_fixincPA[plotmask] - voff) /
                       scale_radprofile,
                       color='cornflowerblue',
                       linewidth=1.5,
                       linestyle='solid',
                       label=linelabel)
        if (DoStellarMass):
            StellarMass.KepMass(axprofile,
                                rrs_fixincPA,
                                cosi,
                                bmaj,
                                v_Phi_prof_mid_fixincPA,
                                sv_Phi_prof_fixincPA,
                                distance,
                                a_min,
                                a_max,
                                linecolor='cornflowerblue',
                                RadialScaling=vnorm)

    if (label != ''):
        axprofile.text(a_min + (a_max - a_min) * 0.05,
                       ymax - (ymax - ymin) * 0.12, label)

    #axprofile.legend(fontsize=16)
    axprofile.legend()  # loc='upper left')

    axprofile.tick_params(axis='both', length=8, direction='in', pad=10)
    axprofile.tick_params(top='on', right='on', direction='in')
    axprofile.tick_params(which='minor',
                          top='on',
                          length=4,
                          right='on',
                          direction='in')

    plt.setp(axprofile.get_xticklabels(), visible=False)  #, fontsize=6)

    if (title != ''):
        plt.title(title)

    return (ymin, ymax, voff, scale_radprofile)


def PlotV_R(axprofile,
            rrs_fixincPA,
            a_min,
            a_max,
            v_R_prof_fixincPA,
            sv_R_prof_fixincPA,
            ContinuumGaps=False,
            label='',
            VisibleXaxis=False,
            RadialScaling=True):

    rmax = np.max(rrs_fixincPA)

    sv_R_prof_fixincPA = np.nan_to_num(sv_R_prof_fixincPA)

    axprofile.set_xlim(a_min, a_max)

    #maskrange=( (rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))
    plotmask = np.where((rrs_fixincPA >= a_min) & (rrs_fixincPA <= a_max))
    maskrange = plotmask

    scale_radprofile = 1.
    voff = 0.
    if isinstance(RadialScaling, np.ndarray):
        scale_radprofile = RadialScaling
        voff = 0.
    elif RadialScaling:
        scale_radprofile = (np.sqrt(rrs_fixincPA[maskrange]) / np.sqrt(a_max))

    #from pprint import pprint

    #pprint( list(zip(rrs_fixincPA, sv_R_prof_fixincPA, maskrange) ))

    ymin = 1.01 * np.min(
        ((v_R_prof_fixincPA[maskrange] - sv_R_prof_fixincPA[maskrange]) - voff)
        / scale_radprofile)
    ymax = 1.01 * np.max(
        ((v_R_prof_fixincPA[maskrange] + sv_R_prof_fixincPA[maskrange]) - voff)
        / scale_radprofile)

    if isinstance(RadialScaling, np.ndarray):
        axprofile.set_ylabel('')
        linelabel = r'$\overline{v}_{R} / v_K$'
    elif RadialScaling:
        axprofile.set_ylabel(r'$\sqrt{R/' + str(a_max) +
                             '} \\times \\overline{v}_{R}$ / km s$^{-1}$')
        linelabel = r'$\overline{v}_R$'
    else:
        axprofile.set_ylabel(r'$\overline{v}_{R}$ / km s$^{-1}$')
        linelabel = r'$\overline{v}_R$'

    axprofile.plot(rrs_fixincPA[plotmask],
                   (v_R_prof_fixincPA[plotmask] - voff) / scale_radprofile,
                   color='C0',
                   linewidth=1.5,
                   linestyle='solid',
                   label=linelabel)

    #print( "wcols v_R")
    #print((np.array(zip(rrs_fixincPA[plotmask],v_R_prof_fixincPA[plotmask]+sv_R_prof_fixincPA[plotmask]))))

    axprofile.fill_between(
        rrs_fixincPA[plotmask],
        ((v_R_prof_fixincPA[plotmask] + sv_R_prof_fixincPA[plotmask]) - voff) /
        scale_radprofile,
        ((v_R_prof_fixincPA[plotmask] - sv_R_prof_fixincPA[plotmask]) - voff) /
        scale_radprofile,
        lw=0.1,
        color='C0',
        alpha=0.2,
        interpolate=True)  #, step='mid'

    axprofile.plot(rrs_fixincPA[plotmask],
                   rrs_fixincPA[plotmask] * 0.,
                   color='grey',
                   linewidth=0.5,
                   linestyle='dotted')

    if (label != ''):
        axprofile.text(a_min + (a_max - a_min) * 0.05,
                       ymax - (ymax - ymin) * 0.12, label)

    axprofile.set_xlim(a_min, a_max)
    axprofile.set_ylim(ymin, ymax)
    #axprofile.legend(fontsize=16)
    axprofile.legend(loc='lower left')

    axprofile.tick_params(axis='both', length=8, direction='in', pad=10)
    axprofile.tick_params(top='on', right='on', direction='in')
    axprofile.tick_params(which='minor',
                          top='on',
                          length=4,
                          right='on',
                          direction='in')

    plt.setp(axprofile.get_xticklabels(), visible=VisibleXaxis)  #, fontsize=6)

    if (ContinuumGaps):
        drawgaps(axprofile, ContinuumGaps, ymin, ymax)

    return (ymin, ymax)


def PlotV_z(axprofile,
            rrs_fixincPA,
            a_min,
            a_max,
            v_z_prof_fixincPA,
            sv_z_prof_fixincPA,
            BackSide=False,
            ContinuumGaps=False,
            label='',
            RadialScaling=False,
            VisibleXaxis=False):

    if BackSide:
        v_z_prof_fixincPA *= -1

    rmax = np.max(rrs_fixincPA)

    sv_z_prof_fixincPA = np.nan_to_num(sv_z_prof_fixincPA)

    axprofile.set_xlim(a_min, a_max)
    #maskrange=( (rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))
    plotmask = np.where((rrs_fixincPA >= a_min) & (rrs_fixincPA <= a_max))
    maskrange = plotmask

    scale_radprofile = 1.
    if RadialScaling:
        scale_radprofile = (np.sqrt(rrs_fixincPA[maskrange]) / np.sqrt(a_max))

    #from pprint import pprint

    #pprint( list(zip(rrs_fixincPA, sv_R_prof_fixincPA, maskrange) ))

    #print(">>>>> median error:",np.median(sv_z_prof_fixincPA[maskrange]))
    #print(">>>>> median values:",np.median(v_z_prof_fixincPA[maskrange]))

    ymin = 1.01 * np.min(
        (v_z_prof_fixincPA[maskrange] - sv_z_prof_fixincPA[maskrange]) *
        scale_radprofile)
    ymax = 1.01 * np.max(
        (v_z_prof_fixincPA[maskrange] + sv_z_prof_fixincPA[maskrange]) *
        scale_radprofile)

    if BackSide:
        prefix = r'$-$'
    else:
        prefix = ''

    if RadialScaling:
        linelabel = r'$\overline{v}_z$'
        axprofile.set_ylabel(prefix + r'$\sqrt{R/' + str(a_max) +
                             '} \\times \\overline{v}_{z}(R)$ / km s$^{-1}$')
    else:
        linelabel = prefix + r'$\overline{v}_{z}$ / km s$^{-1}$'
        axprofile.set_ylabel('')

    axprofile.plot(rrs_fixincPA[plotmask],
                   v_z_prof_fixincPA[plotmask] * scale_radprofile,
                   color='C0',
                   linewidth=1.5,
                   linestyle='solid',
                   label=linelabel)

    #print( "wcols v_z")
    #print((np.array(zip(rrs_fixincPA[plotmask],v_z_prof_fixincPA[plotmask]+sv_z_prof_fixincPA[plotmask]))))

    axprofile.fill_between(
        rrs_fixincPA[plotmask],
        (v_z_prof_fixincPA[plotmask] + sv_z_prof_fixincPA[plotmask]) *
        scale_radprofile,
        (v_z_prof_fixincPA[plotmask] - sv_z_prof_fixincPA[plotmask]) *
        scale_radprofile,
        lw=0.1,
        color='C0',
        alpha=0.2,
        interpolate=True)  #, step='mid'

    axprofile.plot(rrs_fixincPA[plotmask],
                   rrs_fixincPA[plotmask] * 0.,
                   color='grey',
                   linewidth=0.5,
                   linestyle='dotted')

    if (label != ''):
        axprofile.text(a_min + (a_max - a_min) * 0.05,
                       ymax - (ymax - ymin) * 0.12, label)

    axprofile.set_xlim(a_min, a_max)
    axprofile.set_ylim(ymin, ymax)
    #axprofile.legend(fontsize=16)
    axprofile.legend(loc='lower left')

    axprofile.tick_params(axis='both', length=8, direction='in', pad=10)
    axprofile.tick_params(top='on', right='on', direction='in')
    axprofile.tick_params(which='minor',
                          top='on',
                          length=4,
                          right='on',
                          direction='in')

    plt.setp(axprofile.get_xticklabels(), visible=VisibleXaxis)  #, fontsize=6)

    if (ContinuumGaps):
        drawgaps(axprofile, ContinuumGaps, ymin, ymax)

    print(">>>>> v_z ::", ymin, ymax)

    return (ymin, ymax)
