import numpy as np
import matplotlib.pyplot as plt

import ConeRot.RotOrient.StellarMass as StellarMass


def PlotV_phi(axprofile,rrs_fixincPA,a_min,a_max,v_Phi_prof_fixincPA,sv_Phi_prof_fixincPA,v_Phi_prof_mid_fixincPA,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=False,label=''):
         
    ######################################################################
    # Plot rotation curves

 
    rmax=np.max(rrs_fixincPA)

    axprofile.set_xlim(a_min,a_max)
    maskrange=np.where( (rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))
    plotmask = np.where( (rrs_fixincPA >= a_min) & (rrs_fixincPA <= a_max) )


    ymin=0.99*np.min((v_Phi_prof_fixincPA[maskrange]-sv_Phi_prof_fixincPA[maskrange]) * (np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max)))

    dup = (v_Phi_prof_fixincPA[maskrange]+sv_Phi_prof_fixincPA[maskrange]) * (np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max))
    dup_corr = v_Phi_prof_mid_fixincPA[maskrange] * (np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max))

    ymax=1.005*np.max(np.concatenate((dup, dup_corr)))



    axprofile.set_ylim(ymin,ymax)

    axprofile.fill_between(rrs_fixincPA[plotmask], (v_Phi_prof_fixincPA[plotmask]+sv_Phi_prof_fixincPA[plotmask])*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max), (v_Phi_prof_fixincPA[plotmask]-sv_Phi_prof_fixincPA[plotmask])*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max), lw=0.1,color='grey', alpha=0.2, interpolate=True) #, step='mid'

    axprofile.plot(rrs_fixincPA[plotmask],v_Phi_prof_fixincPA[plotmask]*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max),color='grey',linewidth=1.5,linestyle='solid',label='data')


    if (DoStellarMass):
        StellarMass.KepMass(axprofile,rrs_fixincPA,cosi,bmaj,v_Phi_prof_fixincPA,sv_Phi_prof_fixincPA,distance,a_min,a_max)
        StellarMass.KepMass(axprofile,rrs_fixincPA,cosi,bmaj,v_Phi_prof_mid_fixincPA,sv_Phi_prof_fixincPA,distance,a_min,a_max,linecolor='cornflowerblue')

    if (ContinuumGaps):
        for argap in ContinuumGaps:
            axprofile.plot([argap, argap],[ymin,ymax],color='black',linewidth=0.5,linestyle='dotted')


    print( "ymax v_Phi_prof_mid_fixincPA", np.max(v_Phi_prof_mid_fixincPA[plotmask]*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max)))
    print( "ymax v_Phi_prof_mid_fixincPA", np.max(v_Phi_prof_mid_fixincPA[plotmask]*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max)))



    axprofile.plot(rrs_fixincPA[plotmask],v_Phi_prof_mid_fixincPA[plotmask]*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max),color='cornflowerblue',linewidth=1.5,linestyle='solid',label='data mid.')

    if (label != ''):
        axprofile.text(a_min+(a_max-a_min)*0.05,ymax-(ymax-ymin)*0.1,label)
    
    #axprofile.legend(fontsize=16)
    axprofile.legend()


    axprofile.set_ylabel(r'$\sqrt{r/'+str(a_max)+'} \\times \\tilde{v}_{\phi}(r)$ / km s$^{-1}$')


    axprofile.tick_params(axis='both', length = 8,  direction='in', pad=10)
    axprofile.tick_params(top='on',right='on', direction='in')
    axprofile.tick_params(which='minor', top='on', length = 4, right='on', direction='in')

    plt.setp(axprofile.get_xticklabels(),visible=False) #, fontsize=6)

    return


def PlotV_R(axprofile,rrs_fixincPA,a_min,a_max,v_R_prof_fixincPA,sv_R_prof_fixincPA,ContinuumGaps=False,label=''):


    rmax=np.max(rrs_fixincPA)

    sv_R_prof_fixincPA=np.nan_to_num(sv_R_prof_fixincPA)
    
    axprofile.set_xlim(a_min,a_max)
    maskrange=( (rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))
    plotmask = np.where( (rrs_fixincPA >= a_min) & (rrs_fixincPA <= a_max) )

    #from pprint import pprint
    
    #pprint( list(zip(rrs_fixincPA, sv_R_prof_fixincPA, maskrange) ))
    

    
    ymin=1.01*np.min((v_R_prof_fixincPA[maskrange] - sv_R_prof_fixincPA[maskrange]) * (np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max)))
    ymax=1.01*np.max((v_R_prof_fixincPA[maskrange] + sv_R_prof_fixincPA[maskrange]) * (np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max)))



    axprofile.plot(rrs_fixincPA[plotmask],v_R_prof_fixincPA[plotmask]*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max),color='red',linewidth=1.5,linestyle='solid',label=r'data $v_R$')
    axprofile.set_ylabel(r'$\sqrt{r/'+str(a_max)+'} \\times \\tilde{v}_{R}(r)$ / km s$^{-1}$')

    print( "wcols v_R")
    print((np.array(zip(rrs_fixincPA[plotmask],v_R_prof_fixincPA[plotmask]+sv_R_prof_fixincPA[plotmask]))))

    axprofile.fill_between(rrs_fixincPA[plotmask], (v_R_prof_fixincPA[plotmask]+sv_R_prof_fixincPA[plotmask])*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max), (v_R_prof_fixincPA[plotmask]-sv_R_prof_fixincPA[plotmask])*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max), lw=0.1,color='red', alpha=0.2, interpolate=True) #, step='mid'


    if (label != ''):
        axprofile.text(a_min+(a_max-a_min)*0.05,ymax-(ymax-ymin)*0.1,label)

    axprofile.set_xlim(a_min,a_max)
    axprofile.set_ylim(ymin,ymax)
    #axprofile.legend(fontsize=16)
    axprofile.legend(loc='lower left')

    axprofile.tick_params(axis='both', length = 8,  direction='in', pad=10)
    axprofile.tick_params(top='on',right='on', direction='in')
    axprofile.tick_params(which='minor', top='on', length = 4, right='on', direction='in')
    
    plt.setp(axprofile.get_xticklabels(),visible=False) #, fontsize=6)
    
    if (ContinuumGaps):
        for argap in ContinuumGaps:
            axprofile.plot([argap, argap],[ymin,ymax],color='black',linewidth=0.5,linestyle='dotted')


    return
