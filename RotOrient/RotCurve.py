import numpy as np
import matplotlib.pyplot as plt
import sys

import ConeRot.RotOrient.StellarMass as StellarMass


def PlotV_phi(axprofile,rrs_fixincPA,a_min,a_max,v_Phi_prof_fixincPA,sv_Phi_prof_fixincPA,v_Phi_prof_mid_fixincPA,distance,cosi,bmaj, DoStellarMass=True, ContinuumGaps=False,label='',RadialScaling=True,title=''):
         
    ######################################################################
    # Plot rotation curves

 
    rmax=np.max(rrs_fixincPA)

    axprofile.set_xlim(a_min,a_max)
    #maskrange=np.where( (rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))
    plotmask = np.where( (rrs_fixincPA >= a_min) & (rrs_fixincPA <= a_max) )
    maskrange=plotmask
    
    scale_radprofile=1.
    if RadialScaling:
        scale_radprofile=(np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max))
        
    ymin=0.99*np.min((v_Phi_prof_fixincPA[maskrange]-sv_Phi_prof_fixincPA[maskrange]) * scale_radprofile)

    dup = (v_Phi_prof_fixincPA[maskrange]+sv_Phi_prof_fixincPA[maskrange]) * scale_radprofile # (np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max))
    dup_corr = v_Phi_prof_mid_fixincPA[maskrange] * scale_radprofile # (np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max))

    ymax=1.005*np.max(np.concatenate((dup, dup_corr)))
    
    axprofile.set_ylim(ymin,ymax)

    axprofile.fill_between(rrs_fixincPA[plotmask], (v_Phi_prof_fixincPA[plotmask]+sv_Phi_prof_fixincPA[plotmask])*scale_radprofile, (v_Phi_prof_fixincPA[plotmask]-sv_Phi_prof_fixincPA[plotmask])*scale_radprofile, lw=0.1,color='grey', alpha=0.2, interpolate=True) #, step='mid'

    axprofile.plot(rrs_fixincPA[plotmask],v_Phi_prof_fixincPA[plotmask]*scale_radprofile,color='grey',linewidth=1.5,linestyle='solid',label=r'$v_\phi$')


    if (DoStellarMass):
        StellarMass.KepMass(axprofile,rrs_fixincPA,cosi,bmaj,v_Phi_prof_fixincPA,sv_Phi_prof_fixincPA,distance,a_min,a_max,RadialScaling=RadialScaling)

    if (ContinuumGaps):
        for argap in ContinuumGaps:
            axprofile.plot([argap, argap],[ymin,ymax],color='black',linewidth=0.5,linestyle='dotted')


    #print( "ymax v_Phi_prof_mid_fixincPA", np.max(v_Phi_prof_mid_fixincPA[plotmask]*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max)))
    #print( "ymax v_Phi_prof_mid_fixincPA", np.max(v_Phi_prof_mid_fixincPA[plotmask]*np.sqrt(rrs_fixincPA[plotmask])/np.sqrt(a_max)))



    axprofile.plot(rrs_fixincPA[plotmask],v_Phi_prof_mid_fixincPA[plotmask]*scale_radprofile,color='cornflowerblue',linewidth=1.5,linestyle='solid',label=r'$v_\phi$ mid.')
    if (DoStellarMass):
        StellarMass.KepMass(axprofile,rrs_fixincPA,cosi,bmaj,v_Phi_prof_mid_fixincPA,sv_Phi_prof_fixincPA,distance,a_min,a_max,linecolor='cornflowerblue',RadialScaling=RadialScaling)

    if (label != ''):
        axprofile.text(a_min+(a_max-a_min)*0.05,ymax-(ymax-ymin)*0.12,label)
    
    #axprofile.legend(fontsize=16)
    axprofile.legend(loc='upper left')


    if RadialScaling:
        axprofile.set_ylabel(r'$\sqrt{R/'+str(a_max)+'} \\times \\tilde{v}_{\phi}(R)$ / km s$^{-1}$')
    else:
        axprofile.set_ylabel(r'$\tilde{v}_{\phi}(R)$ / km s$^{-1}$')


    axprofile.tick_params(axis='both', length = 8,  direction='in', pad=10)
    axprofile.tick_params(top='on',right='on', direction='in')
    axprofile.tick_params(which='minor', top='on', length = 4, right='on', direction='in')

    plt.setp(axprofile.get_xticklabels(),visible=False) #, fontsize=6)

    if (title != ''):
        plt.title(title)
        
    return (ymin,ymax)


def PlotV_R(axprofile,rrs_fixincPA,a_min,a_max,v_R_prof_fixincPA,sv_R_prof_fixincPA,ContinuumGaps=False,label='',VisibleXaxis=False,RadialScaling=True):


    rmax=np.max(rrs_fixincPA)

    sv_R_prof_fixincPA=np.nan_to_num(sv_R_prof_fixincPA)
                     
    axprofile.set_xlim(a_min,a_max)
    
    #maskrange=( (rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))
    plotmask = np.where( (rrs_fixincPA >= a_min) & (rrs_fixincPA <= a_max) )
    maskrange=plotmask
    
    scale_radprofile=1.
    if RadialScaling:
        scale_radprofile=(np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max))

    #from pprint import pprint
    
    #pprint( list(zip(rrs_fixincPA, sv_R_prof_fixincPA, maskrange) ))
    

    
    ymin=1.01*np.min((v_R_prof_fixincPA[maskrange] - sv_R_prof_fixincPA[maskrange]) * scale_radprofile)
    ymax=1.01*np.max((v_R_prof_fixincPA[maskrange] + sv_R_prof_fixincPA[maskrange]) * scale_radprofile)



    axprofile.plot(rrs_fixincPA[plotmask],v_R_prof_fixincPA[plotmask]*scale_radprofile,color='C0',linewidth=1.5,linestyle='solid',label=r'$v_R$')

    if RadialScaling:
        axprofile.set_ylabel(r'$\sqrt{R/'+str(a_max)+'} \\times \\tilde{v}_{R}(R)$ / km s$^{-1}$')
    else:
        axprofile.set_ylabel(r'$\tilde{v}_{R}(R)$ / km s$^{-1}$')

    #print( "wcols v_R")
    #print((np.array(zip(rrs_fixincPA[plotmask],v_R_prof_fixincPA[plotmask]+sv_R_prof_fixincPA[plotmask]))))

    axprofile.fill_between(rrs_fixincPA[plotmask], (v_R_prof_fixincPA[plotmask]+sv_R_prof_fixincPA[plotmask])*scale_radprofile, (v_R_prof_fixincPA[plotmask]-sv_R_prof_fixincPA[plotmask])*scale_radprofile, lw=0.1,color='C0', alpha=0.2, interpolate=True) #, step='mid'

    axprofile.plot(rrs_fixincPA[plotmask], rrs_fixincPA[plotmask]*0.,color='grey',linewidth=0.5,linestyle='dotted')

    if (label != ''):
        axprofile.text(a_min+(a_max-a_min)*0.05,ymax-(ymax-ymin)*0.12,label)

    axprofile.set_xlim(a_min,a_max)
    axprofile.set_ylim(ymin,ymax)
    #axprofile.legend(fontsize=16)
    axprofile.legend(loc='lower left')

    axprofile.tick_params(axis='both', length = 8,  direction='in', pad=10)
    axprofile.tick_params(top='on',right='on', direction='in')
    axprofile.tick_params(which='minor', top='on', length = 4, right='on', direction='in')

    plt.setp(axprofile.get_xticklabels(),visible=VisibleXaxis) #, fontsize=6)
    
    if (ContinuumGaps):
        for argap in ContinuumGaps:
            axprofile.plot([argap, argap],[ymin,ymax],color='black',linewidth=0.5,linestyle='dotted')


    return (ymin,ymax)




def PlotV_z(axprofile,rrs_fixincPA,a_min,a_max,v_z_prof_fixincPA,sv_z_prof_fixincPA,BackSide=False,ContinuumGaps=False,label='',RadialScaling=False):


    if BackSide:
        v_z_prof_fixincPA *= -1
    
    rmax=np.max(rrs_fixincPA)

    sv_z_prof_fixincPA=np.nan_to_num(sv_z_prof_fixincPA)
    
    axprofile.set_xlim(a_min,a_max)
    #maskrange=( (rrs_fixincPA > a_min) & (rrs_fixincPA < a_max))
    plotmask = np.where( (rrs_fixincPA >= a_min) & (rrs_fixincPA <= a_max) )
    maskrange=plotmask

    scale_radprofile=1.
    if RadialScaling:
        scale_radprofile=(np.sqrt(rrs_fixincPA[maskrange])/np.sqrt(a_max))

    #from pprint import pprint
    
    #pprint( list(zip(rrs_fixincPA, sv_R_prof_fixincPA, maskrange) ))

    #print(">>>>> median error:",np.median(sv_z_prof_fixincPA[maskrange]))
    #print(">>>>> median values:",np.median(v_z_prof_fixincPA[maskrange]))
    

    ymin=1.01*np.min((v_z_prof_fixincPA[maskrange] - sv_z_prof_fixincPA[maskrange]) * scale_radprofile)
    ymax=1.01*np.max((v_z_prof_fixincPA[maskrange] + sv_z_prof_fixincPA[maskrange]) * scale_radprofile)

          
    axprofile.plot(rrs_fixincPA[plotmask],v_z_prof_fixincPA[plotmask]*scale_radprofile,color='C0',linewidth=1.5,linestyle='solid',label=r'$v_z$')
    if BackSide:
        prefix=r'$-$'
    else:
        prefix=''

    if RadialScaling:
        axprofile.set_ylabel(prefix+r'$\sqrt{R/'+str(a_max)+'} \\times \\tilde{v}_{z}(R)$ / km s$^{-1}$')
    else:
        axprofile.set_ylabel(prefix+r'$\tilde{v}_{z}(R)$ / km s$^{-1}$')

    #print( "wcols v_z")
    #print((np.array(zip(rrs_fixincPA[plotmask],v_z_prof_fixincPA[plotmask]+sv_z_prof_fixincPA[plotmask]))))

    axprofile.fill_between(rrs_fixincPA[plotmask], (v_z_prof_fixincPA[plotmask]+sv_z_prof_fixincPA[plotmask])*scale_radprofile, (v_z_prof_fixincPA[plotmask]-sv_z_prof_fixincPA[plotmask])*scale_radprofile, lw=0.1,color='C0', alpha=0.2, interpolate=True) #, step='mid'

    axprofile.plot(rrs_fixincPA[plotmask], rrs_fixincPA[plotmask]*0.,color='grey',linewidth=0.5,linestyle='dotted')


    if (label != ''):
        axprofile.text(a_min+(a_max-a_min)*0.05,ymax-(ymax-ymin)*0.12,label)

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
            


    print(">>>>> v_z ::",ymin,ymax)


    return (ymin,ymax)
