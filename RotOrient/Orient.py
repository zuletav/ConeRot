import numpy as np



def PlotOrientProfile(axprofile,rregions, PAs, allradsPA, PAerrors, incs, allradsinc,incerrors, psis, allradspsi, psierrors):

    dPAs=PAs-allradsPA
    dincs=incs-allradsinc
    #print( "dPAs",dPAs)
    #print( "allradsPA",allradsPA)
    #print( "PAs",PAs)
    #print( "dincs",dincs)
    
    if (allradspsi < 0):
        dpsis=-psis+allradspsi
    else:
        dpsis=psis-allradspsi

    thislabel = r'$\mathrm{PA}-$%.1f$^\circ$'  %  allradsPA 
    axprofile.plot(rregions,(PAs-allradsPA)  ,linewidth=1.,linestyle='solid',color='C0',label=thislabel,marker='o',markersize=4.)
    #print( "PA errors:",PAerrors)
    axprofile.fill_between(rregions, (PAs-allradsPA)+PAerrors[1,:], (PAs-allradsPA)-PAerrors[0,:], lw=0.1,color='C0', alpha=0.2, interpolate=True) #, step='mid'


    #print( "PAerrors[1,:]",PAerrors[1,:])

    thislabel = r'$i-%.1f^\circ$'  %  allradsinc 
    axprofile.plot(rregions,(incs-allradsinc),linewidth=1.,linestyle='solid',color='C2',label=thislabel,marker='o',markersize=4.)
    #print( "inc errors:",incerrors)
    axprofile.fill_between(rregions, (incs-allradsinc)+incerrors[1,:], (incs-allradsinc)-incerrors[0,:],lw=0.1,color='C2', alpha=0.2, interpolate=True) #, step='mid'

    #print(  "incerrors",incerrors[1,:])

    if (allradspsi < 0):
        thislabel = r'$-\psi-%.1f^\circ$'  %  -allradspsi
        axprofile.plot(rregions,-(psis-allradspsi),linewidth=1.,linestyle='solid',color='C1',label=thislabel,marker='o',markersize=4.)
        axprofile.fill_between(rregions, -(psis-allradspsi)-psierrors[1,:], -(psis-allradspsi)+psierrors[0,:],lw=0.1,color='C1', alpha=0.2, interpolate=True) #, step='mid'
        ymax=1.1*np.max(np.concatenate( (dPAs+PAerrors[1,:], dincs+incerrors[1,:], dpsis+psierrors[1,:], -(psis-allradspsi)+psierrors[1,:]   )))
        ymin=1.05*np.min(np.concatenate( (dPAs-PAerrors[0,:], dincs-incerrors[0,:], dpsis-psierrors[0,:], -(psis-allradspsi)-psierrors[0,:] )))
    else:
        thislabel = r'$\psi-%.1f^\circ$'  %  allradspsi

        #print("VARIABLE INC PA")
        #print(np.array(list(zip(rregions,psis,(psis-allradspsi),(psis-allradspsi)+psierrors[1,:], (psis-allradspsi)-psierrors[0,:]))))
        axprofile.plot(rregions,(psis-allradspsi),linewidth=1.,linestyle='solid',color='C1',label=thislabel,marker='o',markersize=4.)
        #axprofile.plot(rregions,(psis-allradspsi),linewidth=4.,linestyle='solid',color='black')
        axprofile.fill_between(rregions, (psis-allradspsi)+psierrors[1,:], (psis-allradspsi)-psierrors[0,:],lw=0.1,color='C1', alpha=0.2, interpolate=True) #, step='mid'
        ymax=1.05*np.max(np.concatenate( (dPAs+PAerrors[1,:], dincs+incerrors[1,:], dpsis+psierrors[1,:], (psis-allradspsi)+psierrors[1,:]   )))
        ymin=1.05*np.min(np.concatenate( (dPAs-PAerrors[0,:], dincs-incerrors[0,:], dpsis-psierrors[0,:], (psis-allradspsi)-psierrors[0,:] )))



    axprofile.set_ylabel(r'deg')
    axprofile.set_xlabel(r'$r$ / arcsec')


    return (ymin,ymax)



def  PlotOrientProfile_fixincPA(axprofile,rregions_fixincPA, psis_fixincPA, allradspsi_fixincPA, psierrors_fixincPA):


    

    #print("psis_fixincPA",psis_fixincPA)
    #print("allradspsi_fixincPA",allradspsi_fixincPA)

    if (allradspsi_fixincPA < 0):
        dpsis_fixincPA=-psis_fixincPA+allradspsi_fixincPA
    else:
        dpsis_fixincPA=psis_fixincPA-allradspsi_fixincPA



    if (allradspsi_fixincPA < 0):
        thislabel = r'$-\psi-%.1f^\circ$ fix $i,\mathrm{PA}$'  %  -allradspsi_fixincPA
        axprofile.plot(rregions_fixincPA,-(psis_fixincPA-allradspsi_fixincPA),linewidth=1.,linestyle='solid',color='grey',label=thislabel,marker='o',markersize=4.)
        #print( "psi errors fixincPA:",psierrors_fixincPA)
        axprofile.fill_between(rregions_fixincPA, -(psis_fixincPA-allradspsi_fixincPA)-psierrors_fixincPA[1,:], -(psis_fixincPA-allradspsi_fixincPA)+psierrors_fixincPA[0,:],lw=0.1,color='grey', alpha=0.2, interpolate=True) #, step='mid'


        ymax=1.1*np.max(-(psis_fixincPA-allradspsi_fixincPA)+psierrors_fixincPA[0,:])
        ymin=1.05*np.min( -(psis_fixincPA-allradspsi_fixincPA)-psierrors_fixincPA[1,:])

        #print( "ymax=",ymax)

    else:
        thislabel = r'$\psi-%.1f^\circ$ fix $i,\mathrm{PA}$'  %  allradspsi_fixincPA
        axprofile.plot(rregions_fixincPA,(psis_fixincPA-allradspsi_fixincPA),linewidth=1.,linestyle='solid',color='grey',label=thislabel,marker='o',markersize=4.)
        #print( "psi errors fixincPA:",psierrors_fixincPA)
        
        #print("FIX INC PA")
        #print(np.array(list(zip(rregions_fixincPA,psis_fixincPA,psis_fixincPA-allradspsi_fixincPA,(psis_fixincPA-allradspsi_fixincPA)+psierrors_fixincPA[1,:], (psis_fixincPA-allradspsi_fixincPA)-psierrors_fixincPA[0,:]))))

        #axprofile.errorbar(rregions_fixincPA,(psis_fixincPA-allradspsi_fixincPA),yerr=psierrors_fixincPA,color='black')
        axprofile.fill_between(rregions_fixincPA, (psis_fixincPA-allradspsi_fixincPA)+psierrors_fixincPA[1,:], (psis_fixincPA-allradspsi_fixincPA)-psierrors_fixincPA[0,:],lw=0.1,color='grey', alpha=0.2, interpolate=True) #, step='mid'

        ymax=1.05*np.max((psis_fixincPA-allradspsi_fixincPA)+psierrors_fixincPA[1,:])
        ymin=1.05*np.min((psis_fixincPA-allradspsi_fixincPA)-psierrors_fixincPA[0,:])

        #print( "ymax=",ymax)


    return (ymin,ymax)




