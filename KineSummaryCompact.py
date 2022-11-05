import os
import re
from astropy.io import fits 
import scipy
import scipy.signal
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from  scipy.signal import medfilt2d

include_path='/home/azuleta/common/python/include/'
sys.path.append(include_path)

import ImUtils.Resamp as Resamp
import ImUtils.Cube2Im as Cube2Im

def addimage(iplotpos,label,atitle,filename_grey,filename_contours=False,VisibleXaxis=False,VisibleYaxis=False,DoBeamEllipse=False,DoGreyCont=False,DoCB=False,Clevs=False,Region=False,vsyst=0.,nplotsx=1,nplotsy=1,Region_Contours=False,SymmetricRange=False,UseScatter=False,cmap='ocean_r',filename_weights='',SubtractVsyst=False,ColorBarScale=1.,cblabel='km/s',Zoom=False,side=3.5,RegionOfInterest=False,a_min=-1,a_max=-1,FilterOutliers=True,MaskRegions=True):

        print("vsyst",vsyst)


        print( "nplotsx ", nplotsx, iplotpos)
        ax = plt.subplot(nplotsy, nplotsx, iplotpos)
        # ax=axes[iplotpos]

        plt.setp(ax.get_xticklabels(),visible=VisibleXaxis)
        plt.setp(ax.get_yticklabels(),visible=VisibleYaxis)

        ax.tick_params(axis='both',length = 5, width=1., color = 'grey',direction='in',left=True, right=True,bottom=True, top=True)

        ax.spines['right'].set_color('grey')
        ax.spines['left'].set_color('grey')
        ax.spines['top'].set_color('grey')
        ax.spines['bottom'].set_color('grey')


        if ((iplotpos % nplotsx) == 1):
                ax.set_ylabel(r'$\delta$  offset / arcsec')
        if (iplotpos > (nplotsx*(nplotsy-1))):
                ax.set_xlabel(r'$\alpha$ offset / arcsec')


        print( "loading filename_grey",filename_grey)
        flog.write("loading filename_grey "+filename_grey+"\n")

        f = fits.open(filename_grey)
        im_grey = f[0].data
        hdr_grey= f[0].header
        cdelt=3600.*hdr_grey['CDELT2']

        side0=hdr_grey['NAXIS2']*cdelt

        if Zoom:
                if (side > side0):
                        sys.exit("side too large")
                        
                nx=np.rint(side/cdelt)
                ny=np.rint(side/cdelt)
                
                i_star = np.rint(( (0.)/ hdr_grey['CDELT1'] )+ (hdr_grey['CRPIX1']-1.))
                j_star = np.rint(( (0.)/ hdr_grey['CDELT2'] )+ (hdr_grey['CRPIX2']-1.))

                j0=int(j_star-(ny-1.)/2.+1)
                j1=int(j_star+(ny-1.)/2.+1)
                i0=int(i_star-(nx-1.)/2.+1)
                i1=int(i_star+(nx-1.)/2.+1)
                subim_grey = im_grey[j0:j1,i0:i1]
                
        else:
                side=side0
                i0=0
                i1=hdr_grey['NAXIS1']
                j0=0
                j1=hdr_grey['NAXIS2']

                subim_grey=im_grey.copy()
 

        a0 = side/2.
        a1 = -side/2.
        d0 = -side/2.
        d1 = side/2.


        subim_grey = im_grey[int(j0):int(j1),int(i0):int(i1)]
        print( "i0 "+str(i0),"hdr_grey['CRPIX2']", hdr_grey['CRPIX2'])


        if SubtractVsyst:
                print( "subtracting vsyst=",vsyst)
                subim_grey = subim_grey - vsyst

        if (('_azim_av_drot' in filename_grey) and not ('diff' in filename_grey)):
                print( "vsyst",vsyst)
                subim_grey = subim_grey - vsyst

        range2=subim_grey.max()

        range1=np.amin(subim_grey)

        #clevs = [np.amin(subim_grey),0.,np.amax(subim_grey)]

        clevs = [np.amin(subim_grey),np.amax(subim_grey)]

        #clabels=['%.0f' % (ColorBarScale*clevs[0]),'','%.0f' % (ColorBarScale*clevs[2])]

        clabels=['%.2f' % (ColorBarScale*clevs[0]), '%.2f' % (ColorBarScale*clevs[1])]

        if (isinstance(Clevs,tuple)):
                print( "inherit Clevs")
                clevs = Clevs[0]
                range1=clevs[0]
                range2=clevs[2]
                #clabels=['%.0f' % (clevs[0]),'0.','%.0f' % (clevs[2])]
                #clabels=['%.1f' % (clevs[0]-clevs[1]),'0.','%.1f' % (clevs[2]-clevs[1])]
                clabels=Clevs[1]

        if (isinstance(Region,str)):
                file_region=Region
                print( "loading ", file_region)
                f = fits.open(file_region)
                im_region = f[0].data
                hdr_region= f[0].header
                subim_region = im_region[int(j0):int(j1),int(i0):int(i1)]

        else:
                
                (ny,nx) = subim_grey.shape 
                x=np.arange(1,nx+1)
                y=np.arange(1,ny+1)
                X, Y = np.meshgrid(x, y)
                
                X0 = np.floor(nx/2)+1
                Y0 = np.floor(ny/2)+1
                dxxs = -cdelt *(X-X0)
                dyys = cdelt *(Y-Y0)
                rrs = np.sqrt( (dxxs)**2 + (dyys)**2)
                if ((a_min > 0) and (a_max > 0.)):
                        mask =  (  (rrs > a_min) & (rrs < a_max) )
                else:
                        mask= (rrs >0.)
                subim_region=np.zeros(rrs.shape)
                subim_region[mask]=1.

                

        dum=subim_region*subim_grey
        print( "dum max:",np.max(dum))
        print( "dum min:",np.min(dum))

        subim_region[np.where(subim_region < 0.)] = 0.
        subim_region[np.where(subim_region > 1.)] = 1.

        if MaskRegions:
                subim_grey[np.where(subim_region < 0.5)] = 0.
                
        print("vsyst",vsyst)


        if (Clevs=='Region'): # and (not 'centered.fits' in filename_grey)):
                if FilterOutliers:
                        subim_grey_filt=scipy.signal.medfilt(subim_grey,3)
                else:
                        subim_grey_filt=subim_grey
                range1=np.min(subim_grey_filt[np.where(subim_region > 0.99)])
                range2=np.max(subim_grey_filt[np.where(subim_region > 0.99)])

                scatter_subim_raw=np.std(subim_grey[np.where(subim_region > 0.9)])
                mean_subim=np.median(subim_grey[np.where(subim_region > 0.9)])
                scatter_subim=np.sqrt(np.median((subim_grey[np.where(subim_region > 0.9)] - mean_subim)**2))

                print("Region: range1, range2",range1,range2)

                if (os.path.exists(filename_weights)):
                        print( "loading ", filename_weights)

                        f = fits.open(filename_weights)
                        im_w = f[0].data
                        hdr_w= f[0].header
                        subim_w = im_w[int(j0):int(j1),int(i0):int(i1)]

                        medianval_w = np.median(subim_w)
                        print("medianval_w",medianval_w)

                        range1=np.min(subim_grey_filt[np.where( (subim_region > 0.99) & (subim_w > (medianval_w)))])
                        range2=np.max(subim_grey_filt[np.where( (subim_region > 0.99) & (subim_w > (medianval_w)))])


                        scatter_subim_raw=np.std(subim_grey[np.where( (subim_region > 0.9) & (subim_w > 1E-3))])
                        mean_subim=np.median(subim_grey[np.where( (subim_region > 0.9) & (subim_w > 1E-3))])
                        scatter_subim=np.sqrt(np.median((subim_grey[np.where( (subim_region > 0.9) & (subim_w > 1E-3))] - mean_subim)**2))

                        print("Region + weights: range1, range2",range1,range2)




                        print( "Region range: ",range1," ",range2, "soft scatter ",scatter_subim,"hard scatter ",scatter_subim_raw)
                        flog.write( "Region range: "+str(range1)+" "+str(range2)+" soft scatter "+str(scatter_subim)+" hard scatter "+str(scatter_subim_raw)+"\n")

                av_range0=(range1+range2)/2.
                if (SymmetricRange):
                        range0=max(np.fabs(range1),np.fabs(range2))
                        if (UseScatter):
                                print( "using scatter for ranges")
                                range0=5.*scatter_subim
                        range0= float('%.1f' % (range0))
                        range1=-range0
                        range2=range0
                        clevs = [range1,0.,range2]
                        clabels=['%.2f' % (clevs[0]),'0.','%.2f' % (clevs[2])]
                else:
                        delta_range0=max(np.fabs(range1-av_range0),np.fabs(range2-av_range0))
                        if (UseScatter):
                                range1=-3.*scatter_subim+av_range0
                                range2=3.*scatter_subim+av_range0

                        clevs = [range1,range2]
                        clabels=['%.2f' % (ColorBarScale*clevs[0]),'%.2f' % (ColorBarScale*clevs[1])]



        # norm=colors.PowerNorm(gamma=gamma)




        print( "IMSHOW range1",range1,"range2",range2)
        range1=clevs[0]

        RegionOfInterestStats=True
        if RegionOfInterest and RegionOfInterestStats:
                                
                (ny,nx) = subim_grey.shape 
                x=np.arange(1,nx+1)
                y=np.arange(1,ny+1)
                X, Y = np.meshgrid(x, y)
                
                X0 = np.floor(nx/2)+1
                Y0 = np.floor(ny/2)+1
                dxxs = -cdelt *(X-X0)
                dyys = cdelt *(Y-Y0)
                [x_planet,y_planet] = RegionOfInterest 

                rrs = np.sqrt( (dxxs-x_planet)**2 + (dyys-y_planet)**2)
                mask_planet =  (rrs < 0.5)
                print( "cdelt = ",cdelt,"\n")
                
                #subim_grey[np.where(mask_planet)] = 0.
                sigmavelo_planet = np.std(subim_grey[np.where(mask_planet)])
                print( "sigmavelo_planet ",sigmavelo_planet)


        
        rasterimage=ax.imshow(subim_grey, origin='lower', cmap=cmap, #norm=norm,
                              extent=[a0,a1,d0,d1], vmin=range1, vmax=range2,interpolation='nearest')







        maxval=np.amax(subim_grey)


        print( "a0",a0)
        print( "d0",a0)

        ax.text(a1*0.9,d0*0.9,label,weight='bold',fontsize=12,ha='right')

        ax.text(a0*0.9,d0*0.9,atitle,weight='bold',fontsize=12)

        #clevs = [float(np.amin(subim_grey)),maxval.astype('float')]


        #        clevs = [float(np.amin(subim_grey)),maxval]

        #        if (filename_grey == 'VLA_ABC_modout_nostar_noclump1_smooth.fits'):



        plt.plot(0.,0.,marker='*',color='yellow',markersize=0.2,markeredgecolor='black')
        plt.plot(0.,0.,marker='*',color='yellow',markersize=1.)


        if RegionOfInterest:
                from matplotlib.patches import Ellipse
                [x_planet,y_planet] = RegionOfInterest 
                e = Ellipse(xy=[x_planet,y_planet], width=0.32, height=0.32, angle=0.,color='yellow',fill=False,lw=1.)
                e.set_clip_box(ax.bbox)
                e.set_facecolor('yellow')
                e.set_alpha(0.7)
                ax.add_artist(e)

        if (DoCB):

                #cmap1 = cmap
                #norm = mpl.colors.Normalize(vmin=range1, vmax=range2)
                #fig=plt.gcf()
                #cbar_ax = fig.add_axes([0.92, 0.62, 0.01, 0.15])
                #if (iplotpos > 2):
                #	cbar_ax = fig.add_axes([0.92, 0.22, 0.01, 0.15])
                #cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=colors.Colormap(cmap1), norm=norm, orientation='vertical', ticks=clevs)
                #cb.ax.set_yticklabels(clabels)
                ## cb.ax.set_ylabel('km/s', rotation=270)                
                #cb.ax.tick_params(labelsize=12) 


                #cmap1 = cmap
                #norm = mpl.colors.Normalize(vmin=range1, vmax=range2)
                #fig=plt.gcf()
                #cbar_ax = fig.add_axes([0.92, 0.32, 0.01, 0.35])
                #cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=colors.Colormap(cmap1), norm=norm, orientation='vertical', ticks=clevs)
                #cb.ax.set_yticklabels(clabels)
                #cb.ax.tick_params(labelsize=12) 

                from mpl_toolkits.axes_grid1 import make_axes_locatable
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="3%", pad=0.2)
                cax.yaxis.set_ticks_position('left')
                cax.xaxis.set_ticks_position('top')
                cax.xaxis.set_tick_params(labelsize=12, direction='in')
                cax.yaxis.set_tick_params(labelsize=12, direction='in')
                fmt='%.2f'
                print("drawing colorbar")
                print("clevs",clevs)
                print("clabels",clabels)
                #cb = plt.colorbar(rasterimage, cax=cax, format=fmt, ticks=clevs)
                cb = plt.colorbar(rasterimage, cax=cax, ticks=clevs)
                cb.ax.set_yticklabels(clabels)        
                cb.ax.set_ylabel(cblabel, rotation=270)                



        if (Region_Contours):
                levels = np.array([0.1])
                CS = ax.contour(subim_region, levels, origin='lower', linewidths=0.2,
                                 linestyles = 'solid', 
                                 extent=[a0,a1,d0,d1], colors='green')



        if (isinstance(DoGreyCont,str)):
                filename_grey_cont=DoGreyCont



                f_cont = fits.open(filename_grey_cont)
                im_grey_cont = f_cont[0].data
                hdr_grey_cont= f_cont[0].header
                subim_grey_cont = im_grey_cont[int(j0):int(j1),int(i0):int(i1)]
                print( "i0 "+str(i0),"hdr_grey['CRPIX2']", hdr_grey['CRPIX2'])
                levels = [-1.65344452058093, -1.6129840879707 , -1.55158301988206, -1.51002707264227   ]


                CS = ax.contour(subim_grey_cont,levels , origin='lower', linewidths=1.0,
                                  linestyles = 'solid', 
                                  extent=[a0,a1,d0,d1], colors='red')

        elif (DoGreyCont):
                levels=np.array((vsyst))
                CS = ax.contour(subim_grey,levels , origin='lower', linewidths=0.5,
                                  linestyles = 'solid', 
                                  extent=[a0,a1,d0,d1], colors='green')


        #ax.plot(0.,0.,marker='*',color='yellow',markersize=8)








        if (filename_contours!=False):

                ######################################################################
                #  contours

                print( "loading filename_contours",filename_contours)
                f = fits.open(filename_contours)
                im_cont = f[0].data
                hdr_cont= f[0].header

                subim_cont = im_cont[int(j0):int(j1),int(i0):int(i1)]

                #subim_cont[np.where(subim_region < 0.5)] = 0.


                if ('m2' in filename_contours):
                        levels = np.array([0.6,0.7,0.8]) * np.amax(subim_cont)
                        print( "m2 levels",levels)

                if ('continuum' in filename_contours):
                        levels = np.array([0.2,0.4,0.6,0.8,0.9]) * np.amax(subim_cont)
                        levels = np.array([0.2]) * np.amax(subim_cont)
                        print( "continuum levels",levels)


                print( "subim_cont.shape",subim_cont.shape)

                #CS = axcb.contour(subim_cont,levels , origin='lower', linewidths=1.0,
                #                 linestyles = 'solid', 
                #                 extent=[a0,a1,d0,d1], colors='white')

                #CS = axcb.contour(subim_cont,levels , origin='lower', linewidths=1.0,
                #                  linestyles = 'solid', 
                #                  extent=[a0,a1,d0,d1], colors='black', alpha = 0.5  )

                #alphas=(1.+np.arange(len(levels)))/(len(levels)) * 0.6
                #linewidths=1.0-((0.+np.arange(len(levels)))/(len(levels)-1)) *0.5

                #print( "linewidths")

                alphas=np.ones(len(levels))*0.5
                linewidths=np.ones(len(levels))*2.

                levels_list = levels.tolist()        
                for ilevel in range(len(levels)):
                        alpha_val = alphas[ilevel]
                        alinew=linewidths[ilevel]
                        #print( "alpha_val",alpha_val,"lw",alinew)
                        alevel = levels[ilevel]
                        #print( "alevel ", alevel)
                        CS = ax.contour(subim_cont, [alevel,], origin='lower', linewidths=alinew,
                                          linestyles = 'solid', 
                                          extent=[a0,a1,d0,d1], colors='black', alpha = alpha_val  )

        if DoBeamEllipse:


                #Bmax/2 0.0579669470623286; Bmin/2 0.038567442164739;
                #PA-51.682370436407deg (South of East);

                bmaj = hdr_grey['BMAJ'] * 3600.
                bmin = hdr_grey['BMIN'] * 3600.
                bpa = hdr_grey['BPA']

                print( "bpa ", bpa,"\n")
                print( "bmaj ", bmaj,"\n")
                print( "bmin ", bmin,"\n")
                e = Ellipse(xy=[a1*0.8,d0*0.8], width=bmin, height=bmaj, angle=-bpa,color='blue')
                e.set_clip_box(axcb.bbox)
                e.set_facecolor('yellow')
                e.set_alpha(0.5)
                ax.add_artist(e)






        return  clevs, clabels




        


                

def exec_summary_allrads(workdir,filename_source,vsyst=0.,basename_errormap=False,file_m0=False,file_m2=False,file_continuum=False,Zoom=False,RegionOfInterest=False,side=3.5,AllRads=True,a_min=-1,a_max=-1,FilterOutliers=True,UseScatter=True):


        
        if (not re.search(r"\/$",workdir)):
            workdir+='/'
            print("added trailing back slach to workdir")


        inbasename=os.path.basename(filename_source)
        inbasename=re.sub('\.fits','',inbasename)
        print("inbasename",inbasename)

        global flog
        flog=open(workdir+inbasename+"_region_scatter_diff.txt","w+")


        usescattag=''
        if UseScatter:
                usescattag='_sat'
        

        
        if AllRads:
                fileout = workdir+'fig_'+inbasename+'_report_allrads_diff'+usescattag+'.pdf'
        else:
                fileout = workdir+'fig_'+inbasename+'_report_diff'+usescattag+'.pdf'

        #file_continuum='/Users/simon/common/ppdisks/HD97048/data/b7_2016/continuum/mod_out.fits'

        print( "inbasename:",inbasename)
        matplotlib.rc('font', family='sans-serif') 
        #matplotlib.rcParams.update({'font.size': 16})
        font = {'family' : 'Arial',
                'weight' : 'normal',
                'size'   : 12}

        matplotlib.rc('font', **font)


        size_marker=10

        #cmaps = ['magma', 'inferno', 'plasma', 'viridis', 'bone', 'afmhot', 'gist_heat', 'CMRmap', 'gnuplot', 'Blues_r', 'Purples_r', 'ocean', 'hot', 'seismic_r']
        gamma=1.0

        if file_continuum:
                figsize=(20., 5)
                nplotsx=3
                nplotsy=1
        else:
                figsize=(13, 5)
                nplotsx=2
                nplotsy=1

                
                
        plt.figure(figsize=figsize)

        iplotpos=0

        #inbasename='work_optim_DConeF_vsyst0/mom_1'
        #inbasename='work_optim_DConeT_vsyst0/mom_1'

        vrange=3.0



        #inc=''
        #label=r'$v_\circ - v^m_\circ$ '
        #atitle='c'
        #filename_contours=False
        ##filename_grey='/Users/simon/common/ppdisks/MWC758/VLA-2016/C-array/MWC758_C-array_clean.fits'
        ##filename_grey='/Users/simon/common/ppdisks/MWC758/VLA-2016/reproc/C-array/clean_MWC758_C-array.16chan.ms.fits'
        #filename_grey=workdir+inbasename+'_allrads_azim_av_drot_diff.fits'
        #filename_region=workdir+inbasename+'_imregions.fits'
        #iplotpos += 1
        #(clevs, clabels)=addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=True)


        #inc=''
        #label=r'cont.'
        #atitle='a'
        #filename_contours=False
        ## filename_region=False
        #import Resamp
        #file_im_continuum=workdir+inbasename+'_fullimcont.fits'
        #Resamp.cube2im(file_continuum,file_im_continuum)
        #Resamp.gridding(file_im_continuum,workdir+inbasename+'_centered.fits',fileout=workdir+inbasename+'_subimcont.fits')                
        #filename_grey=workdir+inbasename+'_subimcont.fits'
        #iplotpos += 1
        #addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=True,DoBeamEllipse=False,Clevs=False,Region=False,nplotsx=nplotsx,nplotsy=nplotsy,cmap='ocean_r')

        #cmap='magma_r'


        
        inc=''
        label=r'$v_\circ$'
        atitle=''
        filename_contours=False
        if (file_continuum):
                file_im_continuum=workdir+inbasename+'_fullim_continuum.fits'
                Cube2Im.slice0(file_continuum,file_im_continuum)
                Resamp.gridding(file_im_continuum,workdir+inbasename+'_centered.fits',fileout=workdir+inbasename+'_subim_continuum.fits')                
                filename_contours=workdir+inbasename+'_subim_continuum.fits'
                #filename_contours=False
                

        label=r'a) $v_\circ$'
        filename_grey=workdir+inbasename+'_centered.fits'

        filename_region=workdir+inbasename+'_imregions.fits'
        if (not os.path.isfile(filename_region)):
                filename_region=False
                
        #filename_weights=workdir+inbasename+'_e_wcentered.fits'
        iplotpos += 1
        print("iplotpos",iplotpos)
        (clevs, clabels)=addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=True,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=True,cmap='RdBu_r',UseScatter=False,DoCB=True, SubtractVsyst=True,vsyst=vsyst, Zoom=Zoom,RegionOfInterest=RegionOfInterest,side=side,a_min=a_min,a_max=a_max,FilterOutliers=FilterOutliers)



        label=r'b) $v_\circ - v^m_\circ$'

        if AllRads:
                filename_grey=workdir+inbasename+'_allrads_azim_av_drot_diff.fits'
        else:
                filename_grey=workdir+inbasename+'_azim_av_drot_diff.fits'
                
        #filename_region=workdir+inbasename+'_imregions.fits'
        filename_weights=workdir+inbasename+'_e_wcentered.fits'
        iplotpos += 1
        print("iplotpos",iplotpos)
        (clevs, clabels)=addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=True,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=True,cmap='RdBu_r',filename_weights=filename_weights,UseScatter=UseScatter,DoCB=True,vsyst=vsyst, Zoom=Zoom,RegionOfInterest=RegionOfInterest,side=side,a_min=a_min,a_max=a_max,FilterOutliers=FilterOutliers)

        if file_continuum:
                label=r'c) 225GHz continuum'

                filename_grey=workdir+inbasename+'_subim_continuum.fits'
                filename_contours=False
                iplotpos += 1
                print("iplotpos",iplotpos)
                addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=False,cmap='ocean_r',DoCB=True,ColorBarScale=1E6,cblabel=r'$\mu{\rm Jy\,beam}^{-1}$',vsyst=vsyst, Zoom=Zoom,RegionOfInterest=RegionOfInterest,side=side,a_min=a_min,a_max=a_max,FilterOutliers=FilterOutliers,MaskRegions=False)






        #
        #
        #
        #inc=''
        #atitle='e'
        #import Resamp
        #file_im_m2=workdir+inbasename+'_fullim_m2.fits'
        #label=r'$v_\circ - v^m_\circ + \sigma$'
        #Resamp.cube2im(file_m2,file_im_m2)
        #Resamp.gridding(file_im_m2,workdir+inbasename+'_centered.fits',fileout=workdir+inbasename+'_subim_m2.fits')                
        #filename_grey=workdir+inbasename+'_allrads_azim_av_drot_diff.fits'
        #filename_region=workdir+inbasename+'_imregions.fits'
        #filename_contours=workdir+inbasename+'_subim_m2.fits'
        #iplotpos += 1
        #addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=True)


        plt.subplots_adjust(hspace=0.15)
        plt.subplots_adjust(wspace=0.)


        print( fileout)
        #plt.tight_layout()

        print( "USED VSYST=",vsyst)
        plt.savefig(fileout, bbox_inches='tight', dpi=600)
        #plt.savefig(fileout, bbox_inches='tight')
        #plt.savefig(fileout)

        flog.close()
        return



def exec_summary_faceon(workdir,filename_source,vsyst=0.,basename_errormap=False,file_m0=False,file_m2=False,file_continuum=False, Zoom=False,RegionOfInterest=False,side=3.5,AllRads=True,a_min=-1,a_max=-1,FilterOutliers=True,UseScatter=True):


        
        if (not re.search(r"\/$",workdir)):
            workdir+='/'
            print("added trailing back slach to workdir")


        inbasename=os.path.basename(filename_source)
        inbasename=re.sub('\.fits','',inbasename)
        print("inbasename",inbasename)

        global flog
        flog=open(workdir+inbasename+"_region_scatter_diff.txt","w+")

        allradsstr=''
        if AllRads:
                allradsstr='_allrads'

        usescattag=''
        if UseScatter:
                usescattag='_sat'
        
        if Zoom:
                fileout = workdir+'fig_'+inbasename+'_report'+allradsstr+'_faceon_zoom'+usescattag+'.pdf'
        else:
                fileout = workdir+'fig_'+inbasename+'_report'+allradsstr+'_faceon'+usescattag+'.pdf'

        #file_continuum='/Users/simon/common/ppdisks/HD97048/data/b7_2016/continuum/mod_out.fits'

        print( "inbasename:",inbasename)
        matplotlib.rc('font', family='sans-serif') 
        #matplotlib.rcParams.update({'font.size': 16})
        font = {'family' : 'Arial',
                'weight' : 'normal',
                'size'   : 12}

        matplotlib.rc('font', **font)


        size_marker=10

        #cmaps = ['magma', 'inferno', 'plasma', 'viridis', 'bone', 'afmhot', 'gist_heat', 'CMRmap', 'gnuplot', 'Blues_r', 'Purples_r', 'ocean', 'hot', 'seismic_r']
        gamma=1.0

        if file_continuum:
                figsize=(20., 5)
                nplotsx=3
                nplotsy=1
        else:
                figsize=(13, 5)
                nplotsx=2
                nplotsy=1
                
                
        plt.figure(figsize=figsize)

        iplotpos=0

        #inbasename='work_optim_DConeF_vsyst0/mom_1'
        #inbasename='work_optim_DConeT_vsyst0/mom_1'

        vrange=3.0



        #inc=''
        #label=r'$v_\circ - v^m_\circ$ '
        #atitle='c'
        #filename_contours=False
        ##filename_grey='/Users/simon/common/ppdisks/MWC758/VLA-2016/C-array/MWC758_C-array_clean.fits'
        ##filename_grey='/Users/simon/common/ppdisks/MWC758/VLA-2016/reproc/C-array/clean_MWC758_C-array.16chan.ms.fits'
        #filename_grey=workdir+inbasename+'_allrads_azim_av_drot_diff.fits'
        #filename_region=workdir+inbasename+'_imregions.fits'
        #iplotpos += 1
        #(clevs, clabels)=addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=True)


        #inc=''
        #label=r'cont.'
        #atitle='a'
        #filename_contours=False
        ## filename_region=False
        #import Resamp
        #file_im_continuum=workdir+inbasename+'_fullimcont.fits'
        #Resamp.cube2im(file_continuum,file_im_continuum)
        #Resamp.gridding(file_im_continuum,workdir+inbasename+'_centered.fits',fileout=workdir+inbasename+'_subimcont.fits')                
        #filename_grey=workdir+inbasename+'_subimcont.fits'
        #iplotpos += 1
        #addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=True,DoBeamEllipse=False,Clevs=False,Region=False,nplotsx=nplotsx,nplotsy=nplotsy,cmap='ocean_r')

        #cmap='magma_r'



        inc=''
        label=r'a) $v_\circ$'
        atitle=''
        filename_contours=False
        if (file_continuum):
                file_im_continuum=workdir+inbasename+'_fullim_continuum.fits'
                Cube2Im.slice0(file_continuum,file_im_continuum)

                if AllRads:
                        fileref=workdir+inbasename+'_allrads_resamp_faceon.fits'
                else:
                        fileref=workdir+inbasename+'_resamp_faceon.fits'
                        
                Resamp.gridding(file_im_continuum,fileref,fileout=workdir+inbasename+'_subim_continuum_faceon.fits')                
                filename_contours=workdir+inbasename+'_subim_continuum_faceon.fits'
                #filename_contours=False

        #work_dgaussmoments_LOSNoise_tclean_briggs0.5_pix2_Merid_fixPAinc/im_g_v0_allrads_diff_faceon.fits
        #work_dgaussmoments_LOSNoise_tclean_briggs0.5_pix2_Merid_fixPAinc/im_g_v0_imregions_faceon.fits
        #work_dgaussmoments_LOSNoise_tclean_briggs0.5_pix2_Merid_fixPAinc/im_g_v0_allrads_resamp_faceon.fits

        if AllRads:
                filename_grey=workdir+inbasename+'_allrads_resamp_faceon.fits'
        else:
                filename_grey=workdir+inbasename+'_resamp_faceon.fits'
                
        filename_region=workdir+inbasename+'_imregions_faceon.fits'
        if (not os.path.isfile(filename_region)):
                filename_region=False

        #filename_weights=workdir+inbasename+'_e_wcentered.fits'
        iplotpos += 1
        print("iplotpos",iplotpos)
        (clevs, clabels)=addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=True,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=False,cmap='RdBu_r',UseScatter=False,DoCB=True, SubtractVsyst=True,vsyst=vsyst,Zoom=Zoom,RegionOfInterest=RegionOfInterest,side=side,a_min=a_min,a_max=a_max,FilterOutliers=FilterOutliers)

        label=r'b) $v_\circ - v^m_\circ$'

        if AllRads:
                filename_grey=workdir+inbasename+'_allrads_diff_faceon.fits'
        else:
                filename_grey=workdir+inbasename+'_diff_faceon.fits'
        #filename_region=workdir+inbasename+'_imregions_faceon.fits'
        #filename_weights=workdir+inbasename+'_e_wcentered.fits'
        iplotpos += 1
        print("iplotpos",iplotpos)
        (clevs, clabels)=addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=True,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=True,cmap='RdBu_r',UseScatter=UseScatter,DoCB=True,Zoom=Zoom,RegionOfInterest=RegionOfInterest,side=side,a_min=a_min,a_max=a_max,FilterOutliers=FilterOutliers) # ,vsyst=vsyst)

        if file_continuum:
                label=r'c) 225GHz continuum'

                filename_grey=workdir+inbasename+'_subim_continuum_faceon.fits'
                #filename_contours=False
                iplotpos += 1
                print("iplotpos",iplotpos)
                addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=True,SymmetricRange=False,cmap='ocean_r',DoCB=True,ColorBarScale=1E6,cblabel=r'$\mu{\rm Jy\,beam}^{-1}$',Zoom=Zoom,RegionOfInterest=RegionOfInterest,side=side,a_min=a_min,a_max=a_max,FilterOutliers=FilterOutliers) #,vsyst=vsyst)






        #
        #
        #
        #inc=''
        #atitle='e'
        #import Resamp
        #file_im_m2=workdir+inbasename+'_fullim_m2.fits'
        #label=r'$v_\circ - v^m_\circ + \sigma$'
        #Resamp.cube2im(file_m2,file_im_m2)
        #Resamp.gridding(file_im_m2,workdir+inbasename+'_centered.fits',fileout=workdir+inbasename+'_subim_m2.fits')                
        #filename_grey=workdir+inbasename+'_allrads_azim_av_drot_diff.fits'
        #filename_region=workdir+inbasename+'_imregions.fits'
        #filename_contours=workdir+inbasename+'_subim_m2.fits'
        #iplotpos += 1
        #addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=True)


        plt.subplots_adjust(hspace=0.15)
        plt.subplots_adjust(wspace=0.)


        print( fileout)
        #plt.tight_layout()

        print( "USED VSYST=",vsyst)
        plt.savefig(fileout, bbox_inches='tight', dpi=600)
        #plt.savefig(fileout, bbox_inches='tight')
        #plt.savefig(fileout)

        flog.close()
        return



