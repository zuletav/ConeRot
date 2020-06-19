from astropy.io import fits 
import scipy
import scipy.signal

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable


def addimage(iplotpos,label,atitle,filename_grey,filename_contours=False,VisibleXaxis=False,VisibleYaxis=False,DoBeamEllipse=False,DoGreyCont=False,DoCB=True,Clevs=False,Region=False,vsyst=0.,nplotsx=1,nplotsy=1,Region_Contours=False,SymmetricRange=False,UseScatter=False):

        #print( "nplotsx ", nplotsx, iplotpos)
        
        ax = plt.subplot(nplotsy, nplotsx, iplotpos)
        
        plt.setp(ax.get_xticklabels(),visible=VisibleXaxis)#, fontsize=6)
        plt.setp(ax.get_yticklabels(),visible=VisibleYaxis)#, fontsize=6)
        #ax.tick_params(axis='both',length = 5, width=1., color = 'grey')
        
        ax.tick_params(axis='both',length = 5, width=1., color = 'grey',direction='in',left=True, right=True,bottom=True, top=True)
        
        #plt.gca().set_xticks([-1.,-0.5,0.,0.5,1.0])
        #plt.gca().set_yticks([-1.,-0.5,0.,0.5,1.0])
        
        ax.spines['right'].set_color('grey')
        ax.spines['left'].set_color('grey')
        ax.spines['top'].set_color('grey')
        ax.spines['bottom'].set_color('grey')
        
        
        if ((atitle == 'a')):
                ax.set_ylabel(r'$\delta$  offset / arcsec')
                ax.set_xlabel(r'$\alpha$ offset / arcsec')
                
        #if ((atitle == 'd')):
        #ax.set_ylabel(r'$\delta$  offset / arcsec')

        ######################################################################
        #  GREY + CONTOURS

        # print( "loading filename_grey",filename_grey)
        flog.write("loading filename_grey "+filename_grey+"\n")

        f = fits.open(filename_grey)
        im_grey = f[0].data
        hdr_grey= f[0].header
        cdelt=3600.*hdr_grey['CDELT2']



        j0=0.
        j1=hdr_grey['NAXIS2']-1
        i0=0.
        i1=hdr_grey['NAXIS1']-1

        side=hdr_grey['NAXIS2']*cdelt
        
        a0 = side/2.
        a1 = -side/2.
        d0 = -side/2.
        d1 = side/2.
        
                
        subim_grey = im_grey[int(j0):int(j1),int(i0):int(i1)]
        # print( "i0 "+str(i0),"hdr_grey['CRPIX2']", hdr_grey['CRPIX2'])

        range2=subim_grey.max()

        range1=np.amin(subim_grey)

        clevs = [np.amin(subim_grey),0.,np.amax(subim_grey)]

        clabels=['%.0f' % (clevs[0]),'%.0f' % (clevs[1])]
        if (isinstance(Region,str)):
                file_region=Region
                #print("loading "+file_region)
                
                f = fits.open(file_region)
                im_region = f[0].data
                hdr_region= f[0].header
                f.close()
                subim_region = im_region[int(j0):int(j1),int(i0):int(i1)]
                
                dum=subim_region*subim_grey

                subim_region[np.where(subim_region < 0.)] = 0.
                subim_region[np.where(subim_region > 1.)] = 1.


                if (Clevs=='Region'):
                        subim_grey_filt=scipy.signal.medfilt(subim_grey,5)
                        range1=np.min(subim_grey_filt[np.where(subim_region > 0.99)])
                        range2=np.max(subim_grey_filt[np.where(subim_region > 0.99)])
                        
                        scatter_subim_raw=np.std(subim_grey[np.where(subim_region > 0.9)])

                        
                        
                        mean_subim = np.median(subim_grey[np.where(subim_region > 0.9)])


                        scatter_subim=np.sqrt(np.median((subim_grey[np.where(subim_region > 0.9)] - mean_subim)**2))
                        # print(("Region range: ",range1," ",range2, "soft scatter ",scatter_subim,"hard scatter ",scatter_subim_raw))
                        flog.write("Region range: "+str(range1)+" "+str(range2)+" soft scatter "+str(scatter_subim)+" hard scatter "+str(scatter_subim_raw)+"\n")

                        av_range0=(range1+range2)/2.

                if (SymmetricRange):
                        range0=max(np.fabs(range1),np.fabs(range2))
                        if (UseScatter):
                                #print( "using scatter for ranges")
                                range0=5.*scatter_subim
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
                        clabels=['%.1f' % (clevs[0]),'%.1f' % (clevs[1])]

        elif (isinstance(Clevs,bool)):
                if ('diff' in filename_grey):
                        range1=-0.1
                        range2=0.1
                        clevs = [range1,0.,range2]
                        clabels=['%.1f' % (clevs[0]),'0.','%.1f' % (clevs[2])]
                elif ('wcentered' in filename_grey):
                        range1=np.amin(subim_grey[np.where(subim_grey > 0.)])
                        range2=np.amax(subim_grey)
                        clevs = [range1,range2]
                        clabels=['%.0f' % (clevs[0]),'%.0f' % (clevs[1])]
                elif ('region' in filename_grey):
                        range1=0.
                        range2=np.amax(subim_grey)
                        clevs = [range1,range2]
                        clabels=['%.0f' % (clevs[0]),'%.0f' % (clevs[1])]

                        
        else:
                range1=Clevs[0]
                range2=Clevs[2]
                clevs = Clevs
                clabels=['%.1f' % (clevs[0]-clevs[1]),'0.','%.1f' % (clevs[2]-clevs[1])]


        cmap='ocean_r'
        cmap='RdBu_r'
        if ('wcentered' in filename_grey):
                cmap='magma_r'


        
        plt.imshow(subim_grey, origin='lower', cmap=cmap, #norm=norm,
                   extent=[a0,a1,d0,d1], vmin=range1, vmax=range2)

        plt.plot(0.,0.,marker='*',color='yellow',markersize=1.)

        maxval=np.amax(subim_grey)


        plt.text(a1*0.5,d1*1.1,label,weight='bold',fontsize=12,ha='center')
        
        plt.text(a0*0.9,d1*1.1,atitle,weight='bold',fontsize=12)


        if (DoCB):

                
                axcb=plt.gca()
                divider = make_axes_locatable(axcb)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cb = plt.colorbar(cax=cax,ticks=clevs) #
                
                cax.set_yticklabels(clabels)





        if (Region_Contours):
                levels = np.array([0.1])
                CS = axcb.contour(subim_region, levels, origin='lower', linewidths=0.2,
                                 linestyles = 'solid', 
                                 extent=[a0,a1,d0,d1], colors='green')


                        
        if (isinstance(DoGreyCont,str)):
                filename_grey_cont=DoGreyCont



                f_cont = fits.open(filename_grey_cont)
                im_grey_cont = f_cont[0].data
                hdr_grey_cont= f_cont[0].header
                subim_grey_cont = im_grey_cont[int(j0):int(j1),int(i0):int(i1)]
                #print( "i0 "+str(i0),"hdr_grey['CRPIX2']", hdr_grey['CRPIX2'])
                levels = [-1.65344452058093, -1.6129840879707 , -1.55158301988206, -1.51002707264227   ]


                CS = axcb.contour(subim_grey_cont,levels , origin='lower', linewidths=1.0,
                                  linestyles = 'solid', 
                                  extent=[a0,a1,d0,d1], colors='red')
                
        elif (DoGreyCont):
                levels=np.array((vsyst,))
                CS = axcb.contour(subim_grey,levels , origin='lower', linewidths=0.5,
                                  linestyles = 'solid', 
                                  extent=[a0,a1,d0,d1], colors='green')
                
        

        


        


                
        if (filename_contours!=False):
                
                
                f = fits.open(filename_contours)
                im_cont = f[0].data
                hdr_cont= f[0].header
                cdelt_cont=3600*hdr_cont['CDELT2']
                        

                j0=0.
                j1=hdr_cont['NAXIS2']-1
                i0=0.
                i1=hdr_cont['NAXIS1']-1
                
                
                if ( np.fabs(i0-np.rint(i0))>0):
                        #print( filename_contours)
                        print( 'WARNING: need to interpolate! i0'+str(i0))
            
                subim_cont = im_cont[int(j0):int(j1),int(i0):int(i1)]


                CS = axcb.contour(subim_cont,levels , origin='lower', linewidths=1.5,
                                 linestyles = 'solid', 
                                 extent=[a0,a1,d0,d1], colors='grey')

                CS = axcb.contour(subim_cont,levels , origin='lower', linewidths=0.5,
                                  linestyles = 'solid', 
                                  extent=[a0,a1,d0,d1], colors='black')
                
        if DoBeamEllipse:
                from matplotlib.patches import Ellipse
                
                
                
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
                axcb.add_artist(e)


                




        

                
def exec_summary(basename,fileout,vsyst=0.,basename_errormap=False,file_continuum=False,nplots=3):

        matplotlib.rc('font', family='sans-serif') 
        matplotlib.rcParams.update({'font.size': 12})
        #font = {'family' : 'Arial',
        #        'weight' : 'normal',
        #        'size'   : 10}
        #matplotlib.rc('font', **font)
        
        global flog
        flog=open(basename+"_region_scatter.txt","w+")

        size_marker=10

        #cmaps = ['magma', 'inferno', 'plasma', 'viridis', 'bone', 'afmhot', 'gist_heat', 'CMRmap', 'gnuplot', 'Blues_r', 'Purples_r', 'ocean', 'hot', 'seismic_r']
        gamma=1.0

        plt.figure(figsize=(17, 8))
        #plt.figure(figsize=(4.13, 4))
        nplotsx=nplots
        nplotsy=1
        iplotpos=0


        vrange=3.0

        inc=''
        label=r'$v_\circ$'
        atitle='a'
        filename_contours=False
        filename_grey=basename+'_centered.fits'
        filename_region=basename+'_region_drot.fits'
        iplotpos += 1
        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=True,DoBeamEllipse=False,DoGreyCont=True,Clevs='Region',Region=filename_region,vsyst=vsyst,nplotsx=nplotsx,nplotsy=nplotsy)



        inc=''
        label=r'$\tilde{v}_\circ$'
        atitle='b'
        filename_contours=False
        filename_grey=basename+'_azim_av_drot.fits'
        filename_region=basename+'_region_drot.fits'
        iplotpos += 1
        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=True)



        inc=''
        label=r'$v_\circ - \tilde{v}_\circ$'
        atitle='c'
        filename_contours=False
        filename_grey=basename+'_azim_av_drot_diff.fits' #'_azim_av_drot_diff_b.fits'
        filename_region=basename+'_region_drot.fits'
        iplotpos += 1
        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=True)

        

                        
        if (isinstance(basename_errormap,str)):
                inc=''
                label=r'$\sigma(v_\circ)$'
                atitle='d'
                filename_contours=False
                filename_region=basename+'_region_drot.fits'
                filename_grey=basename_errormap+'_centered.fits'
                #print( "testing for filename_grey",filename_grey)
                import os
                if  (os.path.isfile(filename_grey)):
                        iplotpos += 1
                        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy)
                else:
                        print( "error map not found")


        if (isinstance(file_continuum,str)):
                inc=''
                label=r'cont.'
                atitle='e'
                filename_contours=False
                import  Resamp
                file_im_continuum=basename+'_fullimcont.fits'
                
                Resamp.cube2im(file_continuum,file_im_continuum)
                Resamp.gridding(file_im_continuum,basename+'_centered.fits',fileout=basename+'_subimcont.fits')                
                filename_grey=basename+'_subimcont.fits'
                import os
                if  (os.path.isfile(filename_grey)):
                        iplotpos += 1
                        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs=False,Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy)
                else:
                        print( "continuum  map not found")

                
        plt.subplots_adjust(hspace=0.15)
        plt.subplots_adjust(wspace=0.2)


        #print( fileout)
        #plt.tight_layout()
        print( "USED VSYST=",vsyst)

        plt.savefig(fileout, bbox_inches='tight')

        flog.close()
        return


def exec_summary_allrads(basename,fileout,vsyst=0.,basename_errormap=False,file_m0=False,file_m2=False,file_continuum=False,nplots=4):

        #global nplotsx
        #global nplotsy

        print( "basename:",basename)

        global flog
        flog=open(basename+"_region_scatter.txt","w+")

        matplotlib.rc('font', family='sans-serif') 
        matplotlib.rcParams.update({'font.size': 12})

        size_marker=10

        #cmaps = ['magma', 'inferno', 'plasma', 'viridis', 'bone', 'afmhot', 'gist_heat', 'CMRmap', 'gnuplot', 'Blues_r', 'Purples_r', 'ocean', 'hot', 'seismic_r']
        gamma=1.0

        #plt.figure(figsize=(4.13, 4))
        figsize=(17, 8)
        if (nplots > 6):
                nplotsx=4
                nplotsy=2
        elif (nplots > 3):
                figsize=(14, 8)
                nplotsx=3
                nplotsy=2
        else: 
                figsize=(14, 5)
                nplotsx=3
                nplotsy=1

        plt.figure(figsize=figsize)

        iplotpos=0

        vrange=3.0

        inc=''
        label=r'$v_\circ$'
        atitle='a'
        filename_contours=False
        filename_grey=basename+'_centered.fits'
        filename_region=basename+'_region_drot.fits'
        iplotpos += 1
        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=True,DoBeamEllipse=False,DoGreyCont=True,Clevs='Region',Region=filename_region,vsyst=vsyst,nplotsx=nplotsx,nplotsy=nplotsy)




        inc=''
        label=r'$\tilde{v}_\circ$'
        atitle='b'
        filename_contours=False
        filename_grey=basename+'_allrads_azim_av_drot.fits'
        filename_region=basename+'_region_drot.fits'
        iplotpos += 1
        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=True)





        inc=''
        label=r'$v_\circ - \tilde{v}_\circ$'
        atitle='c'
        filename_contours=False
        filename_grey=basename+'_allrads_azim_av_drot_diff.fits'
        filename_region=basename+'_imregions.fits'
        iplotpos += 1
        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy,Region_Contours=False,SymmetricRange=True)


        inc=''
        label=r'regions'
        atitle='d'
        filename_contours=False
        filename_grey=basename+'_imregions.fits'
        iplotpos += 1
        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=True,DoBeamEllipse=False,Clevs=False,nplotsx=nplotsx,nplotsy=nplotsy)

        
                        
        if (isinstance(basename_errormap,str)):
                inc=''
                label=r'$\sigma(v_\circ)$'
                atitle='d'
                filename_contours=False
                filename_region=basename+'_region_drot.fits'
                filename_grey=basename_errormap+'_centered.fits'
                import os
                if  (os.path.isfile(filename_grey)):
                        iplotpos += 1
                        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy)
                else:
                        print( "error map not found")


        if (isinstance(file_continuum,str)):
                inc=''
                label=r'cont.'
                atitle='e'
                filename_contours=False
                import Resamp
                file_im_continuum=basename+'_fullimcont.fits'
                
                Resamp.cube2im(file_continuum,file_im_continuum)
                                  
                Resamp.gridding(file_im_continuum,basename+'_centered.fits',fileout=basename+'_subimcont.fits')                
                filename_grey=basename+'_subimcont.fits'
                import os
                if  (os.path.isfile(filename_grey)):
                        iplotpos += 1
                        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs=False,Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy)
                else:
                        print( "continuum  map not found")


        if (isinstance(file_m0,str)):
                inc=''
                label=r'I0'
                atitle='f'
                filename_contours=False
                #filename_region=False
                import Resamp 
                file_im_m0=basename+'_fullim_m0.fits'
                
                Resamp.cube2im(file_m0,file_im_m0)
                                  
                Resamp.gridding(file_im_m0,basename+'_centered.fits',fileout=basename+'_subim_m0.fits')                
                filename_grey=basename+'_subim_m0.fits'
                import os
                if  (os.path.isfile(filename_grey)):
                        iplotpos += 1
                        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs=False,Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy)
                else:
                        print( "m0  map not found")


        if (isinstance(file_m2,str)):
                inc=''
                label=r'$\sigma$'
                atitle='g'
                filename_contours=False
                #filename_region=False
                import Resamp
                file_im_m2=basename+'_fullim_m2.fits'
                
                Resamp.cube2im(file_m2,file_im_m2)
                                  
                Resamp.gridding(file_im_m2,basename+'_centered.fits',fileout=basename+'_subim_m2.fits')                
                filename_grey=basename+'_subim_m2.fits'
                import os
                if  (os.path.isfile(filename_grey)):
                        iplotpos += 1
                        addimage(iplotpos,label,atitle,filename_grey,filename_contours,VisibleXaxis=True,VisibleYaxis=False,DoBeamEllipse=False,Clevs='Region',Region=filename_region,nplotsx=nplotsx,nplotsy=nplotsy)
                else:
                        print( "m2  map not found")

                
        plt.subplots_adjust(hspace=0.2)
        plt.subplots_adjust(wspace=0.1)


        print( fileout)
        #plt.tight_layout()

        print( "USED VSYST=",vsyst)
        plt.savefig(fileout, bbox_inches='tight')
        #plt.savefig(fileout)

        flog.close()
        return

