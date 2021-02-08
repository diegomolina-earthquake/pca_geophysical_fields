#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 14:25:04 2021

@author: diego
"""

'''
Definition of function to plot 1) the explained variance, 
2) the PCs and EOFs, 3) the reconstructed fields taking into
account a number of k modes 
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.ndimage import gaussian_filter
from matplotlib.colors import LinearSegmentedColormap
 
def explained_variance(c,k,fields,path,Lambda,largo):
    figsize = (10,9)
    c = int(c)
    fields = fields
    L = Lambda[0:50]
    S = np.sum(L)
    VE = (L*100)/S
    modoss = np.linspace(1,50,50) 
    folder   = ['Gravity', 'Locking', 'Friction']    
    
    f, axarr = plt.subplots(2, sharex = True,figsize = figsize)
    axarr[0].plot(modoss,VE,'b', label = 'Percentage Variance')
    axarr[0].plot(modoss,VE,'b*')
    axarr[0].set_title('Explained variance')
    axarr[0].legend(loc='upper right').get_frame().set_alpha(0.5)
    axarr[0].set_ylabel('Explained Variance [%]')
    axarr[0].grid()
    axarr[1].plot(modoss,np.cumsum(VE),'k',label = 'Percentage Cumulative Variance')
    axarr[1].plot(modoss,np.cumsum(VE),'k*')
    axarr[1].legend(loc='lower right').get_frame().set_alpha(0.5)
    axarr[1].set_xlabel(r'$\lambda_k$ modes or components')
    axarr[1].set_ylabel('Cumulative explained variance [%]')
    axarr[1].grid()
    print 'The explained variance by the '+str(k)+' mode is %'+str(VE[k-1])+''
    print 'The cumulative explained variance by the '+str(k)+' modes is %'+str(np.cumsum(VE)[k-1])+''
 
    for i in range(len(fields)):
        
        if c == 1: 
            h = fields[0]
	    p = folder[h-1]
            
        if c == 2:

            if sum(fields) == 3:
            
            	p = '/'+folder[0]+'_'+folder[1]+'/'    
            ###################################################################
            ###################################################################
            if sum(fields) == 4:
                                          
                p = '/'+folder[0]+'_'+folder[2]+'/'    
            ###################################################################
            ###################################################################
            if sum(fields) == 5:
              
                p = '/'+folder[1]+'_'+folder[2]+'/'                
            
        if c == 3:
            
            p = '/'+folder[0]+'_'+folder[1]+'_'+folder[2]+'/'  
     
    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/		     var_explicada.png', format='png',dpi=300,transparent = True)
    
    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/		     var_explicada.svg', format='svg',dpi=300,transparent = True)
    
def plot_eofs(c,k,lat,EOFs,EOFsum,fields,path,largo):
    fields = fields
    c = int(c)

    figsize = (10,6)
    folder   = ['Gravity', 'Locking', 'Friction'] 
    label   = ['Gravity eof', 'Locking eof', 'Friction eof']
    label_s = ['Gravity eof sum', 'Locking eof sum', 'Friction eof sum']

    color = ['b','r','g']
    f, axarr = plt.subplots(1, sharex = True,figsize = figsize)
    
    
    for i in range(len(fields)):
        
        if c == 1: 
            h = fields[0]
	    p = folder[h-1]

            
            axarr.plot(lat,EOFs[:,0],color[h-1], 
                       linewidth = 4, label = label[h-1])
            axarr.plot(lat,EOFsum[:,0],color[h-1], 
                       linewidth = 1, label = label_s[h-1],alpha = 0.6)
            
        if c == 2:

            if sum(fields) == 3 and fields[i] == 1:
                           
                axarr.plot(lat,EOFs[:,0],color[0], 
                           linewidth = 4, label = label[0])
                axarr.plot(lat,EOFsum[:,0],color[0], 
                           linewidth = 1, label = label_s[0],alpha = 0.6)
                p = '/'+folder[0]+'_'+folder[1]+'/'    
            if sum(fields) == 3 and fields[i] == 2:
                           
                axarr.plot(lat,EOFs[:,1],color[1], 
                           linewidth = 4, label = label[1])
                axarr.plot(lat,EOFsum[:,1],color[1], 
                           linewidth = 1, label = label_s[1],alpha = 0.6)
                p = '/'+folder[0]+'_'+folder[1]+'/'    
            ###################################################################
            ###################################################################
            if sum(fields) == 4 and fields[i] == 1:
                           
                axarr.plot(lat,EOFs[:,0],color[0], 
                           linewidth = 4, label = label[0])
                axarr.plot(lat,EOFsum[:,0],color[0], 
                           linewidth = 1, label = label_s[0],alpha = 0.6)
                p = '/'+folder[0]+'_'+folder[2]+'/'                    
            if sum(fields) == 4 and fields[i] == 3:
                           
                axarr.plot(lat,EOFs[:,1],color[2], 
                           linewidth = 4, label = label[2])
                axarr.plot(lat,EOFsum[:,1],color[2], 
                           linewidth = 1, label = label_s[2],alpha = 0.6)
                p = '/'+folder[0]+'_'+folder[2]+'/'    
            ###################################################################
            ###################################################################
            if sum(fields) == 5 and fields[i] == 2:
                           
                axarr.plot(lat,EOFs[:,0],color[1], 
                           linewidth = 4, label = label[1])
                axarr.plot(lat,EOFsum[:,0],color[1], 
                           linewidth = 1, label = label_s[1],alpha = 0.6)
                p = '/'+folder[1]+'_'+folder[2]+'/'                    
            if sum(fields) == 5 and fields[i] == 3:
                           
                axarr.plot(lat,EOFs[:,1],color[2], 
                           linewidth = 4, label = label[2])
                axarr.plot(lat,EOFsum[:,1],color[2], 
                           linewidth = 1, label = label_s[2],alpha = 0.6)
                p = '/'+folder[1]+'_'+folder[2]+'/'                
            
        if c == 3:
            
            axarr.plot(lat,EOFs[:,fields[i]-1],color[fields[i]-1], 
                       linewidth = 4, label = label[fields[i]-1])
            axarr.plot(lat,EOFsum[:,fields[i]-1],color[fields[i]-1], 
                       linewidth = 1, label = label_s[fields[i]-1],alpha = 0.6)
            p = '/'+folder[0]+'_'+folder[1]+'_'+folder[2]+'/'                
    axarr.set_title('EOFs')
    axarr.legend(loc='upper left',fontsize = 'x-small').get_frame().set_alpha(0.5)
    axarr.set_ylabel('Principal vectors mode '+str(k)+'')
    axarr.set_xlabel('Latitud [Degrees]')
    #axarr.set_xlim(-44,-18)
    axarr.grid()
    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/EOFs_mode_'+str(k)+'.png', format='png',dpi=300,transparent = True)
    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/EOFs_mode_'+str(k)+'.svg', format='svg',dpi=300,transparent = True)


    
def plot_pcs(c,k,dist,PCs,PCsum,fields,path,largo):
    fields = fields
    c = int(c)
    folder   = ['Gravity', 'Locking', 'Friction'] 
    figsize = (10,7)
    label   = ['Gravity eof', 'Locking eof', 'Friction eof']
    label_s = ['Gravity eof sum', 'Locking eof sum', 'Friction eof sum']

    color = ['b','r','g']
    f, axarr = plt.subplots(1, sharex = True,figsize = figsize)
    
    n = 0
    dist_PCs = []
    for i in range(len(PCs[0,:])) :
        dist_PCs.append(n)
        n = n+dist/1000.  
        
    for i in range(len(fields)):
        
        if c == 1: 
            h = fields[0]
	    p = folder[h-1]

            
            axarr.plot(dist_PCs,PCs[0,:],color[h-1], 
                       linewidth = 4, label = label[h-1])
            axarr.plot(dist_PCs,PCsum[0,:],color[h-1], 
                       linewidth = 1, label = label_s[h-1],alpha = 0.6)
            
        if c == 2:
            
            if sum(fields) == 3 and fields[i] == 1:
                           
                axarr.plot(dist_PCs,PCs[0,:],color[0], 
                           linewidth = 4, label = label[0])
                axarr.plot(dist_PCs,PCsum[0,:],color[0], 
                           linewidth = 1, label = label_s[0],alpha = 0.6)
                p = '/'+folder[0]+'_'+folder[1]+'/'                    
            if sum(fields) == 3 and fields[i] == 2:
                           
                axarr.plot(dist_PCs,PCs[1,:],color[1], 
                           linewidth = 4, label = label[1])
                axarr.plot(dist_PCs,PCsum[1,:],color[1], 
                           linewidth = 1, label = label_s[1],alpha = 0.6)
                p = '/'+folder[0]+'_'+folder[1]+'/'    
            ###################################################################
            ###################################################################
            if sum(fields) == 4 and fields[i] == 1:
                           
                axarr.plot(dist_PCs,PCs[0,:],color[0], 
                           linewidth = 4, label = label[0])
                axarr.plot(dist_PCs,PCsum[0,:],color[0], 
                           linewidth = 1, label = label_s[0],alpha = 0.6)
                p = '/'+folder[0]+'_'+folder[2]+'/'                                    
            if sum(fields) == 4 and fields[i] == 3:
                           
                axarr.plot(dist_PCs,PCs[1,:],color[2], 
                           linewidth = 4, label = label[2])
                axarr.plot(dist_PCs,PCsum[1,:],color[2], 
                           linewidth = 1, label = label_s[2],alpha = 0.6)
                p = '/'+folder[0]+'_'+folder[2]+'/'                    
            ###################################################################
            ###################################################################
            if sum(fields) == 5 and fields[i] == 2:
                           
                axarr.plot(dist_PCs,PCs[0,:],color[1], 
                           linewidth = 4, label = label[1])
                axarr.plot(dist_PCs,PCsum[0,:],color[1], 
                           linewidth = 1, label = label_s[1],alpha = 0.6)
                p = '/'+folder[1]+'_'+folder[2]+'/'                
            if sum(fields) == 5 and fields[i] == 3:
                           
                axarr.plot(dist_PCs,PCs[1,:],color[2], 
                           linewidth = 4, label = label[2])
                axarr.plot(dist_PCs,PCsum[1,:],color[2], 
                           linewidth = 1, label = label_s[2],alpha = 0.6)
                p = '/'+folder[1]+'_'+folder[2]+'/'                
            
        if c == 3:
            
            axarr.plot(dist_PCs,PCs[fields[i]-1,:],color[fields[i]-1], 
                       linewidth = 4, label = label[fields[i]-1])
            axarr.plot(dist_PCs,PCsum[fields[i]-1,:],color[fields[i]-1], 
                       linewidth = 1, label = label_s[fields[i]-1],alpha = 0.6)
            p = '/'+folder[0]+'_'+folder[1]+'_'+folder[2]+'/'                            
    axarr.set_title('PCs')
    axarr.legend(loc='upper left',fontsize = 'x-small').get_frame().set_alpha(0.5)
    axarr.set_ylabel('Principal components mode '+str(k)+'')
    axarr.set_xlabel('Trench distance [km]')
    axarr.grid()
    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/PCs_mode_'+str(k)+'.png', 	     format='png',dpi=300,transparent = True)
    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/PCs_mode_'+str(k)+'.svg', 	     format='svg',dpi=300,transparent = True)
     
def plot_recontruction(var_field,c,k,fields,nx,ny,xi,yi,path,path_g,largo):
    c = int(c)
    fields = fields
    figsize = (12,10)
    fig = plt.figure(4,figsize = figsize)
    folder   = ['Gravity', 'Locking', 'Friction']  

    #definition of map geographycal extension
    latmin_map = -46     #<-Latitud  min (degrees) (South)
    latmax_map = -18     #<-Latitud  max (degrees) (North)
    lonmin_map = -76     #<-Longitud min (degrees) (West)
    lonmax_map = -68     #<-Longitud max (degrees) (East)

    fosa = np.loadtxt(path_g+'/Trench.xy')
    fosa_lon = fosa[:,0]-360 ; fosa_lat = fosa[:,1]            #>> Latitud y longitud de los puntos de la fosa.
    #needed to plot trench geometry
    info_nzsa = fosa
    x = info_nzsa.T[0]; y = info_nzsa.T[1]
    x2 = x[0::12]; y2 = y[0::12]
    dx = np.diff(x2); dy = np.diff(y2)
    length = np.sqrt(dx**2 + dy**2)
    dx /= length; dy /= length




    #3.2) se definen los valores maximos y minimos que tendra la paleta
    mas=15
    mis=-15
    blevels_top =np.linspace(mis,mas,300)  
    mindem = 5000
    rango  = 10000
    
    
    mapu1 = LinearSegmentedColormap.from_list('mycmap', 
    										[( 0., 'midnightblue'),
                                                  ((-4200.+mindem)/rango, 'navy'),
    										((-4000.+mindem)/rango, 'darkblue'),
                                                 ((-2000.+mindem)/rango, 'steelblue'),
    										((-1800.+mindem)/rango, 'c'),
    										((-1100.+mindem)/rango, 'c'),
                                             ((-200.+mindem)/rango, 'w'),
                                                         ((0.+mindem)/rango, 'white'),
                                             ((200.+mindem)/rango,'w'),
    										((1100.+mindem)/rango,'y'),
    										((1800.+mindem)/rango, 'y'), 
                                                 ((2000.+mindem)/rango, 'orangered'), 
    										((4000.+mindem)/rango,'red'),
                                                    ((4200.+mindem)/rango, 'firebrick'),										 
    										(1., 'darkred')])


    mapu2 = LinearSegmentedColormap.from_list('mycmap', 
    										[( 0., 'darkmagenta'),
                                                  ((-4300.+mindem)/rango, 'magenta'),
    										((-4000.+mindem)/rango, 'darkblue'),
                                                 ((-2300.+mindem)/rango, 'c'),
    										((-1800.+mindem)/rango, 'c'),
    										((-1100.+mindem)/rango, 'w'),
                                             ((-200.+mindem)/rango, 'w'),
                                                         ((0.+mindem)/rango, 'white'),
                                             ((200.+mindem)/rango,'w'),
    										((1100.+mindem)/rango,'w'),
    										((1800.+mindem)/rango, 'y'), 
                                                 ((2300.+mindem)/rango, 'y'), 
    										((4000.+mindem)/rango, 'orange'),
                                                    ((4300.+mindem)/rango, 'tomato'),										 
    										(1., 'orangered')])



    palete   = [plt.cm.seismic, mapu1, mapu2]
    mis = np.array([-0.4,-0.4,-0.4]); mas = np.array([0.4,0.4,0.4])

    for i in range(len(fields)):
        
        ax = fig.add_subplot(1,c,i+1)
        plt.subplots_adjust(wspace=-0.7)
        
        map = Basemap(projection='merc', resolution="h", llcrnrlon=lonmin_map, llcrnrlat=latmin_map, 
              urcrnrlon=lonmax_map, urcrnrlat=latmax_map)                        
        map.drawcoastlines(zorder=21,color = 'k')
        map.drawcountries(zorder=21)
        parallels = map.drawparallels(np.linspace(-16, -46,7), labels=[0, 0, 0, 0], fmt="%.0f", #4
                        dashes=[2000, 2000], zorder=24,linewidth = 0.0000001)
        for m in parallels:
            try:
                parallels[m][1][0].set_rotation(90)
            except:
                pass      
        map.drawmeridians([-71,-74], labels=[0, 0, 0, 0], fmt="%.0f",
                       dashes=[5000, 5000],zorder=24,linewidth = 0.0000001)    
        
        map.plot(x, y, color= 'k', linewidth=0.5, zorder=41, latlon=True, alpha=1,clip_on=True)
        map.quiver(x2[0:20], y2[0:20], dy[0:20], -dx[0:20], latlon=True,headaxislength=8, headlength=8,headwidth=18, color= 'k', scale=40, zorder=41, alpha=1, clip_on=True)
            
        
        if c == 1:
            h = fields[0]
	    p = folder[h-1]
            Y = gaussian_filter(var_field[0], sigma=1)
            mis1=mis[h-1];mas1=mas[h-1]
            blevels_top =np.linspace(mis1,mas1,100)  
            z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[h-1],extend='both',zorder=20,latlon = True)
            cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%")                           
        if c == 2:
            
            if sum(fields) == 3 and fields[i] == 1:
                
                Y = gaussian_filter(var_field[0], sigma=1)
                mis1=mis[0];mas1=mas[0]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[0],extend='both',zorder=20,latlon = True)
           	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%") 
                p = '/'+folder[0]+'_'+folder[1]+'/'                                    
            if sum(fields) == 3 and fields[i] == 2:
                
                Y = gaussian_filter(var_field[1], sigma=1)
                mis1=mis[1];mas1=mas[1]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[1],extend='both',zorder=20,latlon = True)
            	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%")                           
                p = '/'+folder[0]+'_'+folder[1]+'/'                                   
            ###################################################################
            ###################################################################
            if sum(fields) == 4 and fields[i] == 1:
                
                Y = gaussian_filter(var_field[0], sigma=1)
                mis1=mis[0];mas1=mas[0]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[0],extend='both',zorder=20,latlon = True)
            	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%")  
                p = '/'+folder[0]+'_'+folder[2]+'/'                                                   
            if sum(fields) == 4 and fields[i] == 3:
                           
                Y = gaussian_filter(var_field[1], sigma=1)
                mis1=mis[2];mas1=mas[2]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[2],extend='both',zorder=20,latlon = True)
            	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%") 
                p = '/'+folder[0]+'_'+folder[2]+'/'                                                  
            ###################################################################
            ###################################################################
            if sum(fields) == 5 and fields[i] == 2:
                           
                Y = gaussian_filter(var_field[0], sigma=1)
                mis1=mis[1];mas1=mas[1]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[1],extend='both',zorder=20,latlon = True)
            	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%")                 
                p = '/'+folder[1]+'_'+folder[2]+'/'                                                               
            if sum(fields) == 5 and fields[i] == 3:
                           
                Y = gaussian_filter(var_field[1], sigma=1)
                mis1=mis[2];mas1=mas[2]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[2],extend='both',zorder=20,latlon = True)
		cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%") 
           	p = '/'+folder[1]+'_'+folder[2]+'/'                                                               
            
        if c == 3:
            
            
            Y = gaussian_filter(var_field[fields[i]-1], sigma=1)
            mis1=mis[fields[i]-1];mas1=mas[fields[i]-1]
            blevels_top =np.linspace(mis1,mas1,100)  
            z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[fields[i]-1],extend='both',zorder=20,latlon = True ,ax = ax)
            cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%") 
            p = '/'+folder[0]+'_'+folder[1]+'_'+folder[2]+'/'                            
	    
    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/fields_mode_'+str(k)+'.png', format='png',dpi=300,transparent=True)
    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/fields_mode_'+str(k)+'.svg', format='svg',dpi=300,transparent=True)

    
                
                
def plot_recontruction_sum(var_field,c,k,fields,nx,ny,xi,yi,path,path_g,largo):

    c = int(c)
    fields = fields
    figsize = (12,10)
    fig = plt.figure(6,figsize = figsize)
    folder   = ['Gravity', 'Locking', 'Friction']  

    #definition of map geographycal extension
    latmin_map = -46     #<-Latitud  min (degrees) (South)
    latmax_map = -18     #<-Latitud  max (degrees) (North)
    lonmin_map = -76     #<-Longitud min (degrees) (West)
    lonmax_map = -68     #<-Longitud max (degrees) (East)

    fosa = np.loadtxt(path_g+'/Trench.xy')
    fosa_lon = fosa[:,0]-360 ; fosa_lat = fosa[:,1]            #>> Latitud y longitud de los puntos de la fosa.
    #needed to plot trench geometry
    info_nzsa = fosa
    x = info_nzsa.T[0]; y = info_nzsa.T[1]
    x2 = x[0::12]; y2 = y[0::12]
    dx = np.diff(x2); dy = np.diff(y2)
    length = np.sqrt(dx**2 + dy**2)
    dx /= length; dy /= length

    #3.2) se definen los valores maximos y minimos que tendra la paleta
    mas=15
    mis=-15
    blevels_top =np.linspace(mis,mas,300)  
    mindem = 5000
    rango  = 10000
    
    mapu1 = LinearSegmentedColormap.from_list('mycmap', 
    										[( 0., 'midnightblue'),
                                                  ((-4200.+mindem)/rango, 'navy'),
    										((-4000.+mindem)/rango, 'darkblue'),
                                                 ((-2000.+mindem)/rango, 'steelblue'),
    										((-1800.+mindem)/rango, 'c'),
    										((-1100.+mindem)/rango, 'c'),
                                             ((-200.+mindem)/rango, 'w'),
                                                         ((0.+mindem)/rango, 'white'),
                                             ((200.+mindem)/rango,'w'),
    										((1100.+mindem)/rango,'y'),
    										((1800.+mindem)/rango, 'y'), 
                                                 ((2000.+mindem)/rango, 'orangered'), 
    										((4000.+mindem)/rango,'red'),
                                                    ((4200.+mindem)/rango, 'firebrick'),										 
    										(1., 'darkred')])
    mapu2 = LinearSegmentedColormap.from_list('mycmap', 
    										[( 0., 'darkmagenta'),
                                                  ((-4300.+mindem)/rango, 'magenta'),
    										((-4000.+mindem)/rango, 'darkblue'),
                                                 ((-2300.+mindem)/rango, 'c'),
    										((-1800.+mindem)/rango, 'c'),
    										((-1100.+mindem)/rango, 'w'),
                                             ((-200.+mindem)/rango, 'w'),
                                                         ((0.+mindem)/rango, 'white'),
                                             ((200.+mindem)/rango,'w'),
    										((1100.+mindem)/rango,'w'),
    										((1800.+mindem)/rango, 'y'), 
                                                 ((2300.+mindem)/rango, 'y'), 
    										((4000.+mindem)/rango, 'orange'),
                                                    ((4300.+mindem)/rango, 'tomato'),										 
    										(1., 'orangered')])


    palete   = [plt.cm.seismic, mapu1, mapu2]
    
    mis = np.array([-0.4,-0.4,-0.4]); mas = np.array([0.4,0.4,0.4])

    for i in range(len(fields)):
        
        ax = fig.add_subplot(1,c,i+1)
        plt.subplots_adjust(wspace=-0.7)
        
        map = Basemap(projection='merc', resolution="h", llcrnrlon=lonmin_map, llcrnrlat=latmin_map, 
              urcrnrlon=lonmax_map, urcrnrlat=latmax_map)                        
        map.drawcoastlines(zorder=21,color = 'k')
        map.drawcountries(zorder=21)
        parallels = map.drawparallels(np.linspace(-16, -46,7), labels=[0, 0, 0, 0], fmt="%.0f", #4
                        dashes=[2000, 2000], zorder=24,linewidth = 0.000001)
        for m in parallels:
            try:
                parallels[m][1][0].set_rotation(90)
            except:
                pass      
        map.drawmeridians([-71,-74], labels=[0, 0, 0, 0], fmt="%.0f",
                       dashes=[5000, 5000],zorder=24,linewidth = 0.000001)    
        
        map.plot(x, y, color= 'k', linewidth=0.5, zorder=41, latlon=True, alpha=1,clip_on=True)
        map.quiver(x2[0:20], y2[0:20], dy[0:20], -dx[0:20], latlon=True,headaxislength=8, headlength=8,headwidth=18, color= 'k', scale=40, zorder=41, alpha=1, clip_on=True)
            
        
        if c == 1:
            h = fields[0]
	    p = folder[h-1]
            Y = gaussian_filter(var_field[0], sigma=1)
            mis1=mis[h-1];mas1=mas[h-1]
            blevels_top =np.linspace(mis1,mas1,100)  
            z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[h-1],extend='both',zorder=20,latlon = True)
            cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%")                           
        if c == 2:
            
            if sum(fields) == 3 and fields[i] == 1:
                
                Y = gaussian_filter(var_field[0], sigma=1)
                mis1=mis[0];mas1=mas[0]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[0],extend='both',zorder=20,latlon = True)
            	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%") 
                p = '/'+folder[0]+'_'+folder[1]+'/'                                                    
            if sum(fields) == 3 and fields[i] == 2:
                
                Y = gaussian_filter(var_field[1], sigma=1)
                mis1=mis[1];mas1=mas[1]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[1],extend='both',zorder=20,latlon = True)
            	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%")                           
                p = '/'+folder[0]+'_'+folder[1]+'/'                                                   
            ###################################################################
            ###################################################################
            if sum(fields) == 4 and fields[i] == 1:
                
                Y = gaussian_filter(var_field[0], sigma=1)
                mis1=mis[0];mas1=mas[0]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[0],extend='both',zorder=20,latlon = True)
            	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%") 
                p = '/'+folder[0]+'_'+folder[2]+'/'                                                                   
            if sum(fields) == 4 and fields[i] == 3:
                           
                Y = gaussian_filter(var_field[1], sigma=1)
                mis1=mis[2];mas1=mas[2]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[2],extend='both',zorder=20,latlon = True)
            	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%")   
                p = '/'+folder[0]+'_'+folder[2]+'/'                                                               
            ###################################################################
            ###################################################################
            if sum(fields) == 5 and fields[i] == 2:
                           
                Y = gaussian_filter(var_field[0], sigma=1)
                mis1=mis[1];mas1=mas[1]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[1],extend='both',zorder=20,latlon = True)
            	cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%")                 
		p = '/'+folder[1]+'_'+folder[2]+'/'                  
            if sum(fields) == 5 and fields[i] == 3:
                           
                Y = gaussian_filter(var_field[1], sigma=1)
                mis1=mis[2];mas1=mas[2]
                blevels_top =np.linspace(mis1,mas1,100)  
                z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[2],extend='both',zorder=20,latlon = True)
		cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%") 
		p = '/'+folder[1]+'_'+folder[2]+'/'                                
            
        if c == 3:
            
            print i
            
            Y = gaussian_filter(var_field[fields[i]-1], sigma=1)
            mis1=mis[fields[i]-1];mas1=mas[fields[i]-1]
            blevels_top =np.linspace(mis1,mas1,100)  
            z  = map.contourf(xi,yi,Y,vmin=mis1,vmax=mas1,levels=blevels_top,cmap=palete[fields[i]-1],extend='both',zorder=20,latlon = True ,ax = ax)
            cb = map.colorbar(z,ticks= [-0.4,0,0.4],location='bottom', extend = 'both',pad="4%") 
            p = '/'+folder[0]+'_'+folder[1]+'_'+folder[2]+'/'                            
	    

    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/fields_mode_sum_'+str(k)+'.png', format='png',dpi=300,transparent=True)
    plt.savefig(path+'/PCA_plots/ppf_'+str(largo/1000.)+'km/'+str(c)+'/'+p+'/fields_mode_sum_'+str(k)+'.svg', format='svg',dpi=300,transparent=True)

            
                
    
     
        
         
     
        
     
