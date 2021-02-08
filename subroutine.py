
#contribution of Felipe Veras. PhD studendt
#felipe.orlando.vera.sanhueza@gfz-potsdam.de
from math import *
from numpy import *
import scipy as scp
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from fractions import Fraction
import numpy as np




def subfalla(dip,W,L,nx,ny,lon0,lat0,strike):
    #Se obtienen las coordenadas centrales asumuendo una geometrica de falla rectangular.
    dx = L / nx; dy = W * cos(dip) / ny   #largo y ancho de las subfallas
     
    xs = []; ys = []

    #Se obtiene las coordenadas del centro de cada subfalla (latitud y longitud)
    #Esto permmite asociarlo con el deslizamiento obtenido y graficarlo en un mapa.
    alpha = degrees(strike) ; beta = alpha + 270.0
    
    lat1,  lon1, alpha21 = vinc_pt( lat0, lon0, alpha , dx/2 )
    lat2,  lon2, alpha21 = vinc_pt( lat1, lon1,  beta , dy/2 )  #Coordenadas centrales de la primera subfalla.
    for j in range(ny):
        for i in range(nx):
            lat3 = lat2 ; lon3 = lon2
            if i > 0:
             lat3,  lon3, alpha21 = vinc_pt( lat2, lon2, alpha , i*dx )
            lat4 = lat3 ; lon4 = lon3
            if j > 0 :
             lat4 , lon4, alpha21 = vinc_pt( lat3, lon3,  beta , j*dy )
            xs.append(lon4)
            ys.append(lat4)

    lon_slip = array([xs])[0]
    lat_slip = array([ys])[0]

    return  lon_slip, lat_slip

def coordenadas_subfallas(ny,nx,dip,W,L,lon0,lat0,strike,fosa_lon,fosa_lat,delta_lat):

    #La funcion obtiene las coordenadas centrales de las subfallas de acuerdo a la
    #configuracion del slab model, junto con las coordendas de los vertices de esta.    
    
    #El problema comienza generando una falla rectangular dividida en ny*nx subfallas,
    #para posteriormente desplazarlas en direccion hacia la fosa y obtener asi una configuracion
    #de subfallas de acuerdo al slab model.
    #Felipe Vera S.
    #fveras@udec.cl    
    
    #1. Se determina latitud y longitud de cada subfalla asumiendo una destribucion rectangular.
    lon_sl, lat_sl = subfalla(dip,W,L,nx,ny,lon0,lat0,strike) 
    dx = L / nx; dy = W * cos(dip) / ny
    
    #Variables para almacenar la posicion de central de las subfallas de acuerso al Slab Model..
    lon_central=[]
    lat_central=[]

    #1.1 Se obtiene posicion central de las subfallas de acuerdo a la geometria rectangular.
    k=0
    for i_nx in range(nx):
        for i_ny in range(ny):
            lon_central.append(lon_sl[k])
            lat_central.append(lat_sl[k])               
            k=k+1; 
        
    #1.2 Se generan vertices de las subfallas de acuerdo a la geometria rectangular.
    FIL=nx*ny  #Se evaluaran todas las subfallas.
    COL=5      #Cada subfalla tiene cinco vertices a dibujar (el ultimo coincide con el primero)

    #Se decalran variables para almacenar vertices de todas las subfallas.
    vertices_lon = [] 
    vertices_lat = []
    
    #Se crean matrices de ceros para rellenar con los vertices de todas las subfallas.
    for i in range(FIL):                #Numero de filas de la matriz
        vertices_lon.append([0]*COL)    #Numero de columnas de la matriz, se inicializa con ceros
        vertices_lat.append([0]*COL)
    
    #1.3 Se declaran vertices claves: Son los vertices superior e inferior de las subfallas que estan
    #Junto a la fosa. Sirven para buscar los valores de referencia con el modelo slab1.0    
    vertice_claveinf_lon=[]
    vertice_claveinf_lat=[]
    vertice_clavesup_lon=[]
    vertice_clavesup_lat=[]
    
    for i in range(nx):                       #Numero de filas de la matriz
        vertice_claveinf_lon.append([0]*1)    #Numero de columnas de la matriz, se inicializa con ceros
        vertice_claveinf_lat.append([0]*1)
        vertice_clavesup_lon.append([0]*1)    
        vertice_clavesup_lat.append([0]*1)

    #1.4 Se generan los vertices de todas las subfallas y se almacenan los vertices claves (todavia una geometria rectangular).
    k=0
    for i_ny in range(ny):
        fila=0
        for j_nx in range(nx):
                frac=Fraction((ny-i_ny),ny)        

                #Latitud,longitud inicial
                lataux,lonaux,n = vinc_pt(lat_central[k], lon_central[k], np.degrees(strike)+90, dy/2 )
                lat_ini,lon_ini,n = vinc_pt(lataux, lonaux, np.degrees(strike)+180, dx/2 )
        
                #Segundo vertice
                lataux,lonaux,n = vinc_pt(lat_central[k], lon_central[k], np.degrees(strike)+90, dy/2 )
                lat1,lon1,n = vinc_pt(lataux, lonaux, np.degrees(strike), dx/2 )
            
                #Tercer vertice
                lat2,lon2,aux2 = vinc_pt(lat1, lon1, np.degrees(strike)+270, frac*W*np.cos(dip) )
              
                #Cuarto vertice
                lat3,lon3,aux2 = vinc_pt(lat2, lon2, np.degrees(strike)+180, dx)
        
                #Quinto Vertice
                lat4 = lat_ini  
                lon4 = lon_ini

                vertices_lon[k][0]=lon_ini
                vertices_lon[k][1]=lon1
                vertices_lon[k][2]=lon2
                vertices_lon[k][3]=lon3
                vertices_lon[k][4]=lon4

                vertices_lat[k][0]=lat_ini
                vertices_lat[k][1]=lat1
                vertices_lat[k][2]=lat2
                vertices_lat[k][3]=lat3
                vertices_lat[k][4]=lat4
            
               #Vertices claves.
                if (i_ny==ny-1):  #si se esta en la ultima fila (junto a la fosa)
                    vertice_claveinf_lon[j_nx][0]=lon3
                    vertice_claveinf_lat[j_nx][0]=lat3
                    vertice_clavesup_lon[j_nx][0]=lon2
                    vertice_clavesup_lat[j_nx][0]=lat2

                k=k+1
        fila=fila+1;
        
        


    #1.5 A partir de los vertices claves se obtiene la distancia de esto con la fosa (modelo slab 1.0)
    dl=[]
 
    for i_nx in range(nx):
        #Se ingresa a una subfalla y se busca la maxima longitud proyectrada a la fosa
        minlon= vertice_claveinf_lon[i_nx][0]
        for j in range(len(fosa_lat)):
            if (fosa_lat[j]>=vertice_claveinf_lat[i_nx][0]+delta_lat and fosa_lat[j]<=vertice_clavesup_lat[i_nx][0]+delta_lat):
                if fosa_lon[j]<minlon:
                    minlon=fosa_lon[j]
        #Se calcula distancia a desplazar para todos los puntos centrales.
        #Se esta probando delta LAT, TEST
        L, aux, aux = vinc_dist(  vertice_claveinf_lat[i_nx][0],  vertice_claveinf_lon[i_nx][0],   vertice_claveinf_lat[i_nx][0],  minlon ) 
        dl.append(L)


    #1.6 Se obtienen posiciones finales de las subfallas asociadas a una geometria del slab model.                
            
    lon_central_new=[]
    lat_central_new=[]

    k=0
    for i_ny in range(ny):
        ll=0
        for i_nx in range(nx):
            dist=dl[ll]
            lat_aux,  lon_aux,  alpha21  = vinc_pt( lat_central[k], lon_central[k], np.degrees(strike)+270,dist) 
            lon_central_new.append(lon_aux)
            lat_central_new.append(lat_aux)               
            k=k+1; 
            ll=ll+1
    
    #1.7 Se reesciben los vertices desde una distribucion rectangular a una adaptada al slab        
    k=0
    for i_ny in range(ny):
        fila=0
        for j_nx in range(nx):
                #frac=Fraction((ny-i_ny),ny)   
                frac_W=Fraction(W,ny)

                #Latitud,longitud inicial
                lataux,lonaux,n = vinc_pt(lat_central_new[k], lon_central_new[k], np.degrees(strike)+90, dy/2 )
                lat_ini,lon_ini,n = vinc_pt(lataux, lonaux, np.degrees(strike)+180, dx/2 )
        
                #Segundo vertice
                lataux,lonaux,n = vinc_pt(lat_central_new[k], lon_central_new[k], np.degrees(strike)+90, dy/2 )
                lat1,lon1,n = vinc_pt(lataux, lonaux, np.degrees(strike), dx/2 )
            
                #Tercer vertice
                lat2,lon2,aux2 = vinc_pt(lat1, lon1, np.degrees(strike)+270, frac_W*np.cos(dip) )
              
                #Cuarto vertice
                lat3,lon3,aux2 = vinc_pt(lat2, lon2, np.degrees(strike)+180, dx)
        
                #Quinto Vertice
                lat4 = lat_ini  
                lon4 = lon_ini

                vertices_lon[k][0]=lon_ini
                vertices_lon[k][1]=lon1
                vertices_lon[k][2]=lon2
                vertices_lon[k][3]=lon3
                vertices_lon[k][4]=lon4

                vertices_lat[k][0]=lat_ini
                vertices_lat[k][1]=lat1
                vertices_lat[k][2]=lat2
                vertices_lat[k][3]=lat3
                vertices_lat[k][4]=lat4
        
                k=k+1
        fila=fila+1;    
        
    return lon_central_new, lat_central_new,vertices_lon, vertices_lat

###############################################################################

def coordenadas_subfallas_EW(ny,nx,dip,W,L,lon0,lat0,strike,fosa_lon,fosa_lat,delta_lat):
    #La funcion obtiene las coordenadas centrales de las subfallas de acuerdo a la
    #configuracion del slab model, junto con las coordendas de los vertices de esta.    
    
    #El problema comienza generando una falla rectangular dividida en ny*nx subfallas,
    #para posteriormente desplazarlas en direccion hacia la fosa y obtener asi una configuracion
    #de subfallas de acuerdo al slab model.    
    
    #1. Se determina latitud y longitud de cada subfalla asumiendo una destribucion rectangular.
    lon_sl, lat_sl = subfalla(dip,W,L,nx,ny,lon0,lat0,strike) 
    dx = L / nx; dy = W * cos(dip) / ny
    
    #Variables para almacenar la posicion de central de las subfallas de acuerso al Slab Model..
    lon_central=[]
    lat_central=[]

    #1.1 Se obtiene posicion central de las subfallas de acuerdo a la geometria rectangular.
    k=0
    for i_nx in range(nx):
        for i_ny in range(ny):
            lon_central.append(lon_sl[k])
            lat_central.append(lat_sl[k])               
            k=k+1; 
        
    #1.2 Se generan vertices de las subfallas de acuerdo a la geometria rectangular.
    FIL=nx*ny  #Se evaluaran todas las subfallas.
    COL=5      #Cada subfalla tiene cinco vertices a dibujar (el ultimo coincide con el primero)

    #Se decalran variables para almacenar vertices de todas las subfallas.
    vertices_lon = [] 
    vertices_lat = []
    
    #Se crean matrices de ceros para rellenar con los vertices de todas las subfallas.
    for i in range(FIL):                #Numero de filas de la matriz
        vertices_lon.append([0]*COL)    #Numero de columnas de la matriz, se inicializa con ceros
        vertices_lat.append([0]*COL)
    
    #1.3 Se declaran vertices claves: Son los vertices superior e inferior de las subfallas que estan
    #Junto a la fosa. Sirven para buscar los valores de referencia con el modelo slab1.0    
    vertice_claveinf_lon=[]
    vertice_claveinf_lat=[]
    vertice_clavesup_lon=[]
    vertice_clavesup_lat=[]
    
    for i in range(nx):                       #Numero de filas de la matriz
        vertice_claveinf_lon.append([0]*1)    #Numero de columnas de la matriz, se inicializa con ceros
        vertice_claveinf_lat.append([0]*1)
        vertice_clavesup_lon.append([0]*1)    
        vertice_clavesup_lat.append([0]*1)

    #1.4 Se generan los vertices de todas las subfallas y se almacenan los vertices claves (todavia una geometria rectangular).
    k=0
    for i_ny in range(ny):
        fila=0
        for j_nx in range(nx):
                frac=Fraction((ny-i_ny),ny)        

                #Latitud,longitud inicial
                lataux,lonaux,n = vinc_pt(lat_central[k], lon_central[k], np.degrees(strike)+90, dy/2 )
                lat_ini,lon_ini,n = vinc_pt(lataux, lonaux, np.degrees(strike)+180, dx/2 )
        
                #Segundo vertice
                lataux,lonaux,n = vinc_pt(lat_central[k], lon_central[k], np.degrees(strike)+90, dy/2 )
                lat1,lon1,n = vinc_pt(lataux, lonaux, np.degrees(strike), dx/2 )
            
                #Tercer vertice
                lat2,lon2,aux2 = vinc_pt(lat1, lon1, np.degrees(strike)+270, frac*W*np.cos(dip) )
              
                #Cuarto vertice
                lat3,lon3,aux2 = vinc_pt(lat2, lon2, np.degrees(strike)+180, dx)
        
                #Quinto Vertice
                lat4 = lat_ini  
                lon4 = lon_ini

                vertices_lon[k][0]=lon_ini
                vertices_lon[k][1]=lon1
                vertices_lon[k][2]=lon2
                vertices_lon[k][3]=lon3
                vertices_lon[k][4]=lon4

                vertices_lat[k][0]=lat_ini
                vertices_lat[k][1]=lat1
                vertices_lat[k][2]=lat2
                vertices_lat[k][3]=lat3
                vertices_lat[k][4]=lat4
            
               #Vertices claves.
                if (i_ny==ny-1):  #si se esta en la ultima fila (junto a la fosa)
                    vertice_claveinf_lon[j_nx][0]=lon3
                    vertice_claveinf_lat[j_nx][0]=lat3
                    vertice_clavesup_lon[j_nx][0]=lon2
                    vertice_clavesup_lat[j_nx][0]=lat2

                k=k+1
        fila=fila+1;

    #1.5 A partir de los vertices claves se obtiene la distancia de esto con la fosa (modelo slab 1.0)
    dl=[]
 
    for i_nx in range(nx):
        #Se ingresa a una subfalla y se busca la maxima longitud proyectrada a la fosa
        minlon= vertice_claveinf_lon[i_nx][0]
        for j in range(len(fosa_lat)):
            if (fosa_lat[j]<=vertice_claveinf_lat[i_nx][0]+delta_lat and fosa_lat[j]>=vertice_clavesup_lat[i_nx][0]+delta_lat):
                if fosa_lon[j]>minlon:
                    minlon=fosa_lon[j]
        #Se calcula distancia a desplazar para todos los puntos centrales.
        L, aux, aux = vinc_dist(  vertice_claveinf_lat[i_nx][0],  vertice_claveinf_lon[i_nx][0],   vertice_claveinf_lat[i_nx][0],  minlon ) 
        dl.append(L)


    #1.6 Se obtienen posiciones finales de las subfallas asociadas a una geometria del slab model.                
            
    lon_central_new=[]
    lat_central_new=[]

    k=0
    for i_ny in range(ny):
        ll=0
        for i_nx in range(nx):
            dist=dl[ll]
            lat_aux,  lon_aux,  alpha21  = vinc_pt( lat_central[k], lon_central[k], np.degrees(strike)-90,dist) 
            lon_central_new.append(lon_aux)
            lat_central_new.append(lat_aux)               
            k=k+1; 
            ll=ll+1
    
    #1.7 Se reesciben los vertices desde una distribucion rectangular a una adaptada al slab        
    k=0
    for i_ny in range(ny):
        fila=0
        for j_nx in range(nx):
                #frac=Fraction((ny-i_ny),ny)   
                frac_W=Fraction(W,ny)

                #Latitud,longitud inicial
                lataux,lonaux,n = vinc_pt(lat_central_new[k], lon_central_new[k], np.degrees(strike)+90, dy/2 )
                lat_ini,lon_ini,n = vinc_pt(lataux, lonaux, np.degrees(strike)+180, dx/2 )
        
                #Segundo vertice
                lataux,lonaux,n = vinc_pt(lat_central_new[k], lon_central_new[k], np.degrees(strike)+90, dy/2 )
                lat1,lon1,n = vinc_pt(lataux, lonaux, np.degrees(strike), dx/2 )
            
                #Tercer vertice
                lat2,lon2,aux2 = vinc_pt(lat1, lon1, np.degrees(strike)+270, frac_W*np.cos(dip) )
              
                #Cuarto vertice
                lat3,lon3,aux2 = vinc_pt(lat2, lon2, np.degrees(strike)+180, dx)
        
                #Quinto Vertice
                lat4 = lat_ini  
                lon4 = lon_ini

                vertices_lon[k][0]=lon_ini
                vertices_lon[k][1]=lon1
                vertices_lon[k][2]=lon2
                vertices_lon[k][3]=lon3
                vertices_lon[k][4]=lon4

                vertices_lat[k][0]=lat_ini
                vertices_lat[k][1]=lat1
                vertices_lat[k][2]=lat2
                vertices_lat[k][3]=lat3
                vertices_lat[k][4]=lat4
        
                k=k+1
        fila=fila+1;    

    return lon_central_new, lat_central_new,vertices_lon, vertices_lat


################################################################################

     



#
# --------------------------------------------------------------------- 
# |                                                                     |
# |	geodetic.cc -  a collection of geodetic functions                   |
# |	Jim Leven  - Dec 99                                                 |
# |                                                                     |
# | originally from:                                                    |
# | http://wegener.mechanik.tu-darmstadt.de/GMT-Help/Archiv/att-8710/Geodetic_py |                                                                   |
# |                                                                     |
# --------------------------------------------------------------------- 
# 
# 
# ----------------------------------------------------------------------
# | Algrothims from Geocentric Datum of Australia Technical Manual	    |
# | 								                                    |
# | http://www.anzlic.org.au/icsm/gdatum/chapter4.html	        		|
# | 								                                    |
# | This page last updated 11 May 1999 	                				|
# | 								                                    |
# | Computations on the Ellipsoid	                    				|
# | 								                                    |
# | There are a number of formulae that are available           		|
# | to calculate accurate geodetic positions, 		            		|
# | azimuths and distances on the ellipsoid.			                |
# | 								                                    |
# | Vincenty's formulae (Vincenty, 1975) may be used 		            |
# | for lines ranging from a few cm to nearly 20,000 km, 	            |
# | with millimetre accuracy. 					                        |
# | The formulae have been extensively tested 		                    |
# | for the Australian region, by comparison with results       		|
# | from other formulae (Rainsford, 1955 & Sodano, 1965). 	            |
# |								                                        |
# | * Inverse problem: azimuth and distance from known 	        		|
# |			latitudes and longitudes 			                        |
# | * Direct problem: Latitude and longitude from known 	            |
# |			position, azimuth and distance. 		                    |
# | * Sample data 						                                |
# | * Excel spreadsheet 			                            		|
# | 								                                    |
# | Vincenty's Inverse formulae				                    		|
# | Given: latitude and longitude of two points                 		|
# |			(phi1, lembda1 and phi2, lembda2), 	                        |
# | Calculate: the ellipsoidal distance (s) and 	            		|
# | forward and reverse azimuths between the points (alpha12, alpha21).	|
# |									                                    |
# ---------------------------------------------------------------------- 

def vinc_dist(  phi1,  lembda1,  phi2,  lembda2 ) :
        """ 
        Returns the distance between two geographic points on the ellipsoid
        and the forward and reverse azimuths between these points.
        lats, longs and azimuths are in decimal degrees, distance in metres 

        Returns ( s, alpha12,  alpha21 ) as a tuple
        """

        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        if (abs( phi2 - phi1 ) < 1e-8) and ( abs( lembda2 - lembda1) < 1e-8 ) :
                return 0.0, 0.0, 0.0

        piD4   = math.atan( 1.0 )
        two_pi = piD4 * 8.0

        phi1    = phi1 * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0		# unfortunately lambda is a key word!
        phi2    = phi2 * piD4 / 45.0
        lembda2 = lembda2 * piD4 / 45.0

        b = a * (1.0 - f)

        TanU1 = (1-f) * tan( phi1 )
        TanU2 = (1-f) * tan( phi2 )

        U1 = atan(TanU1)
        U2 = atan(TanU2)

        lembda = lembda2 - lembda1
        last_lembda = -4000000.0		# an impossibe value
        omega = lembda

        # Iterate the following equations, 
        #  until there is no significant change in lembda 

        while ( last_lembda < -3000000.0 or lembda != 0 and abs( (last_lembda - lembda)/lembda) > 1.0e-9 ) :

                sqr_sin_sigma = pow( cos(U2) * sin(lembda), 2) + \
                        pow( (cos(U1) * sin(U2) - \
                        sin(U1) *  cos(U2) * cos(lembda) ), 2 )

                Sin_sigma = sqrt( sqr_sin_sigma )

                Cos_sigma = sin(U1) * sin(U2) + cos(U1) * cos(U2) * cos(lembda)
        
                sigma = atan2( Sin_sigma, Cos_sigma )

                Sin_alpha = cos(U1) * cos(U2) * sin(lembda) / sin(sigma)
                alpha = asin( Sin_alpha )

                Cos2sigma_m = cos(sigma) - (2 * sin(U1) * sin(U2) / pow(cos(alpha), 2) )

                C = (f/16) * pow(cos(alpha), 2) * (4 + f * (4 - 3 * pow(cos(alpha), 2)))

                last_lembda = lembda

                lembda = omega + (1-C) * f * sin(alpha) * (sigma + C * sin(sigma) * \
                        (Cos2sigma_m + C * cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2) )))

        u2 = pow(cos(alpha),2) * (a*a-b*b) / (b*b)

        A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))

        B = (u2/1024) * (256 + u2 * (-128+ u2 * (74 - 47 * u2)))

        delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B/4) * \
                (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2) ) - \
                (B/6) * Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) * \
                (-3 + 4 * pow(Cos2sigma_m,2 ) )))

        s = b * A * (sigma - delta_sigma)

        alpha12 = atan2( (cos(U2) * sin(lembda)), \
                (cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(lembda)))

        alpha21 = atan2( (cos(U1) * sin(lembda)), \
                (-sin(U1) * cos(U2) + cos(U1) * sin(U2) * cos(lembda)))

        if ( alpha12 < 0.0 ) : 
                alpha12 =  alpha12 + two_pi
        if ( alpha12 > two_pi ) : 
                alpha12 = alpha12 - two_pi

        alpha21 = alpha21 + two_pi / 2.0
        if ( alpha21 < 0.0 ) : 
                alpha21 = alpha21 + two_pi
        if ( alpha21 > two_pi ) : 
                alpha21 = alpha21 - two_pi

        alpha12    = alpha12    * 45.0 / piD4
        alpha21    = alpha21    * 45.0 / piD4
        return s, alpha12,  alpha21 

   # END of Vincenty's Inverse formulae 


#-------------------------------------------------------------------------------
# Vincenty's Direct formulae							|
# Given: latitude and longitude of a point (phi1, lembda1) and 			|
# the geodetic azimuth (alpha12) 						|
# and ellipsoidal distance in metres (s) to a second point,			|
# 										|
# Calculate: the latitude and longitude of the second point (phi2, lembda2) 	|
# and the reverse azimuth (alpha21).						|
# 										|
#-------------------------------------------------------------------------------

def  vinc_pt( phi1, lembda1, alpha12, s ) :
        """

        Returns the lat and long of projected point and reverse azimuth
        given a reference point and a distance and azimuth to project.
        lats, longs and azimuths are passed in decimal degrees

        Returns ( phi2,  lambda2,  alpha21 ) as a tuple 

        """
 
        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        piD4 = atan( 1.0 )
        two_pi = piD4 * 8.0

        phi1    = phi1    * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0
        alpha12 = alpha12 * piD4 / 45.0
        if ( alpha12 < 0.0 ) : 
                alpha12 = alpha12 + two_pi
        if ( alpha12 > two_pi ) : 
                alpha12 = alpha12 - two_pi

        b = a * (1.0 - f)

        TanU1 = (1-f) * tan(phi1)
        U1 = atan( TanU1 )
        sigma1 = atan2( TanU1, cos(alpha12) )
        Sinalpha = cos(U1) * sin(alpha12)
        cosalpha_sq = 1.0 - Sinalpha * Sinalpha

        u2 = cosalpha_sq * (a * a - b * b ) / (b * b)
        A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * \
                (320 - 175 * u2) ) )
        B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2) ) )

        # Starting with the approximation
        sigma = (s / (b * A))

        last_sigma = 2.0 * sigma + 2.0	# something impossible

        # Iterate the following three equations 
        #  until there is no significant change in sigma 

        # two_sigma_m , delta_sigma
        while ( abs( (last_sigma - sigma) / sigma) > 1.0e-9 ) :
                two_sigma_m = 2 * sigma1 + sigma

                delta_sigma = B * sin(sigma) * ( cos(two_sigma_m) \
                        + (B/4) * (cos(sigma) * \
                        (-1 + 2 * pow( cos(two_sigma_m), 2 ) -  \
                        (B/6) * cos(two_sigma_m) * \
                        (-3 + 4 * pow(sin(sigma), 2 )) *  \
                        (-3 + 4 * pow( cos (two_sigma_m), 2 ))))) \

                last_sigma = sigma
                sigma = (s / (b * A)) + delta_sigma

        phi2 = atan2 ( (sin(U1) * cos(sigma) + cos(U1) * sin(sigma) * cos(alpha12) ), \
                ((1-f) * sqrt( pow(Sinalpha, 2) +  \
                pow(sin(U1) * sin(sigma) - cos(U1) * cos(sigma) * cos(alpha12), 2))))

        lembda = atan2( (sin(sigma) * sin(alpha12 )), (cos(U1) * cos(sigma) -  \
                sin(U1) *  sin(sigma) * cos(alpha12)))

        C = (f/16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq ))

        omega = lembda - (1-C) * f * Sinalpha *  \
                (sigma + C * sin(sigma) * (cos(two_sigma_m) + \
                C * cos(sigma) * (-1 + 2 * pow(cos(two_sigma_m),2) )))

        lembda2 = lembda1 + omega

        alpha21 = atan2 ( Sinalpha, (-sin(U1) * sin(sigma) +  \
                cos(U1) * cos(sigma) * cos(alpha12)))

        alpha21 = alpha21 + two_pi / 2.0
        if ( alpha21 < 0.0 ) :
                alpha21 = alpha21 + two_pi
        if ( alpha21 > two_pi ) :
                alpha21 = alpha21 - two_pi

        phi2       = phi2       * 45.0 / piD4
        lembda2    = lembda2    * 45.0 / piD4
        alpha21    = alpha21    * 45.0 / piD4

        return phi2,  lembda2,  alpha21 

  # END of Vincenty's Direct formulae

#--------------------------------------------------------------------------
# Notes: 
# 
# * "The inverse formulae may give no solution over a line 
# 	between two nearly antipodal points. This will occur when 
# 	lembda ... is greater than pi in absolute value". (Vincenty, 1975)
#  
# * In Vincenty (1975) L is used for the difference in longitude, 
# 	however for consistency with other formulae in this Manual, 
# 	omega is used here. 
# 
# * Variables specific to Vincenty's formulae are shown below, 
# 	others common throughout the manual are shown in the Glossary. 
# 
# 
# alpha = Azimuth of the geodesic at the equator
# U = Reduced latitude
# lembda = Difference in longitude on an auxiliary sphere (lembda1 & lembda2 
# 		are the geodetic longitudes of points 1 & 2)
# sigma = Angular distance on a sphere, from point 1 to point 2
# sigma1 = Angular distance on a sphere, from the equator to point 1
# sigma2 = Angular distance on a sphere, from the equator to point 2
# sigma_m = Angular distance on a sphere, from the equator to the 
# 		midpoint of the line from point 1 to point 2
# u, A, B, C = Internal variables
# 
# 

#*******************************************************************

def est_dist( phi1,  lembda1,  phi2,  lembda2 ) :
        """ 

        Returns an estimate of the distance between two geographic points
        This is a quick and dirty vinc_dist 
        which will generally estimate the distance to within 1%
        Returns distance in metres

        """

        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        piD4   = 0.785398163397 

        phi1    = phi1 * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0
        phi2    = phi2 * piD4 / 45.0
        lembda2 = lembda2 * piD4 / 45.0

        c = cos((phi2+phi1)/2.0)

        return sqrt( pow(fabs(phi2-phi1), 2) + \
                pow(fabs(lembda2-lembda1)*c, 2) ) * a * ( 1.0 - f + f * c )

   # END of rough estimate of the distance.

 

      
