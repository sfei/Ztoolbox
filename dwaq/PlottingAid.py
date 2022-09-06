# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 10:48:50 2018

@author: zhenlinz
"""
import matplotlib.pyplot as plt
import numpy as np
#import utm
from sys import platform



#%%
def ll2utm(lonlat):
    
    from shapely.geometry import Point
    from functools import partial
    import pyproj
    from shapely.ops import transform
    
    point1 = Point(lonlat[0],lonlat[1])

    project = partial(
    pyproj.transform,
    pyproj.Proj(init='epsg:4269'), #NAD 83
    pyproj.Proj(init='epsg:32610')) #NAD83 / UTM zone 10N

    point2 = transform(project, point1)
    
    utm_x = point2.xy[0][0]
    utm_y = point2.xy[1][0]
    
    return np.asarray([utm_x,utm_y])

#%% 

def PlotGridBoundaries(ax=None,domain='Delta',**kwargs):
    """
    It is now only working for Suisun Bay bathymetry, but can certianly be extended to other domain
    """
    
    import geopandas
    
    if platform=='win32':         
        if domain=='Delta':
            geodf = geopandas.read_file('C:/Users/zhenlinz/Google Drive/1_Nutrient_Share/1_Projects_NUTRIENTS/Modeling/NOTES_ZZ/SuisunBay/Data/DeltaSubregions_v6/bathymetry_boundary_dissolve.shp')        
    else:
        if domain=='Delta':
            geodf = geopandas.read_file('Y:/zhenlin/dwaq/cascade/Data/bathymetry_boundary_dissolve.shp')
    
    # plot the outline of the bathymetry
    boundary = geodf['geometry'].boundary[0]
    
    
    for i in range(len(boundary)):
        boundary_xy = np.asarray( list(boundary[i].coords)).T
        utmvalue = [ll2utm([x,y]) for x, y in zip(boundary_xy[0,:],boundary_xy[1,:])]
        utmvalue = np.asarray(utmvalue).T
        if ax:
            ax.plot(utmvalue[0,:].astype(float),utmvalue[1,:].astype(float),**kwargs,zorder=-1)  #need to change to float, otherwise the figure messes up
        else:
            plt.plot(utmvalue[0,:].astype(float),utmvalue[1,:].astype(float),**kwargs,zorder=-1) 