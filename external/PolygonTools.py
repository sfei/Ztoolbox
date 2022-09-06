# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 10:02:48 2018

@author: zhenlinz
"""
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import numpy as np

def PolyRings(poly):
    """ This function finds the coordinates of inner and outer rings of the 
    supplied polygon obj
    To see the coordinates of inner or outer rings, do Outerring_coords[i]
    """
    Outerring_coords = poly.exterior.coords
    Innerring_coords = []
    Innerrings = poly.interiors
    for ring in Innerrings:
        Innerring_coords.append(ring.coords)
    return Outerring_coords, Innerring_coords

def PolyHollow(poly):
    """ quick algorithm to check if the polygon is hollow
    """
    if poly.exterior.length==poly.length:
        return False
    else:
        return True
    

def DissAggPoly(poly):
    """ dis-aggregate multi-polygons into multiple polygons
    """
    polyv = []
    if type(poly)==Polygon: # single polygon
        polyv.append(poly)
    else:
        for pi in poly:
            polyv.append(pi)
    return polyv

def FindNearbyPoly(poly1,poly_grid):
    """A quick algorithm that finds nearby polygons from poly_grid (geopandas 
    dataframe) that may be interecepting with poly1
    """
    si = poly_grid.sindex
    ni = list(si.intersection(poly1.bounds)) # finding nearby polygon indices
    return ni

def plotPolygon(poly):
    xy = poly.boundary.xy
    x = xy[0]
    y = xy[1]
    plt.plot(x,y)
    
def Polylen(poly):
    # find the number of polygons in the geometry
    if isinstance(poly,Polygon):
        poly_len = 1
    elif isinstance(poly, MultiPolygon):
        poly_len = len(poly)  
    return poly_len

def FindMultiPoly(sub_grid):
    # find the indices for multipolygons
    plens =[]
    for poly in sub_grid.geometry:
        poly_len = Polylen(poly)
        plens.append(poly_len)
        
    ind = np.where(np.array(plens)>1)[0]
    return ind
            