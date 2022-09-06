# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 11:47:27 2017
Check if the input nutrient sources were assigned correctly to the boundaries
@author: zhenlinz
"""
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import pandas as pd
import logging

def CreateMapfromFile(gridfile,datafile,varname='FlowElem_bl',layer=0,time=-1,plotnan=True):
    """Create spatial map of input varialbe (varname) 
    If no varname is supplied, the default varname is bottom level ('FlowElem_bl')
    The default layer is surface layer=0 and the default timestep is the last one.
    """


    mapfile = xr.open_dataset(datafile)    
    
    try:
        var = mapfile[varname].isel(time=time,layer=layer).values
    except NameError:
        var = mapfile[varname].values   
        print("The plotted variable is not time dependent")
    
    coll = CreateMapFromData(gridfile,var,plotnan)    

    return coll


def CreateMapFromData(gridfile,var,plotnan=False,color=False,fig=None, 
                      ax=None,colorbar='on',region=None, **kwargs):
    """Create spatial map of var from data directly
    var needs to be 1D
    plotnan = True: plot nan values as the minimum; False: do not show nan valu polygons
    color=True: pass on color tuples to plot instead of array values
    color tuples could be 3 or 4 numbers, the first 3 corresponds to RGB color and
    the last one means alpha value
    fig: figure name; when none is passed, generate a new one
    region: None means the entire region; otherwise a list of indices are supplied. 
    """    
    bath = xr.open_dataset(gridfile)    
    xc = bath['FlowElemContour_x'].values
    yc = bath['FlowElemContour_y'].values
    
    
    if region is not None:
        xc = xc[region]
        yc = yc[region]
        var = var[region]
    
    nodes = np.asanyarray([xc,yc]).transpose([1,2,0])
    ind = np.where(nodes<=0)
    if len(ind)>0:
        nodes[ind] =np.nan # assign nan values to the empty points
    
    if plotnan == False:
        ind = np.where(np.isnan(var))
        if len(ind)>0:
            nodes[ind,:,:] = np.nan
    
    if not ax:        
        if not fig:    
            fig, ax = plt.subplots()        
        ax = fig.gca()

    
    if color!=False:
        coll = PolyCollection(nodes,facecolors=color,edgecolors=color,**kwargs)
    else:
        coll = PolyCollection(nodes,array=var,edgecolors="face",**kwargs)
    col = ax.add_collection(coll)
    ax.axis('equal') # This is important so that the image will be centered. 
    #coll.set_clim(0,0.1)
    if not color:
        if colorbar=='on':
            plt.colorbar(col)
        
    return coll    

def PlotGrid(gridfile,fig=None,plotz=True,**kwargs):
    bath = xr.open_dataset(gridfile)    
    xc = bath['FlowElemContour_x'].values
    yc = bath['FlowElemContour_y'].values
    nodes = np.asanyarray([xc,yc]).transpose([1,2,0])
    nodes[np.where(nodes<=0)] =np.nan # assign nan values to the empty points
    
    if plotz:
        var = bath['FlowElem_bl'].values
        coll = PolyCollection(nodes,array=var,**kwargs)
    else:
        coll = PolyCollection(nodes,**kwargs)

    if not fig:    
        fig, ax = plt.subplots()        
    ax = fig.gca()

        
    ax.add_collection(coll)
    ax.axis('equal') # This is important so that the image will be centered. 
    #coll.set_clim(0,0.1)

        
    return coll        


def Poly3DCollection_color(gridfile,var, ax=None, fig=None,cvar=None,
                           inpoly=None,cmap='jet',LogNorm=False, 
                           minc=None,maxc=None,**kwargs):
    """ var is the z level to be plotted in the 3D polygon
        cvar is the color to be plotted in the 3D polygon
    """
    
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib import colors

    if ax is None:        
        if fig is None:    
            fig, ax = plt.subplots()        
        ax = fig.gca(projection='3d')       
    
    if cvar is None:
        cvar = var 
        
        
    xe,ye,ze,f = Interp3d(gridfile,var,inpoly)
    verts = np.asarray([xe,ye,ze])
    verts = verts.transpose([1,2,0])
    col = Poly3DCollection(verts, **kwargs)    
    
    # set colors
    if not minc:
            minc = min(cvar)
    
    if not maxc:
            maxc = min(cvar)
        
    if LogNorm:
        norm = colors.LogNorm(vmin=minc,vmax=maxc)
    else:
        norm = colors.Normalize(vmin=minc,vmax=maxc)
        
    plt.set_cmap(cmap)
    cm1 = plt.get_cmap()
    colorvar = cm1(norm(cvar))
    col.set_facecolor(colorvar)
    col.set_edgecolor(colorvar)

    try:
        pla = ax.add_collection3d(col)
    except AttributeError as err:
        print("create input ax by ax=fig.gca(projection='3d')" )
        
    ax.axis('equal')
    ax.set_xlim([np.nanmin(xe),np.nanmax(xe)])
    ax.set_ylim([np.nanmin(ye),np.nanmax(ye)])
    ax.set_zlim([np.nanmin(ze),np.nanmax(ze)])
    
    return col
    
def Interp3d(gridfile,var,inpoly=None):
    from scipy import interpolate
    g = dwaqGrid(gridfile,inpoly=inpoly)
    xc = g.xc
    yc = g.yc
    xe = g.xe
    ye = g.ye
    
#    gdata = xr.open_dataset(gridfile)
#    xe =  gdata['FlowElemContour_x'].values
#    ye = gdata['FlowElemContour_y'].values
#    xe[np.where(xe<=0)]=np.nan
#    ye[np.where(ye<=0)]=np.nan
#    
#    xe = np.squeeze(xe[inpoly,:])
#    ye = np.squeeze(ye[inpoly,:])
#    
#    xc = gdata['FlowElem_xcc'].values
#    yc = gdata['FlowElem_ycc'].values
#    xc = np.squeeze(xc[inpoly])
#    yc = np.squeeze(yc[inpoly])
        
    x = np.r_[xc[None,:],yc[None,:]]
    f = interpolate.NearestNDInterpolator(x.T,var)
    
    # 3D interpolation from cell centers to nodes. 
    
    ze = []
    for i in range(np.shape(xe)[1]):
        #z = [f(xe[j,i],ye[j,i]) for j in range(np.shape(xe)[0])]
        z = np.nan*np.ones(np.shape(xe[:,i]))
        
        z[~np.isnan(xe[:,i])] = f(xe[~np.isnan(xe[:,i]),i],ye[~np.isnan(xe[:,i]),i])
        #z =  [griddata(xc,yc,depth, xe[j,i],ye[j,i]) for j in range(np.shape(xe)[0])]
        ze.append(z)
    
    ze = np.asarray(ze).T

    return xe, ye, ze,f  

class dwaqGrid(object):
    
    def __init__(self, gridfile,inpoly=None):
        self.gridfile = gridfile
        self.gdata = xr.open_dataset(gridfile)
        self.inpoly = inpoly

    def readgrid(self,var,nonzero=True):
        x = self.gdata[var].values
        if nonzero:
            x[np.where(x<=0)]=np.nan
        x = np.squeeze(x[self.inpoly])
        return x        
    
    @property
    def xe(self): 
        return self.readgrid('FlowElemContour_x')
    
    @property
    def ye(self):
        return self.readgrid('FlowElemContour_y')  

    @property
    def xc(self):
        return self.readgrid('FlowElem_xcc')
    
    @property
    def yc(self):
        return self.readgrid('FlowElem_ycc')    
     
    @property        
    def bl(self):
        return self.readgrid('FlowElem_bl',nonzero=False)
    
    def LocationIndex(self,utm_x,utm_y):
        # Find the corresponding grid index at input utm_x,utm_y
        xc = self.xc
        yc = self.yc
        dist = (xc-utm_x)**2 + (yc-utm_y)**2
        ind = np.where(dist==min(dist))[0][0]
        return ind
    
    def boundary(self):
        # merge all cells to generate one polygon
        from shapely.geometry import Polygon
        from shapely.ops import cascaded_union
        Polyi = []
        xe = self.xe
        ye = self.ye
        for xei,yei in zip(xe,ye):
            indv = np.where(~np.isnan(xei))[0]
            if len(indv)>0:        
                Polyi.append( Polygon( (np.asarray([xei[indv],yei[indv]])).T ) )
            else:
                Polyi.append( Polygon( (np.asarray([xei,yei])).T) )        
        Polyunion = cascaded_union(Polyi)  # merge polygons to create the boundary of the sectioned grid
        return Polyunion
    
    def toGeopandasPoly(self,createshpfile=False):
        # create geopandas polygon based on grid cell boundary coordinates
        from shapely.geometry import Polygon
        import geopandas as gpd
        xe = self.xe
        ye = self.ye
        poly = []      
        ncells = self.gdata.nFlowElem.values
        for i in ncells:
            polyi_x = xe[i]
            polyi_y = ye[i]
            polyi = Polygon([px,py] for px, py in zip(polyi_x,polyi_y) if (px>=0) & (py>=0))
            poly.append(polyi)
            
        df = pd.DataFrame()
        df['geometry'] = poly
        gdf = gpd.GeoDataFrame(df,geometry='geometry')
        gdf.crs = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"
        if createshpfile:            
            out_shp_fn = self.gridfile+'.shp'
            gdf.to_file(out_shp_fn)
            logging.info("%s generated"%out_shp_fn)            
        return gdf
    
    def cell2cell(self, celli):
        # find the cells connected with celli
        # the cell numbers are 1 based. 
        FlowLink= self.readgrid('FlowLink',nonzero=False)
        indi = np.where(FlowLink[:,1]==celli)[0] 
        celllink = FlowLink[indi,0]
        maxi = np.max(FlowLink[:,1])
        return celllink[celllink<=maxi]
        
        
    
def selectPointsInShpfile(gridfile,shpfile):
    """Find the indices of the grid within a shapefile polygon
    """
    import shapefile
    import warnings
    import utm
    from shapely.geometry import Polygon, Point

    bath = xr.open_dataset(gridfile)    
    node_x = bath.FlowElem_xcc
    node_y = bath.FlowElem_ycc
    
    shapeobj = shapefile.Reader(shpfile).shapes()
    
    if len(shapeobj)>1:
        warnings.warn('Only the first shape file object is used')
        
    poly = shapeobj[0].points

    polyutm = [utm.from_latlon(lat,lon)[0:2] for lon, lat in poly]
    polysh = Polygon(polyutm)
    inpoly = []
    for xi, yi in zip(node_x,node_y):
        pointi = Point((xi,yi))
        inpoly.append(int(pointi.within(polysh)))
    inpoly = np.asarray(inpoly)
    inpoly = np.where(inpoly==1)[0]
    return inpoly
    
    
    
if __name__ == "__main__":
    model_dir = 'cascade_nutrients_v01_003'
    gridfile = os.path.join(model_dir,"flowgeom.nc")
    datafile = os.path.join(model_dir,"dwaq_map.nc")
    
    coll = CreateMapfromFile(gridfile,datafile,"OXY",layer=-1,time=-1)
    coll.set_clim(0,8)


    
    