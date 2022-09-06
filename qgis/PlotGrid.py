# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 16:12:31 2017

@author: zhenlinz
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

fn='C:\Workspace\Projects\SuisunBay\ModelData\sal_temp_waqgeom.nc'
fn=r'C:\Users\zhenlinz\Google Drive\1_Nutrient_Share\1_Projects_NUTRIENTS\Modeling\NOTES_ZZ\Projects\ModelAggregation\SFBay\flowgeom_Rusty_agg.nc'
dnc = xr.open_dataset(fn)

xc = dnc['FlowElemContour_x'].values
yc = dnc['FlowElemContour_y'].values
FlowElem_bl = dnc['FlowElem_bl'].values*-1 #bl here means bottom level. 

nodes = np.asanyarray([xc,yc]).transpose([1,2,0])
nodes[np.where(nodes<0)] =np.nan # assign nan values to the empty points

# First make grid plot by matplotlib. 
ax = plt.gca()
coll = PolyCollection(nodes,array=FlowElem_bl)
col = ax.add_collection(coll)
ax.axis('equal') # This is important so that the image will be centered. 



#%% The second half the script needs to link with qgis

def CreateQgisGrid(polys,depth):
    """
    This function creates a bathymetry plot in Qgis for San Francisco Bay area
    """
    # Create a new polygon memory layer
    bathLayer = QgsVectorLayer("Polygon?crs=EPSG:32610", "bathymetry","memory") #32610 corresponds to WGS4-UTMzone 10N
     #lineLayer = QgsVectorLayer("LineString", 'test layer', "memory")
     
    pr = bathLayer.dataProvider()
     
      # Add attributes fieldname and type
    res = pr.addAttributes([QgsField("Index", QVariant.Int), 
    QgsField("Depth", QVariant.Double)])
    
    bathLayer.updateFields() #Now updating the layer fields. 
        
    for i in range(len(polys)):
        polysi = polys[i,:,:]
        
            #ax = plt.gca()
            #coll = PolyCollection(polys)
            #col = ax.add_collection(coll)
            #plt.show()
            #ax.axis('equal') 
              
        
        
            #create layer
            
            
        
       
        
        # create a feature
        
        poly = QgsFeature()
        points = [QgsPoint(polysi[j,0],polysi[j,1]) for j in range(len(polysi)) if ~np.isnan(polysi[j,0])]
        poly.setGeometry(QgsGeometry.fromPolygon([points]))        
        di = float(depth[i]) # need to convert to python float first, otherwise it won't work        
        poly. setAttributes([i, di]) #set attribute index for poly. 
        
        pr.addFeatures([poly])
        bathLayer.updateExtents()
    # register the layer
    QgsMapLayerRegistry.instance().addMapLayers([bathLayer])


if __name__ =="__console__":
    import qgis
    from PyQt4.QtCore import QVariant
    CreateQgisGrid(nodes,FlowElem_bl)

