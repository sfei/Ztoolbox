# -*- coding: utf-8 -*-
"""
Created on Fri May  3 14:43:12 2019
Check the contiguity of the spatial cluster on dwaq grid and make corrections if not
@author: zhenlinz
"""

import dwaq.PostProcessing as dpp
import geopandas as gpd
import logging
import numpy as np

# high resolutino dwaq grid file
hr_grid_fn = r"C:\Users\zhenlinz\Google Drive\1_Nutrient_Share\1_Projects_NUTRIENTS\Modeling\NOTES_ZZ\Projects\ModelAggregation\SFBay\wy2013c_waqgeom.nc"  
# shapefile polygons
shp_fn = r"C:\Users\zhenlinz\Google Drive\1_Nutrient_Share\1_Projects_NUTRIENTS\Modeling\NOTES_ZZ\Projects\ModelAggregation\SFBay\Rusty_agg_mod\Agg_mod.shp" 
outshp_fn = r"C:\Users\zhenlinz\Google Drive\1_Nutrient_Share\1_Projects_NUTRIENTS\Modeling\NOTES_ZZ\Projects\ModelAggregation\SFBay\Rusty_agg_mod\Agg_mod_continuous.shp" 

grid = dpp.dwaqGrid(hr_grid_fn)
grid_gpd = grid.toGeopandasPoly()

poly_gpd = gpd.read_file(shp_fn)
poly_gpd.Index = poly_gpd.index

grid_cluster = gpd.sjoin(grid_gpd, poly_gpd, how="left", op='within')

if np.any(np.isnan(grid_cluster.Index.values)):
    contiguous=False
    logging.warning("unclassfied cells identified: correction needed")
    
indnan = np.where(np.isnan(grid_cluster.Index.values))[0]

cindex = []
for i in indnan:
    area = 0
    for pi in np.arange(len(poly_gpd)):
        poly_1 = grid_gpd.geometry[i]
        poly_2 = poly_gpd.geometry[pi]
        if poly_1.overlaps(poly_2):
            area_new = poly_1.intersection(poly_2).area
            if area_new>area:
                area=area_new
                ci = pi
            else:
                ci = np.nan
                logging.warning("cell %d does not intersect with any polygons"%i)
    cindex.append(ci)


# Need a better algorithm than this. 
# For those cells not within any polygons, they are classified as in a polygon if they touch it. 
grid_cluster.Index[indnan] = cindex
grid_cluster_new = grid_cluster.dissolve(by='Index')

grid_cluster_new.drop(['index_right','Depth'],axis=1,inplace=True)
grid_cluster_new.to_file(outshp_fn)
    

    

