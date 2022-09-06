# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 09:39:41 2018
Perform grid aggregation with supporting scripts
@author: zhenlinz
"""

import geopandas as gpd
import xarray as xr 
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import dwaq.PostProcessing as dpp
from shapely.geometry import Polygon
import pandas as pd
import copy


workpath = r'C:\Users\zhenlinz\Google Drive\1_Nutrient_Share\1_Projects_NUTRIENTS\Modeling\NOTES_ZZ\Projects\ModelAggregation\SFBay'
hgridFile = 'wy2013c_waqgeom.nc' # DFM or DWAQ grid file. 
edge_weighted = False    # if we want to assign edge weight

if edge_weighted:
    mapfile = r'Y:\zhenlin\dwaq\cascade\suisun_nutrient_cycling\ModelResults\dwaq_map6.nc'
    dset= xr.open_dataset(mapfile)
    varname = ['NO3']   # or None


#%% load the grid file and find edge to cell
# Note that I should not use NetLink and nodes supplied by the grid file; for
# some reason, they may not match those in the polygons. So start only from the polygons
nc = xr.open_dataset(os.path.join(workpath,hgridFile))
poly_x = nc.FlowElemContour_x.values
poly_y = nc.FlowElemContour_y.values
cell_x = nc.FlowElem_xcc.values
cell_y = nc.FlowElem_ycc.values
icells = nc.nFlowElem.values

#%%  First find all the nodes
def nodes():
    """Function used to find all the nodes that connecct the polygons
    """
    node_xy = []
    for i in icells:
        poly_nodes = [ (x, y) for x,y in zip(poly_x[i],poly_y[i]) if (x>=0) & (y>=0)]            
        [node_xy.append((nx,ny)) for nx, ny in poly_nodes]
    node_xy = np.asarray(node_xy)
    node_xy = np.unique(node_xy,axis=0)
    return node_xy

# Then find the corresponding node indices (this may take a minute to run)
def cell_to_nodes(node_xy):
    """ Function used to map cell to nodes
    This may take a minute or two. 
    """
    poly2nodes = np.ones_like(poly_x)*-99
    for i in icells:
        poly_nodes = [ (x, y) for x,y in zip(poly_x[i],poly_y[i]) if (x>=0) & (y>=0)] 
        for p, pos in enumerate(poly_nodes):
            ind  = np.argwhere( (pos[0]==node_xy[:,0]) & (pos[1]==node_xy[:,1]) )
            poly2nodes[i,p] = ind
    return poly2nodes

#
def cell_to_cell(poly2nodes):
    print("""map cells to other cells;
          This may take a couple of min to run. 
          """)
    cell2cell=[]
    for i in icells:
        polyi = poly2nodes[i,poly2nodes[i,:]>0]
        max_edges = len(polyi)
        poly_edges = [(polyi[ni-1],polyi[ni]) for ni in range(max_edges)]
        poly_edges = np.sort(poly_edges,axis=1) #sort the edge indices
        celli = []
        for p in np.arange(len(poly_edges)):   
            ind = np.argwhere(poly_edges[p][0]==poly2nodes)[:,0]
            ind = ind[ind!=i]
            cellip = [indi for indi in ind if poly_edges[p][1] in poly2nodes[indi]]
            if cellip:
                celli.append(cellip[0])
                if len(cellip)>1:
                    print("more than one cell share boundaries with current cell:"+str(i))
        cell2cell.append(celli)
    return cell2cell  

def createDepthWeighing(varname, cell2cell,scale=None):
    """ This script creates weighing based on depth output from DWAQ. 
    varname: depth variable name
    scale: normal or 'log'
    """
    var = dset[varname]
    var = -1*var.values #all depth are nagative 
    
    if scale:
        if scale == 'log':
            var = var + var.mean()/1e5 # adding a very small value to avoide log(0)
            var = np.log(var)    
    var_norm = (var-var.mean())/var.std() # standardization the variable using z-transform
            
    cellgrad = []
    for i in icells:
        ci = cell2cell[i]            
        ci = cell2cell[i]
        dist = np.sqrt( (cell_x[i]-cell_x[ci])**2 + \
                      (cell_y[i]-cell_y[ci])**2 ) 
        grad = np.abs(var_norm[i]-var_norm[ci])/dist  # this is dC/dx
        cellgrad.append(grad)   
    return cellgrad, var

def createScalarWeighing(varname, cell2cell,surface=True):
    """ This script creates weighing based on 4D scalar output from DWAQ. 
    varname: lisf of scalar variable name; for one scalar, list ['var']
    surface: use surface value
    """
    multi_cellgrad = []
    multi_var = []
    for a in varname: 
        var = dset[a]
    
        if surface:
            var = var.sel(layer=0).values
        else:
            var = var.mean(dim='layer').values 
    
        var_norm = (var-np.nanmean(var))/np.nanstd(var) # standardization the variable using z-transform
        
        cellgrad = []
        for i in icells:
            ci = cell2cell[i]            
            ci = cell2cell[i]
            dist = np.sqrt( (cell_x[i]-cell_x[ci])**2 + \
                          (cell_y[i]-cell_y[ci])**2 )  
            grad = np.sqrt(np.nanmean(np.square([ np.abs( (var[:,i]-var[:,c]) )/dist[cj] for cj, c in enumerate(ci) ]),axis=1))
            cellgrad.append(grad) 
        multi_cellgrad.append(cellgrad)
        multi_var.append(var)
    return multi_cellgrad,multi_var


def combineWeighing(Grad):
    """combine the gradients from different var to a single weighing for each edge
    """
    Grad_sum = []
    for i in icells:
        Gradi = np.stack(Grad[:,i])
        Gradi_sum = np.sqrt(np.square(Gradi).mean(axis=0))
        Grad_sum.append(Gradi_sum)
    return Grad_sum


def Calculate_R2(var,labels):
    """calculate R2 for each indiviual variable using method given by:
    http://pro.arcgis.com/en/pro-app/tool-reference/spatial-statistics/how-spatially-constrained-multivariate-clustering-works.htm
    """
    TSS = np.nansum( (var -np.nanmean(var))**2) 
    ESS = 0
    for i in np.arange(max(labels)+1):
        celli = np.argwhere(labels==i)
        cluster_mean = np.nanmean(var[celli])
        cluster_ESS = np.nansum( (var[celli]-cluster_mean)**2 )
        if cluster_mean:
            ESS += cluster_ESS
    R2 = (TSS-ESS)/TSS
    return R2,TSS, ESS
        
def Calculate_PseudoF(var_all,labels):
    """ calculate PseudoF statistics using the method given by
    http://pro.arcgis.com/en/pro-app/tool-reference/spatial-statistics/how-spatially-constrained-multivariate-clustering-works.htm
    Fstat: a ratio reflecting within-group similarity and between-group difference
    """
    SST = 0 
    SSE = 0 
    nc = max(labels)+1
    n = max(icells)+1
    for i in np.arange(np.shape(var_all)[0]):
        var = var_all[i]
        SST += np.nansum( (var -np.nanmean(var))**2) 
        for i in np.arange(max(labels)+1):
            celli = np.argwhere(labels==i)
            cluster_mean = np.nanmean(var[celli])
            cluster_SSE = np.nansum( (var[celli]-cluster_mean)**2 )
            if cluster_mean:
                SSE += cluster_SSE
    R2 = (SST-SSE)/SST
    
    Fstat = (R2/(nc-1))/((1-R2)/(n-nc)) 
    return R2, Fstat           
 
#%% calculating cell to cell information from only the polygon files (this can be slow)
node_xy = nodes()
poly2nodes = cell_to_nodes(node_xy)
cell2cell = cell_to_cell(poly2nodes)

#%% Find the number of shared edges
alledges = []
for i in icells:
    celli = cell2cell[i]
    [alledges.append([i, j]) for j in celli]
alledges = np.sort(alledges,axis=1)
alledges = np.unique(alledges,axis=0)
N_internal_edge = len(alledges)

#%% using the attributes from results file and generate weighing. 
if edge_weighted:    
    depthGrad,depthvar = createDepthWeighing('bedlevel',cell2cell,scale='log')
    depthGrad = np.array(depthGrad)[None,:]
    if varname: # extra atrributes    
        scalarGrad,scalar = createScalarWeighing(varname,cell2cell)
        Grad = np.concatenate((depthGrad,scalarGrad))
        Grad_sum = combineWeighing(Grad)
    else:
        Grad_sum = depthGrad[0]

#%% Getting the maximum gradient
        
#Grad_max = max([max(Gradi) for Gradi in Grad_sum])
#Grad_small = Grad_max*1e-5
        
Grad_sum_2series = np.array([])
for Gradi in Grad_sum:
    Grad_sum_2series = np.append( Grad_sum_2series,Gradi)        
Grad_max = np.percentile(Grad_sum_2series,95)
#Grad_small = Grad_max/100


#Grad_small = np.nanpercentile(Grad_sum_2series,5) # use the 5 percentile Gradient as the minimum    
#%% Finding the scaling factors: make sure that (max_Grad+B)/Z makes 100 and (median_Grad+B)/Z makes 50
#Z = (max_Grad-median_Grad)/50
#B = max_Grad-2*median_Grad

Grad_sum_scaled = copy.deepcopy(Grad_sum)
for i, Gradi in enumerate(Grad_sum):
    #Grad_c = (Gradi+B)/Z
    Grad_c = Grad_max/Gradi
    if np.any(np.isnan(Grad_c)):
        Grad_c[np.isnan(Grad_c)] = 1 
    if np.any(Grad_c<1):
        Grad_c[ Grad_c<1 ] = 1 #The weights have to start from one

    Grad_sum_scaled[i] = Grad_c.astype(int)  

#%% Plot the gradients
#from matplotlib import cm
#
#sb = np.argwhere(cell_y<4182717.9959078836)[:,0]
#for i in icells[sb]:
#    
#    print(i)
#    nbrs = cell2cell[i]
#    Gradi = Grad_sum_scaled[i]
#    for g,c in enumerate(nbrs):
#        ints = np.intersect1d(poly2nodes[i],poly2nodes[c]).astype(int)
#        ints = ints[ints>0]
#        color = cm.jet( int(Gradi[g]/(Grad_max/Grad_small)*255) )
#        dist = np.sqrt(np.diff(node_xy[ints][:,0])**2 + np.diff(node_xy[ints][:,1])**2)
#        plt.plot(node_xy[ints][:,0],node_xy[ints][:,1],color = color)     

#%% Generate input file for METIS

with open("grid.graph","wt") as fp:
    # This line could additionally including fmt and ncon, needed
    # for variable weights, and something about multiple constraints.
    if edge_weighted:
        fp.write("%d %d %s\n"%(icells[-1]+1, N_internal_edge,'001'))
    else:
        fp.write("%d %d\n"%(icells[-1]+1, N_internal_edge))

    if edge_weighted:
        for i in icells:
            nbrs=cell2cell[i]      
            Gradi = Grad_sum_scaled[i]
            nbrs_weighted = []
            for nbr,g in zip(nbrs,Gradi):
                nbrs_weighted.append(str(nbr+1))
                nbrs_weighted.append(str(g))         
            fp.write(" ".join(nbrs_weighted) + "\n")             
    else:
        for i in icells:
            nbrs=cell2cell[i]       
            nbrs=[str(nbr+1) for nbr in nbrs if nbr>=0]        
            fp.write(" ".join(nbrs) + "\n")     

  
#%% run metis code
import subprocess
subprocess.call("gpmetis -contig grid.graph 5000 >metis.out",
                shell=True)

#%% Read metis output file and make sure that there is no error or warning. 


#%% load the generated file and plot
labels=np.loadtxt("grid.graph.part.50",dtype=np.int32)

# randomly permute to get more color contrast:
n=labels.max()+1
perm=np.argsort(np.random.random(n))
## 
labels=perm[labels]

fig,ax = plt.subplots()
dpp.CreateMapFromData(hgridFile,labels,fig=fig,cmap='jet')
dpp.PlotGrid(hgridFile,fig=fig,plotz=False,facecolors='none',edgecolors='w',linewidth=0.05)

#%% calculate F statitiscs
  
if edge_weighted:
    # calculate rms for the scalars
    R2 = [Calculate_R2(depthvar,labels)]
    var_all = depthvar[None,:]
    if varname:
        for i in range(len(varname)):
            var_i=np.sqrt(np.nanmean(np.square(scalar[i]),axis=0))
            R2.append( Calculate_R2(var_i,labels)   )
            var_all = np.concatenate( (var_all,var_i[None,:]))
    R2_all, Fstat = Calculate_PseudoF(var_all,labels)
            

#%% generate shape file

poly = []
indinv= [] # keep track and remove the invalid polygons
for i in icells:
    polyi_x = poly_x[i]
    polyi_y = poly_y[i]
    if (np.any(polyi_x<1)) | (np.any(polyi_y<1)):
        indinv.append(i)
    polyi = Polygon([px,py] for px, py in zip(polyi_x,polyi_y) if (px>=0) & (py>=0))
    poly.append(polyi)
    
df = pd.DataFrame()
df['geometry'] = poly
gdf = gpd.GeoDataFrame(df,geometry='geometry')
gdf['CLUSTER_ID'] = labels 

gdf.geometry=gdf.buffer(0.1) # create a buffer to avoid small self-intersection
agg_grid = gdf.dissolve(by='CLUSTER_ID')
agg_grid.crs = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"
agg_grid.to_file("aggGrid100.shp")

#%%



def AddShpfileProperties(gridshp,outgridshp, datafile,varname,ts=1,pos='top'):
    """ 
    Not checked... 
    Add one or a list of variabiles from datafile to the 
    shapefile    
    gridshp: normally means the dwaq grid shape file
    outgridshp: output grid shape file
    datafile: normally means model results file
    varname: one or list of variables in the datafile;
    ts: time step
    pos: either 'top' or 'avg' over the depth. 
    """
        
    grid = gpd.read_file(gridshp)
    dset = xr.open_dataset(datafile)
    
    if pos == 'top':
        var =  dset[varname][ts,0,:].values
    elif pos == 'avg':
        var = dset[varname][ts,:,:].values.mean(axis='depth')
    else:
        raise ValueError(" 'pos' can either be 'top' or 'avg' ")
           
    grid[varname] = var
    grid.to_file(outgridshp)
    

#%% Find all the edges based on nodes number
def edges_link(poly2nodes):  
    """
    Find all the edges based on nodes number
    """
    edge2nodes = [] #the corresonding node number for edges 
    max_edges = np.shape(poly2nodes)[1]
    for i in icells:        
        polyi = poly2nodes[i,poly2nodes[i,:]>0]
        n_edges = len(polyi)
        poly_edges = [(polyi[ni-1],polyi[ni]) for ni in range(n_edges)]
        if n_edges<max_edges:
            poly_edges = poly_edges + [(-99,-99)]*(max_edges-n_edges)
        edge2nodes.append(poly_edges)   
    edge2nodes = np.reshape(edge2nodes,[np.shape(edge2nodes)[0]*
                                        np.shape(edge2nodes)[1],
                                        np.shape(edge2nodes)[2]])
    
    
    edge2nodes = np.sort(edge2nodes,axis=1) # order the edge indices
    edge2nodes = np.unique(edge2nodes,axis=0) # remove duplicated edges
    edge2nodes = edge2nodes[edge2nodes.min(axis=1)>=0,:] #remove rows with negative values
    edge2nodes = edge2nodes.astype(int)
    return edge2nodes

def cells_to_edges(poly2nodes,edge2nodes):
    """
    map cells to edges; this may take a minute or two
    """
    cells2edges = np.ones_like(poly2nodes)*-99
    for i in icells:
        polyi = poly2nodes[i,poly2nodes[i,:]>0]
        max_edges = len(polyi)
        poly_edges = [(polyi[ni-1],polyi[ni]) for ni in range(max_edges)]
        poly_edges = np.sort(poly_edges,axis=1) #sort the edge indices
        for p, pos in enumerate(poly_edges):
            ind = np.argwhere( (pos[0]==edge2nodes[:,0]) & (pos[1]==edge2nodes[:,1]))
            if ind:
                cells2edges[i,p]=ind

    cells2edges = cells2edges.astype(int)    
    return cells2edges
    
#%%
def edges_to_cell(edgeNumber):           
    ind = np.argwhere(cells2edges==edgeNumber)
    if len(ind)==1:
        edges2cell = np.asarray([ind[:,0][0],-99])
    else:
        edges2cell = ind[:,0]
    return edges2cell 
#%%
    
def cells_to_cells(cells2edges):
    """ map cells to cells; this may take a couple of minutes
    """       
    cell2cell = []
    for i in icells:
        edgei = cells2edges[i]
        celli_cell = np.asarray([edges_to_cell(e) for e in edgei if e>=0])
        celli_cell = celli_cell[(celli_cell!=i) & (celli_cell>0)]
        cell2cell.append(celli_cell)
    return cell2cell

            
#%% This is a fairly slow algorithm: find edge_cells.  
#Edge_cells=np.ones_like(Edges_x)*-999
#for e, (ex,ey) in enumerate(zip(Edges_x,Edges_y)):
#    ind = np.argwhere((poly_x==ex[0]) & (poly_y==ey[0])) 
#    indc = ind[:,0] #  the indices of the cell
#    indi = ind[:,1] # the indies of the node on the cell
#    if isinstance(indc,np.int64):       
#        maxnode = sum(~np.isnan(poly_x[indc]))
#        ind2 = indi+1
#        ind1 = indi-1
#        if ind2>=maxnode:
#            ind2= ind2-maxnode
#        if ([ex[1],ey[1]] ==[poly_x[indc,ind1],poly_y[indc,ind1]] ):
#            Edge_cells[e,1]=indc+1 #to cell
#        elif ([ex[1],ey[1]] ==[poly_x[indc,ind2],poly_y[indc,ind2]] ):
#            Edge_cells[e,0]=indc+1 #from cell 
#    else:
#        for i,j in zip(indc,indi):
#            maxnode = sum(~np.isnan(poly_x[i]))
#            ind2 = j+1
#            ind1 = j-1
#            if ind2>=maxnode:
#                ind2= ind2-maxnode
#            if ([ex[1],ey[1]] ==[poly_x[i,ind1],poly_y[i,ind1]] ):
#                Edge_cells[e,1]=i #to cell
#            elif ([ex[1],ey[1]] ==[poly_x[i,ind2],poly_y[i,ind2]] ):
#                Edge_cells[e,0]=i #from cell 
#    

