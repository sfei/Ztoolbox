# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 15:42:36 2017
Access the database and plot the station map of chlorophyll for 2011. 

@author: zhenlinz
"""

import sqlite3
import pandas as pd
import os

#%%
databasepath = "C:/Users/zhenlinz/Google Drive/1_Nutrient_Share/1_Projects_NUTRIENTS/Modeling/NOTES_ZZ/SuisunBay/Data/database"
databasefile = os.path.join(databasepath,"DeltaSuisun2.sqlite")
coon = sqlite3.connect(databasefile)
df = pd.read_sql_query("""
                       	SELECT DISTINCT
		                   "STATION NUMBER",                          
		                   StationTable.Latitude AS Lat,
		                   StationTable.Longitude AS Lon,
                          StationTable."Site Name" AS SiteName
	                    FROM 
		                   DiscreteTable
	                    INNER JOIN StationTable on DiscreteTable."Station Number" = StationTable."Site Number"
	                    WHERE Analyte="Chlorophyll a"
                        AND julianday("Collection Date")>=julianday('2010-08-01') 
		                   AND julianday("Collection Date")<julianday('2011-10-01');                       
                       """,
                       coon)

df2 = pd.read_sql_query("""
                       	SELECT DISTINCT
		                   "STATION NUMBER",                          
		                   StationTable.Latitude AS Lat,
		                   StationTable.Longitude AS Lon,
                          StationTable."Site Name" AS SiteName
	                    FROM 
		                   DiscreteTable
	                    INNER JOIN StationTable on DiscreteTable."Station Number" = StationTable."Site Number"
	                    WHERE Analyte="Dissolved Nitrite + Nitrate" 
                        AND julianday("Collection Date")>=julianday('2010-08-01') 
		                   AND julianday("Collection Date")<julianday('2011-10-01');                       
                       """,
                       coon)


coon.close()

df.columns = ['Station','Lat','Lon','SiteName'] # Nutrient sites
df2.columns = ['Station','Lat','Lon','SiteName'] # chlorophyll sites




#%% 

# Plot on Bokeh in spyder
if __name__ == "__main__":  
    
    from bokeh.io import output_file, show
    from bokeh.models import (
        GMapPlot, GMapOptions,Circle, Diamond, HoverTool,DataRange1d, ColumnDataSource, PanTool, WheelZoomTool,SaveTool)

    map_options = GMapOptions(lat=38.18, lng=-121.8, map_type='terrain',zoom=9) 
    gplot = GMapPlot(x_range=DataRange1d(),y_range=DataRange1d(),map_options = map_options)
    gplot.title.text = 'Station Locations'
    
    # For GMaps to function, Google requires you obtain and enable an API key:
    #
    #     https://developers.google.com/maps/documentation/javascript/get-api-key
    #
    # Replace the value below with your personal API key:
    GOOGLE_API_KEY = "AIzaSyAX1D9G6UVjaZrr63qbZ_aoMJpixUyIHEc"
    gplot.api_key = GOOGLE_API_KEY
    
    psource =  ColumnDataSource(df)
    psource2 = ColumnDataSource(df2)
    
    circle = Circle(x="Lon",y="Lat", fill_color= 'red',size=10,fill_alpha =0.8)
    diamond= Diamond(x="Lon",y="Lat",fill_color = 'green',size=8)
    #gplot.add_glyph(psource,circle)
    gplot.add_glyph(psource,diamond)
    my_hover = HoverTool()
    my_hover.tooltips = [('Station',"@SiteName")]
    
    gplot.add_tools(PanTool(),WheelZoomTool(),my_hover,SaveTool())
    output_file("gmap_example_P.html")
    show(gplot)

#%%
def CreateQgisPoints(df,latname,lonname):
    """
    This function plots the data points from a pandas dataframe
    """
    # Create a new polygon memory layer
    bathLayer = QgsVectorLayer("Point?crs=EPSG:4269", "PointLayer","memory") #4269 corresponds to NAD 83   
     
    pr = bathLayer.dataProvider()
    
    
    # Add attributes fieldname and type
    FieldColumns = [column for column in df.columns if column not in [latname,lonname]]
    datatype = [type(df[Field].values[0]) for Field in FieldColumns]
    # Qgis only accepts three types at attributes: double, str, and int, so need to make conversion.
    try:
        datatype = ['Double' for typelist in datatype if 'float' in typelist]
    except TypeError:
        print('Converting float to double type field')

    try:
        datatype = ['Int' for typelist in datatype if 'int' in typelist]
    except TypeError:
        print('Int type field')
        
    try: 
        datatype = ['String' for typelist in datatype if typelist not in ['double','int']] # all other datatypes will be converted to string
    except TypeError:
        print('String type field')    
        
       
    AttributeList = [QgsField(column,eval('QVariant.'+typelist)) for column, typelist in zip(FieldColumns, datatype)]
    res = pr.addAttributes(AttributeList)
    
    bathLayer.updateFields() #Now updating the layer fields. 
    
      
    for i in range(len(df)):       
        
        # create a point feature
        pointi = QgsPoint(df[lonname].values[i],df[latname].values[i])  
        
        point = QgsFeature()        
        point.setGeometry(QgsGeometry.fromPoint(pointi))   
           
        point.setAttributes(list(df[FieldColumns].values[i])) #set attribute index for poly. 
        
        pr.addFeatures([point])
        bathLayer.updateExtents()
    # register the layer
    QgsMapLayerRegistry.instance().addMapLayers([bathLayer])
    

if __name__ == "__console__":
    import qgis
    from PyQt4.QtCore import QVariant
    CreateQgisPoints(df,latname = 'Lat',lonname = 'Lon')








