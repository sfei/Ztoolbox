
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 16:57:29 2018

@author: zhenlinz
"""

import json
import psycopg2
import pandas as pd
import os

dbpath = r'C:\Users\zhenlinz\Google Drive\1_Nutrient_Share\1_Projects_NUTRIENTS\Modeling\NOTES_ZZ\SuisunBay\Data\database'

with open(os.path.join(dbpath,'postgreConfig.json')) as f:
    conf = json.load(f)
    
conn_str = "host={} dbname={} user={} password={}".format(\
                 conf['host'],conf['database'],conf['user'],conf['passw'])


conn = psycopg2.connect(conn_str)


text = """
        SELECT distinct          
          allresults.stationcode,
          allresults.stationname,
          stationlookup.targetlatitude,
          stationlookup.targetlongitude
        FROM 
          allresults 
        LEFT JOIN stationlookup on stationlookup.stationcode=allresults.stationcode
        WHERE 
          parametercode = '80154'         
        """
#AND sampleeventdatetime< '2011-10-01' :: date    
df = pd.read_sql(text,con=conn)

# parametercode:
# ssc: 80154
# chl: 32283

conn.close()

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
    CreateQgisPoints(df,latname = 'targetlatitude',lonname = 'targetlongitude')