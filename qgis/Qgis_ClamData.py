import pandas as pd
import os

databasepath = "C:/Users/zhenlinz/Google Drive/1_Nutrient_Share/1_Projects_NUTRIENTS/Modeling/NOTES_ZZ/SuisunBay/Data/database/csv"
datafile= os.path.join(databasepath,"Corbicula.csv")
df = pd.read_csv(datafile)

#%% get data for 2011
#df['Datetime'] = pd.to_datetime(df['Date'])
#df = df[(df['Datetime']>=pd.to_datetime('2010-08-01')) & (df['Datetime']<pd.to_datetime('2011-10-01'))] 

#%%
df_site = pd.DataFrame()
for site in df.Site.unique():
    df_site = df_site.append(df[df.Site==site].iloc[0]) 

df_site = df_site[['Lat','Lon','Site']] 


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
    CreateQgisPoints(df_site,latname = 'Lat',lonname = 'Lon')

