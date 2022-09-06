# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 17:13:35 2017

@author: zhenlinz
"""
import pandas as pd
import sqlite3
import sys

dbfile = r'C:\Users\zhenlinz\Google Drive\1_Nutrient_Share\1_Projects_NUTRIENTS\Modeling\NOTES_ZZ\SuisunBay\Data\database\DeltaSuisun2.sqlite'
def AccessSQLITE(text):        
    
    #    with open('postgreConfig.json') as f:
    #        conf = json.load(f)
            
    conn = sqlite3.connect(dbfile)        
    
    try:
        df = pd.read_sql(text,con=conn)
        conn.close()
    except:
        print(sys.exc_info())
        conn.close() 
        print("data base closed")        
        df={}
    
    return df


def readsql_analyte(dbfile,start_date,end_date,Analyte):
    """
    Read all data related to the analyte or analytes for the time period
    """
    coon = sqlite3.connect(dbfile)     
    if (type(Analyte)== list) & (len(Analyte)>1): # more than one analyte
        AnalyteList = ["Analyte=\"""" + Analytei + "\"" for Analytei in Analyte] 
        AnalyteText = '('+ ' OR '.join(AnalyteList) + ')'
    else:
        AnalyteText = "Analyte=\"""" + Analyte + "\""
            
    text = """
            SELECT 
		       "STATION NUMBER",   
              "Collection Date",
              "Analyte",
              "Result",
              "Depth",
		       StationTable.Latitude AS Lat,
		       StationTable.Longitude AS Lon,
              StationTable."Site Name" AS SiteName
	         FROM 
		        DiscreteTable
	         INNER JOIN StationTable on DiscreteTable."Station Number" = StationTable."Site Number"
	         WHERE """ + AnalyteText + \
               """ AND julianday("Collection Date")>=julianday(\"""" + start_date + "\")" + \
               """ AND julianday("Collection Date")<julianday(\"""" + end_date + "\")" + \
		      """ AND typeof(Result)="real"
            ORDER BY Lat;
           """
    df = pd.read_sql_query(text,coon)
    coon.close()
    return df

def readsql_site(dbfile,start_date,end_date,SiteNumber):
    """ 
    Read all data from a specific station
    """
    coon = sqlite3.connect(dbfile) 
    text = """
            SELECT 
		       "STATION NUMBER",   
              "Collection Date",
              "Analyte",
              "Result",
              "Depth",
		       StationTable.Latitude AS Lat,
		       StationTable.Longitude AS Lon,
              StationTable."Site Name" AS SiteName
	         FROM 
		        DiscreteTable
	         INNER JOIN StationTable on DiscreteTable."Station Number" = StationTable."Site Number"
	         WHERE "Station Number" = \"""" + SiteNumber + "\"" \
               """ AND julianday("Collection Date")>=julianday(\"""" + start_date + "\")" + \
               """ AND julianday("Collection Date")<julianday(\"""" + end_date + "\")" + \
		      """ AND typeof(Result)="real"
            ORDER BY Lat;
           """
    df = pd.read_sql_query(text,coon)
    coon.close()
    return df

def readsql(dbfile,start_date,end_date,SiteNumber,Analyte):
    """ 
    Read specific Analytes from specific stations
    SiteNumber could be station number or site name. 
    """
    coon = sqlite3.connect(dbfile)     
    if (type(Analyte)== list) & (len(Analyte)>1): # more than one analyte
        AnalyteList = ["Analyte=\"""" + Analytei + "\"" for Analytei in Analyte] 
        AnalyteText = '('+ ' OR '.join(AnalyteList) + ')'
    else:
        AnalyteText = "Analyte=\"""" + Analyte + "\""
    
    if (type(SiteNumber)== list) & (len(SiteNumber)>1): # more than one analyte
        SiteList = [""" "Station Number" =\"""" + Sitei + "\"" for Sitei in SiteNumber] 
        SiteList = SiteList + [""" "SiteName" =\"""" + Sitei + "\"" for Sitei in SiteNumber] 
        SiteText = '('+ ' OR '.join(SiteList) + ')'
    else:
        SiteText = """ ("Station Number"  =\"""" + SiteNumber + "\""  + \
                    ' OR ' + \
                    """ "SiteName"  =\"""" + SiteNumber + "\")"

    
    text = """
            SELECT 
		       "STATION NUMBER",   
              "Collection Date",
              "Analyte",
              "Result",
              "Depth",
              "Rpt Limit",
		       StationTable.Latitude AS Lat,
		       StationTable.Longitude AS Lon,
              StationTable."Site Name" AS SiteName
	         FROM 
		        DiscreteTable
	         INNER JOIN StationTable on DiscreteTable."Station Number" = StationTable."Site Number"
	         WHERE """ + AnalyteText + \
               """ AND """ + SiteText + \
               """ AND julianday("Collection Date")>=julianday(\"""" + start_date + "\")" + \
               """ AND julianday("Collection Date")<julianday(\"""" + end_date + "\")" + \
		      """ AND typeof(Result)="real"
            ORDER BY Lat;
           """

    
    df = pd.read_sql_query(text,coon)
    coon.close()
    
    return df    
    

def dftoGeodf(df):
    """convert pandas dataframe to geopandas dataframe using 'Lat' and 'Lon'
       values provided in the database
    """
    import geopandas
    from shapely.geometry import Point
    geometry = [Point(xy) for xy in zip(df.Lon,df.Lat)]
    df = df.drop(['Lon','Lat'],axis=1)
    crs = {'init':'epsg:4326'}
    gdf = geopandas.GeoDataFrame(df,crs=crs,geometry=geometry)
    return gdf
    
    

if __name__ == "__main__": 
    start_date = '2011-09-01'
    end_date =  '2012-06-01'
    dbfile = "../database/DeltaSuisun.sqlite"
    df = readsql(dbfile,start_date,end_date,["Dissolved Ammonia","Dissolved Nitrate + Nitrite" ]) #Dissolved Ammonia "Dissolved Nitrate