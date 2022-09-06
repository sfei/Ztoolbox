# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 10:32:39 2018

@author: zhenlinz

updated by alliek in september 2022

note from alliek: ask Shira or Gemma for the NUTVIZ password! she will share it in last pass, along with the
server and the username. if password doesn't work, check if you need to change the server (currently
toyon-dev.sfei.org) or the username (currently nutviz_ro) specified in AccessPGSQL function below...

another note from alliek: i can't get this to work on richmond or chicago, so i'm running it on my laptop for the 
time being, may try to fix later

"""


import psycopg2
import pandas as pd
import sys

def AccessPGSQL(nutviz_password,text):        
    
    #    with open('postgreConfig.json') as f:
    #        conf = json.load(f)
            
    # conn_str = "host={} dbname={} user={} password={}".format(\
    #                   "redbud_dev.sfei.org","nutviz","nutviz_ro","crispywolf")

    conn_str = "host={} dbname={} user={} password={}".format(\
                      "toyon-dev.sfei.org","nutviz","nutviz_ro",nutviz_password)
    
    # edited by SW in 2020 to have new hostname, password, etc 

    conn = psycopg2.connect(conn_str)
    
    try:
        df = pd.read_sql(text,con=conn)
        conn.close()
    except:
        print(sys.exc_info())
        conn.close() 
        print("data base closed")        
        df={}
    
    return df


def GetParameters(nutviz_password,*param):
    text = """
            SELECT 
                parametercode, 
                parameternameshort, 
                parameternamelong, 
                unit
            FROM
                parameterlookup"""
    
    textadd = []         

    if len(param)>0:
        text = text + " where "
        
    for p in param:
        textadd.append("""UPPER(parameternamelong) like UPPER(""" + "'%"+ p +"%')"+ """ """)
    
    if len(textadd)>1:
        textadd = " OR ".join(textadd)
    
    if len(textadd)==1:
        text = text+textadd[0]
    elif len(textadd)>1:
        text = text+textadd
            
    df = AccessPGSQL(nutviz_password,text)
    
    return df
        
def GetStations(nutviz_password,*stationname):
    text = """
            SELECT *            
            FROM
                stationlookup"""
    
    textadd = []         

    if len(stationname)>0:
        text = text + " where "
        
    for p in stationname:
        textadd.append("""UPPER(stationname) like UPPER(""" + "'%"+ p +"%')"+ """ """)
    
    if len(textadd)>1:
        textadd = " OR ".join(textadd)
    
    if len(textadd)==1:
        text = text+textadd[0]
    elif len(textadd)>1:
        text = text+textadd
            
    df = AccessPGSQL(nutviz_password,text)
    
    return df


def GetAnalyte(nutviz_password,parametercode,start_date,end_date,*stationcode):  
    """
    parametercode could be either one text code or a list of the codes. 
    """
    if (type(parametercode)== list) & (len(parametercode)>1): # more than one analyte
        AnalyteList = ["parametercode=\'""" + codei + "\'" for codei in parametercode] 
        AnalyteText = '('+ ' OR '.join(AnalyteList) + ')'
    else:
        AnalyteText = "parametercode='" + parametercode + "'"
        
        
    text = """
            SELECT 
                parameternameshort,
                parametercode,
                sitecode, 
                stationcode,
                sampleeventdatetime,
                result,
                unit
            FROM
                allresults
            WHERE """ + AnalyteText + \
                """ AND sampleeventdatetime<'""" + end_date +"""' :: date
                    AND sampleeventdatetime>'""" + start_date +"""' :: date
            ORDER BY sampleeventdatetime """ 
    
    df = AccessPGSQL(nutviz_password,text)
    
    return df

def GetStation(nutviz_password,stationcode,start_date,end_date):  
    """
    stationcode could be either one text code or a list of the codes. 
    """
    if (type(stationcode)== list) & (len(stationcode)>1): # more than one analyte
        StationList = ["stationcode=\'""" + codei + "\'" for codei in stationcode] 
        StationText = '('+ ' OR '.join(StationList) + ')'
    else:
        StationText = "stationcode='" + stationcode + "'"
        
        
    text = """
            SELECT 
                parameternameshort,
                parametercode,
                sitecode, 
                sampleeventdatetime,
                result,
                unit
            FROM
                allresults
            WHERE """ + StationText + \
                """ AND sampleeventdatetime<'""" + end_date +"""' :: date
                    AND sampleeventdatetime>'""" + start_date +"""' :: date
            ORDER BY sampleeventdatetime """ 
    
    df = AccessPGSQL(nutviz_password,text)
    
    return df    


def GetAnalyteStation(nutviz_password,stationcode,parametercode,start_date,end_date):  
    """
    stationcode could be either one text code or a list of the codes. 
    """
    if (type(stationcode)== list) & (len(stationcode)>1): # more than one analyte
        StationList = ["stationcode=\'""" + codei + "\'" for codei in stationcode] 
        StationText = '('+ ' OR '.join(StationList) + ')'
    else:
        StationText = "stationcode='" + stationcode + "'"
        
    if (type(parametercode)== list) & (len(parametercode)>1): # more than one analyte
        AnalyteList = ["parametercode=\'""" + codei + "\'" for codei in parametercode] 
        AnalyteText = '('+ ' OR '.join(AnalyteList) + ')'
    else:
        AnalyteText = "parametercode='" + parametercode + "'"
        
        
    text = """
            SELECT 
                parameternameshort,
                parametercode,
                sitecode, 
                sampleeventdatetime,
                result,
                unit
            FROM
                allresults
            WHERE """ + StationText + \
                " AND " + AnalyteText + \
                """ AND sampleeventdatetime<'""" + end_date +"""' :: date
                    AND sampleeventdatetime>'""" + start_date +"""' :: date
            ORDER BY sampleeventdatetime """ 
    
    df = AccessPGSQL(nutviz_password,text)
    
    return df    
        
if __name__ == "__main__":
    from database import io_postgre
    HFdf = io_postgre.GetAnalyteStation('DMB',['32315','32283'],'2013-01-01','2013-05-01')
            
            
    
    
