# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 12:23:17 2017
TimeUtility.py
manipulate and make conversion between time objects
@author: zhenlinz
"""

import datetime
import numpy as np
import dateutil

def TimeConverter(TimeIn,Informat,Outformat):
    """
    The inputs must all be string    
    strptime converts a formatted date string to datetime; strftime converts datatime to formatted date string 
    Informat e.g., "%d/%m/%y %H:%M"   '21/11/06 16:30'
             e.g., "%Y-%m-%d %H:%M"   '2002-01-01 00:06'
    """    
    if np.size(TimeIn)==1:
        TimeOut = datetime.datetime.strptime(TimeIn, Informat).strftime(Outformat)
    else:
        TimeOut = [datetime.datetime.strptime(Timei, Informat).strftime(Outformat) for Timei  in TimeIn]
    
    return TimeOut


def Str2Datetime(Timestr,Informat):
    """
    Informat e.g., "%d/%m/%y %H:%M"   '21/11/06 16:30'
             e.g., "%Y-%m-%d %H:%M"   '2002-01-01 00:06'
    """
    if np.size(Timestr)==1:
        TimeOut = datetime.datetime.strptime(Timestr, Informat)
    else:
        TimeOut = [datetime.datetime.strptime(Timei, Informat) for Timei  in Timestr]
    
    return TimeOut

def Datetime2Str(Time):
    """
    This is a quick converstion from datetime to a default string format
    """
    Outformat = "%Y-%m-%d %H:%M"
    if np.size(Time)==1:
        Timeout = Time.strftime(Outformat)
    else:
        Timeout = [Timei.strftime(Outformat) for Timei in Time]
    
    return Timeout

                

def StrTimeDiff(Timestr1,Timestr2, Informat):
    """
    Finding the time difference in seconds between Timestr1 and Timestr2
    Timestr1 and Timestr2 need to be the same size
    Informat e.g., "%d/%m/%y %H:%M"   '21/11/06 16:30'
             e.g., "%Y-%m-%d %H:%M"   '2002-01-01 00:06'
    """
    
    if np.size(Timestr1)==1:
        Timedate1 = datetime.datetime.strptime(Timestr1,Informat)
        Timedate2 = datetime.datetime.strptime(Timestr2,Informat)
    else:
        Timedate1 = [datetime.datetime.strptime(Timei, Informat) for Timei  in Timestr1]
        Timedate2 = [datetime.datetime.strptime(Timei, Informat) for Timei  in Timestr2]
    
    DiffTime = (Timedate2-Timedate1).total_seconds()
    
    return DiffTime

def GMT2Local(timein):
    #Timezone conversion between GMT and local time (note that loal time takes into account of daylight savings)
    #timein:datetime object
    
    from_zone = dateutil.tz.gettz('GMT')
    to_zone =  dateutil.tz.tzlocal() #dateutil.tz.gettz('PDT')
    
    # Tell the datetime object that it's in UTC time zone since 
    # datetime objects are 'naive' by default
    gmt = timein.replace(tzinfo=from_zone)

    # Convert time zone
    timeout = gmt.astimezone(to_zone)
    
    return timeout



        
    
    