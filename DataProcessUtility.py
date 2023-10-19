## -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 10:21:13 2017
Various functions to process the data
@author: zhenlinz
"""

import numpy as np
from scipy.signal import butter, lfilter
import datetime
from matplotlib.dates import date2num 
from matplotlib.dates import num2date
import warnings

def remove_duplicates(var):
    """ This function returns duplicated indices and the ordered unique vars
    input var could be either string or number 
    output:
        uorder:ordered var with duplicates removed
        indices: the indices where var is unique
        ind: the indices of output indices where duplicates were removed; This is removed  
    """
    u, indices = np.unique(var, return_index=True)
    indices = np.sort(indices)
    uorder = var[indices]
#    if len(u)<len(var):
#        ind = np.argwhere(np.diff(indices)>1)+1
#        if indices[-1]!=len(var)-1:  
#            ind = np.append(ind,np.r_[indices[-1]+1:len(var):1])
#    else:
#        ind = []
            
    return(uorder,indices)


def butter_bandpass(lowcut, highcut, fs, order=5):
    """Calculate frequency response for different orders
    when lowcut = 0: pass any freuqncy lower than highcut (lowpass): eg. moving avargesave 
    when highcut = 0: pass any frequency greater than lowcut (highpass): e.g., high frequency signal
    btype : {‘lowpass’, ‘highpass’, ‘bandpass’, ‘bandstop’}, optional
        The type of filter. Default is ‘lowpass’. """    
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    if lowcut == 0:
        btype = 'lowpass'
        b,a = butter(order,high,btype)
    elif highcut == 0:
        btype = 'highpass'
        b,a = butter(order,low,btype)
    else:
        btype = 'bandpass'                
        b, a = butter(order, [low, high], btype)
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order)
    y = lfilter(b, a, data)
    return y


def autoFillnan(time, interval,*var):
    """ This function finds the gap in time series and fill in . 
    time is datetime object in UTC (or else, convert to UTC) and interval means the time interval between the objects in days. 
    timeout is a datetime object in UTC. 
    outvar of both output time and time series.   
    """
    
    # convert datetime to number for easy manipulation
    timein= date2num(time)
    #varin = var[:] # using varin=var will cause varin to be associated with var. 
    varin = np.array(var,dtype=np.float64).squeeze()
#%%    
    if (np.size(varin)!=0) & (np.size(varin)!=np.size(timein)):
        warnings.warn('the shapes of the input var and time do not match')
    
    difftime = np.diff(timein)
    # For some reason, after converting to number, the interval was not exactly the same even for the same time interval, so round to integer. 
    ind = np.squeeze(np.argwhere(np.round((difftime-interval)/interval)>=1))
    inserttimes = np.squeeze(np.round((difftime[ind]-interval)/interval).astype(int)) 
    timeout = np.zeros(len(timein)+np.sum(inserttimes))*np.nan
    varout = np.zeros(len(varin)+np.sum(inserttimes))*np.nan
    
    if np.size(ind)==0:
        print('No auto fill required')
        if np.size(varin)==np.size(timein):
            outvar = np.vstack((time,varin))
        else:
            outvar = time[:]
    else:
        if np.size(ind)==1:
            timeout[0:ind+1]=timein[0:ind+1]
            timeout[ind+1:ind+1+inserttimes] = [timein[ind]+interval*(j+1) for j in range(inserttimes)]
            timeout[ind+1+inserttimes:]=timein[ind+1:]
            if np.size(varin)==np.size(timein):
                varout[0:ind+1]=varin[0:ind+1]
                varout[ind+1+inserttimes:]=varin[ind+1:]
        else:        
            timeout[0:ind[0]+1]=timein[0:ind[0]+1]
            if np.size(varin)==np.size(timein):
                varout[0:ind[0]+1]=varin[0:ind[0]+1]
            for i in np.arange(np.size(ind)):            
                timeout[ind[i]+1+np.sum(inserttimes[:i]):ind[i]+1+np.sum(inserttimes[:i])+inserttimes[i]] = [timein[ind[i]]+interval*(j+1) for j in range(inserttimes[i])]
                if i != np.size(ind)-1:
                    timeout[ind[i]+1+np.sum(inserttimes[:i+1]):ind[i+1]+1+np.sum(inserttimes[:i+1])] = timein[ind[i]+1:ind[i+1]+1]
                    if np.size(varin)==np.size(timein):
                        varout[ind[i]+1+np.sum(inserttimes[:i+1]):ind[i+1]+1+np.sum(inserttimes[:i+1])] = varin[ind[i]+1:ind[i+1]+1]
                else:
                    timeout[ind[i]+1+np.sum(inserttimes[:i+1]):] = timein[ind[i]+1:]  
                    if np.size(varin)==np.size(timein):
                        varout[ind[i]+1+np.sum(inserttimes[:i+1]):] = varin[ind[i]+1:]
        
        if np.size(varin)==np.size(timein):
            outvar = np.vstack((num2date(timeout),varout.astype('float')))
        else:
            outvar = num2date(timeout)
 #%%   
    return(outvar)

class TimeSeries:
    def __init__(self,time,var):
        self.time = time
        self.var = var

    @staticmethod
    def reorgnize(var,nx,ind_start=0,*ny):
        """ 
        Reshape the input time series to a matrix according to input rules        
        """
        varin = var     
        
        ny = ny or np.size(var)//nx
        
        if len(varin)>nx*ny+ind_start+1:
            varin = varin[:nx*ny+ind_start]  
        
        varout = np.reshape(varin,(ny,nx))
        return varout
        
    
    def Climatology(self):
        """ 
        Calculate climatologic condition based on the provided time series 
        interval: of the output climatological series; default is the time sereis interval in days
        Output:
            Cvar: climatology value of the var
            Cint: climatologic time series
            Yvar: Matrix of the var used to calculate Cvar
            Year: the year corresponds to Yvar
        """
        
        time = date2num(self.time)
        var = self.var
        
        month = np.array([timei.month for timei in self.time])
        day = np.array([timei.day for timei in self.time])
        
        #check if the time series is unevenly spaced, if so, make correction
        difftime = np.diff(time) 
        interval=  np.min(difftime)
        ind = np.argwhere(np.round((difftime-interval)/interval)>=1).squeeze() 
         
        
        if ind:            
            autoFillnan(time, interval,var) 
            
        # I will remove all the leap year days (29 Feb) so that the number of days is always 365. 
        indleapday = np.squeeze(np.argwhere((month==2) & (day==29) ))
        if len(indleapday)>0:
            time = np.delete(time,indleapday)
            var = np.delete(var,indleapday)
            
        # Now we can reorganize the time series into matrix
        nx = np.round(365/interval)
        Yvar = self.reorgnize(var,nx=nx)
        Ytime = self.reorgnize(time,nx=nx)
        Ydate = num2date(Ytime)
        Year =[Ydate[i][0].year for i in range(np.shape(Ydate)[0])]
        Cint = date2num(Ydate[0][:]) - date2num(datetime.datetime(Year[0],1,1))
            
        Cint[np.argwhere(Cint>365)] = Cint[np.argwhere(Cint>365)]-365
        Cvar = np.nanmean(Yvar,axis=0)        
              
        return(Cvar,Cint,Yvar,Year)   
    
    def dailyAve(self):
        """
        Calculate daily average
        """
        time= date2num(self.time)
        
        var = self.var
        #check if the time series is unevenly spaced, if not, make correction
        difftime = np.diff(time) 
        interval=  np.min(difftime) 
        
        ind = np.argwhere(np.round((difftime-interval)/interval)>=1).squeeze()         
        
        if ind:            
            autoFillnan(time, interval,var)    
        
        nx = np.round(1/interval) #number of input       
        
       
        var_matrix = self.reorgnize(var,nx)
        time_matrix = self.reorgnize(time,nx)
        var_daily = np.nanmean(var_matrix,axis=1)
        time_daily = np.nanmean(time_matrix,axis =1)
        date_daily = num2date(time_daily)
        
        return(date_daily,var_daily)
        
        
        
        
  

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import freqz

    # Sample rate and desired cutoff frequencies (in Hz).
    fs = 5000.0
    lowcut = 500.0
    highcut = 1250.0

    # Plot the frequency response for a few different orders.
    plt.figure(1)
    plt.clf()
    for order in [3, 6, 9]:
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        w, h = freqz(b, a, worN=2000)
        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
             '--', label='sqrt(0.5)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Gain')
    plt.grid(True)
    plt.legend(loc='best')

    # Filter a noisy signal.
    T = 0.05
    nsamples = T * fs
    t = np.linspace(0, T, nsamples, endpoint=False)
    a = 0.02
    f0 = 600.0
    x = 0.1 * np.sin(2 * np.pi * 1.2 * np.sqrt(t))
    x += 0.01 * np.cos(2 * np.pi * 312 * t + 0.1)
    x += a * np.cos(2 * np.pi * f0 * t + .11)
    x += 0.03 * np.cos(2 * np.pi * 2000 * t)
    plt.figure(2)
    plt.clf()
    plt.plot(t, x, label='Noisy signal')

    y = butter_bandpass_filter(x, lowcut, highcut, fs, order=6)
    plt.plot(t, y, label='Filtered signal (%g Hz)' % f0)
    plt.xlabel('time (seconds)')
    plt.hlines([-a, a], 0, T, linestyles='--')
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='upper left')

    plt.show()
