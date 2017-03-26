import csv
import sys
import numpy as np
#import scipy as sp
import math
from lininterp import interp
import interpbelow as ib

#Applies the virtual temperature correction to a pressure (in mb) , temperature and dew point (in Kelvin)
def virtualtempcorr(p,t,q):
    
    eps=0.608
    tv=t*(1+eps*q)
    return tv
    
#Takes a single pressure and temperature (in Kelvin) measurement and computes the (saturation) mixing ratio
def w_s(p,t):
    e_s=6.1121*math.exp(17.502*(t-273.15)/(t+240.97-273.15))
    w_s=0.622*e_s/(p-e_s)
    return w_s

#Takes lists of pressure and virtual temperature and calculates thickness
def hypsometric(p,tv,sfcp):
    Rd = 287.058 # J/Kg*k
    g = 9.81 ; #m/s^2
    
    thickness = 0 - 0.5*(tv[0]+tv[0]*Rd/g*np.log(p[0]/sfcp))
    thickness = [thickness]
    
    for i in range(1,len(p),1):
        thickness.append(thickness[i-1]-0.5*(tv[i-1]+tv[i])*Rd/g*np.log(p[i]/p[i-1]))
        
    return thickness
#Take a single pressure (mb) and mixing ratio (kg/kg) measurement and converts it to dew point temperature
def Tdcalc(p,w):
    e = w*p/(0.622+w)
    Td = 240.97*np.log(e/6.1365)/(17.502-np.log(e/6.1365))
    return Td

def interpolater(htlist, plist,vars,sfcp,resolution,folder,tod):
    
    
    temp = []
    td = []
    qlist=[]
    
    pinterp = []
    hinterp = []
    tinterp = []
    qinterp = []
    uinterp = []
    vinterp = []
    MODthick=[]
    model_hinterp = []
    
    tbar = []
    pbar = []
    tv = []
    model_tv = []
    press=[]
    
    g=9.81      #m/s^2
    R=287.058   #J/Kg*K
    
    
    #Create variables

    if folder == 'LTT' or folder == 'LLT' or folder == 'NARR' or folder == 'F1F2' or folder == 'F3Plus':

        temp = vars[0]
        q = vars[1]
        #Remove 0 points
        for i in range(0,len(plist),1):
            if q[i] == 0:
                q[i] = 0.00000001
        u = vars[2]
        v = vars[3]
        #for i in range(0,len(qlist),1):
        #    td.append(Tdcalc(plist[i],qlist[i]/1000)+273.15)
    
    belowsfc = 0

    #Modify lists for observed surface pressure
    for i in range(0,len(plist),1):
        if plist[i] > sfcp:
            belowsfc += 1
            
    if belowsfc >= 2 :
        #Remove values below surface aside from first measurement below surface to allow for interpolation to surface
        for i in range(belowsfc-2,0,-1):
            plist.pop(i)
            temp.pop(i)
            q.pop(i)
            u.pop(i)
            v.pop(i)
    if belowsfc >= 1 :
        #Interpolate to Surface
        temp.insert(1,interp(temp,plist,sfcp))
        del temp[0]
        q.insert(1,interp(q,plist,sfcp))
        del td[0]
        u.insert(1,interp(u,plist,sfcp))
        del u[0]
        v.insert(1,interp(v,plist,sfcp))
        del v[0]
        plist.insert(1,sfcp)
        del plist[0]
    elif belowsfc == 0 and np.float(plist[0]) != sfcp:
        temp.insert(0,ib.interp(temp,plist,sfcp))
        q.insert(0,ib.interp(q,plist,sfcp))
        u.insert(0,ib.interp(u,plist,sfcp))
        v.insert(0,ib.interp(v,plist,sfcp))
        plist.insert(0,sfcp)

        
    
    #Create virtual temperature
    for i in range(0,len(temp),1):
        tv.append(virtualtempcorr(plist[i],temp[i],q[i]))
    
    #Calculate un-interpolated thickness values
    thick = htlist
    basealtitude = thick[0]
    #Change thickness values to agl
    for i in range(0,len(thick),1):
        thick[i] -= basealtitude

    #Interpolate thickness values
    hinterp.append(thick[0])
    
    
    maxht = thick[-1]
    j = 1
    while hinterp[j-1] <= maxht:
        hinterp.append(hinterp[j-1]+resolution)
        j+=1
        
    #Interpolate p, T, Td, u & v
    pinterp = np.interp(hinterp,thick, plist)
    tinterp = np.interp(hinterp,thick, temp)
    qinterp = np.interp(hinterp,thick, q)
    uinterp = np.interp(hinterp,thick, u)
    vinterp = np.interp(hinterp,thick, v)

    
    return pinterp, hinterp, tinterp, qinterp, uinterp, vinterp
    
    
    
    
        
    
        
        
    
    
    
        

