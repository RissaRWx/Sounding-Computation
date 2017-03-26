#!/usr/bin/env python

import sys
import math
import csv
import numpy as np
import Bisection
import thetae
import lininterp
import operator
from svrparameters import bulkshear


#Takes dewpoint (in Kelvin) and pressure (in mb) and converts to vapor pressure
def e_s(t):
    
    e_s=6.112*math.exp(17.5*(t-273.15)/(t+237.3-273.15))
    return e_s
    

#Applies the virtual temperature correction to a pressure (in mb) , temperature and dew point (in Kelvin)
def virtualtempcorr_env(p,t,q):
    
    eps=0.608
    tv=t*(1+eps*q)
    return tv
    
#Applies the virtual temperature correction to a pressure (in mb) , temperature and dew point (in Kelvin)
def virtualtempcorr_path(p,t,td):
    
    eps=0.608
    tv=t*(1+eps*w_s(p,td))
    return tv

#Equation to calculate the enthalpy of vaporization from temperature in Kelvin
def Lv(T):
    Lv = 1.91846 *math.pow(10.0,6.0)*math.pow(T/(T-33.91),2.0)
    #Lv = 2500300
    return Lv

#Turns the list of dew points (in Kelvin) and pressures in to a list of absolute humidity
def q(p,t):
    qlist = []
    e_s = []
    for i in range(0,len(p), 1):
        e_s.append(6.112*math.exp(17.5*(t[i]-273.15)/(t[i]+237.3-273.15)))
        qlist.append(0.622*e_s[i]/(p[i]-e_s[i]))
    return qlist

#Takes a single pressure and temperature (in Kelvin) measurement and computes the (saturation) mixing ratio
def w_s(p,t):
    if (t+237.3-273.15) < 0:
        w_s = 0.0001
        return w_s
    e_s=6.112*math.exp(17.5*(t-273.15)/(t+237.3-273.15))
    #e_s = 6.1121*math.exp((2500800/461.51)*(t-273.15)/(t*273.15)+(-.0009477)*math.log(t/273.15)-(-.0009477)*(t-273.15)/t)
    w_s=0.622*e_s/(p-e_s)
    return w_s
    
def mlcape(hghtlist,plist,tlist,qlist):

    eps = 0.622
    Cp=1004.0
    Rd=287.04
    g=9.81
    po=1000.0
    L = 2501000
    hincrement = 1
    p_env = []
    t_env = []
    td_env = []
    q_env = []
    h_env = []
    tpath_corr = []
    t_env_corr = []
    tpath = []
    thlist = []
    mldepth = 100
    tempsum = 0
    tdsum = 0
    psum = 0
    qsum = 0
    
    mlcin = 0
    mlcape = 0
    mllfc = 0
    mllcl = 0
    mlel = 0
    mltotalcape = 0
    mlct = 0
    cape3km = 0
    
    onekmtheta = 0
    sfctheta = 0
    lapse = 0
    
    
    #qlist=q(plist,tdlist)

    #Determine lowest 1km lapse rate
    for i in range(0,len(hghtlist)-1,1):
        if hghtlist[i] == 1000 :
            onekmtheta = tlist[i]*math.pow(po/plist[i],Rd/Cp)
        if hghtlist[i] == 500:
            theta500 = tlist[i]*math.pow(po/plist[i],Rd/Cp)
        if hghtlist[i] == 100:
            theta100 = tlist[i]*math.pow(po/plist[i],Rd/Cp)
    
    sfctheta = tlist[0]*math.pow(po/plist[0],Rd/Cp)
    lapse1000 = (onekmtheta - sfctheta)
    #lapse500 = (theta500 - sfctheta)*2
    #lapse100 = (theta100 - sfctheta)*10

    #Determine the characteristics of the 1km - mixed layer parcel
    k = 0
    i = 0
    for x in range(0,len(hghtlist)-1,1):
        thlist.append(tlist[x] * math.pow(po/plist[x],Rd/Cp))
    while plist[0]-plist[i] <= mldepth:
      tempsum += thlist[i]
      qsum += qlist[i]
      k +=1
      i +=1
        
    theta_i = (tempsum/k)
    p_i= plist[0]

    q_i= qsum/k
    
    #Find theta of starting parcel
    #theta_i = T_i * math.pow(po/(p_i), Rd/Cp)
    T_i = theta_i * math.pow(p_i/po,Rd/Cp)
    #Find the height of saturation for the start parcel
    press = p_i
    temp = T_i
    satht = 0

    i = 1
    while w_s(press,temp) > q_i :
        temp = temp - 0.0098*hincrement
        press = po * math.pow(temp/theta_i,Cp/Rd)
        satht += hincrement
        #press = plist[i]
        #temp = theta_i * math.pow(plist[i]/po,Rd/Cp)
        #satht = hghtlist[i]
        #i += 1
    P_sat = press
    T_sat = temp
    h_sat = satht
    Th_sat = T_sat * math.pow(po/P_sat,Rd/Cp)

    #Find the theta_e of the starting parcel
    Theta_e = thetae.compute(Th_sat,q_i,T_sat,P_sat)    

    #Equation for determining parcel temp from thetae
    def fcn(q, p, T_env,thetae_p,Tl):
        Cp=1004.0
        L0 = 2563130
        L1 = 1754
        K2 = 1137000
        po = 1000.0
        T0 = 273.15
        po = 1000.0
        theta = T_env * math.pow( po / p , Rd/Cp)
        x = thetae_p - theta * math.exp((L0 - L1 * (Tl - T0) + K2 * q) * q / Cp / Tl)
        #x = thetae_p - theta * math.pow(po/p, 0.2854*(1.0-0.28*q)) * math.exp((3.376/T-0.00254)*q*math.pow(10.0,3.0)*(1.0+0.81*q))
        return x
        
    #Create a list of p, t, td, and q above the saturation point
    for i in range(0,len(plist), 1):
        if plist[i] <= P_sat:
            p_env.append(plist[i])
            t_env.append(tlist[i])
            q_env.append(qlist[i])
            h_env.append(hghtlist[i])
    #Calculate the temperature of the parcel path upwards 
    guess1 = T_sat - 1.0
    guess2 = T_sat
    j=0
    adiabat = math.log(Th_sat) + L * q_i / (Cp * T_sat)
    for i in range(0,len(p_env),1):
        tpath.append(Bisection.method(fcn,guess1, guess2,Theta_e, p_env[i],w_s,T_sat))
        guess1 = tpath[j] - 1.0
        guess2 = tpath[j] 
        j+=1

    #Apply virtual temperature correction to unsaturated parcel path and environment
    for i in range(0,len(t_env),1):
        tpath_corr.append(virtualtempcorr_path(p_env[i],tpath[i],tpath[i]))
        t_env_corr.append(virtualtempcorr_env(p_env[i],t_env[i],q_env[i]))

    #for i in range(0,len(t_env),1):
    #  tpath_corr[i] = tpath_corr[i] * math.pow(po/p_env[i],Rd/Cp)
    #  t_env_corr[i] = t_env_corr[i] * math.pow(po/p_env[i],Rd/Cp)
    
    #Calculate CAPE 
    for i in range(1,len(p_env),1):
        #CAPE
        if t_env_corr[i] <= tpath_corr[i]: 
            mlcape += g*((tpath_corr[i]-t_env_corr[i])/t_env_corr[i])*(h_env[i]-h_env[i-1])
            mlel = h_env[i]
        mltotalcape += g*((tpath_corr[i]-t_env_corr[i])/t_env_corr[i])*(h_env[i]-h_env[i-1])
        if p_env[i] < mlel and mltotalcape > 0 :
            mlct= h_env[i]/1000.0
        #CIN
        if t_env_corr[i] > tpath_corr[i] and mlcape == 0.0:
            mlcin += g*((tpath_corr[i]-t_env_corr[i])/t_env_corr[i])*(h_env[i]-h_env[i-1])
            mllfc = h_env[i]
            t_lfc = t_env_corr[i]
            
        #Calculate 0-3km CAPE
        if t_env_corr[i] <= tpath_corr[i] and h_env[i] <= (3000.0) :
            cape3km += g*((tpath_corr[i]-t_env_corr[i])/t_env_corr[i])*(h_env[i]-h_env[i-1])
    
    if mllfc == 0:
        mllfc = h_sat
        t_lfc = T_sat

    print 'mlcin', mlcin
    print 'mllcl', h_sat
    print 'mllfc', mllfc
    print 'mlel', mlel
    print 'mlcape', mlcape

    return mlcape, mlcin, h_sat, mllfc, mlel, mlct, lapse1000, T_sat, t_lfc, cape3km, q_i, tpath_corr,p_env,h_env,qlist[0]
    
def sbcape(hghtlist,presslist,tlist,qlist):

    eps = 0.622
    Cp=1004.0
    Rd=287.04
    g=9.81
    po=1000.0
    L = 2501000
    hincrement = 1
    p_unsat = []
    t_unsat = []
    td_unsat = []
    q_unsat = []
    h_unsat = []
    tpath_corr = []
    t_unsat_corr = []
    tpath = []
    
    sbcin = 0
    sbcape = 0
    sblfc = 0
    sblcl = 0
    sbel = 0
    
    #print presslist[0:5]
    #print tdlist[0:5]
    
    #qlist=q(presslist,tdlist)
    
    #Assign attributes of the starting surface-based parcel
    p_i = presslist[0]
    h_i = hghtlist[0]
    T_i = tlist[0]
    #Td_i = tdlist[0]
    q_i = qlist[0]
    
    #Find theta of starting parcel
    theta_i = T_i * math.pow(po/p_i, Rd/Cp)

    #Find the height of saturation for the start parcel
    press=p_i
    temp = T_i
    satht = h_i
    
    while w_s(press,temp) > q_i :
        satht += hincrement
        temp -= 0.0098*hincrement
        press = 1000.0*math.pow(temp/theta_i,Cp/Rd)
    P_sat = press
    T_sat = temp

    # Find the theta_e of the starting parcel
    Theta_e = thetae.compute(theta_i,q_i,T_sat,P_sat)

    #Equation for determining parcel temp from thetae
    def fcn(q, p, T_env,thetae_p,Tl):
        Cp=1004.0
        L0 = 2563130
        L1 = 1754
        K2 = 1137000
        po = 1000.0
        T0 = 273.15
        po = 1000.0
        theta = T_env * math.pow( po / p , Rd/Cp)
        x = thetae_p - theta * math.exp((L0 - L1 * (Tl - T0) + K2 * q) * q / Cp / Tl)
        return x
        
    #Create a list of p, t, td, and q above the saturation point
    for i in range(0,len(presslist), 1):
        if presslist[i] <= P_sat:
            p_unsat.append(presslist[i])
            t_unsat.append(tlist[i])
            #td_unsat.append(tdlist[i])
            q_unsat.append(qlist[i])
            h_unsat.append(hghtlist[i])

    #Calculate the temperature of the parcel path upwards 
    guess1 = T_sat - 0.25
    guess2 = T_sat + 0.25
    j=0
    adiabat = math.log(theta_i) + L * q_i / (Cp * T_sat)

    for i in range(0,len(p_unsat),1):
        tpath.append(Bisection.method(fcn,guess1, guess2,Theta_e, p_unsat[i],w_s,T_sat))
        guess1 = tpath[j] - 0.5
        guess2 = tpath[j] + 0.5 
        j+=1

    #Apply virtual temperature correction to unsaturated parcel path and environment
    for i in range(0,len(t_unsat),1):
        tpath_corr.append(virtualtempcorr_path(p_unsat[i],tpath[i],tpath[i]))
        t_unsat_corr.append(virtualtempcorr_env(p_unsat[i],t_unsat[i],q_unsat[i]))

    
     #Calculate sbCAPE 
    for i in range(0,len(p_unsat)-1,1):
        #CAPE
        if t_unsat_corr[i] <= tpath_corr[i]: 
            sbcape += g*((tpath_corr[i]-t_unsat_corr[i])/t_unsat_corr[i])*(h_unsat[i+1]-h_unsat[i])
            sbel = p_unsat[i]
        #CIN
        if t_unsat_corr[i] > tpath_corr[i] and mlcape == 0.0:
            sbcin += g*((tpath_corr[i]-t_unsat_corr[i])/t_unsat_corr[i])*(h_unsat[i+1]-h_unsat[i])
            sblfc = p_unsat[i]
    
    if sblfc == 0:
        sblfc = P_sat
    #print tpath_corr
    print 'sbcin', sbcin
    print 'sblcl', satht
    print 'sblfc', sblfc
    print 'sbel',sbel
    print 'sbcape', sbcape
    
    return sbcape, sbcin, P_sat, sblfc, sbel, tpath_corr
    
def dcape(hghtlist,presslist,tlist,qlist,startp) :
    
    
    eps = 0.622
    Cp=1004.0
    Rd=287.04
    g=9.81
    po=1000.0
    hincrement = 10
    L = 2501000
    p_unsat = []
    t_unsat = []
    td_unsat = []
    q_unsat = []
    h_unsat = []
    tpath_corr = []
    t_unsat_corr = []
    tpath = []
    
    DCAPE = 0
    
    #qlist=q(presslist,tdlist)
    
    T_i = 99999
    
    #Set up variables for when the user asks for the lowest thetae in a certain depth
    if startp < 0:
        P_sat = []
        T_sat = []
        theta = []
        
        p_sat = []
        t_sat = []
        Theta_e =[]
        
        #Calculate theta for all temperatures and pressures
        for i in range(0,len(presslist),1):
            theta.append(tlist[i] * math.pow(po/presslist[i], Rd/Cp))
            
        #Calculate saturation temperature, height, and pressure for all points
        for i in range(0,len(tlist),1):
            press=presslist[i]
            temp = tlist[i]
            satht = hghtlist[i]
            theta_i = theta[i]
            q_i = qlist[i]
            
            while w_s(press,temp) > q_i :
                temp -= 0.0098*hincrement
                press = 1000.0*math.pow(temp/theta_i,Cp/Rd)
                satht += hincrement
            p_sat.append(press)
            t_sat.append(temp)
        
        #Calculate thetae for all temperatures and pressures
        for i in range(0,len(presslist),1):
            theta_temp = tlist[i] * math.pow(po/presslist[i],Rd/Cp)
            Theta_e.append(thetae.compute(theta_temp,qlist[i],t_sat[i],p_sat[i]))
        
        #Find the lowest thetae in the specified depth
        i=1
        thetaelow = Theta_e[0]
        pressstart = presslist[0]
        while presslist[i] >= presslist[0]-np.absolute(startp):
            if Theta_e[i] < thetaelow:
                thetaelow = Theta_e[i]
                T_i = tlist[i]
                q_i = qlist[i]
                h_i = hghtlist[i]
                P_sat = p_sat[i]
                T_sat = t_sat[i]
                theta_i = theta[i]
                pressstart = presslist[i]
            i += 1
        Theta_e = thetaelow
        startp = pressstart

            
    #Find variablesfor when the user asks for DCAPE from a specific pressure level
    if startp > 0:
        #Find the initial temperature and dewpoint to start at
        
        for i in range(0,len(hghtlist),1):
            if presslist[i] == startp:
                T_i = tlist[i]
                #Td_i = tdlist[i]
                q_i = qlist[i]
                h_i =  hghtlist[i]
       
        if T_i == 99999:
            T_i = lininterp.interp(tlist,presslist,startp)
            #Td_i = lininterp.interp(tdlist,presslist,startp)
            q_i = lininterp.interp(qlist,presslist,startp)
            h_i = lininterp.interp(hghtlist,presslist,startp)
            
        #Find theta of the parcel at the start height
        theta_i = T_i * math.pow(po/startp, Rd/Cp)
       
        #Find the height of saturation for the start parcel
        press=startp
        temp = T_i
        satht = h_i
        
        while w_s(press,temp) > q_i :
            temp -= 0.0098*hincrement
            press = 1000.0*math.pow(temp/theta_i,Cp/Rd)
            satht += hincrement
        P_sat = press
        T_sat = temp

    # Find the theta_e of the starting parcel
    Theta_e = thetae.compute(theta_i,q_i,T_sat,P_sat)
    

    #Equation for determining parcel temp from thetae
    def fcn(q, p, T_env,thetae_p,Tl):
        Cp=1004.0
        L0 = 2563130
        L1 = 1754
        K2 = 1137000
        po = 1000.0
        T0 = 273.15
        po = 1000.0
        theta = T_env * math.pow( po / p , Rd/Cp)
        x = thetae_p - theta * math.exp((L0 - L1 * (Tl - T0) + K2 * q) * q / Cp / Tl)
        return x
        
    #Create a list of p, t, td, and q below the saturation point
    for i in range(len(presslist)-1,-1, -1):
        if presslist[i] >= P_sat:
            p_unsat.append(presslist[i])
            t_unsat.append(tlist[i])
            #td_unsat.append(tdlist[i])
            q_unsat.append(qlist[i])
            h_unsat.append(hghtlist[i])

    #Calculate the temperature of the parcel path downwards 
    guess1 = T_sat + 50.0
    guess2 = T_sat - 50.0
    j=0
    adiabat = math.log(theta_i) + L * q_i / (Cp * T_sat)
    for i in range(0,len(p_unsat),1):
        tpath.append(Bisection.method(fcn,guess1, guess2,Theta_e, p_unsat[i],w_s,T_sat))
        guess1 = tpath[j] + 50.0
        guess2 = guess1 - 50.0 
        j+=1
    
    #Apply virtual temperature correction to unsaturated parcel path and environment
    for i in range(0,len(t_unsat),1):
        tpath_corr.append(virtualtempcorr_path(p_unsat[i],tpath[i],tpath[i]))
        t_unsat_corr.append(virtualtempcorr_env(p_unsat[i],t_unsat[i],q_unsat[i]))
     
    t_dcape = tpath_corr[-1]
    
    rfd_diff = tlist[0] - t_dcape
    #Calculate DCAPE 
    for i in range(0,len(p_unsat)-1,1):
        if p_unsat[i] >= startp: 
            DCAPE += (-g*((tpath_corr[i]-t_unsat_corr[i])/t_unsat_corr[i])*(h_unsat[i+1]-h_unsat[i]))
    
    
    return DCAPE, rfd_diff
    
def mucape(hghtlist,presslist,tlist,qlist):

    eps = 0.622
    Cp=1004.0
    Rd=287.04
    g=9.81
    po=1000.0
    L = 2501000
    hincrement = 5

    cin = []
    cape = []
    lfc = []
    lcl = []
    el = []
    P_sat = []

    #print presslist[0:5]
    #print tdlist[0:5]
    i = 0
    while presslist[i] > (presslist[0]-300) :
        mutop = i
        i +=1
    
    #qlist=q(presslist,tdlist)
    for starth in range(0,len(presslist[0:mutop]),1):
        p_unsat = []
        t_unsat = []
        td_unsat = []
        q_unsat = []
        h_unsat = []
        tpath_corr = []
        t_unsat_corr = []
        tpath = []
        cin_pre = 0
        cape_pre = 0
        lfc_pre = 0
        lcl_pre = 0
        el_pre = 0
        
        #Assign attributes of the starting surface-based parcel
        p_i = presslist[starth]
        h_i = hghtlist[starth]
        T_i = tlist[starth]
        #Td_i = tdlist[0]
        q_i = qlist[starth]
    
        #Find theta of starting parcel
        theta_i = T_i * math.pow(po/p_i, Rd/Cp)

        #Find the height of saturation for the start parcel
        press=p_i
        temp = T_i
        satht = h_i
    
        while w_s(press,temp) > q_i :
            satht += hincrement
            temp -= 0.0098*hincrement
            press = 1000.0*math.pow(temp/theta_i,Cp/Rd)
        P_sat = press
        T_sat = temp

    # Find the theta_e of the starting parcel
        Theta_e = thetae.compute(theta_i,q_i,T_sat,P_sat)

        #Equation for determining parcel temp from thetae
        def fcn(q, p, T_env,thetae_p,Tl):
            Cp=1004.0
            L0 = 2563130
            L1 = 1754
            K2 = 1137000
            po = 1000.0
            T0 = 273.15
            po = 1000.0
            theta = T_env * math.pow( po / p , Rd/Cp)
            x = thetae_p - theta * math.exp((L0 - L1 * (Tl - T0) + K2 * q) * q / Cp / Tl)
            return x
        
        #Create a list of p, t, td, and q above the saturation point
        for i in range(starth,len(presslist), 1):
            if presslist[i] <= P_sat:
                p_unsat.append(presslist[i])
                t_unsat.append(tlist[i])
                #td_unsat.append(tdlist[i])
                q_unsat.append(qlist[i])
                h_unsat.append(hghtlist[i])

        #Calculate the temperature of the parcel path upwards 
        guess1 = T_sat - 0.25
        guess2 = T_sat + 0.25
        j=0
        adiabat = math.log(theta_i) + L * q_i / (Cp * T_sat)

        for i in range(0,len(p_unsat),1):
            tpath.append(Bisection.method(fcn,guess1, guess2,Theta_e, p_unsat[i],w_s,T_sat))
            guess1 = tpath[j] - 0.5
            guess2 = tpath[j] + 0.5 
            j+=1

        #Apply virtual temperature correction to unsaturated parcel path and environment
        for i in range(0,len(t_unsat),1):
            tpath_corr.append(virtualtempcorr_path(p_unsat[i],tpath[i],tpath[i]))
            t_unsat_corr.append(virtualtempcorr_env(p_unsat[i],t_unsat[i],q_unsat[i]))

    
         #Calculate sbCAPE 
        for i in range(0,len(p_unsat)-1,1):
            #CAPE
            if t_unsat_corr[i] <= tpath_corr[i]: 
                cape_pre += g*((tpath_corr[i]-t_unsat_corr[i])/t_unsat_corr[i])*(h_unsat[i+1]-h_unsat[i])
                el_pre = h_unsat[i]
            #CIN
            if t_unsat_corr[i] > tpath_corr[i] and mlcape == 0.0:
                cin_pre += g*((tpath_corr[i]-t_unsat_corr[i])/t_unsat_corr[i])*(h_unsat[i+1]-h_unsat[i])
                lfc_pre = h_unsat[i]
    
        if lfc_pre == 0:
            lfc_pre = P_sat
        cape.append(cape_pre)
        el.append(el_pre)
        cin.append(cin_pre)
        lfc.append(lfc_pre)
        if max(lfc) == 0:
            mulfc = P_sat[starth]
        #print tpath_corr

    mu_loc, mucape = max(enumerate(cape), key=operator.itemgetter(1))
    mucin = cin[mu_loc]
    mulfc = lfc[mu_loc]
    muel = el[mu_loc]
    mulcl = h_unsat[mu_loc]
    print 'mucape', mucape, 'muel', muel
        
    return mucape, mucin, mulcl, mulfc, muel, tpath_corr, hghtlist[mu_loc]
    
def efflayer(hghtlist,presslist,tlist,qlist,uwind,vwind):

    eps = 0.622
    Cp=1004.0
    Rd=287.04
    g=9.81
    po=1000.0
    L = 2501000
    hincrement = 5
    
    cin = []
    cape = []
    lfc = []
    lcl = []
    el = []
    P_sat = []
    loc_in_snd = 0
    eff_bot = 0
    eff_top = 0
    #print presslist[0:5]
    #print tdlist[0:5]
    cape_pre = 0
    cin_pre = -300
    i = 0
    # Find the lowest parcel with >100 CAPE and >=-250 CIN
    starth = 0
    loc_in_snd = 0
    
    while cape_pre <= 100 or cin_pre <= -250 :
        p_unsat = []
        t_unsat = []
        td_unsat = []
        q_unsat = []
        h_unsat = []
        tpath_corr = []
        t_unsat_corr = []
        tpath = []
        cin_pre = 0
        cape_pre = 0
        lfc_pre = 0
        lcl_pre = 0
        el_pre = 0
        
        
        #Assign attributes of the starting parcel
        p_i = presslist[starth]
        h_i = hghtlist[starth]
        T_i = tlist[starth]
        #Td_i = tdlist[starth]
        q_i = qlist[starth]
        
        
    
        #Find theta of starting parcel
        theta_i = T_i * math.pow(po/p_i, Rd/Cp)

        #Find the height of saturation for the start parcel
        press=p_i
        temp = T_i
        satht = h_i
    
        while w_s(press,temp) > q_i :
            temp -= 0.0098*hincrement
            press = 1000.0*math.pow(temp/theta_i,Cp/Rd)
            satht += hincrement
        P_sat.append(press)
        P_sat_temp = press
        T_sat = temp
    
        # Find the theta_e of the starting parcel
        Theta_e = thetae.compute(theta_i,q_i,T_sat,P_sat_temp)
    

        #Equation for determining parcel temp from thetae
        def fcn(q, p, T_env,thetae_p,Tl):
          Cp=1004.0
          L0 = 2563130
          L1 = 1754
          K2 = 1137000
          po = 1000.0
          T0 = 273.15
          po = 1000.0
          theta = T_env * math.pow( po / p , Rd/Cp)
          x = thetae_p - theta * math.exp((L0 - L1 * (Tl - T0) + K2 * q) * q / Cp / Tl)
          return x
        
        #Create a list of p, t, td, and q above the saturation point
        for i in range(0,len(presslist), 1):
            if presslist[i] <= P_sat[starth]:
                p_unsat.append(presslist[i])
                t_unsat.append(tlist[i])
                #td_unsat.append(tdlist[i])
                q_unsat.append(qlist[i])
                h_unsat.append(hghtlist[i])
                
            
        #Calculate the temperature of the parcel path upwards 
        guess1 = T_sat - 20.0
        guess2 = T_sat + 20.0
        j=0
        adiabat = math.log(theta_i) + L * q_i / (Cp * T_sat)
        for i in range(0,len(p_unsat),1):
            tpath.append(Bisection.method(fcn,guess1, guess2,Theta_e, p_unsat[i],w_s,T_sat))
            guess1 = tpath[j] - 20.0
            guess2 = guess1 + 20.0 
            j+=1
    
        #Apply virtual temperature correction to unsaturated parcel path and environment
        for i in range(0,len(t_unsat),1):
            tpath_corr.append(virtualtempcorr_path(p_unsat[i],tpath[i],tpath[i]))
            t_unsat_corr.append(virtualtempcorr_env(p_unsat[i],t_unsat[i],q_unsat[i]))
    
         #Calculate CAPE 
        for i in range(0,len(p_unsat)-1,1):
            #CAPE
            if t_unsat_corr[i] <= tpath_corr[i]: 
                cape_pre += g*((tpath_corr[i]-t_unsat_corr[i])/t_unsat_corr[i])*(h_unsat[i+1]-h_unsat[i])
                
                el_pre = p_unsat[i]
            #CIN
            if t_unsat_corr[i] > tpath_corr[i] and cape_pre == 0:
                cin_pre += g*((tpath_corr[i]-t_unsat_corr[i])/t_unsat_corr[i])*(h_unsat[i+1]-h_unsat[i])
                lfc_pre = p_unsat[i]
        print cape_pre, cin_pre
        cape.append(cape_pre)
        el.append(el_pre)
        cin.append(cin_pre)
        lfc.append(lfc_pre)
        eff_bot = hghtlist[starth]
        eff_bot_loc = starth
        starth += 1
        if starth > len(presslist)-1 :
          break
        
        
    #Find highest parcel with >=100 CAPE and >= 250 CIN
    while cape_pre >= 100 and cin_pre >= -250 :
        p_unsat = []
        t_unsat = []
        td_unsat = []
        q_unsat = []
        h_unsat = []
        tpath_corr = []
        t_unsat_corr = []
        tpath = []
        cin_pre = 0
        cape_pre = 0
        lfc_pre = 0
        lcl_pre = 0
        el_pre = 0

        
        #Assign attributes of the starting parcel
        p_i = presslist[starth]
        h_i = hghtlist[starth]
        T_i = tlist[starth]
        #Td_i = tdlist[starth]
        q_i = qlist[starth]
        
        
    
        #Find theta of starting parcel
        theta_i = T_i * math.pow(po/p_i, Rd/Cp)

        #Find the height of saturation for the start parcel
        press=p_i
        temp = T_i
        satht = h_i
    
        while w_s(press,temp) > q_i :
            temp -= 0.0098*hincrement
            press = 1000.0*math.pow(temp/theta_i,Cp/Rd)
            satht += hincrement
        P_sat.append(press)
        P_sat_temp = press
        T_sat = temp
    
        # Find the theta_e of the starting parcel
        Theta_e = thetae.compute(theta_i,q_i,T_sat,P_sat_temp)
    

        #Equation for determining parcel temp from thetae
        def fcn(q, p, T_env,thetae_p,Tl):
          Cp=1004.0
          L0 = 2563130
          L1 = 1754
          K2 = 1137000
          po = 1000.0
          T0 = 273.15
          po = 1000.0
          theta = T_env * math.pow( po / p , Rd/Cp)
          x = thetae_p - theta * math.exp((L0 - L1 * (Tl - T0) + K2 * q) * q / Cp / Tl)
          return x
        
        #Create a list of p, t, td, and q above the saturation point
        for i in range(0,len(presslist), 1):
            if presslist[i] <= P_sat[starth]:
                p_unsat.append(presslist[i])
                t_unsat.append(tlist[i])
                #td_unsat.append(tdlist[i])
                q_unsat.append(qlist[i])
                h_unsat.append(hghtlist[i])
                
            
        #Calculate the temperature of the parcel path upwards 
        guess1 = T_sat - 20.0
        guess2 = T_sat + 20.0
        j=0
        adiabat = math.log(theta_i) + L * q_i / (Cp * T_sat)
        for i in range(0,len(p_unsat),1):
            tpath.append(Bisection.method(fcn,guess1, guess2,Theta_e, p_unsat[i],w_s,T_sat))
            guess1 = tpath[j] - 20.0
            guess2 = guess1 + 20.0 
            j+=1
    
        #Apply virtual temperature correction to unsaturated parcel path and environment
        for i in range(0,len(t_unsat),1):
            tpath_corr.append(virtualtempcorr_path(p_unsat[i],tpath[i],tpath[i]))
            t_unsat_corr.append(virtualtempcorr_env(p_unsat[i],t_unsat[i],q_unsat[i]))
    
         #Calculate CAPE 
        for i in range(0,len(p_unsat)-1,1):
            #CAPE
            if t_unsat_corr[i] <= tpath_corr[i]: 
                cape_pre += g*((tpath_corr[i]-t_unsat_corr[i])/t_unsat_corr[i])*(h_unsat[i+1]-h_unsat[i])
                el_pre = p_unsat[i]
            #CIN
            if t_unsat_corr[i] > tpath_corr[i] and cape_pre == 0:
                cin_pre += g*((tpath_corr[i]-t_unsat_corr[i])/t_unsat_corr[i])*(h_unsat[i+1]-h_unsat[i])
                lfc_pre = p_unsat[i]

        cape.append(cape_pre)
        el.append(el_pre)
        cin.append(cin_pre)
        lfc.append(lfc_pre)
        eff_top = hghtlist[starth-1]
        eff_top_loc = starth
        starth += 1
        #print tpath_corr
    if eff_top == 0 :
      eff_bot = -999
      eff_top = -999
      eff_bot_loc = 0
      eff_top_loc = 0

    
    print 'eff_bot', eff_bot, 'eff_top', eff_top
    return eff_bot, eff_top, eff_bot_loc,eff_top_loc
    

    
    
    
        
    
    
    
    
    
    
    
    
    
    
    
    