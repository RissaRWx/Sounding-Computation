#! /usr/bin/env python

import math

def compute(theta, q, T,p):
    Cp=1004.0
    L0 = 2563130
    L1 = 1754
    K2 = 1137000
    po = 1000.0
    T0 = 273.15
    #po = 1000.0
    #thetae = theta * math.pow(po/p, 0.2854*(1.0-0.28*q)) * math.exp((3.376/T-0.00254)*q*math.pow(10.0,3.0)*(1.0+0.81*q))
    #thetae = theta * math.exp((3036/T-1.78)*q*(1+0.448 * q))
    #thetae = theta * math.exp(Lv*q/Cp/T)
    thetae = theta * math.exp((L0 - L1 * (T - T0) + K2 * q) * q / Cp / T)
    return thetae
