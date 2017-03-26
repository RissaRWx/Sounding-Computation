#!/usr/bin/env python

import decimal
import math
        
def method(fcn,guess1, guess2,thetae, p, w_s,T_sat):
    Cp=1004.0
    Rd=287.04
    po=1000.0
    while 0.5*math.fabs(guess1-guess2) > 3.0*math.pow(10.0,-7.0):
        x3 = (guess1+guess2)/2.0
        q1 = w_s(p,guess1)
        q2 = w_s(p,guess2)
        if fcn(q1, p, guess1,thetae,T_sat)*fcn(q2,p,guess2,thetae,T_sat) < 0:
            guess1=x3
        else:
            guess2=x3
    return x3


    