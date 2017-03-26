#! /usr/bin/env python

def interp(varlist, hlist, targetlvl):
    i = 0
    
    while hlist[i] > targetlvl:
        i+=1
    dist = hlist[i]-hlist[i-1]
    
    interpvar = (varlist[i]-varlist[i-1])/dist*(targetlvl-hlist[i-1])+varlist[i-1]
    
            
    return interpvar
            