import numpy as np

def uvconvert(u,v):
    direction=57.29578*(np.atan(u,v))+180 
    speed=np.sqrt(u*u+v*v) 
    return speed, direction