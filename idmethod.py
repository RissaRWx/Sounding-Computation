#!/usr/bin/env python
#These functions use another python class call "matfun.py" which has matrix and vector operations.  Available at http://users.rcn.com/python/download/python.htm

#mean wind for a layer AGL
#this algorith uses lists of wind direction and speed, heights and the top and
#bottom of the layer
#the key to finding the mean wind of a layer to move the coordinate system
#where the bottom of the layer is the new origin and the x-axis for the new
#coordinate system is the shear vector of the layer
#While this looks very nice graphically, it's pretty nasty mathematically 

import math
import matfunc
from lininterp_wind import lininterp_wind

def meanwind(ulist,vlist,hghtlist,blayer,tlayer):
  vector = (ulist,vlist)
  tailvect=[]
  headvect=[]
  i = 0
  
  #find the bottom and the top of the layer in the lists
  for i in range(0,len(hghtlist),1):
    if hghtlist[i] == blayer:
      tailvect = [vector[0][i],vector[1][i]]
    elif hghtlist[i] == tlayer:
      headvect = [vector[0][i],vector[1][i]]
    elif hghtlist[i] > tlayer: break
    else: pass
  
  if tailvect==[]:
    tailvect=lininterp_wind(ulist,vlist,hghtlist,blayer)
   
  if headvect==[]:
    headvect=lininterp_wind(ulist,vlist,hghtlist,tlayer)
  
  
 

  #move coordinate system to the tail of the vector [tail becomes origin]
  modheadvect = [headvect[0] - tailvect[0],headvect[1] - tailvect[1]]
  #find the angle between shear vector and x-axis
  angle = math.acos(modheadvect[0] / math.sqrt(modheadvect[0] ** 2.0 + modheadvect[1] ** 2.0))  
  if modheadvect[1] < 0:
    angle = -angle
  rotation = matfunc.Mat([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])
  counter = matfunc.Mat([[math.cos(angle),math.sin(angle)],[-math.sin(angle),math.cos(angle)]])
  i = 0
  j = 0
  vmod = 0
  #rotate the hodograph so that the new x-axis is the shear vector
  #between the top and bottom level of the layer for the mean wind
  #see the comet presentation on hodographs for graphical explanation
  #Section 5.5 at http://deved.comet.ucar.edu/mesoprim/hodograf/print.htm
  for i in range(0,len(hghtlist),1):
    if hghtlist[i] >= blayer and hghtlist[i] <= tlayer:
      #modifty each wind ob to set the tail as origin
      vect = [vector[0][i] - tailvect[0],vector[1][i] - tailvect[1]]
      vectormat = matfunc.Mat([[vect[0]],[vect[1]]])
      #rotate to complete transition to new coordinate system 
      a = rotation.mmul(vectormat)
      vmod += a[1][0]
      j += 1
    elif hghtlist[i] > tlayer: break
    else: pass
    
  #make the head vector a matrix for upcoming math
  headmat = matfunc.Mat([[modheadvect[0]],[modheadvect[1]]])
  #rotate the head vector to make it the next x axis
  modheadmat = rotation.mmul(headmat)
  #to do the next part, it is assumed the wind measurements are evenly
  #distributed---this assumption MUST be met
  #the mean u wind component is just half of the difference of the bottom
  #level u-component and top level u-component 
  umodmean = modheadmat[0][0] / 2.0
  #the mean v wind component is an average of the v obs in coordinate system
  vmodmean = vmod / j
  #make mean wind matrix for upcoming math
  modmat = matfunc.Mat([[umodmean],[vmodmean]])
  #rotate back and move to original origin
  meanmat = counter.mmul(modmat)
  meanwind = [meanmat[0][0] + tailvect[0],meanmat[1][0] + tailvect[1]]
  return meanwind

#The method takes in a wind direction and speed list and a height list to use 
#the ID method described in Bunkers et al 2000 for supercell motion estimate
def bunkersmotion(ulist,vlist,hghtlist,blayer,tlayer):
  #find the mean 0-6km wind
  mean = meanwind(ulist,vlist,hghtlist,blayer,6000)
  #find the upper 500 m wind; this is the head of the shear vector
  head = meanwind(ulist,vlist,hghtlist,5500,6000)
  #find the lower 500 m wind; this is the tail of the shear vector
  tail = meanwind(ulist,vlist,hghtlist,blayer,500)
  #deviation (as described in Bunkers et al 2000)
  D = 7.5
  #u-component of the shear vector
  shearu = head[0] - tail[0]
  #v-component of shear vector
  shearv = head[1] - tail[1]
# shearmag = sqrt(shearu ** 2 + shearv ** 2)
# make shear vector a vector in python
  shear = matfunc.Vec([shearu,shearv,0])
#  rotate = Mat([[cos(-90. * degree),-sin(-90. * degree)],[sin(-90. * degree),cos(90. * degree)]])
#  shearterm = rotate.mmul(unitshear)
  #z-component unit matrix for use
  unitmat = matfunc.Vec([0,0,1])
  #find k cross shear; Bunkers et al 2000, Eqns 1 & 2
  shearterm = shear.cross(unitmat)
  #Bunkers et al 2000, Eqn 1, u-component right mover
  #Basically what is done is the mean wind in the 0-6 km layer has a deviation,
  #that is weighted by the shear vector magnitude and direction, added to it

  meanu = mean[0] + D * (shearterm[0] / shear.norm())
  #v-component of right mover
  meanv = mean[1] + D * (shearterm[1] / shear.norm())
  #Bunkers et al 2000, Eqn 2, u-compnent left mover
  meanleftu = mean[0] - D * (shearterm[0] / shear.norm())
  #Bunkers et al 2000, Eqn 2, v-component left mover
  meanleftv = mean[1] - D * (shearterm[1] / shear.norm())
  #return vector motion components to user
  motion = [meanu,meanv,meanleftu,meanleftv]
  return meanu, meanv, meanleftu, meanleftv