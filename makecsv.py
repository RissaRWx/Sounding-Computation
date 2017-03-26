#from interpolate import interpolater
from interpolate2 import interpolater
import csv
import sys
import os

def make(zlist, presslist,tlist,tdlist,ulist,vlist, filename, sfcp, resolution, foldername,tod):

    newpress=[]
    newtemp=[]
    newtd=[]
    newulist=[]
    newvlist=[]
    finalname=[]
    capedistname=[]
    
    

    newname=list(filename)
    newname.pop(-4)
    newname.pop(-3)
    newname.pop(-2)
    newname.pop(-1)
    newname.insert(0,'N')
    finalname.extend(newname)
    capedistname.extend(newname)
    newname.append('.')
    newname.append('c')
    newname.append('s')
    newname.append('v')

    finalname.append('v')
    finalname.append('a')
    finalname.append('l')
    finalname.append('u')
    finalname.append('e')
    finalname.append('s')
    finalname.append('.')
    finalname.append('c')
    finalname.append('s')
    finalname.append('v')
    
    capedistname.append('c')
    capedistname.append('a')
    capedistname.append('p')
    capedistname.append('e')
    capedistname.append('.')
    capedistname.append('c')
    capedistname.append('s')
    capedistname.append('v')    

    valuesfile="".join(finalname)
    interpfile="".join(newname)
    capefile="".join(capedistname)
    
    root = ''.join(list(os.environ['HOME']))

    finalfile = root + '/Dropbox/'+ foldername + '/Interpolated_Soundings/' + interpfile
    #finalfile = root + '/' + foldername + '/Interpolated_Soundings/' + interpfile
    
    interpresults = interpolater(zlist,presslist,[tlist, tdlist, ulist, vlist],sfcp,resolution,foldername, tod)

    
    pinterp = interpresults[0]
    hinterp = interpresults[1]
    tinterp = interpresults[2]
    tdinterp = interpresults[3]
    uinterp = interpresults[4]
    vinterp = interpresults[5]

    
    #newpress.extend(interpolater(hghtlist,tlist, elevation,resolution,'temperature')[0])
    #newtemp.extend(interpolater(hghtlist,tlist, elevation,resolution,'temperature')[2])
    #newpress.extend(interpolater(hghtlist,tlist, elevation,resolution,'temperature')[1])
    #newtd.extend(interpolater(hghtlist,tdlist, elevation,resolution,'dewpoint')[1])
    #newulist.extend(interpolater(hghtlist,ulist, elevation,resolution,'wind')[1])
    #newvlist.extend(interpolater(hghtlist,vlist, elevation,resolution,'wind')[1])

    
    result=open(finalfile,'w')
    write=csv.writer(result, delimiter=',')

    #write.writerow(newheight)
    #write.writerow(newpress)
    #write.writerow(newtemp)
    #write.writerow(newtd)
    #write.writerow(newulist)
    #write.writerow(newvlist)
     
    write.writerow(hinterp)
    write.writerow(pinterp)
    write.writerow(tinterp)
    write.writerow(tdinterp)
    write.writerow(uinterp)
    write.writerow(vinterp)

    result.close()

    return interpfile, valuesfile, capefile








