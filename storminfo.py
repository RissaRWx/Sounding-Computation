#!/usr/bin/env python

import sys
import math
import csv
from idmethod import meanwind, bunkersmotion
from svrparameters import helicity, ehi, bulkshear, superhelicity,srwind, meanwind, lapserate, mean_RH, srwind_mean
import Capevals
import makecsv
import os
import numpy as np
import numpy.ma as ma
import fnmatch
from uvtostheta import uvconvert
from operator import itemgetter
import codecs


folder = sys.argv[1]

root = ''.join(list(os.environ['HOME']))


directory = root + '/Dropbox/'  + folder
soundingpath= directory + '/RawFiles'  # insert the path to the directory of interest here

dirList=os.listdir(soundingpath)
csvpath = fnmatch.filter(dirList,'*.csv')
#if folder == 'MC3E':
#	csvpath = fnmatch.filter(dirList,'201104151658_S02.txt')

resolution = float(sys.argv[2])

for file in csvpath:
	print file
	
	lcl=[]
	blayer=[]
	tlayer=[]
	altitude=[]
	info=[]
	hght=[]
	press=[]
	temp=[]
	dewpt=[]
	uwind=[]
	vwind=[]
	umotion = 0
	vmotion = 0
	h=[]
	p=[]
	t=[]
	td=[]
	z = []
	u=[]
	v=[]
	ws = []
	wd = []
	stdata=[]
	q = []
	TOD = []
	speeds = []
	sr_u = []
	sr_v = []
	mixr = []

	fullpath = soundingpath + '/' + file

	if folder == 'F1F2' or folder == 'F3Plus':
		data = csv.reader(fullpath, delimiter=',')
		z,p,t,q,u,v = np.genfromtxt(fullpath,dtype = float,delimiter=',',skip_header=0,unpack=True,missing_values = '********',filling_values = -99999)
		z = z.tolist()
		p = p.tolist()
		t = t.tolist()
		q=q.tolist()
		u=u.tolist()
		v = v.tolist()	


	if folder == 'F1F2' or folder =='F3Plus': 
		sfcpress = p[0]
		sfct = t[0] 
		tod = []
		path = []

	#if folder == 'Nighttime_Soundings' or folder == 'LTTRUC' or folder == 'New_Sound' or folder=='F1F2' or folder == 'Test':
	#	print folder
	#	csvresults = makecsv.make(p, t, td, u, v, file, sfcpress, resolution, folder, tod)
	#	csvname = directory + '/Interpolated_Soundings/' + csvresults[0]
	
	if folder == 'F1F2' or folder =='F3Plus':   #### RYAN EDIT THESE NAMES
		csvresults = makecsv.make(z,p, t, q, u, v, file, sfcpress, resolution, folder,tod)
		csvname = directory + '/Interpolated_Soundings/' + csvresults[0]
  
	interpsounding=open(csvname,'rU')
	data=csv.reader(interpsounding, delimiter=',')
	for row in data:
		info.append(row)
		

	prepress=info[1]
	prehght=info[0]
	pretemp=info[2]
	preq=info[3]
	preuwind=info[4]
	prevwind=info[5]
	

	for i in range (0,len(info[0]),1):

		hght.append(float(prehght[i]))
		press.append(float(prepress[i]))
		temp.append(float(pretemp[i]))
		mixr.append(float(preq[i]))
		uwind.append(float(preuwind[i]))
		vwind.append(float(prevwind[i]))
		

	mixr[:] = [x / 1000 for x in mixr]
	basealt=altitude
    
	i = 0


	tlayer=hght[-1]

	umotionR =  bunkersmotion(uwind,vwind,hght,0,tlayer)[0]

	vmotionR =  bunkersmotion(uwind,vwind,hght,0,tlayer)[1]

	umotionL = bunkersmotion(uwind,vwind,hght,0,tlayer)[2]

	vmotionL = bunkersmotion(uwind,vwind,hght,0,tlayer)[3]

	umotion = umotionR
	vmotion = vmotionR
	
	def uvconvert(umotion,vmotion):
		if umotion>0 and vmotion>0:
			degrees = 180 + math.degrees(math.atan(umotion/vmotion))
		elif umotion>0 and vmotion<0:
			degrees = 360 - math.degrees(math.atan(umotion/math.fabs(vmotion)))
		elif umotion<0 and vmotion<0:
			degrees = 90 - math.degrees(math.atan(math.fabs(umotion)/math.fabs(vmotion)))
		else:
			degrees = 180 - math.degrees(math.atan(math.fabs(umotion)/vmotion))
		speed = math.sqrt(math.pow(umotion,2) + math.pow(vmotion,2))
		speed *= 2.23693629
		return speed, degrees

	[stspd,stdir] = uvconvert(umotion,vmotion)

    #dur = path/stspd*60
    
	for i in range(0,len(uwind),1):
		sr_u.append(uwind[i]-umotion)
		sr_v.append(vwind[i]-vmotion)

	for i in range(0,len(uwind),1):
		if hght[i] >= 2000.0 and hght[i] <= 9000.0:
			speeds.append(math.sqrt(math.pow(sr_u[i],2) + math.pow(sr_v[i],2)))
	minwind = min(speeds)

	mlcaperesults = Capevals.mlcape(hght, press,temp,mixr)

	RFD_level = sfcpress - 200
	
	sr_comps_50m = srwind(uwind,vwind,hght,umotion,vmotion,50.0)
	wind_50m= np.sqrt(math.pow(sr_comps_50m[0]+umotion,2.0)+math.pow(sr_comps_50m[1]+vmotion,2.0))
	srwind_50m = np.sqrt(math.pow(sr_comps_50m[0],2.0)+math.pow(sr_comps_50m[1],2.0))
    
	sr_comps_100m = srwind(uwind,vwind,hght,umotion,vmotion,100.0)
	wind_100m= np.sqrt(math.pow(sr_comps_100m[0]+umotion,2.0)+math.pow(sr_comps_100m[1]+vmotion,2.0))
	srwind_100m = np.sqrt(math.pow(sr_comps_100m[0],2.0)+math.pow(sr_comps_100m[1],2.0))
    
	sr_comps_150m = srwind(uwind,vwind,hght,umotion,vmotion,150.0)
	wind_150m= np.sqrt(math.pow(sr_comps_150m[0]+umotion,2.0)+math.pow(sr_comps_150m[1]+vmotion,2.0))
	srwind_150m = np.sqrt(math.pow(sr_comps_150m[0],2.0)+math.pow(sr_comps_150m[1],2.0))
    
	sr_comps_6km = srwind(uwind,vwind,hght,umotion,vmotion,6000.0)
	srwind_6km = np.sqrt(math.pow(sr_comps_6km[0],2.0)+math.pow(sr_comps_6km[1],2.0))
    
	sr_comps_9km = srwind(uwind,vwind,hght,umotion,vmotion,9000.0)
	[srwind_9km,srdir_9km] = uvconvert(sr_comps_9km[0],sr_comps_9km[1])
	
	## 700-500 mb lapse rate calculation

	lapse = lapserate(temp,press,hght,700.,500.)
    
    
	DCAPE_200 = Capevals.dcape(hght,press,temp,mixr,RFD_level)
	DCAPE_400 = Capevals.dcape(hght, press, temp, mixr, -400)
	sbcaperes = Capevals.sbcape(hght,press,temp,mixr)

	mucaperes = Capevals.mucape(hght,press,temp,mixr)

    
	if folder == 'F3Plus':
		addname = 'S'
	elif folder == 'F1F2':
		addname = 'W'
	distcape = 0
	
	[eff_bot,eff_top,eff_bot_loc,eff_top_loc] = Capevals.efflayer(hght,press,temp,mixr,uwind,vwind)
	
	if eff_bot >= 0 :
		effbs_top = (mucaperes[4]-eff_bot)/2
		effbwd = bulkshear(uwind,vwind,hght,eff_bot,effbs_top)
		effbs = effbwd / (effbs_top-eff_bot)
	else:
		effbwd = -999
		effbs = -999
	
	print 'effbwd', effbwd, 'effbs_top', effbs_top
	if distcape==0:
      
		valuesname=directory + '/values_files/' + csvresults[1] 
		results=open(valuesname,'w')
		writer=csv.writer(results, delimiter=',')
		#writer.writerow(['Storm Speed (mph)', uvconvert(umotion,vmotion)[0]])
		#writer.writerow(['Storm motion direction (degrees)', uvconvert(umotion,vmotion)[1]])
		#writer.writerow(['Tornado Duration (mins)',dur])
		writer.writerow(['MLCAPE', mlcaperesults[0]])
		writer.writerow(['MLCIN',mlcaperesults[1]])
		writer.writerow(['MLLCL', mlcaperesults[2]])
		writer.writerow(['MLLFC', mlcaperesults[3]])
		writer.writerow(['MLLCLT', mlcaperesults[7]])
		writer.writerow(['MLLFCT', mlcaperesults[8]])
		writer.writerow(['DCAPE', DCAPE_400[0]])
		writer.writerow(['DCAPET', DCAPE_400[1]])
		writer.writerow(['DCAPE200', DCAPE_200[0]])
		writer.writerow(['DCAPE200T', DCAPE_200[1]])
		writer.writerow(['CAPE03',mlcaperesults[9]])
		writer.writerow(['SBCAPE', sbcaperes[0]])
		
		writer.writerow(['MUCAPE', mucaperes[0]])
		writer.writerow(['MUH', mucaperes[6]])
		writer.writerow(['MUCIN', mucaperes[1]])
        
		writer.writerow(['SRH1', helicity(uwind,vwind,hght,umotion,vmotion,0,1000)[0]])
		writer.writerow(['SRH3', helicity(uwind,vwind,hght,umotion,vmotion,0,3000)[0]])
		writer.writerow(['SRH6', helicity(uwind,vwind,hght,umotion,vmotion,0,6000)[0]])
        
		writer.writerow(['VS01', bulkshear(uwind,vwind,hght,0,1000)])
		writer.writerow(['VS03', bulkshear(uwind,vwind,hght,0,3000)])
		writer.writerow(['VS06', bulkshear(uwind,vwind,hght,0,6000)])
		writer.writerow(['VS09', bulkshear(uwind,vwind,hght,0,9000)])
		writer.writerow(['VS13', bulkshear(uwind,vwind,hght,1000,3000)])
		writer.writerow(['VS16', bulkshear(uwind,vwind,hght,1000,6000)])
		writer.writerow(['VS19', bulkshear(uwind,vwind,hght,1000,9000)])
		writer.writerow(['VS36', bulkshear(uwind,vwind,hght,3000,6000)])
		writer.writerow(['VS39', bulkshear(uwind,vwind,hght,3000,9000)])
		
		writer.writerow(['EFFSRH', helicity(uwind,vwind,hght,umotion,vmotion,eff_bot,eff_top)[1]])
		writer.writerow(['EFFBWD', effbwd])
		writer.writerow(['EFFBS', effbs])
		writer.writerow(['EFFBOTP', eff_bot])
		writer.writerow(['EFFTOPP', eff_top])
        
		writer.writerow(['EHI1', ehi(mlcaperesults[0],uwind,vwind,hght,umotion,vmotion,0,1000)[0]])
		writer.writerow(['EHI3', ehi(mlcaperesults[0],uwind,vwind,hght,umotion,vmotion,0,3000)[0]])
		writer.writerow(['SFCPR', sfcpress])
		
		writer.writerow(['SFCT', sfct])
        
        
        #if folder == 'Nighttime_Soundings' or folder == 'LTTRUC' or folder == 'Hurricane' or folder == 'New_Sound' :
			#writer.writerow(['TOD',tod])

		results.close()
        
	if distcape==1:
		capedist=open(makecsv.make(h, t, td, u, v, file, altitude, resolution, folder)[2],'w')
		writer=csv.writer(capedist,delimiter=',')
		writer.writerow(caperesults[0])
		writer.writerow(caperesults[1])
		writer.writerow([caperesults[2]])
        
		capedist.close()
        
	interpsounding.close()
