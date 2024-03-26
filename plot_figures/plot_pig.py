# script to process and plot sea ice data

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

# read netcdf file

ds = xr.open_dataset('hol_green/surf/green_meltrate.nc')
mlt_green = ds.surfDiag.values

ds = xr.open_dataset('hol_blue/surf/blue_meltrate.nc')
mlt_blue = ds.surfDiag.values

Q = np.fromfile('data/shelfice_bedmach.bin',dtype='float64').byteswap().reshape((384,600))
ice_topo=-Q
ice_topo[ice_topo==0]=np.nan
print(np.nanmin(Q))
shal_shelf = np.zeros((384,600))
shal_shelf[ice_topo<50] = 1
shal_shelf = np.tile(shal_shelf,(84,1,1))
deep_shelf = np.zeros((384,600))
deep_shelf[ice_topo>50] = 1
deep_shelf = np.tile(deep_shelf,(84,1,1))

#ds = xr.open_dataset('hol_green/phys/green_ohc.nc')
#ohc = ds.physDiag.values
#ohc = np.squeeze(ohc[:,0,:,:])

#ice_blue[ohc==0] = np.nan
#ice_green[ohc==0] = np.nan
#opn_blue[ohc==0] = np.nan
#opn_green[ohc==0] = np.nan

import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

rac = mds.rdmds('data/RAC')
rac=np.tile(rac,(84,1,1))
#print(np.nanmax(rac))
mlt_blue = -np.multiply(mlt_blue,rac)*3.154e-5
mlt_green = -np.multiply(mlt_green,rac)*3.154e-5
#plt.figure()
#pcol=plt.pcolormesh(ice_topo)
#plt.show()
#ice_topo = np.tile(ice_topo,(84,1,1))
#mlt_blue[np.isnan(ice_topo)] = np.nan

#shal

#mlt_blue = np.multiply(mlt_blue,shal_shelf)
#mlt_green = np.multiply(mlt_green,shal_shelf) 

#sdeep

#mlt_blue = np.multiply(mlt_blue,deep_shelf)
#mlt_green = np.multiply(mlt_green,deep_shelf) 

print(np.nansum(mlt_green[0,:,:]))
print(np.nansum(mlt_green[14,:,:]))
print('surface melt anomaly: {}'.format(100*(np.nansum(mlt_green)-np.nansum(mlt_blue))/np.nansum(mlt_blue)))
#mlt_green = np.reshape(mlt_green, (7,12,384,600))
#mlt_blue = np.reshape(mlt_blue, (7,12,384,600))

#print(np.nansum(mlt_green[0,0,:,:]))
#print(np.nansum(mlt_green[1,2,:,:]))

### --- plot figures --- ###

# All ice shelves

fig=plt.figure(figsize=(40,10))
#pbar=plt.bar(np.linspace(0,7,7),np.nansum(1e-6*np.nansum(mlt_blue[:,0:100,300:400],2),1))
#pbar=plt.bar(np.linspace(2007.8,2013.8,7),np.nansum(np.nansum(np.nansum(mlt_blue[:,:,0:60,370:420],3),2),1),color='g',width=0.4)
#pbar=plt.bar(np.linspace(2008.2,2014.2,7),np.nansum(np.nansum(np.nansum(mlt_green[:,:,0:60,370:420],3),2),1),color='b',width=0.4)
fig=plt.fill_between(np.linspace(60,84,37),0,620,color=(0.9,0.9,0.9))
plt.ylim(210,620)
plt.xlim(0,84)
plin1,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,:,:],2),1),color=(51/255,160/255,44/255),linewidth=10,label='with chlorophyll')
plin2,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,:,:],2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle='--',label='without chlorophyll')
plt.yticks(fontsize=60)
#plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=60)
#plt.xticks
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.grid()

plt.ylabel('Melt (Gt/yr)',fontsize=60)
#plin = plt.plot(np.nansum(np.nansum(mlt_blue[:,60:93,370:420],2),1),'b:',linewidth=5)
#plin = plt.plot(np.nansum(np.nansum(mlt_green[:,60:93,370:420],2),1),'g',linewidth=5)
#scat = plt.scatter([12,24,48],[
plt.yticks(fontsize=40)
#plt.xticks(fontsize=40)
plt.legend(handles=[plin1,plin2],fontsize=40)
plt.savefig('cdw_melt.png')

# with depth...
print(np.nanmax(ice_topo))
green_flat = np.ravel(np.nanmean(mlt_green[:60,:,:],0))
blue_flat  = np.ravel(np.nanmean(mlt_blue[:60,:,:],0))
#ice_topo   = np.tile(ice_topo,(84,1,1))
topo_flat  = np.ravel(ice_topo)
print(np.nanmax(topo_flat))

fig = plt.figure(figsize=(40,20))
hst_blue=np.histogram(topo_flat,range=(0,1600),bins=32,weights=blue_flat)
hst_green=np.histogram(topo_flat,range=(0,1600),bins=32,weights=green_flat)
print(topo_flat)
print(blue_flat)
plt.grid()
plt.xlim(0,1600)
plt.ylim(0,72)
plin1, = plt.plot(np.linspace(25,1575,32),hst_green[0],color=(51/255,160/255,44/255),linewidth=10,label='with chlorophyll')
plin2, = plt.plot(np.linspace(25,1575,32),hst_blue[0],color=(31/255,120/255,180/255),linewidth=10,linestyle='--',label='without chlorophyll')
plt.yticks([0,10,20,30,40,50,60,70],['','10','20','30','40','50','60','70',],fontsize=60)
plt.xticks(fontsize=60)
plt.legend(handles=[plin1,plin2],fontsize=40)
plt.ylabel('Melt (Gt/yr)',fontsize=60)
plt.xlabel('Depth (m)', fontsize=60)
plt.savefig('melt_wrt_depth.png')

print('total melt from green:{}'.format(np.nansum(green_flat)))
print('total melt from blue:{}'.format(np.nansum(blue_flat)))
print('total melt anomaly:{}'.format(np.nansum(blue_flat-green_flat)))
print('total melt anomaly as %:{}'.format(100*(np.nansum(blue_flat)-np.nansum(green_flat))/np.nansum(blue_flat)))

mlt_diff = np.ravel(np.nanmean(mlt_green-mlt_blue,0))
mlt_neg  = np.copy(mlt_diff)
mlt_pos  = np.copy(mlt_diff)
ice_drft = np.ravel(Q)
mlt_neg[ice_drft>-50]=0
mlt_pos[ice_drft<-50]=0

fig = plt.figure(figsize=(40,25))
plt.hist(-ice_drft,weights=mlt_neg,range=(0,1600),color=(166/255,206/255,227/255),bins=32)
plt.hist(-ice_drft,weights=mlt_pos,range=(0,1600),color=(178/255,223/255,138/255),bins=32)
#plt.gca().invert_xaxis()
plt.yticks(fontsize=60)
plt.xticks(fontsize=60)
plt.ylabel('Melt rate anomaly (Gt/yr)',fontsize=80)
#plin = plt.plot(np.nansum(np.nansum(mlt_blue[:,60:93,370:420],2),1),'b:',linewidth=5)
#plin = plt.plot(np.nansum(np.nansum(mlt_green[:,60:93,370:420],2),1),'g',linewidth=5)
#scat = plt.scatter([12,24,48],[
plt.xlabel('Depth (m)',fontsize=80)
plt.savefig('melt_prof.png')
print(np.nansum(mlt_diff))

#fig=plt.figure()
#plt.scatter(Q,np.nansum(mlt_green-mlt_blue,0))
#plt.show()

# Individually

fig=plt.figure(figsize=(40,20))
#pbar=plt.bar(np.linspace(0,7,7),np.nansum(1e-6*np.nansum(mlt_blue[:,0:100,300:400],2),1))
#pbar=plt.bar(np.linspace(2007.8,2013.8,7),np.nansum(np.nansum(np.nansum(mlt_blue[:,:,0:60,370:420],3),2),1),color='g',width=0.4)
#pbar=plt.bar(np.linspace(2008.2,2014.2,7),np.nansum(np.nansum(np.nansum(mlt_green[:,:,0:60,370:420],3),2),1),color='b',width=0.4)
plin=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,:,255:295],2),1),color=(178/255,223/255,138/255),linewidth=20)
plin=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,:,255:295],2),1),color=(31/255,120/255,180/255),linewidth=20,linestyle=':')
plt.yticks(fontsize=80)
plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=80)
plt.xticks
plt.ylabel('Melt (Gt/yr)',fontsize=100)
#plin = plt.plot(np.nansum(np.nansum(mlt_blue[:,60:93,370:420],2),1),'b:',linewidth=5)
#plin = plt.plot(np.nansum(np.nansum(mlt_green[:,60:93,370:420],2),1),'g',linewidth=5)
#scat = plt.scatter([12,24,48],[
plt.yticks((20,40,60),fontsize=80)
#plt.xticks(fontsize=40)
plt.savefig('melt_dot.png')

fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),10,75,color=(0.9,0.9,0.9))
#pbar=plt.bar(np.linspace(0,7,7),np.nansum(1e-6*np.nansum(mlt_blue[:,0:100,300:400],2),1))
#pbar=plt.bar(np.linspace(2007.8,2013.8,7),np.nansum(np.nansum(np.nansum(mlt_blue[:,:,0:60,370:420],3),2),1),color='g',width=0.4)
#pbar=plt.bar(np.linspace(2008.2,2014.2,7),np.nansum(np.nansum(np.nansum(mlt_green[:,:,0:60,370:420],3),2),1),color='b',width=0.4)
plin1,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,0:60,370:420],2),1),color=(51/255,160/255,44/255),linewidth=10,label='with chlorophyll')
plin2,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,0:60,370:420],2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle='--',label='without chlorophyll')
plt.yticks(fontsize=80)
#plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=80)
#plt.xticks
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.ylabel('Melt (Gt/yr)',fontsize=60)
#plin = plt.plot(np.nansum(np.nansum(mlt_blue[:,60:93,370:420],2),1),'b:',linewidth=5)
#plin = plt.plot(np.nansum(np.nansum(mlt_green[:,60:93,370:420],2),1),'g',linewidth=5)
#scat = plt.scatter([12,24,48],[
plt.ylim(10,75)
plt.xlim(0,84)
plt.grid()
plt.yticks(fontsize=40)
plt.legend(handles=[plin1,plin2],fontsize=40)
#plt.xticks(fontsize=40)
plt.savefig('melt_pig.png')

fig=plt.figure(figsize=(40,20))
#pbar=plt.bar(np.linspace(0,7,7),np.nansum(1e-6*np.nansum(mlt_blue[:,0:100,300:400],2),1))
#pbar=plt.bar(np.linspace(2007.8,2013.8,7),np.nansum(np.nansum(np.nansum(mlt_blue[:,:,0:60,370:420],3),2),1),color='g',width=0.4)
#pbar=plt.bar(np.linspace(2008.2,2014.2,7),np.nansum(np.nansum(np.nansum(mlt_green[:,:,0:60,370:420],3),2),1),color='b',width=0.4)
plin=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,65:95,370:420],2),1),color=(51/255,160/255,44/255),linewidth=20,label='green')
plin=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,65:95,370:420],2),1),color=(31/255,120/255,180/255),linewidth=20,linestyle=':',label='blue')
plt.yticks(fontsize=80)
plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=80)
plt.xticks
plt.ylabel('Melt rate (Gt/yr)',fontsize=100)
#plin = plt.plot(np.nansum(np.nansum(mlt_blue[:,60:93,370:420],2),1),'b:',linewidth=5)
#plin = plt.plot(np.nansum(np.nansum(mlt_green[:,60:93,370:420],2),1),'g',linewidth=5)
#scat = plt.scatter([12,24,48],[
plt.yticks((5,10,15,20),fontsize=80)
#plt.xticks(fontsize=40)
plt.savefig('melt_cos.png')


