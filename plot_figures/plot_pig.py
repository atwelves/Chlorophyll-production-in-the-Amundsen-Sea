# script to process and plot sea ice data
# Depending on the lines commented out, this script can be used to plot total melt,
# melt from deeper than 50m or melt from shallower than 50m

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

### read netcdf files ---

# ice shelf melt rates from GREEN
ds = xr.open_dataset('hol_green/surf/green_meltrate.nc')
mlt_green = ds.surfDiag.values

# ice shelf melt rates from BLUE
ds = xr.open_dataset('hol_blue/surf/blue_meltrate.nc')
mlt_blue = ds.surfDiag.values

# read ice shelf topography
Q = np.fromfile('data/shelfice_bedmach.bin',dtype='float64').byteswap().reshape((384,600))
ice_topo=-Q
ice_topo[ice_topo==0]=np.nan
shal_shelf = np.zeros((384,600))
shal_shelf[ice_topo<50] = 1
shal_shelf = np.tile(shal_shelf,(84,1,1))
deep_shelf = np.zeros((384,600))
deep_shelf[ice_topo>50] = 1
deep_shelf = np.tile(deep_shelf,(84,1,1))

import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

# read model grid dimensions
rac = mds.rdmds('data/RAC')
rac=np.tile(rac,(84,1,1))
mlt_blue = -np.multiply(mlt_blue,rac)*3.154e-5
mlt_green = -np.multiply(mlt_green,rac)*3.154e-5

shal_rac = np.multiply(shal_shelf,rac)
print('area of shallow shelf is {}'.format(np.nansum(shal_rac[1,:,:])*1e-9))
deep_rac = np.multiply(deep_shelf,rac)
print('area of deep shelf is {}'.format(np.nansum(deep_rac[1,:,:])*1e-9))
shal_frac = 100*np.nansum(shal_rac)/np.nansum(shal_rac+deep_rac)
print('shallow fraction of total is {}'.format(shal_frac))

# !!!!!!!!!!!!!!!!
# comment this out to select shallow component of ice shelf melt

#mlt_blue = np.multiply(mlt_blue,shal_shelf)
#mlt_green = np.multiply(mlt_green,shal_shelf) 

# !!!!!!!!!!!!!
# comment this out to select deep  component of ice shelf melt

#mlt_blue = np.multiply(mlt_blue,deep_shelf)
#mlt_green = np.multiply(mlt_green,deep_shelf) 

print('surface melt anomaly: {}'.format(100*(np.nansum(mlt_green)-np.nansum(mlt_blue))/np.nansum(mlt_blue)))
#mlt_green = np.reshape(mlt_green, (7,12,384,600))
#mlt_blue = np.reshape(mlt_blue, (7,12,384,600))

### --- plot figures --- ###

# All ice shelves

# Figure 9b
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),0,620,color=(0.9,0.9,0.9))
plt.ylim(210,620)
plt.xlim(0,84)
plin1,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,:,:],2),1),color=(51/255,160/255,44/255),linewidth=10,label='with chlorophyll')
plin2,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,:,:],2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle='--',label='without chlorophyll')
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.grid()
plt.ylabel('Melt (Gt yr$⁻¹$)',fontsize=80)
plt.yticks(np.linspace(200,600,5),[' 200',' 300',' 400',' 500',' 600'],fontsize=60)
plt.legend(handles=[plin1,plin2],fontsize=60)
#plt.savefig('cdw_melt.png')

pct_anom=100*np.divide(np.nansum(np.nansum(mlt_green[:,:,:],2),1)-np.nansum(np.nansum(mlt_blue[:,:,:],2),1),np.nansum(np.nansum(mlt_blue[:,:,:],2),1))
print('min anom is {}'.format(np.nanmin(-pct_anom[:60])))
print('max anom is {}'.format(np.nanmax(-pct_anom[:60])))

# Figure 8b
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),0,29,color=(0.9,0.9,0.9))
plt.ylim(0,29)
plt.xlim(0,84)
plin1,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,:,:],2),1),color=(51/255,160/255,44/255),linewidth=10,label='with chlorophyll')
plin2,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,:,:],2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle='--',label='without chlorophyll')
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.grid()
plt.ylabel('Melt (Gt yr$⁻¹$)',fontsize=80)
plt.yticks(np.linspace(0,25,6),['  0','  5',' 10',' 15',' 20',' 25'],fontsize=60)
plt.legend(handles=[plin1,plin2],fontsize=60)
plt.savefig('aasw_melt.png')

# with depth...
green_flat = np.ravel(np.nanmean(mlt_green[:60,:,:],0))
blue_flat  = np.ravel(np.nanmean(mlt_blue[:60,:,:],0))
#ice_topo   = np.tile(ice_topo,(84,1,1))
topo_flat  = np.ravel(ice_topo)

# Figure 9c

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
plt.yticks([0,10,20,30,40,50,60,70],['',' 10',' 20',' 30',' 40',' 50',' 60',' 70'],fontsize=60)
plt.xticks(fontsize=60)
plt.legend(handles=[plin1,plin2],fontsize=60)
plt.ylabel('Melt (Gt yr$⁻¹$)',fontsize=80)
plt.xlabel('Depth (m)', fontsize=80)
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

# Individually

# Pine Island Ice Shelf
# Figure 10c
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),0,75,color=(0.9,0.9,0.9))
plin1,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,0:60,370:420],2),1),color=(51/255,160/255,44/255),linewidth=10,label='with chlorophyll')
plin2,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,0:60,370:420],2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle='--',label='without chlorophyll')
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Melt (Gt yr$⁻¹$)',fontsize=80)
plt.ylim(10,75)
plt.xlim(0,84)
plt.grid()
plt.yticks(np.linspace(0,60,4),['  0',' 20',' 40',' 60'],fontsize=60)
plt.legend(handles=[plin1,plin2],fontsize=60)
plt.savefig('melt_pig.png')

# Thwaites Glacier ice shelf
fig=plt.figure(figsize=(40,10))
plin1,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,:,325:370],2),1),color=(51/255,160/255,44/255),linewidth=10,label='with chlorophyll')
plin2,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,:,325:370],2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle='--',label='without chlorophyll')
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Melt (Gt yr$⁻¹$)',fontsize=80)
#plt.ylim(10,75)
plt.xlim(0,84)
plt.grid()
#plt.yticks(np.linspace(0,60,4),['  0',' 20',' 40',' 60'],fontsize=60)
plt.legend(handles=[plin1,plin2],fontsize=60)
plt.savefig('melt_thw.png')

# Getz ice shelf
fig=plt.figure(figsize=(40,10))
plin1,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,:,50:255],2),1),color=(51/255,160/255,44/255),linewidth=10,label='with chlorophyll')
plin2,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,:,50:255],2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle='--',label='without chlorophyll')
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Melt (Gt yr$⁻¹$)',fontsize=80)
#plt.ylim(10,75)
plt.xlim(0,84)
plt.grid()
#plt.yticks(np.linspace(0,60,4),['  0',' 20',' 40',' 60'],fontsize=60)
plt.legend(handles=[plin1,plin2],fontsize=60)
plt.savefig('melt_getz.png')

# Dotson Ice Shelf
fig=plt.figure(figsize=(40,10))
plin1,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_green[:,:,255:310],2),1),color=(51/255,160/255,44/255),linewidth=10,label='with chlorophyll')
plin2,=plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(mlt_blue[:,:,255:310],2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle='--',label='without chlorophyll')
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Melt (Gt yr$⁻¹$)',fontsize=80)
#plt.ylim(10,75)
plt.xlim(0,84)
plt.grid()
#plt.yticks(np.linspace(0,60,4),['  0',' 20',' 40',' 60'],fontsize=60)
plt.legend(handles=[plin1,plin2],fontsize=60)
plt.savefig('melt_dot.png')
