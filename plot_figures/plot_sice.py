# script to process and plot sea ice data

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import cmocean

### --- read netcdf files --- ###

# Sea ice thickness, GREEN
ds = xr.open_dataset('hol_green/seaIce/sit_green.nc')
lon       = ds.x.values
lat       = ds.y.values
sit_green = ds.seaIceDiag.values
xr.Dataset.close(ds)

# Sea ice concentration, GREEN
ds = xr.open_dataset('hol_green/seaIce/ice_green.nc')
sia_green = ds.seaIceDiag.values
xr.Dataset.close(ds)

# Sea ice thickness, BLUE
ds = xr.open_dataset('hol_blue/seaIce/sit_blue.nc')
sit_blue = ds.seaIceDiag.values
xr.Dataset.close(ds)

# Sea ice concentration, BLUE
ds = xr.open_dataset('hol_blue/seaIce/ice_blue.nc')
sia_blue = ds.seaIceDiag.values
xr.Dataset.close(ds)

# Ocean heat content (only for mask)
ds = xr.open_dataset('hol_green/phys/green_ohc.nc')
ohc = ds.physDiag.values
ohc = np.squeeze(ohc[:,0,:,:])
xr.Dataset.close(ds)

### ------ ###

# Apply mask
sit_blue[ohc==0] = np.nan
sit_green[ohc==0] = np.nan
sia_blue[ohc==0] = np.nan
sia_green[ohc==0] = np.nan

#del ohc
import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

# Read in model grid dimensions
rac = mds.rdmds('data/RAC')
#rac=np.tile(rac,(7,1,1))
rac=np.tile(rac,(84,1,1))

#rac[ohc==0] = np.nan

#del rac

# Read in bathymetry
bathy=mds.rdmds('data/Depth');
off_shelf=np.copy(bathy)
off_shelf[220:,:]=np.nan
off_shelf[:50,:]=0
off_shelf[:95,100:]=0
off_shelf[:150,500:]=0
off_shelf[off_shelf>1000]=np.nan
bathy[np.isnan(off_shelf)]=np.nan

bathy = np.tile(bathy,(84,1,1))

del off_shelf

# Apply mask
open_green = np.zeros((84,384,600))
open_green[sia_green<0.1] = 1
open_blue  = np.zeros((84,384,600))
open_blue[sia_blue<0.10] = 1
open_green = np.multiply(open_green,rac[0:84,:,:])*1e-6
open_blue  = np.multiply(open_blue,rac[0:84,:,:])*1e-6
open_green[bathy==0]=np.nan
open_green[np.isnan(bathy)]=np.nan
open_blue[bathy==0]=np.nan
open_blue[np.isnan(bathy)]=np.nan
open_pct = np.divide(np.nansum(np.nansum(open_green-open_blue,2),1),np.nansum(np.nansum(open_blue,2),1))*100
open_pct[np.nansum(np.nansum(open_green,2),1) < 1000] = np.nan

### --- Plot figures --- ###

# Figure 4a
fig, ax=plt.subplots(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),0,490,color=(0.9,0.9,0.9))
plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(open_green,2),1)/1000,color=(51/255,160/255,44/255),linewidth=10)
plt.yticks(fontsize=60,color=(51/255,160/255,44/255))
plt.xlim(0,84)
plt.grid()
plt.ylim(0,490)
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Open water (10³ km²)',fontsize=60,color=(51/255,160/255,44/255))
plt.yticks(fontsize=60)
plt.savefig('open_water.png')

# Figure 4b
fig, ax=plt.subplots(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),-20,80,color=(0.9,0.9,0.9))
plt.plot(np.linspace(0,84,84),np.nansum(np.nansum(open_green-open_blue,2),1)/1000,color=(51/255,160/255,44/255),linewidth=10)
plt.bar(np.linspace(0,84,84),open_pct,width=1.05,color=(178/255,223/255,138/255))
plt.yticks(fontsize=60,color=(51/255,160/255,44/255))
plt.xlim(0,84)
plt.grid()
plt.ylim(-12,70)
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
secax = ax.secondary_yaxis('right')
secax.set_yticks([-40,-20,0,20,40,60],['-40','-20','0','20','40','60'],fontsize=60,color=(178/255,223/255,138/255))
secax.set_ylabel('    %',fontsize=70,rotation=0,color=(178/255,223/255,138/255))
plt.ylabel('Anomaly (10³ km²)',fontsize=60,color=(51/255,160/255,44/255))
plt.yticks([0,20,40,60],['0','20','40','60'],fontsize=60)
plt.savefig('open_diff.png')

sit_green[np.isnan(bathy)]=np.nan
sit_blue[np.isnan(bathy)]=np.nan
sia_green[np.isnan(bathy)] = np.nan

# select only summer months
new_sia_green = np.reshape(sia_green,(7,12,384,600))
new_sia_green = (new_sia_green[:,0,:,:]+new_sia_green[:,1,:,:]+new_sia_green[:,11,:,:])/3
new_sit_green = np.reshape(sit_green,(7,12,384,600))
new_sit_green = (new_sit_green[:,0,:,:]+new_sit_green[:,1,:,:]+new_sit_green[:,11,:,:])/3
new_sit_blue  = np.reshape(sit_blue,(7,12,384,600))
new_sit_blue = (new_sit_blue[:,0,:,:]+new_sit_blue[:,1,:,:]+new_sit_blue[:,11,:,:])/3
summer_sit_blue = np.copy(new_sit_blue)
summer_sit_green = np.copy(new_sit_green)

fig_mask = np.zeros((384,600))
fig_mask[np.squeeze(bathy[0,:,:])==0] = 1.5
fig_mask[np.isnan(np.squeeze(bathy[0,:,:]))] = 2.5
ice_cmap = ((118/255,42/255,131/255),(153/255,112/255,171/255),(194/255,165/255,207/255),(231/255,212/255,232/255),(217/255,240/255,211/255),(166/255,219/255,160/255),(90/255,174/255,97/255),(27/255,120/255,55/255))

mask_cmap = ((1,1,1),(0,0,0),(0.7,0.7,0.7))

# Figure 4c
plt.figure(figsize=(40,10)); 
pcol=plt.contourf(lon,lat,100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:],0),np.linspace(-8,8,9),colors=ice_cmap,extend='both');
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=20)
plt.ylim(-75.5,-69)
pcol=plt.contourf(lon,lat,fig_mask,[0,1,2,3],colors=mask_cmap,alpha=0.1)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
# 10% ice cover threshold
plt.contour(lon,lat,np.nanmean(new_sia_green[:5,:,:],0),[0.1],linewidths=8,linestyles='--',colors='black')
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=60)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=60)
#plt.ylabel('$\degree$S',fontsize=60)
#plt.xlabel('$\degree$W',fontsize=60)
plt.savefig('sit_summer.png')

rac_div = np.squeeze(rac[:7,:,:])
rac_div[np.isnan(new_sit_green)] = np.nan

# integrate over domain
new_sit_green = np.multiply(new_sit_green,rac_div)
new_sit_blue  = np.multiply(new_sit_blue,rac_div)

print('Avg change in summertime sit: {} cm'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(rac_div[:5,:,:])))
print('Avg change in summertime sit: {} %'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(new_sit_green[:5,:,:])))

# select only winter months
new_sia_green = np.reshape(sia_green,(7,12,384,600))
new_sia_green = (new_sia_green[:,5,:,:]+new_sia_green[:,6,:,:]+new_sia_green[:,7,:,:])/3
new_sit_green = np.reshape(sit_green,(7,12,384,600))
new_sit_green = (new_sit_green[:,5,:,:]+new_sit_green[:,6,:,:]+new_sit_green[:,7,:,:])/3
new_sit_blue  = np.reshape(sit_blue,(7,12,384,600))
new_sit_blue = (new_sit_blue[:,5,:,:]+new_sit_blue[:,6,:,:]+new_sit_blue[:,7,:,:])/3

# Figure 4d
plt.figure(figsize=(40,10));
pcol=plt.contourf(lon,lat,100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:],0),np.linspace(-8,8,9),colors=ice_cmap,extend='both');
cbar=plt.colorbar(extend='both')
cbar.ax.set_yticklabels(['-8','','-4','','0','','4','','8'])
cbar.ax.tick_params(labelsize=30)
plt.ylim(-75.5,-69)
pcol=plt.contourf(lon,lat,fig_mask,[0,1,2,3],colors=mask_cmap,alpha=0.1)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
# 10% ice cover threshold
#plt.contour(lon,lat,np.nanmean(new_sia_green[:5,:,:],0),[0.1],linewidths=5,linestyles='--',colors='black')
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=60)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=60)
#plt.ylabel('$\degree$S',fontsize=60)
#plt.xlabel('$\degree$W',fontsize=60)
plt.savefig('sit_winter.png')

# integrate over domain
new_sit_green = np.multiply(new_sit_green,rac_div)
new_sit_blue  = np.multiply(new_sit_blue,rac_div)

print('Avg change in wintertime sit: {} cm'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(rac_div[:5,:,:])))
print('Avg change in wintertime sit: {} %'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(new_sit_green[:5,:,:])))
