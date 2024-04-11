# script to process and plot sea ice outputs

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import cmocean

### --- Extract model outputs --- ###

# specify file path here
ds = xr.open_dataset('hol_green/seaIce/sit_green.nc')
# read in longitudes
lon       = ds.x.values
# read in latitudes
lat       = ds.y.values
# read in effective thickness (sea ice volume/unit area) from GREEN experiment
sit_green = ds.seaIceDiag.values
xr.Dataset.close(ds)

# specify file path here
ds = xr.open_dataset('hol_blue/seaIce/sit_blue.nc')
# read in effective thickness (sea ice volume/unit area) from BLUE experiment
sit_blue = ds.seaIceDiag.values
xr.Dataset.close(ds)

# specify file path here
ds = xr.open_dataset('hol_green/phys/green_ohc.nc')
# read in ocean heat content to use as mask
ohc = ds.physDiag.values
ohc = np.squeeze(ohc[:,0,:,:])
xr.Dataset.close(ds)

# mask out land, ice shelves
sit_blue[ohc==0] = np.nan
sit_green[ohc==0] = np.nan

import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

# read in cell areas for multiplication
rac = mds.rdmds('data/RAC')
rac=np.tile(rac,(7,1,1))

# read in bathymetry file
bathy=mds.rdmds('data/Depth');
off_shelf=np.copy(bathy)
# set mask to include only continental shelf
off_shelf[220:,:]=np.nan
off_shelf[:50,:]=0
off_shelf[:95,100:]=0
off_shelf[:150,500:]=0
off_shelf[off_shelf>1000]=np.nan
bathy[np.isnan(off_shelf)]=np.nan

bathy = np.tile(bathy,(84,1,1))

del off_shelf

# apply continental shelf mask to effective thicknesses
sit_green[np.isnan(bathy)]=np.nan
sit_blue[np.isnan(bathy)]=np.nan

### ------- ###

### --- Plot figures --- ###

# convert from months to years, months
new_sit_green = np.reshape(sit_green,(7,12,384,600))
# extract only December, January and February
new_sit_green = (new_sit_green[:,0,:,:]+new_sit_green[:,1,:,:]+new_sit_green[:,11,:,:])/3
# convert from months to years, months
new_sit_blue  = np.reshape(sit_blue,(7,12,384,600))
# extract only December, January and February
new_sit_blue = (new_sit_blue[:,0,:,:]+new_sit_blue[:,1,:,:]+new_sit_blue[:,11,:,:])/3

# apply masks for land and ice shelves
fig_mask = np.zeros((384,600))
fig_mask[np.squeeze(bathy[0,:,:])==0] = 1.5
fig_mask[np.isnan(np.squeeze(bathy[0,:,:]))] = 2.5
ice_cmap = ((118/255,42/255,131/255),(153/255,112/255,171/255),(194/255,165/255,207/255),(231/255,212/255,232/255),(217/255,240/255,211/255),(166/255,219/255,160/255),(90/255,174/255,97/255),(27/255,120/255,55/255))

mask_cmap = ((1,1,1),(0,0,0),(0.7,0.7,0.7))

plt.figure(figsize=(40,10)); 
pcol=plt.contourf(lon,lat,100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:],0),np.linspace(-8,8,9),colors=ice_cmap,extend='both');
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=20)
plt.ylim(-75.5,-69)
pcol=plt.contourf(lon,lat,fig_mask,[0,1,2,3],colors=mask_cmap,alpha=0.1)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=40)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=40)
plt.ylabel('$\degree$S',fontsize=40)
plt.xlabel('$\degree$W',fontsize=40)
plt.savefig('sit_summer.png')

# apply mask to cell areas also
rac[np.isnan(new_sit_green)] = np.nan

# multiply effective thickness by cell area to get volume
new_sit_green = np.multiply(new_sit_green,rac)
new_sit_blue  = np.multiply(new_sit_blue,rac)

print('Avg change in summertime sit: {} cm'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(rac[:5,:,:])))
print('Avg change in summertime sit: {} %'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(new_sit_green[:5,:,:])))

# convert from months to years, months
new_sit_green = np.reshape(sit_green,(7,12,384,600))
# extract only June, July and August
new_sit_green = (new_sit_green[:,5,:,:]+new_sit_green[:,6,:,:]+new_sit_green[:,7,:,:])/3
# convert from months to years, months
new_sit_blue  = np.reshape(sit_blue,(7,12,384,600))
# extract only June, July and August
new_sit_blue = (new_sit_blue[:,5,:,:]+new_sit_blue[:,6,:,:]+new_sit_blue[:,7,:,:])/3

plt.figure(figsize=(40,10));
pcol=plt.contourf(lon,lat,100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:],0),np.linspace(-8,8,9),colors=ice_cmap,extend='both');
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=20)
plt.ylim(-75.5,-69)
pcol=plt.contourf(lon,lat,fig_mask,[0,1,2,3],colors=mask_cmap,alpha=0.1)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=40)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=40)
plt.ylabel('$\degree$S',fontsize=40)
plt.xlabel('$\degree$W',fontsize=40)
plt.savefig('sit_winter.png')

# multiply effective thickness by cell area to get volume
new_sit_green = np.multiply(new_sit_green,rac)
new_sit_blue  = np.multiply(new_sit_blue,rac)

print('Avg change in wintertime sit: {} cm'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(rac[:5,:,:])))
print('Avg change in wintertime sit: {} %'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(new_sit_green[:5,:,:])))
