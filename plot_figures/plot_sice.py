# script to process and plot sea ice data

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import cmocean

# read netcdf file

ds = xr.open_dataset('hol_green/seaIce/ice_green.nc')
ice_green = ds.seaIceDiag.values
opn_green = np.copy(ice_green)
opn_green[:,:,:]=1
opn_green[ice_green<0.15] = 0
xr.Dataset.close(ds)

ds = xr.open_dataset('hol_blue/seaIce/ice_blue.nc')
ice_blue = ds.seaIceDiag.values
opn_blue = np.copy(ice_blue)
opn_blue[:,:,:]=1
opn_blue[ice_blue<0.15] = 0
lon = ds.x.values
lat = ds.y.values
xr.Dataset.close(ds)

ds = xr.open_dataset('hol_green/seaIce/sit_green.nc')
sit_green = ds.seaIceDiag.values
#sit_green = np.multiply(sit_green,ice_green)
xr.Dataset.close(ds)

ds = xr.open_dataset('hol_blue/seaIce/sit_blue.nc')
sit_blue = ds.seaIceDiag.values
#sit_blue = np.multiply(sit_blue,ice_blue)
xr.Dataset.close(ds)

#ice_diff = ice_green-ice_blue

ds = xr.open_dataset('hol_green/phys/green_ohc.nc')
ohc = ds.physDiag.values
ohc = np.squeeze(ohc[:,0,:,:])
xr.Dataset.close(ds)

sit_blue[ohc==0] = np.nan
sit_green[ohc==0] = np.nan
opn_blue[ohc==0] = np.nan
opn_green[ohc==0] = np.nan
ice_blue[ohc==0] = np.nan
ice_green[ohc==0] = np.nan
#rac[ohc==0] = np.nan

#del ohc
print('got here')
import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

rac = mds.rdmds('data/RAC')
rac=np.tile(rac,(84,1,1))
print(np.nanmax(rac))
#sit_blue = np.multiply(sit_blue,rac)
#sit_green = np.multiply(sit_green,rac)
opn_blue = np.multiply(opn_blue,rac)
opn_green = np.multiply(opn_green,rac)
#ice_diff = np.multiply(ice_diff,rac)

rac[ohc==0] = np.nan

#del rac

print('got here')

bathy=mds.rdmds('data/Depth');
off_shelf=np.copy(bathy)
off_shelf[220:,:]=np.nan
off_shelf[:50,:]=0
off_shelf[:95,100:]=0
off_shelf[:150,500:]=0
off_shelf[off_shelf>1000]=np.nan
bathy[np.isnan(off_shelf)]=np.nan
#plt.pcolormesh(bathy); plt.show()
#shelf_mask=np.ones((384,600))
#shelf_mask[np.isnan(bathy)]=0
#shelf_mask = np.tile(shelf_mask,(84,1,1,1))

bathy = np.tile(bathy,(84,1,1))

del off_shelf

print('got here')

sit_green[np.isnan(bathy)]=np.nan
sit_blue[np.isnan(bathy)]=np.nan
opn_green[np.isnan(bathy)]=np.nan
opn_blue[np.isnan(bathy)]=np.nan
ice_green[np.isnan(bathy)]=np.nan
ice_blue[np.isnan(bathy)]=np.nan
rac[np.isnan(bathy)]=np.nan

#sit_blue = np.multiply(sit_blue,rac)/np.nansum(rac)
#sit_green = np.multiply(sit_green,rac)/np.nansum(rac)

rac_wt = rac/np.nanmean(rac)

vol_blue=np.multiply(ice_blue,sit_blue)
vol_green=np.multiply(ice_green,sit_green)

new_sit_green = np.reshape(sit_green,(7,12,384,600))
#print(np.shape(new_sit_green))
new_sit_green = (new_sit_green[:,0,:,:]+new_sit_green[:,1,:,:]+new_sit_green[:,11,:,:])/3
new_sit_blue  = np.reshape(sit_blue,(7,12,384,600))
new_sit_blue = (new_sit_blue[:,0,:,:]+new_sit_blue[:,1,:,:]+new_sit_blue[:,11,:,:])/3
summer_sit_blue = np.copy(new_sit_blue)
summer_sit_green = np.copy(new_sit_green)

fig_mask = np.zeros((384,600))
fig_mask[np.squeeze(bathy[0,:,:])==0] = 1.5
fig_mask[np.isnan(np.squeeze(bathy[0,:,:]))] = 2.5
ice_cmap = ((118/255,42/255,131/255),(153/255,112/255,171/255),(194/255,165/255,207/255),(231/255,212/255,232/255),(217/255,240/255,211/255),(166/255,219/255,160/255),(90/255,174/255,97/255),(27/255,120/255,55/255))

#mask_cmap = ((1,1,1),(0,0,0),(33/255,102/255,172/255))
mask_cmap = ((1,1,1),(0,0,0),(0.7,0.7,0.7))

plt.figure(figsize=(40,10)); 
#pcol=plt.contourf(100*np.nanmean(sit_green-sit_blue,0)/np.nanmean(sit_green,0),np.linspace(-10,10,9),cmap='cmo.delta'); 
pcol=plt.contourf(lon,lat,100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:],0),np.linspace(-8,8,9),colors=ice_cmap,extend='both');
#pcol=plt.contourf(100*np.nanmean(new_sit_green-new_sit_blue,0),cmap='cmo.delta');
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
print('Avg change in summertime sit: {}'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])))
print('Avg change in summertime sit: {}'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(new_sit_green[:5,:,:])))

new_sit_green = np.reshape(sit_green,(7,12,384,600))
#print(np.shape(new_sit_green))
new_sit_green = (new_sit_green[:,5,:,:]+new_sit_green[:,6,:,:]+new_sit_green[:,7,:,:])/3
new_sit_blue  = np.reshape(sit_blue,(7,12,384,600))
new_sit_blue = (new_sit_blue[:,5,:,:]+new_sit_blue[:,6,:,:]+new_sit_blue[:,7,:,:])/3

winter_sit_blue = np.copy(new_sit_blue)
winter_sit_green = np.copy(new_sit_green)

plt.figure(figsize=(40,10));
#pcol=plt.contourf(100*np.nanmean(sit_green-sit_blue,0)/np.nanmean(sit_green,0),np.linspace(-10,10,9),cmap='cmo.delta');
pcol=plt.contourf(lon,lat,100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:],0),np.linspace(-8,8,9),colors=ice_cmap,extend='both');
#pcol=plt.contourf(100*np.nanmean(new_sit_green-new_sit_blue,0),cmap='cmo.delta');
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
print('Avg change in wintertime sit: {}'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])))
print('Avg change in wintertime sit: {}'.format(100*np.nanmean(new_sit_green[:5,:,:]-new_sit_blue[:5,:,:])/np.nanmean(new_sit_green[:5,:,:])))

new_ice_green = np.reshape(ice_green,(7,12,384,600))
#print(np.shape(new_sit_green))
new_ice_green = (new_ice_green[:,0,:,:]+new_ice_green[:,1,:,:]+new_ice_green[:,11,:,:])/3
new_ice_blue  = np.reshape(ice_blue,(7,12,384,600))
new_ice_blue = (new_ice_blue[:,0,:,:]+new_ice_blue[:,1,:,:]+new_ice_blue[:,11,:,:])/3
summer_ice_blue = np.copy(new_ice_blue)
summer_ice_green = np.copy(new_ice_green)

plt.figure(figsize=(40,10));
#pcol=plt.contourf(100*np.nanmean(sit_green-sit_blue,0)/np.nanmean(sit_green,0),np.linspace(-10,10,9),cmap='cmo.delta'); 
pcol=plt.contourf(lon,lat,100*np.nanmean(new_ice_green-new_ice_blue,0),np.linspace(-4,4,9),colors=ice_cmap,extend='both');
#pcol=plt.contourf(100*np.nanmean(new_sit_green-new_sit_blue,0),cmap='cmo.delta');
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=20)
plt.ylim(-75.5,-69)
pcol=plt.contourf(lon,lat,fig_mask,[0,1,2,3],colors=mask_cmap,alpha=0.1)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=40)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=40)
plt.ylabel('$\degree$S',fontsize=40)
plt.xlabel('$\degree$W',fontsize=40)
plt.savefig('sic_summer.png')

new_ice_green = np.reshape(ice_green,(7,12,384,600))
#print(np.shape(new_sit_green))
new_ice_green = (new_ice_green[:,5,:,:]+new_ice_green[:,6,:,:]+new_ice_green[:,7,:,:])/3
new_ice_blue  = np.reshape(ice_blue,(7,12,384,600))
new_ice_blue = (new_ice_blue[:,5,:,:]+new_ice_blue[:,6,:,:]+new_ice_blue[:,7,:,:])/3
winter_ice_blue = np.copy(new_ice_blue)
winter_ice_green = np.copy(new_ice_green)


plt.figure(figsize=(40,10));
#pcol=plt.contourf(100*np.nanmean(sit_green-sit_blue,0)/np.nanmean(sit_green,0),np.linspace(-10,10,9),cmap='cmo.delta');
pcol=plt.contourf(lon,lat,100*np.nanmean(new_ice_green-new_ice_blue,0),np.linspace(-4,4,9),colors=ice_cmap,extend='both');
#pcol=plt.contourf(100*np.nanmean(new_sit_green-new_sit_blue,0),cmap='cmo.delta');
cbar=plt.colorbar(extend='both');
cbar.ax.tick_params(labelsize=20)
plt.ylim(-75.5,-69)
pcol=plt.contourf(lon,lat,fig_mask,[0,1,2,3],colors=mask_cmap,alpha=0.1)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=40)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=40)
plt.ylabel('$\degree$S',fontsize=40)
plt.xlabel('$\degree$W',fontsize=40)
plt.savefig('sic_winter.png')

#new_ice_green = np.reshape(ice_green,(7,12,384,600))
#print(np.shape(new_sit_green))
#new_ice_green = (new_ice_green[:,0,:,:]+new_ice_green[:,1,:,:]+new_ice_green[:,11,:,:])/3
#new_ice_blue  = np.reshape(ice_blue,(7,12,384,600))
#new_ice_blue = (new_ice_blue[:,0,:,:]+new_ice_blue[:,1,:,:]+new_ice_blue[:,11,:,:])/3

summer_vol_green = np.multiply(summer_ice_green,summer_sit_green)#,np.squeeze(rac[:7,:,:]))
summer_vol_blue  = np.multiply(summer_ice_blue,summer_sit_blue)#,np.squeeze(rac[:7,:,:]))
winter_vol_green = np.multiply(winter_ice_green,winter_sit_green)#,np.squeeze(rac[:7,:,:]))
winter_vol_blue  = np.multiply(winter_ice_blue,winter_sit_blue)#,np.squeeze(rac[:7,:,:]))

plt.figure()
pcol = plt.pcolormesh(np.nanmean(rac,0))
rac = np.nanmean(rac,0)

plt.figure(figsize=(40,10));
#pnncol=plt.contourf(100*np.nanmean(sit_green-sit_blue,0)/np.nanmean(sit_green,0),np.linspace(-10,10,9),cmap='cmo.delta'); 
#pcol=plt.contourf(lon,lat,100*np.nanmean(winter_vol_green-winter_vol_blue,0),np.linspace(-4,4,9),cmap='cmo.delta',extend='both');
pcol=plt.contourf(lon,lat,1e-6*np.multiply(rac,np.nanmean(summer_vol_green-summer_vol_blue,0)),np.linspace(-1,1,9),colors=ice_cmap,extend='both');
#pcol=plt.contourf(100*np.nanmean(new_sit_green-new_sit_blue,0),cmap='cmo.delta');
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=20)
plt.ylim(-75.5,-69)
pcol=plt.contourf(lon,lat,fig_mask,[0,1,2,3],colors=mask_cmap,alpha=0.1)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=40)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=40)
plt.ylabel('$\degree$S',fontsize=40)
plt.xlabel('$\degree$W',fontsize=40)
plt.savefig('siv_summer.png')

plt.figure(figsize=(40,10));
#pcol=plt.contourf(100*np.nanmean(sit_green-sit_blue,0)/np.nanmean(sit_green,0),np.linspace(-10,10,9),cmap='cmo.delta');
#pcol=plt.contourf(lon,lat,100*np.nanmean(winter_vol_green-winter_vol_blue,0),np.linspace(-4,4,9),cmap='cmo.delta',extend='both');
pcol=plt.contourf(lon,lat,1e-6*np.nanmean(winter_vol_green-winter_vol_blue,0),np.linspace(-1,1,9),colors=ice_cmap,extend='both');
#pcol=plt.contourf(100*np.nanmean(new_sit_green-new_sit_blue,0),cmap='cmo.delta');
cbar=plt.colorbar(extend='both');
cbar.ax.tick_params(labelsize=20)
plt.ylim(-75.5,-69)
pcol=plt.contourf(lon,lat,fig_mask,[0,1,2,3],colors=mask_cmap,alpha=0.1)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=40)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=40)
plt.ylabel('$\degree$S',fontsize=40)
plt.xlabel('$\degree$W',fontsize=40)
plt.savefig('siv_winter.png')

new_ice_green = np.reshape(ice_green,(7,12,384,600))
print(np.shape(new_ice_green))
new_ice_green = np.nanmean(new_ice_green[:,6:8,:,:],1)
new_ice_blue  = np.reshape(ice_blue,(7,12,384,600))
new_ice_blue  = np.nanmean(new_ice_blue[:,6:8,:,:],1)

plt.figure(figsize=(8,5)); 
#pcol=plt.contourf(100*np.nanmean(new_ice_green-new_ice_blue,0)/np.nanmean(new_ice_green,0),np.linspace(-10,10,9),cmap='cmo.delta'); 
pcol=plt.contourf(100*np.nanmean(new_ice_green-new_ice_blue,0),np.linspace(-10,10,9),cmap='cmo.delta');
cbar=plt.colorbar(pcol); 
#cbar.set_label(fontsize=80)
plt.yticks([])
plt.xticks([])
plt.savefig('sic.png')

plt.figure(figsize=(8,5)); 
#pcol=plt.contourf(100*np.nanmean(vol_green-vol_blue,0)/np.nanmean(vol_green,0),np.linspace(-10,10,9),cmap='cmo.delta'); 
pcol=plt.contourf(100*np.nanmean(vol_green-vol_blue,0),np.linspace(-10,10,9),cmap='cmo.delta');
cbar=plt.colorbar(pcol);
#cbar.set_label(fontsize=80)
plt.yticks([])
plt.xticks([])
plt.savefig('siv.png')

ice_green_bin = np.array([np.nanmean(ice_green[10:11,:,:],0),np.nanmean(ice_green[22:23,:,:],0),np.nanmean(ice_green[34:35,:,:],0),np.nanmean(ice_green[46:47,:,:],0),np.nanmean(ice_green[58:59,:,:],0),np.nanmean(ice_green[70:71,:,:],0),np.nanmean(ice_green[82:83,:,:],0)])
#ice_green_bin[ice_green>0.15] = 1
#ice_green_bin[ice_green<0.15] = 0

plt.figure(figsize=(8,5));
pcol=plt.contourf(np.nanmean(ice_green_bin,0),np.linspace(0,0.15,2),cmap='cmo.delta');
cbar=plt.colorbar(pcol);
#cbar.set_label(fontsize=80)
plt.yticks([])
plt.xticks([])
plt.savefig('clima.png')

#del shelf_mask

### --- plot figures --- ###

# sea ice coverage

fig=plt.figure(figsize=(40,10))
#fig=plt.fill_between(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(ice_blue[:,0:100,300:400],2),1),60000)
#plin=plt.plot(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(ice_green[:,0:100,300:400],2),1),'g',linewidth=5)
plin=plt.plot(np.linspace(0,84,84),np.nanmean(np.nanmean(np.multiply(sit_green,rac_wt),2),1),color=(178/255,223/255,138/255),linewidth=10)
plin=plt.plot(np.linspace(0,84,84),np.nanmean(np.nanmean(np.multiply(sit_blue,rac_wt),2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle=':')
plt.yticks(fontsize=40)
plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=40)
plt.xticks
plt.ylabel('SIT (m)',fontsize=60)
#plt.xticks(fontsize=0)
plt.savefig('ice_thickness.png')

# open water area
fig=plt.figure(figsize=(40,10))
#plin=plt.plot(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(opn_blue[:,0:100,300:400],2),1),'b',linewidth=5)
#fig=plt.fill_between(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(opn_blue[:,0:100,300:400],2),1))
#plin=plt.plot(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(opn_green[:,0:100,300:400],2),1),'g',linewidth=5)
plin=plt.plot(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(opn_green,2),1),color=(178/255,223/255,138/255),linewidth=10)
plin=plt.plot(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(opn_blue,2),1),color=(31/255,120/255,180/255),linewidth=10,linestyle=':')
plt.yticks(fontsize=40)
plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=40)
plt.xticks
plt.ylabel('SIE (km²)',fontsize=60)
plt.savefig('open_water.png')
