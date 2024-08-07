# script to plot heat flux climatology

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

# convert temperature to energy density

spec = 4.186 * 1.02*1e6

# seconds in day

day = 86400

# read netcdf file

ds = xr.open_dataset('hol_green/surf/tflux_diff.nc')
flx = ds.surfDiag.values

ds = xr.open_dataset('hol_green/phys/diff_ohc.nc')
ohc = ds.physDiag.values

ds = xr.open_dataset('diff_iceheat.nc')
shf = ds.surfHeatDiag.values

flx[np.squeeze(ohc[:,0,:,:])==0] = np.nan
ohc[ohc==0] = np.nan
ice_mask = np.tile(flx[0,:,:],(50,1,1))
ice_mask = np.tile(ice_mask,(84,1,1,1))
ohc[np.isnan(ice_mask)]=np.nan

import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

rac = mds.rdmds('data/RAC')
rac=np.tile(rac,(50,1,1))
rac=np.tile(rac,(84,1,1,1))
flx = np.multiply(flx,np.squeeze(rac[:,0,:,:]))
shf = np.multiply(shf,np.squeeze(rac[:,0,:,:]))

RC = np.squeeze(mds.rdmds('data/RC'))
RC = -RC
depth = np.zeros((50))
depth[1:50] = np.diff(RC)
depth = np.tile(depth,(384,1))
depth = np.tile(depth,(600,1,1))
depth = np.transpose(depth)
depth = np.tile(depth,(84,1,1,1))
vol   = np.multiply(depth,rac)
ohc   = ohc*vol*spec/day

bathy=mds.rdmds('data/Depth');
off_shelf=np.copy(bathy)
off_shelf[220:,:]=np.nan
off_shelf[:50,:]=0
off_shelf[:95,100:]=0
off_shelf[:150,500:]=0
off_shelf[off_shelf>1000]=np.nan
bathy[np.isnan(off_shelf)]=np.nan
#plt.pcolormesh(bathy); plt.show()
shelf_mask=np.ones((384,600))
shelf_mask[np.isnan(bathy)]=0
shelf_mask = np.tile(shelf_mask,(50,1,1))
shelf_mask = np.tile(shelf_mask,(84,1,1,1))

flx = np.multiply(flx,shelf_mask[:,0,:,:])
ohc = np.multiply(ohc,shelf_mask)

flx = flx*3.154e-11
ohc = ohc*3.154e-11

flx_shf = np.nansum(np.nansum(flx,2),1)
ohc_shf = np.nansum(np.nansum(np.nansum(ohc,3),2),1)
shf_shf = np.nansum(np.nansum(shf,2),1)

adv_shf = ohc_shf-flx_shf

fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(0,84,84),ohc_shf,color=(255/255,237/255,111/255),alpha=0.7)
plin=plt.plot(np.linspace(0,84,84),flx_shf,linewidth=10,color=(253/255,180/255,98/255))
plin=plt.plot(np.linspace(0,84,84),adv_shf,linestyle=':',linewidth=10,color=(253/255,180/255,98/255))
plt.yticks(fontsize=40)
plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=40)
plt.ylabel('Heat anomaly (EJ/yr)',fontsize=60)
plt.savefig('shelf_heat_budget.png')

fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(0,84,84),np.cumsum(ohc_shf)/12,color=(255/255,237/255,111/255),alpha=0.7)
plin=plt.plot(np.linspace(0,84,84),np.cumsum(flx_shf)/12,linewidth=10,color=(253/255,180/255,98/255))
plin=plt.plot(np.linspace(0,84,84),np.cumsum(adv_shf)/12,linestyle=':',linewidth=10,color=(253/255,180/255,98/255))
#plin=plt.plot(np.linspace(0,84,84),np.cumsum(shf_shf),'r',linewidth=5)
plt.yticks(fontsize=40)
plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=40)
plt.ylabel('Cumalitive anomaly (EJ)',fontsize=60)
plt.savefig('shelf_heat_budget_cumalitive.png')

# restrict to pine island bay area

print(np.nanmean(flx_shf[:]))
print(np.nanmean(ohc_shf[:]))
print(np.nanmean(adv_shf[:]))

flx_shf = np.reshape(flx_shf,(7,12))
ohc_shf = np.reshape(ohc_shf,(7,12))
adv_shf = np.reshape(adv_shf,(7,12))

flx_shf_min = np.nanmin(flx_shf[0:5,:],0)
flx_shf_max = np.nanmax(flx_shf[0:5,:],0)
ohc_shf_min = np.nanmin(ohc_shf[0:5,:],0)
ohc_shf_max = np.nanmax(ohc_shf[0:5,:],0)
adv_shf_min = np.nanmin(adv_shf[0:5,:],0)
adv_shf_max = np.nanmax(adv_shf[0:5,:],0)

flx_shf = np.nanmean(flx_shf[0:5,:],0)
ohc_shf = np.nanmean(ohc_shf[0:5,:],0)
adv_shf = np.nanmean(adv_shf[0:5,:],0)

plt.figure(figsize=(40,10))
plt.plot(np.linspace(0,12,13),np.linspace(0,0,13),linewidth=5,color='k')
plt.plot(np.linspace(0.5,11.5,12),flx_shf,linewidth=10,color=(27/255,158/255,119/255))
plt.fill_between(np.linspace(0.5,11.5,12),flx_shf_min,flx_shf_max,color=(27/255,158/255,119/255),alpha=0.2)
plt.plot(np.linspace(0.5,11.5,12),ohc_shf,linewidth=10,color=(117/255,112/255,179/255))
plt.fill_between(np.linspace(0.5,11.5,12),ohc_shf_min,ohc_shf_max,color=(177/255,112/255,179/255),alpha=0.2)
plt.plot(np.linspace(0.5,11.5,12),adv_shf,linewidth=10,color=(217/255,95/255,2/255),linestyle='--')
plt.fill_between(np.linspace(0.5,11.5,12),adv_shf_min,adv_shf_max,color=(217/255,95/255,2/255),alpha=0.2)
#plt.plot(np.repeat(np.nanmean(flx_shf),12),'g',linestyle=':')
#plt.plot(np.repeat(np.nanmean(ohc_shf),12),'k',linestyle=':')
#plt.plot(np.repeat(np.nanmean(adv_shf),12),'b',linestyle=':')
plt.xlim(0,12)
plt.ylim(-140,140)
plt.grid()
plt.xticks(np.linspace(0,12,13),['            Jan  ',
                                 '            Feb  ',
                                 '            Mar  ',
                                 '            Apr  ',
                                 '            May  ',
                                 '            Jun  ',
                                 '            Jul  ',
                                 '            Aug  ',
                                 '            Sep  ',
                                 '            Oct  ',
                                 '            Nov  ',
                                 '            Dec  ',''],fontsize=60)

plt.legend()
plt.ylabel('Heat anom. (EJ yr⁻¹)',fontsize=70)
plt.yticks(np.linspace(-100,100,5),[' -100','  -50','   0','  50','   100'],fontsize=60)
plt.savefig('seasonal_heat.png')

print('surf={}'.format(np.nanmean(flx_shf)))
print('ohc={}'.format(np.nanmean(ohc_shf)))
print('adv={}'.format(np.nanmean(adv_shf)))

flx = np.nansum(np.nansum(flx,2),1)
shf = np.nansum(np.nansum(shf,2),1)
ohc = np.nansum(np.nansum(np.nansum(ohc,3),2),1)
#flx = np.nansum(np.nansum(flx,2),1)
#ohc = np.nansum(np.nansum(np.nansum(ohc,3),2),1)

#del adv
adv = ohc - flx

### --- plot figures --- ###

# sea ice coverage

fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(0,84,84),ohc,color=(255/255,237/255,111/255),alpha=0.7)
plin=plt.plot(np.linspace(0,84,84),flx,linewidth=10,color=(253/255,180/255,98/255))
plin=plt.plot(np.linspace(0,84,84),adv,linestyle=':',linewidth=10,color=(253/255,180/255,98/255))
plt.yticks(fontsize=40)
plt.savefig('pib_heat_budget.png')

fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(0,84,84),np.cumsum(ohc),color='y',alpha=0.3)
plin=plt.plot(np.linspace(0,84,84),np.cumsum(flx),'k',linewidth=5)
plin=plt.plot(np.linspace(0,84,84),np.cumsum(adv),'k:',linewidth=5)
#plin=plt.plot(np.linspace(0,84,84),np.cumsum(shf),'k',linewidth=5)
plt.yticks(fontsize=40)
plt.savefig('pib_heat_budget_cumalitive.png')
