# script to plot hovmoller diagrams of temperature, meltwater

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import cmocean

import MITgcmutils as mitgcm
import MITgcmutils.mds as mds


rac = mds.rdmds('data/RAC')
rac=np.tile(rac,(50,1,1))
#rac=np.tile(rac,(84,1,1,1))

RC = np.squeeze(mds.rdmds('data/RC'))
RC = -RC
depth = np.zeros((50))
depth[1:50] = np.diff(RC)
depth = np.tile(depth,(384,1))
depth = np.tile(depth,(600,1,1))
depth = np.transpose(depth)
#depth = np.tile(depth,(84,1,1,1))
vol   = np.multiply(depth,rac)

del depth
del RC
del rac

# plot volume of water
ds = xr.open_dataset('/media/twelves/My Passport/edinburgh_work/blue_temp.nc')
blue_temp = ds.physDiag.values
xr.Dataset.close(ds)

ds = xr.open_dataset('/media/twelves/My Passport/edinburgh_work/blue_salt.nc')
blue_salt = ds.physDiag.values
xr.Dataset.close(ds)

#blue_temp = blue_temp*vol
blue_temp[blue_salt>34] = 0
#blue_temp[blue_salt<34.5] = 0

ds = xr.open_dataset('/media/twelves/My Passport/edinburgh_work/green_temp.nc')
green_temp = ds.physDiag.values
xr.Dataset.close(ds)

ds = xr.open_dataset('/media/twelves/My Passport/edinburgh_work/green_salt.nc')
green_salt = ds.physDiag.values
xr.Dataset.close(ds)

#green_temp = green_temp*vol
green_temp[green_salt>34] = 0
#green_temp[green_salt<34.5] = 0

# continental shelf only
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
bathy = np.tile(bathy,(50,1,1))
bathy = np.tile(bathy,(84,1,1,1))

#mw[np.isnan(bathy)]=0
#temp[np.isnan(bathy)]=np.nan
blue_temp[np.isnan(bathy)] = np.nan
green_temp[np.isnan(bathy)] = np.nan

#rac = mds.rdmds('data/RAC')
#rac = np.tile(rac,(50,1,1))
#rac = np.tile(rac,(84,50,1,1))

Q = np.fromfile('data/shelfice_bedmach.bin',dtype='float64').byteswap().reshape((384,600))
ice_topo=-Q
ice_topo=np.tile(ice_topo,(50,1,1))
ice_topo=np.tile(ice_topo,(84,1,1,1))

blue_temp[ice_topo>0]=np.nan
green_temp[ice_topo>0]=np.nan

del ice_topo

temp=green_temp-blue_temp
temp[np.isnan(bathy)]=np.nan
del bathy
pos_sst=np.nanmean(np.nanmean(np.squeeze(temp[:,0,:,:]),2),1)
pos_sst[pos_sst<0]=0

neg_sst=np.nanmean(np.nanmean(np.squeeze(temp[:,0,:,:]),2),1)
neg_sst[neg_sst>0]=0

pos_sst_pib=np.nanmean(np.nanmean(np.squeeze(temp[:,0,0:100,300:400]),2),1)
pos_sst_pib[pos_sst_pib<0]=0

neg_sst_pib=np.nanmean(np.nanmean(np.squeeze(temp[:,0,0:100,300:400]),2),1)
neg_sst_pib[neg_sst_pib>0]=0

fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),-0.03,0.42,color=(0.9,0.9,0.9))
#fig=plt.fill_between(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(ice_blue[:,0:100,300:400],2),1),60000)
#plin=plt.plot(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(ice_green[:,0:100,300:400],2),1),'g',linewidth=5)
#pcol=plt.pcolormesh(month,depth,100*np.nanmean(np.nanmean(temp[:,:,0:100,300:400],3),2),cmap='cmo.delta')
pbar=plt.bar(np.linspace(0.5,83.5,84),pos_sst,width=1.05,color=(178/255,223/255,138/255))
pbar=plt.bar(np.linspace(0.5,83.5,84),neg_sst,width=1.05,color=(166/255,206/255,227/255))
plt.grid(axis='x')
plt.ylim(-0.03,0.42)
#pbar=plt.bar(np.linspace(0.5,83.5,84),pos_sst,width=1.05,color=(178/255,223/255,138/255))
#pbar=plt.bar(np.linspace(0.5,83.5,84),neg_sst,width=1.05,color=(166/255,206/255,227/255))
#plt.xlim(-0.5,84.5)
plt.margins(x=0)
plt.yticks([0,0.1,0.2,0.3,0.4])
#plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=40)
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)


plt.ylabel('$\Delta$ SST (deg. C)',fontsize=60)
plt.yticks(fontsize=40)
plt.savefig('sst_anom.png')

blue_above = np.zeros((84))
green_above = np.zeros((84))

for t in range(0,84):
    blue_month = np.squeeze(blue_temp[t,:,:,:])
    blue_vol = np.copy(vol)
    blue_vol[blue_month<0] = 0
    blue_vol[blue_month==0] = 0
    blue_vol[np.isnan(blue_month)] = 0
    blue_above[t] = np.nansum(blue_vol[:,:,:])
    green_month = np.squeeze(green_temp[t,:,:,:])
    green_vol = np.copy(vol)
    green_vol[green_month<0] = 0
    green_vol[green_month==0] = 0
    green_vol[np.isnan(green_month)] = 0
    green_above[t] = np.nansum(green_vol[:,:,:])

plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),0,8900,color=(0.9,0.9,0.9))
plin = plt.plot(np.linspace(0,84,84),green_above*1e-9,color=(51/255,160/255,44/255),linewidth=10)
plin = plt.plot(np.linspace(0,84,84),blue_above*1e-9,color=(31/255,120/255,180/255),linewidth=10,linestyle='--')
#plt.ylim(2500,17500)
#plt.ylim(99000,137000)
plt.ylim(0,8900)
#fig=plt.fill_between(np.linspace(0,84,84),blue_above,color='blue',alpha=0.2)
#plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=60)
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.xlim(0,84)
plt.grid()
plt.ylabel('Volume (10³ km³)',fontsize=60)
#plt.yticks([5000,10000,15000],['5','10','15'],fontsize=40)
#plt.yticks([100000,110000,120000,130000],['100','110','120','130'],fontsize=40)
plt.yticks([0,2000,4000,6000,8000],['0','2','4','6','8'],fontsize=40)
plt.savefig('aasw.png')

#ice_blue = np.multiply(ice_blue,rac)

depth = mds.rdmds('data/RC')
print(np.shape(depth))
depth = np.squeeze(depth[:,0,0])
depth = np.tile(depth,(84,1))

month = np.tile(np.linspace(0,84,84),(50,1))
month = np.transpose(month)

### --- plot figures --- ###

temp = green_temp - blue_temp

pos_sst=np.nanmean(np.nanmean(np.squeeze(temp[:,0,:,:]),2),1)
print(pos_sst)
pos_sst[pos_sst<0]=0
#print(np.nanmean(np.nanmean(pos_sst,2),1))

neg_sst=np.nanmean(np.nanmean(np.squeeze(temp[:,0,:,:]),2),1)
neg_sst[neg_sst>0]=0
#print(np.nanmean(np.nanmean(neg_sst,2),1))

pos_sst_pib=np.nanmean(np.nanmean(np.squeeze(temp[:,0,0:100,300:400]),2),1)
pos_sst_pib[pos_sst_pib<0]=0
#print(np.nanmean(np.nanmean(pos_sst,2),1))

neg_sst_pib=np.nanmean(np.nanmean(np.squeeze(temp[:,0,0:100,300:400]),2),1)
neg_sst_pib[neg_sst_pib>0]=0
#print(np.nanmean(np.nanmean(neg_sst,2),1))

fig=plt.figure(figsize=(40,10))
#fig=plt.fill_between(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(ice_blue[:,0:100,300:400],2),1),60000)
#plin=plt.plot(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(ice_green[:,0:100,300:400],2),1),'g',linewidth=5)
#pcol=plt.pcolormesh(month,depth,100*np.nanmean(np.nanmean(temp[:,:,0:100,300:400],3),2),cmap='cmo.delta')
#pbar=plt.bar(np.linspace(0,84,84),pos_sst_pib,width=1.05,color=(178/255,223/255,138/255))
#pbar=plt.bar(np.linspace(0,84,84),neg_sst_pib,width=1.05,color=(166/255,206/255,227/255))
pbar=plt.bar(np.linspace(0,84,84),pos_sst,width=1.05,color=(51/255,160/255,44/255))
pbar=plt.bar(np.linspace(0,84,84),neg_sst,width=1.05,color=(31/255,120/255,180/255))
plt.xlim(-0.5,84.5)
plt.margins(x=0)
#plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=40)
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)


plt.ylabel('$\Delta$ SST (deg. C)',fontsize=60)
plt.yticks(fontsize=40)
#plt.savefig('sst_anom.png')

# sea ice coverage

fig=plt.figure(figsize=(40,15))
#fig=plt.fill_between(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(ice_blue[:,0:100,300:400],2),1),60000)
#plin=plt.plot(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(ice_green[:,0:100,300:400],2),1),'g',linewidth=5)
#pcol=plt.pcolormesh(month,depth,100*np.nanmean(np.nanmean(temp[:,:,0:100,300:400],3),2),cmap='cmo.delta')
#pcol=plt.contourf(month,-depth,np.nanmean(np.nanmean(temp[:,:,:100,250:300],3),2),np.array([-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25]),cmap='cmo.delta',extend='both')
pcol=plt.contourf(month,-depth,np.nanmean(np.nanmean(green_temp[:,:,:100,250:300],3),2),np.array([-1.5,-1,-0.5,0,0.5,1,1.5]),cmap='cmo.balance',extend='both')
#pcol=plt.contourf(month,depth,np.nanmean(np.nanmean(temp[:,:,0:100,300:400],3),2),np.linspace(-0.8,0.8,11),cmap='cmo.delta')
plt.ylim(0,1500)
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.yticks(fontsize=40)
plt.xlim(0,84)
plt.ylim(0,1000)
plt.margins(x=0)
plt.gca().invert_yaxis()
plt.ylabel('Depth (m)',fontsize=60)
cbar=plt.colorbar(pcol,orientation='horizontal',extend='both')
cbar.ax.tick_params(labelsize=40)
plt.grid(axis='y')
plt.savefig('temp_dot_hov.png')
