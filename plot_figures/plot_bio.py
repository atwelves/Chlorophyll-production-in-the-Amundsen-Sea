# script to process and plot bio data

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

# read netcdf file

ds = xr.open_dataset('hol_green/phyto/chlorophyll.nc')
chl = ds.phytoDiag.values

ds = xr.open_dataset('hol_green/lightSurf/euphotic_depth.nc')
print(ds)
lon = ds.x.values
lat = ds.y.values
zeu = ds.lightSurfDiag.values

print(np.shape(lat))

ds = xr.open_dataset('hol_green/seaIce/ice_green.nc')
ice = ds.seaIceDiag.values

zeu_cp = np.copy(zeu)

chl[zeu==0]=np.nan
zeu[zeu==0]=np.nan

import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

#rac = mds.rdmds('data/RAC')
#print(np.tile(rac,(84,1,1)))
#print(np.nanmax(rac))
#chl = np.multiply(chl,rac)
#zeu = np.multiply(zeu,rac)

### --- plot figures --- ###

#restrict to continental shelf

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

bathy = np.tile(bathy,(144,1,1))

fig_mask = np.zeros((384,600))
fig_mask[np.squeeze(bathy[0,:,:])==0] = 1.5
fig_mask[np.isnan(np.squeeze(bathy[0,:,:]))] = 2.5
ice_cmap = ((118/255,42/255,131/255),(153/255,112/255,171/255),(194/255,165/255,207/255),(231/255,212/255,232/255),(217/255,240/255,211/255),(166/255,219/255,160/255),(90/255,174/255,97/255),(27/255,120/255,55/255))

mask_cmap = ((1,1,1),(0,0,0),(0.7,0.7,0.7))

del off_shelf

print('got here')

zeu_cp = np.copy(zeu)
ice_cp = np.copy(ice)

chl[np.isnan(bathy)]=np.nan
zeu[np.isnan(bathy)]=np.nan
bathy = np.squeeze(bathy[60:,:,:])
ice[np.isnan(bathy)]=np.nan

### -- Read in GlobColour data for same region

import calendar
from calendar import monthrange
import os.path as pth

# nested loops

zeu_avg = np.zeros((7,12))

for yr in range(2008,2015):
    year = '{}'.format(yr)
#    print(year)
    for m in range(1,13):
        month = '{:02d}'.format(m)
#        print(month)
        ym = year + month
        end_day=monthrange(yr,m)[1]
        fname = '../../globColour/L3m_{}01-{}{}__455288665_4_AVW-MODVIR_ZEU_MO_00.nc'.format(ym,ym,end_day)
        if pth.isfile(fname):
            print(fname)
            ds = xr.open_dataset(fname)
            zeu_obs = ds.ZEU_mean.values
            zeu_obs[zeu_obs==-999] = np.nan
            zeu_avg[yr-2008,m-1] = np.nanmean(zeu_obs[263:320,866:1106])
            del zeu_obs
            xr.Dataset.close(ds)
        fname = '../../globColour/L3m_{}01-{}{}__455288665_4_AVW-MERMODVIR_ZEU_MO_00.nc'.format(ym,ym,end_day)
        if pth.isfile(fname):
            print(fname)
            ds = xr.open_dataset(fname)
            zeu_obs = ds.ZEU_mean.values
            zeu_obs[zeu_obs==-999] = np.nan
            zeu_avg[yr-2008,m-1] = np.nanmean(zeu_obs[263:320,866:1106])
            del zeu_obs
            xr.Dataset.close(ds)
        fname = '../../globColour/L3m_{}01-{}{}__455288665_4_AVW-MERMOD_ZEU_MO_00.nc'.format(ym,ym,end_day)
        if pth.isfile(fname):
            print(fname)
            ds = xr.open_dataset(fname)
            zeu_obs = ds.ZEU_mean.values
            zeu_obs[zeu_obs==-999] = np.nan
            zeu_avg[yr-2008,m-1] = np.nanmean(zeu_obs[263:320,866:1106])
            del zeu_obs
            xr.Dataset.close(ds)
        fname = '../../globColour/L3m_{}01-{}{}__455288665_4_AVW-MERMODSWF_ZEU_MO_00.nc'.format(ym,ym,end_day)
        if pth.isfile(fname):
            print(fname)
            ds = xr.open_dataset(fname)
            zeu_obs = ds.ZEU_mean.values
            zeu_obs[zeu_obs==-999] = np.nan
            zeu_avg[yr-2008,m-1] = np.nanmean(zeu_obs[263:320,866:1106])
            del zeu_obs
            xr.Dataset.close(ds)
                   
zeu_ts=np.zeros((85))   
zeu_avg[zeu_avg==0] = np.nan
zeu_ts[1:85]=np.ravel(zeu_avg)
zeu_ts[0]=zeu_ts[1]

### ------ ###

fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),0,1.3,color=(0.9,0.9,0.9))
#fig=plt.fill_between(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(ice_blue[:,0:100,300:400],2),1),60000)
pbar=plt.bar(np.linspace(0.5,83.5,84),np.nanmean(np.nanmean(chl[60:,0:100,300:400],2),1),color=(178/255,223/255,138/255),width=1.05)
#pbar=plt.bar(np.linspace(0,84,84),np.nanmean(np.nanmean(chl[60:,:,:],2),1),color=(51/255,160/255,44/255),width=1)
#pbar=plt.plot(np.linspace(0,84,84),np.nanmedian(np.nanmedian(chl[60:,0:100,300:400],2),1))
#plt.xlim(-0.5,84.5)
#fig=plt.fill_between(np.linspace(48,84,37),0.6,color=(0.9,0.9,0.9))
plt.margins(x=0)
plt.ylim(0,1.3)
plt.grid(axis='x')
#plt.xticks(np.linspace(0,84,15),['','08','','09','','10','','11','','12','','13','','14',''],fontsize=40)
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.ylabel('Chlorophyll (mg/m$^3$)',fontsize=60)
plt.yticks(fontsize=40)
plt.savefig('chlorophyll_pib.png')

import cmocean

ds = xr.open_dataset('hol_green/surf/green_meltrate.nc')
mlt_green = ds.surfDiag.values

new_map=((1,1,1),(0.9,0.9,0.9),(0.7,0.7,0.7),(0,104/255,55/255),(49/255,163/255,84/255),(120/255,198/255,121/255),(194/255,230/255,153/255),(255/255,255/255,204/255))
#new_map=((0,104/255,55/255),(49/255,163/255,84/255),(120/255,198/255,121/255),(194/255,230/255,153/255),(255/255,255/255,204/255))
ice_dec = np.array([np.nanmean(ice_cp[0:1,:,:],0),np.nanmean(ice_cp[12:13,:,:],0),np.nanmean(ice_cp[24:25,:,:],0),np.nanmean(ice_cp[36:37,:,:],0),np.nanmean(ice_cp[48:49,:,:],0),np.nanmean(ice_cp[60:61,:,:],0),np.nanmean(ice_cp[72:73,:,:],0)])
chl_dec = np.array([np.nanmean(zeu_cp[60:61,:,:],0),np.nanmean(zeu_cp[72:73,:,:],0),np.nanmean(zeu_cp[84:85,:,:],0),np.nanmean(zeu_cp[96:97,:,:],0),np.nanmean(zeu_cp[108:109,:,:],0),np.nanmean(zeu_cp[120:121,:,:],0),np.nanmean(zeu_cp[132:133,:,:],0)])
chl_dec[chl_dec<40]=40
chl_dec[chl_dec>60]=60
chl_dec[ice_dec>0.15]=-999
zeu_avg = np.nanmean(chl_dec,0)
zeu_avg[np.nansum(zeu_cp,0)==0]=37
zeu_avg[np.nansum(mlt_green,0)<0]=33
fig=plt.figure(figsize=(20,10))
pcol=plt.contourf(lon,lat,zeu_avg,np.linspace(28,60,9),colors=new_map,extend='both')
#pcol=plt.pcolormesh(np.squeeze(bathy[0,:,:]))
print('GOT HERE NOW')
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=20)
#pcon=plt.contour(lon,lat,zeu_avg,[0,32,60],colors=[166/255,97/255,26/255])
plt.clim(18,60)
plt.gca().yaxis.tick_right()
plt.yticks([-74,-72,-70,-68,-66,-64],['74','72','70','68','66','64'],fontsize=15)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=15)
plt.savefig('zeu_dec.png')

### --- Here is the Supplementary figure 1 --- ###

ice_dec = np.squeeze(ice_cp[49,:,:])
chl_dec = np.squeeze(zeu[109,:,:])
#chl_dec[chl_dec==0]=100
#chl_dec[chl_dec>60]=np.nan
#chl_dec[ice_dec>0.15]=-999
zeu_feb = chl_dec
zeu_ist_lon = 360-np.array([116.78,117.67,117.72,116.5,116.5,115,113.5,113,114,116.37,117.58,113.07,105,102.1])
zeu_ist_lat = [-71.66,-71.79,-72.39,-72.85,-73.5,-73.25,-72.99,-73.5,-73.5,-72.46,-72.93,-73.82,-74.37,-74.86]
zeu_ist = [27.9,25.6,24.7,22.8,20.9,14.8,18.1,18.1,17.6,30.7,15.3,14.8,37.7,40]
#zeu_feb[np.nansum(zeu_cp,0)==0]=37
#zeu_feb[np.nansum(mlt_green,0)<0]=33
fig=plt.figure(figsize=(40,20))
pcol=plt.contourf(lon,lat,zeu_feb,[0,10,20,30,40,50,60,70,80],cmap='cmo.matter_r',extend='max')
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=50)
plt.clim(10,80)
#pcon=plt.contour(lon,lat,zeu_feb,[0,32,60],colors=[166/255,97/255,26/255])
#plt.yticks([-74,-72,-70,-68,-66,-64],['74','72','70','68','66','64'],fontsize=15)
#plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=15)
print(new_map)
plt.scatter(zeu_ist_lon,zeu_ist_lat,c=zeu_ist,s=700,linewidths=5,edgecolor=[1,1,1],cmap='cmo.matter_r')
plt.clim(10,80)
plt.ylim(-75.5,-69)
plt.xlim(235,265)
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=60)
plt.xticks([240,250,260],['120','110','100'],fontsize=60)
plt.ylabel('$\degree$S',fontsize=60)
plt.xlabel('$\degree$W',fontsize=60)
#pcol=plt.contourf(lon,lat,fig_mask,[0,1,2,3],colors=mask_cmap,alpha=0.5)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
plt.title('Euphotic depth (m), February 2012',fontsize=60)
plt.savefig('zeu_insitu.png')

print('size={}'.format(np.size(zeu_ts)))
# open water area
fig=plt.figure(figsize=(40,10))
#plin=plt.plot(np.linspace(0,84,84),np.nansum(1e-6*np.nansum(opn_blue[:,0:100,300:400],2),1),'b',linewidth=5)
fig=plt.fill_between(np.linspace(60,84,37),0,170,color=(0.9,0.9,0.9))#np.nanmean(np.nanmean(zeu[108:,0:100,300:400],2),1),160,color=(0.9,0.9,0.9))
#pbar=plt.bar(np.linspace(0,84,84),np.nanmean(np.nanmean(zeu[60:,:,:],2),1),color=(51/255,160/255,44/255),width=1)
fig=plt.bar(np.linspace(0.5,83.5,84),np.nanmean(np.nanmean(zeu[60:144,0:100,300:400],2),1),color=(178/255,223/255,138/255),width=1.05)
#pbar=plt.bar(np.linspace(0.5,84.5,84),np.nanmean(np.nanmean(zeu[60:144,0:100,300:400],2),1),color=(178/255,223/255,138/255),width=1.05)
#pbar=plt.plot(np.linspace(0,84,84),np.nanmedian(np.nanmedian(zeu[60:,0:100,300:400],2),1))
scat=plt.plot(np.linspace(0.5,83.5,84),zeu_ts[1:],linewidth=5,color='k',alpha=0.7)
plt.ylim(0,170)
plt.grid(axis='x')
plt.gca().invert_yaxis()
#plt.xlim(0,84)
plt.margins(x=0)
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.ylabel('Euphotic depth (m)',fontsize=60)
plt.yticks(fontsize=40)
plt.savefig('euphotic_pib.png')
print(np.linspace(0,83,84))
