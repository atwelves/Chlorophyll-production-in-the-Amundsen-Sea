# script to process and plot bio data

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

# specify threshold for sea ice concentration that indicates coverage
threshold = 0.1

### ---  read netCDFs --- ###

# (surface) chlorophyll concentration from GREEN
ds = xr.open_dataset('hol_green/phyto/chlorophyll.nc')
chl = ds.phytoDiag.values

# euphotic depth from GREEN
ds = xr.open_dataset('hol_green/lightSurf/euphotic_depth.nc')
lon = ds.x.values
lat = ds.y.values
zeu = ds.lightSurfDiag.values

# sea ice concentration from GREEN
ds = xr.open_dataset('hol_green/seaIce/ice_green.nc')
ice = ds.seaIceDiag.values

# sea ice concentration 
ds = xr.open_dataset('hol_blue/seaIce/ice_blue.nc')
ice_blue = ds.seaIceDiag.values

# set land mask
chl[zeu==0]=np.nan
zeu[zeu==0]=np.nan

### ------ ###

import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

# read in model grid dimensions and tile over months
rac = mds.rdmds('data/RAC')
rac = np.tile(rac,(144,1,1))

### --- plot figures --- ###

#restrict to continental shelf

# read in bathymetry and consider only continental shelf
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

# tile bathymetry over months in model run
bathy = np.tile(bathy,(144,1,1))

# create maps and colormaps for plotting
fig_mask = np.zeros((384,600))
fig_mask[np.squeeze(bathy[0,:,:])==0] = 1.5
fig_mask[np.isnan(np.squeeze(bathy[0,:,:]))] = 2.5
ice_cmap = ((118/255,42/255,131/255),(153/255,112/255,171/255),(194/255,165/255,207/255),(231/255,212/255,232/255),(217/255,240/255,211/255),(166/255,219/255,160/255),(90/255,174/255,97/255),(27/255,120/255,55/255))
mask_cmap = ((1,1,1),(0,0,0),(0.7,0.7,0.7))

del off_shelf

zeu_cp = np.copy(zeu)
ice_cp = np.copy(ice)

chl[np.isnan(bathy)]=np.nan
zeu[np.isnan(bathy)]=np.nan
bathy = np.squeeze(bathy[60:,:,:])
ice[np.isnan(bathy)]=np.nan
ice_blue[np.isnan(bathy)] = np.nan

# now for obs***

import calendar
from calendar import monthrange
import os.path as pth

# nested loops

zeu_avg = np.zeros((7,12))

### --- Read in satellite observations of chlorophyll in Amundsen Sea from GlobColour database --- ###

for yr in range(2008,2015):
    year = '{}'.format(yr)
    for m in range(1,13):
        month = '{:02d}'.format(m)
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

# declare arrays for time series
zeu_ts=np.zeros((85))   
zeu_avg[zeu_avg==0] = np.nan
zeu_ts[1:85]=np.ravel(zeu_avg)
zeu_ts[0]=zeu_ts[1]

import cmocean

ds = xr.open_dataset('hol_green/surf/green_meltrate.nc')
mlt_green = ds.surfDiag.values

# Figure 1a
new_map=((1,1,1),(0.9,0.9,0.9),(0.7,0.7,0.7),(0,104/255,55/255),(49/255,163/255,84/255),(120/255,198/255,121/255),(194/255,230/255,153/255),(255/255,255/255,204/255))
ice_dec = np.array([np.nanmean(ice_cp[0:1,:,:],0),np.nanmean(ice_cp[12:13,:,:],0),np.nanmean(ice_cp[24:25,:,:],0),np.nanmean(ice_cp[36:37,:,:],0),np.nanmean(ice_cp[48:49,:,:],0),np.nanmean(ice_cp[60:61,:,:],0),np.nanmean(ice_cp[72:73,:,:],0)])
chl_dec = np.array([np.nanmean(zeu_cp[60:61,:,:],0),np.nanmean(zeu_cp[72:73,:,:],0),np.nanmean(zeu_cp[84:85,:,:],0),np.nanmean(zeu_cp[96:97,:,:],0),np.nanmean(zeu_cp[108:109,:,:],0),np.nanmean(zeu_cp[120:121,:,:],0),np.nanmean(zeu_cp[132:133,:,:],0)])
chl_dec[chl_dec<40]=40
chl_dec[chl_dec>60]=60
chl_dec[ice_dec>0.1]=-999
zeu_avg = np.nanmean(chl_dec,0)
zeu_avg[np.nansum(zeu_cp,0)==0]=37
zeu_avg[np.nansum(mlt_green,0)<0]=33
fig=plt.figure(figsize=(20,10))
pcol=plt.contourf(lon,lat,zeu_avg,np.linspace(28,60,9),colors=new_map,extend='both')
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=20)
plt.clim(18,60)
plt.gca().yaxis.tick_right()
plt.yticks([-74,-72,-70,-68,-66,-64],['74','72','70','68','66','64'],fontsize=15)
plt.xticks([230,240,250,260,270],['130','120','110','100','90'],fontsize=15)
plt.savefig('zeu_dec.png')

# Supplementary Figure 2b
ice_dec = np.squeeze(ice_cp[49,:,:])
zeu_feb = np.squeeze(zeu[109,:,:])
zeu_ist_lon = 360-np.array([116.78,117.67,117.72,116.5,116.5,115,113.5,113,114,116.37,117.58,113.07,105,102.1])
zeu_ist_lat = [-71.66,-71.79,-72.39,-72.85,-73.5,-73.25,-72.99,-73.5,-73.5,-72.46,-72.93,-73.82,-74.37,-74.86]
zeu_ist = [27.9,25.6,24.7,22.8,20.9,14.8,18.1,18.1,17.6,30.7,15.3,14.8,37.7,40]
fig=plt.figure(figsize=(40,20))
pcol=plt.contourf(lon,lat,zeu_feb,[0,10,20,30,40,50,60,70,80],cmap='cmo.matter_r',extend='max')
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=50)
plt.clim(10,80)
plt.scatter(zeu_ist_lon,zeu_ist_lat,c=zeu_ist,s=700,linewidths=5,edgecolor=[1,1,1],cmap='cmo.matter_r')
plt.clim(10,80)
plt.ylim(-75.5,-69)
plt.xlim(235,265)
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=60)
plt.xticks([240,250,260],['120','110','100'],fontsize=60)
plt.ylabel('$\degree$S',fontsize=60)
plt.xlabel('$\degree$W',fontsize=60)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
plt.title('Euphotic depth (m), February 2012',fontsize=60)
plt.savefig('zeu_insitu.png')

# Supplementary Figure 2a
ice_dec = np.squeeze(ice_cp[49,:,:])
zeu_feb = np.squeeze(chl[109,:,:])
zeu_ist_lon = 360-np.array([116.78,117.67,117.72,116.5,116.5,115,113.5,113,114,116.37,117.58,113.07,105,102.1])
zeu_ist_lat = [-71.66,-71.79,-72.39,-72.85,-73.5,-73.25,-72.99,-73.5,-73.5,-72.46,-72.93,-73.82,-74.37,-74.86]
zeu_ist = [3.21,3.03,3.29,3.41,3.84,2.95,4.46,6.20,5.72,2.82,2.00,4.89,2.33,1.56]
fig=plt.figure(figsize=(40,20))
pcol=plt.contourf(lon,lat,zeu_feb,[0,0.5,1,1.5,2.0,2.5,3,3.5,4],cmap='cmo.matter_r',extend='max')
cbar=plt.colorbar(extend='both')
cbar.ax.tick_params(labelsize=50)
plt.clim(0,4)
plt.scatter(zeu_ist_lon,zeu_ist_lat,c=zeu_ist,s=700,linewidths=5,edgecolor=[1,1,1],cmap='cmo.matter_r')
plt.clim(0,4)
plt.ylim(-75.5,-69)
plt.xlim(235,265)
plt.yticks([-74,-72,-70],['74','72','70'],fontsize=60)
plt.xticks([240,250,260],['120','110','100'],fontsize=60)
plt.ylabel('$\degree$S',fontsize=60)
plt.xlabel('$\degree$W',fontsize=60)
pcol=plt.contour(lon,lat,fig_mask,[0,1,2,3],linewidths=3,colors='black')
plt.title('Surface chlorophyll (mg m⁻³), February 2012',fontsize=60)
plt.savefig('chl_insitu.png')

rac[np.isnan(chl)] = np.nan
chl = np.multiply(chl,rac)

# Figure 2a
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),0,1.6,color=(0.9,0.9,0.9))
pbar=plt.bar(np.linspace(0.5,83.5,84),np.nanmean(np.nanmean(chl[60:,0:100,300:400],2),1)/np.nanmean(rac[0,0:100,300:400]),color=(178/255,223/255,138/255),width=1.05)
plt.margins(x=0)
plt.ylim(0,1.6)
plt.yticks(np.linspace(0,1.5,4),['  0.0','  0.5','  1.0','  1.5'])
plt.grid(axis='x')
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Chlorophyll (mg m$^{-3}$)',fontsize=60)
plt.yticks(fontsize=60)
plt.savefig('chlorophyll_pib.png')

rac[np.isnan(zeu)] = np.nan
zeu = np.multiply(zeu,rac)

# Figure 2b
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),0,170,color=(0.9,0.9,0.9))
fig=plt.bar(np.linspace(0.5,83.5,84),np.nanmean(np.nanmean(zeu[60:144,0:100,300:400],2),1)/np.nanmean(rac[0,0:100,300:400]),color=(178/255,223/255,138/255),width=1.05)
scat=plt.plot(np.linspace(0.5,83.5,84),zeu_ts[1:],linewidth=5,color='k',alpha=0.7)
plt.ylim(0,170)
plt.yticks(np.linspace(0,150,4),['  0','  50','  100','  150'])
plt.grid(axis='x')
plt.gca().invert_yaxis()
plt.margins(x=0)
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Euphotic depth (m)',fontsize=60)
plt.yticks(fontsize=60)
plt.savefig('euphotic_pib.png')
