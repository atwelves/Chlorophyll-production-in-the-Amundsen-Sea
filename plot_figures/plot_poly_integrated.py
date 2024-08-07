# script to process and plot sea ice data

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

# convert temperature to energy density

spec = 4.186 * 1.02*1e6
 
thresh = 0.1
# seconds in day

day = 86400

# read netcdf file

ds = xr.open_dataset('hol_green/surf/tflux_green.nc')
green_flx = ds.surfDiag.values
green_flx[green_flx==0] = np.nan

ds = xr.open_dataset('hol_blue/surf/tflux_blue.nc')
blue_flx = ds.surfDiag.values
blue_flx[blue_flx==0] = np.nan

ds = xr.open_dataset('hol_green/seaIce/ice_green.nc')
green_ice = ds.seaIceDiag.values
green_ice[green_flx==0] = np.nan

ds = xr.open_dataset('hol_blue/seaIce/ice_blue.nc')
blue_ice = ds.seaIceDiag.values
blue_ice[blue_flx==0] = np.nan

green_opn = np.copy(green_flx)
blue_opn = np.copy(blue_flx)
green_cov = np.copy(green_flx)
blue_cov = np.copy(blue_flx)

ds = xr.open_dataset('hol_green/surfHeat/latent_green.nc')
green_lat = ds.surfHeatDiag.values
green_lat[green_flx==0] = np.nan

ds = xr.open_dataset('hol_blue/surfHeat/latent_blue.nc')
blue_lat = ds.surfHeatDiag.values
blue_lat[blue_flx==0] = np.nan

ds = xr.open_dataset('hol_green/surfHeat/sensible_green.nc')
green_sen = ds.surfHeatDiag.values
green_sen[green_flx==0] = np.nan

ds = xr.open_dataset('hol_blue/surfHeat/sensible_blue.nc')
blue_sen = ds.surfHeatDiag.values
blue_sen[blue_flx==0] = np.nan

ds = xr.open_dataset('hol_green/surfHeat/longwave_green.nc')
green_lw = -ds.surfHeatDiag.values
green_lw[green_flx==0] = np.nan

ds = xr.open_dataset('hol_blue/surfHeat/longwave_blue.nc')
blue_lw = -ds.surfHeatDiag.values
blue_lw[blue_flx==0] = np.nan

ds = xr.open_dataset('hol_green/surfHeat/shortwave_green.nc')
green_sw = -ds.surfHeatDiag.values
green_sw[green_flx==0] = np.nan

ds = xr.open_dataset('hol_blue/surfHeat/shortwave_blue.nc')
blue_sw = -ds.surfHeatDiag.values
blue_sw[blue_flx==0] = np.nan

poly_blue  = np.zeros((84,384,600))
poly_green = np.zeros((84,384,600))
miz_blue  = np.zeros((84,384,600))
miz_green = np.zeros((84,384,600))

# seperate open and ice covered components
green_opn[green_ice>thresh]=0
green_cov[green_ice<thresh]=0
blue_opn[blue_ice>thresh]=0
blue_cov[blue_ice<thresh]=0

poly_blue[blue_ice<thresh] = 1
poly_green[green_ice<thresh] = 1
poly_common = np.multiply(poly_blue,poly_green)
miz_blue[blue_ice>thresh] = 1
miz_green[green_ice>thresh] = 1

green_lat[green_ice>thresh]=0
green_sen[green_ice>thresh]=0
blue_lat[blue_ice>thresh]=0
blue_sen[blue_ice>thresh]=0
green_lw[green_ice>thresh]=0
green_sw[green_ice>thresh]=0
blue_lw[blue_ice>thresh]=0
blue_sw[blue_ice>thresh]=0

diff_opn = green_opn-blue_opn
diff_cov = green_cov-blue_cov

diff_lat = green_lat-blue_lat
diff_sen = green_sen-blue_sen
diff_lw = green_lw-blue_lw
diff_sw = green_sw-blue_sw

import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

rac = mds.rdmds('data/RAC')
rac=np.tile(rac,(84,1,1))
#rac = rac/np.nansum(rac[0,:,:])

diff_opn = np.multiply(diff_opn,rac)
diff_cov = np.multiply(diff_cov,rac)
diff_lat = np.multiply(diff_lat,rac)
diff_sen = np.multiply(diff_sen,rac)
diff_lw = np.multiply(diff_lw,rac)
diff_sw = np.multiply(diff_sw,rac)

open_green = np.multiply(green_opn,rac)
open_blue = np.multiply(blue_opn,rac)

cov_green = np.multiply(green_cov,rac)
cov_blue = np.multiply(blue_cov,rac)

lw_green = np.multiply(green_lw,rac)
lw_blue = np.multiply(blue_lw,rac)
sw_green = np.multiply(green_sw,rac)
sw_blue = np.multiply(blue_sw,rac)

lat_green = np.multiply(green_lat,rac)
lat_blue = np.multiply(blue_lat,rac)
sen_green = np.multiply(green_sen,rac)
sen_blue = np.multiply(blue_sen,rac)

poly_blue = np.multiply(poly_blue,rac)
poly_green = np.multiply(poly_green,rac)
poly_common = np.multiply(poly_common,rac)
miz_blue = np.multiply(miz_blue,rac)
miz_green = np.multiply(miz_green,rac)

### --- here also want component of air-sea flux, averaged annualy
#diff_lat = np.reshape(diff_lat, (7,12,384,600))
#diff_sen = np.reshape(diff_sen, (7,12,384,600))
#diff_lw = np.reshape(diff_lw, (7,12,384,600))
#diff_sw = np.reshape(diff_sw, (7,12,384,600))

#diff_lat = np.nanmean(diff_lat,0)
#diff_sen = np.nanmean(diff_sen,0)
#diff_lw = np.nanmean(diff_lw,0)
#diff_sw = np.nanmean(diff_sw,0)

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
bathy[bathy==0] = np.nan
print('got here')

# do masking
diff_lat[np.isnan(bathy)]=np.nan
diff_sen[np.isnan(bathy)]=np.nan
diff_lw[np.isnan(bathy)]=np.nan
diff_sw[np.isnan(bathy)]=np.nan
diff_opn[np.isnan(bathy)]=np.nan
diff_cov[np.isnan(bathy)]=np.nan

open_blue[np.isnan(bathy)] = np.nan
open_green[np.isnan(bathy)] = np.nan

cov_blue[np.isnan(bathy)] = np.nan
cov_green[np.isnan(bathy)] = np.nan

lw_blue[np.isnan(bathy)] = np.nan
lw_green[np.isnan(bathy)] = np.nan
sw_blue[np.isnan(bathy)] = np.nan
sw_green[np.isnan(bathy)] = np.nan

lat_blue[np.isnan(bathy)] = np.nan
lat_green[np.isnan(bathy)] = np.nan
sen_blue[np.isnan(bathy)] = np.nan
sen_green[np.isnan(bathy)] = np.nan

poly_blue[np.isnan(bathy)] = np.nan
poly_green[np.isnan(bathy)] = np.nan
poly_common[np.isnan(bathy)] = np.nan
miz_blue[np.isnan(bathy)] = np.nan
miz_green[np.isnan(bathy)] = np.nan

lw_green_uni = np.copy(lw_green)
lw_blue_uni  = np.copy(lw_blue )
sen_green_uni = np.copy(sen_green)
sen_blue_uni  = np.copy(sen_blue )
lat_green_uni = np.copy(lat_green)
lat_blue_uni  = np.copy(lat_blue )
sw_green_uni = np.copy(sw_green)
sw_blue_uni  = np.copy(sw_blue )


lw_green[poly_common==0]=np.nan
lw_blue[poly_common==0]=np.nan
sen_green[poly_common==0]=np.nan
sen_blue[poly_common==0]=np.nan
lat_green[poly_common==0]=np.nan
lat_blue[poly_common==0]=np.nan

#lw_green_uni = np.copy(lw_green)
#lw_blue_uni  = np.copy(lw_blue )

lw_green_uni[poly_common>0]=np.nan
lw_blue_uni[poly_common>0]=np.nan
sen_green_uni[poly_common>0]=np.nan
sen_blue_uni[poly_common>0]=np.nan
lat_green_uni[poly_common>0]=np.nan
lat_blue_uni[poly_common>0]=np.nan
sw_green_uni[poly_common>0]=np.nan
sw_blue_uni[poly_common>0]=np.nan


diff_lat = np.nansum(np.nansum(diff_lat,2),1)
diff_sen = np.nansum(np.nansum(diff_sen,2),1)
diff_lw = np.nansum(np.nansum(diff_lw,2),1)
diff_sw = np.nansum(np.nansum(diff_sw,2),1)
diff_opn = np.nansum(np.nansum(diff_opn,2),1)
diff_cov = np.nansum(np.nansum(diff_cov,2),1)

open_blue = np.nansum(np.nansum(open_blue,2),1)
open_green = np.nansum(np.nansum(open_green,2),1)

cov_blue = np.nansum(np.nansum(cov_blue,2),1)
cov_green = np.nansum(np.nansum(cov_green,2),1)

lw_blue = np.nansum(np.nansum(lw_blue,2),1)
lw_green = np.nansum(np.nansum(lw_green,2),1)
sw_blue = np.nanmean(np.nanmean(sw_blue,2),1)
sw_green = np.nanmean(np.nanmean(sw_green,2),1)

lw_blue_uni  = np.nansum(np.nansum(lw_blue_uni,2),1)
lw_green_uni = np.nansum(np.nansum(lw_green_uni,2),1)
sen_blue_uni  = np.nansum(np.nansum(sen_blue_uni,2),1)
sen_green_uni = np.nansum(np.nansum(sen_green_uni,2),1)
lat_blue_uni  = np.nansum(np.nansum(lat_blue_uni,2),1)
lat_green_uni = np.nansum(np.nansum(lat_green_uni,2),1)
sw_blue_uni  = np.nansum(np.nansum(sw_blue_uni,2),1)
sw_green_uni = np.nansum(np.nansum(sw_green_uni,2),1)

lat_blue = np.nansum(np.nansum(lat_blue,2),1)
lat_green = np.nansum(np.nansum(lat_green,2),1)
sen_blue = np.nansum(np.nansum(sen_blue,2),1)
sen_green = np.nansum(np.nansum(sen_green,2),1)

poly_blue = np.nansum(np.nansum(poly_blue,2),1)
poly_green = np.nansum(np.nansum(poly_green,2),1)
poly_common = np.nansum(np.nansum(poly_common,2),1)
miz_blue = np.nansum(np.nansum(miz_blue,2),1)
miz_green = np.nansum(np.nansum(miz_green,2),1)

# convert units

diff_opn = diff_opn*3.154e-11
diff_cov = diff_cov*3.154e-11

open_green = open_green*3.154e-11
open_blue  = open_blue*3.154e-11

cov_green = cov_green*3.154e-11
cov_blue  = cov_blue*3.154e-11

lw_green = lw_green*3.154e-11
lw_blue  = lw_blue*3.154e-11
sw_green = sw_green*3.154e-11
sw_blue  = sw_blue*3.154e-11
lat_green = lat_green*3.154e-11
lat_blue  = lat_blue*3.154e-11
sen_green = sen_green*3.154e-11
sen_blue  = sen_blue*3.154e-11

lw_green_uni = lw_green_uni*3.154e-11
lw_blue_uni  = lw_blue_uni*3.154e-11
sw_green_uni = sw_green_uni*3.154e-11
sw_blue_uni  = sw_blue_uni*3.154e-11
lat_green_uni = lat_green_uni*3.154e-11
lat_blue_uni  = lat_blue_uni*3.154e-11
sen_green_uni = sen_green_uni*3.154e-11
sen_blue_uni  = sen_blue_uni*3.154e-11

#poly_green = poly_green*3.154e-11
#poly_blue  = poly_blue*3.154e-11

### --- plot figures --- ###

# Figure 6c
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),-118,4,color=(0.9,0.9,0.9))
plt.fill_between(np.linspace(0,84,84),lw_green_uni-lw_blue_uni,color=(217/255,95/255,2/255),alpha=1)
plt.fill_between(np.linspace(0,84,84),lw_green_uni-lw_blue_uni,lw_green+lw_green_uni-lw_blue-lw_blue_uni,color=(117/255,112/255,179/255),alpha=1)
plt.plot(np.linspace(0,84,84),lw_green+lw_green_uni-lw_blue-lw_blue_uni,color='k',linewidth=5)
plt.yticks(fontsize=60)
plt.xlim(0,84)
plt.ylim(-118,4)
plt.grid(axis='x')
plt.grid(axis='y')
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Anomaly (EJ yr⁻¹)',fontsize=70)
plt.savefig('lw.png')

print('lw com area: {}'.format(np.nansum(lw_green-lw_blue)))
print('lw new area: {}'.format(np.nansum(lw_green_uni-lw_blue_uni)))

# Figure 6a
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),-118,4,color=(0.9,0.9,0.9))
plt.fill_between(np.linspace(0,84,84),sen_green_uni-sen_blue_uni,color=(217/255,95/255,2/255),alpha=1)
plt.fill_between(np.linspace(0,84,84),sen_green_uni-sen_blue_uni,sen_green+sen_green_uni-sen_blue-sen_blue_uni,color=(117/255,112/255,179/255),alpha=1)
plt.plot(np.linspace(0,84,84),sen_green+sen_green_uni-sen_blue-sen_blue_uni,color='k',linewidth=5)
plt.xlim(0,84)
plt.ylim(-118,4)
plt.yticks(fontsize=60)
plt.grid()
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Anomaly (EJ yr⁻¹)',fontsize=70)
plt.savefig('sen.png')
print('sen com area: {}'.format(np.nansum(sen_green-sen_blue)))
print('sen new area: {}'.format(np.nansum(sen_green_uni-sen_blue_uni)))

# Figure 6b
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),-118,4,color=(0.9,0.9,0.9))
plt.fill_between(np.linspace(0,84,84),lat_green_uni-lat_blue_uni,color=(217/255,95/255,2/255),alpha=1)
plt.fill_between(np.linspace(0,84,84),lat_green_uni-lat_blue_uni,lat_green_uni+lat_green-lat_blue_uni-lat_blue,color=(117/255,112/255,179/255),alpha=1)
plt.plot(np.linspace(0,84,84),lat_green+lat_green_uni-lat_blue-lat_blue_uni,color='k',linewidth=5)
plt.xlim(0,84)
plt.ylim(-118,4)
plt.yticks(fontsize=60)
plt.grid()
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Anomaly (EJ yr⁻¹)',fontsize=70)
plt.savefig('lat.png')

print('lat com area: {}'.format(np.nansum(lat_green-lat_blue)))
print('lat new area: {}'.format(np.nansum(lat_green_uni-lat_blue_uni)))

# Figure 6e
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),-10,300,color=(0.9,0.9,0.9))
plt.fill_between(np.linspace(0,84,84),sw_green_uni-sw_blue_uni,color=(217/255,95/255,2/255),alpha=1)
plt.fill_between(np.linspace(0,84,84),sw_green_uni-sw_blue_uni,sw_green_uni+sw_green-sw_blue_uni-sw_blue,color=(117/255,112/255,179/255),alpha=1)
plt.plot(np.linspace(0,84,84),sw_green+sw_green_uni-sw_blue-sw_blue_uni,color='k',linewidth=5)
plt.xlim(0,84)
plt.ylim(-10,280)
plt.yticks(fontsize=60)
plt.grid()
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Anomaly (EJ yr⁻¹)',fontsize=70)
plt.savefig('sw.png')

print('sw com area: {}'.format(np.nansum(sw_green-sw_blue)))
print('sw new area: {}'.format(np.nansum(sw_green_uni-sw_blue_uni)))

total_com = lw_green + sen_green + lat_green - lw_blue - sen_blue - lat_blue
total_uni = lw_green_uni + sen_green_uni + lat_green_uni - lw_blue_uni - sen_blue_uni - lat_blue_uni

# Figure 6f
fig=plt.figure(figsize=(40,10))
fig=plt.fill_between(np.linspace(60,84,37),-300,10,color=(0.9,0.9,0.9))
plt.fill_between(np.linspace(0,84,84),total_uni,color=(217/255,95/255,2/255),alpha=1)
plt.fill_between(np.linspace(0,84,84),total_uni,total_com+total_uni,color=(117/255,112/255,179/255),alpha=1)
plt.plot(np.linspace(0,84,84),total_com + total_uni,color='k',linewidth=5)
plt.xlim(0,84)
plt.ylim(-280,10)
plt.yticks(fontsize=60)
plt.grid()
plt.xticks(np.linspace(0,84,29),['                   2008  ','','','',
                                 '              |    2009  ','','','',
                                 '              |    2010  ','','','',
                                 '              |    2011  ','','','',
                                 '              |    2012  ','','','',
                                 '              |    2013  ','','','',
                                 '              |    2014  ','','','',''],fontsize=60)
plt.ylabel('Anomaly (EJ yr⁻¹)',fontsize=70)
plt.savefig('total_outflux.png')

print('sen')
print((np.nanmean(np.divide(sen_green,poly_green))-np.nanmean(np.divide(sen_blue,poly_blue))))#*100/np.nanmean(np.divide(sen_blue,poly_blue)))
print('lat')
print((np.nanmean(np.divide(lat_green,poly_green))-np.nanmean(np.divide(lat_blue,poly_blue))))#*100/np.nanmean(np.divide(lat_blue,poly_blue)))
print('lw')
print((np.nanmean(np.divide(lw_green,poly_green))-np.nanmean(np.divide(lw_blue,poly_blue))))#*100/np.nanmean(np.divide(lw_blue,poly_blue)))
print('open water')
print((np.nanmean(poly_green)-np.nanmean(poly_blue))*100/np.nanmean(poly_blue))

ct=60
print('longwave:{}'.format(100*np.nanmean(lw_green[:ct]-lw_blue[:ct])/np.nanmean(lw_blue[:ct])))
print('shortwave:{}'.format(100*np.nanmean(sw_green[:ct]-sw_blue[:ct])/np.nanmean(sw_blue[:ct])))
print('sensible:{}'.format(100*np.nanmean(sen_green[:ct]-sen_blue[:ct])/np.nanmean(sen_blue[:ct])))
print('latent:{}'.format(100*np.nanmean(lat_green[:ct]-lat_blue[:ct])/np.nanmean(lat_blue[:ct])))
print('sum:{}'.format(np.nanmean(lw_green[:ct]+sen_green[:ct]+lat_green[:ct]-lw_blue[:ct]-sen_blue[:ct]-lat_blue[:ct])))
print('open water:{}'.format(np.nanmean(open_green[:ct]-open_blue[:ct])))
print('ice covered:{}'.format(np.nanmean(cov_green[:ct]-cov_blue[:ct])))
print('total:{}'.format(np.nanmean(open_green[:ct]+cov_green[:ct]-cov_blue[:ct]-open_blue[:ct])))

### --- pie charts --- ###

# Figure 6d
fig, ax = plt.subplots()
sizes = [-np.nansum(lat_green-lat_blue),-np.nansum(lat_green_uni-lat_blue_uni)]
#labels = "$F_{com}$","$F_{miz}$"
rds = (-np.nansum(lat_green-lat_blue)-np.nansum(lat_green_uni-lat_blue_uni))*1e-3
ax.pie(sizes,radius=rds,colors=((117/255,112/255,179/255),(217/255,95/255,2/255)),textprops={'fontsize': 30})
plt.title('Latent',fontsize=35)
plt.savefig('lat_com_v_uni.png')

fig, ax = plt.subplots()
sizes = [-np.nansum(sen_green-sen_blue),-np.nansum(sen_green_uni-sen_blue_uni)]
#labels = "common","unique"
rds = (-np.nansum(sen_green-sen_blue)-np.nansum(sen_green_uni-sen_blue_uni))*1e-3
ax.pie(sizes,radius=rds,colors=((117/255,112/255,179/255),(217/255,95/255,2/255)),textprops={'fontsize': 30})
plt.title('Sensible',fontsize=35)
plt.savefig('sen_com_v_uni.png')

fig, ax = plt.subplots()
sizes = [-np.nansum(lw_green-lw_blue),-np.nansum(lw_green_uni-lw_blue_uni)]
labels = "common","unique"
rds = (-np.nansum(lw_green-lw_blue)-np.nansum(lw_green_uni-lw_blue_uni))*1e-3 
#labels = "{:.1f}%".format(100*net/srf), "{:.1f}%".format(100*res/srf)
ax.pie(sizes,radius=rds,colors=((117/255,112/255,179/255),(217/255,95/255,2/255)),textprops={'fontsize': 25})
plt.title('Longwave',fontsize=35)
plt.savefig('lw_com_v_uni.png')
