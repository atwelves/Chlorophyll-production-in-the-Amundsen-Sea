# script to process and plot surface heat flux outputs

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

### --- Set parameters --- ###

# convert temperature to energy density
spec = 4.186 * 1.02*1e6
 
# set threshold for ice cover: 10% for main figures, vary from 0% to 25% to get supplementary figures 
thresh = 0.10

# seconds in day
day = 86400

### ------ ###

### --- Read in model outputs --- ###

# read in total GREEN surface heat flux
ds = xr.open_dataset('hol_green/surf/tflux_green.nc')
green_flx = ds.surfDiag.values
green_flx[green_flx==0] = np.nan

# read in total BLUE surface heat flux
ds = xr.open_dataset('hol_blue/surf/tflux_blue.nc')
blue_flx = ds.surfDiag.values
blue_flx[blue_flx==0] = np.nan

# read in GREEN ice cover
ds = xr.open_dataset('hol_green/seaIce/ice_green.nc')
green_ice = ds.seaIceDiag.values
green_ice[green_flx==0] = np.nan

# read in BLUE ice cover
ds = xr.open_dataset('hol_blue/seaIce/ice_blue.nc')
blue_ice = ds.seaIceDiag.values
blue_ice[blue_flx==0] = np.nan

green_opn = np.copy(green_flx)
blue_opn = np.copy(blue_flx)
green_cov = np.copy(green_flx)
blue_cov = np.copy(blue_flx)

# read in GREEN latent heat flux
ds = xr.open_dataset('hol_green/surfHeat/latent_green.nc')
green_lat = ds.surfHeatDiag.values
green_lat[green_flx==0] = np.nan

# read in BLUE latent heat flux
ds = xr.open_dataset('hol_blue/surfHeat/latent_blue.nc')
blue_lat = ds.surfHeatDiag.values
blue_lat[blue_flx==0] = np.nan

# read in GREEN sensible heat flux
ds = xr.open_dataset('hol_green/surfHeat/sensible_green.nc')
green_sen = ds.surfHeatDiag.values
green_sen[green_flx==0] = np.nan

# read in BLUE sensible heat flux
ds = xr.open_dataset('hol_blue/surfHeat/sensible_blue.nc')
blue_sen = ds.surfHeatDiag.values
blue_sen[blue_flx==0] = np.nan

# read in GREEN longwave heat flux
ds = xr.open_dataset('hol_green/surfHeat/longwave_green.nc')
green_lw = -ds.surfHeatDiag.values
green_lw[green_flx==0] = np.nan

# read in BLUE longwave heat flux
ds = xr.open_dataset('hol_blue/surfHeat/longwave_blue.nc')
blue_lw = -ds.surfHeatDiag.values
blue_lw[blue_flx==0] = np.nan

# read in GREEN shortwave heat flux
ds = xr.open_dataset('hol_green/surfHeat/shortwave_green.nc')
green_sw = -ds.surfHeatDiag.values
green_sw[green_flx==0] = np.nan

# read in BLUE shortwave heat flux
ds = xr.open_dataset('hol_blue/surfHeat/shortwave_blue.nc')
blue_sw = -ds.surfHeatDiag.values
blue_sw[blue_flx==0] = np.nan

### ------ ###

### --- apply ice cover threshold --- ###

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

### ------ ###

### --- Integrate over open water area --- ###

import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

rac = mds.rdmds('data/RAC')
rac=np.tile(rac,(84,1,1))

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
miz_blue = np.multiply(miz_blue,rac)
miz_green = np.multiply(miz_green,rac)

### ------ ###

### --- Apply mask --- ###

bathy=mds.rdmds('data/Depth');
off_shelf=np.copy(bathy)
off_shelf[220:,:]=np.nan
off_shelf[:50,:]=0
off_shelf[:95,100:]=0
off_shelf[:150,500:]=0
off_shelf[off_shelf>1000]=np.nan
bathy[np.isnan(off_shelf)]=np.nan

bathy = np.tile(bathy,(84,1,1))

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
miz_blue[np.isnan(bathy)] = np.nan
miz_green[np.isnan(bathy)] = np.nan

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
sw_blue = np.nansum(np.nansum(sw_blue,2),1)
sw_green = np.nansum(np.nansum(sw_green,2),1)

lat_blue = np.nansum(np.nansum(lat_blue,2),1)
lat_green = np.nansum(np.nansum(lat_green,2),1)
sen_blue = np.nansum(np.nansum(sen_blue,2),1)
sen_green = np.nansum(np.nansum(sen_green,2),1)

poly_blue = np.nansum(np.nansum(poly_blue,2),1)
poly_green = np.nansum(np.nansum(poly_green,2),1)
miz_blue = np.nansum(np.nansum(miz_blue,2),1)
miz_green = np.nansum(np.nansum(miz_green,2),1)

### ------ ###

### --- convert units --- ###

rad = 3.154e-11*(diff_sw+diff_lw)
bulk = 3.154e-11*(diff_lat+diff_sen)
tot = rad+bulk

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

poly_green = poly_green*3.154e-11
poly_blue  = poly_blue*3.154e-11

### ------ ###

### --- LONGWAVE --- ###

fig=plt.figure(figsize=(40,10))
plin=plt.plot(np.linspace(0.5,83.5,84),lw_green-lw_blue,color=(51/255,160/255,44/255),linewidth=10)
fig=plt.fill_between(np.linspace(60,84,37),-320,320,color=(0.9,0.9,0.9))
fig=plt.bar(np.linspace(0.5,83.5,84),(lw_blue)*(poly_green-poly_blue)/poly_blue,color=(166/255,206/255,227/255),width=1.05)
plt.yticks(fontsize=40)
plt.ylim(-320,320)
plt.xlim(0,84)
plt.grid(axis='x')
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.ylabel('Heat anomaly (EJ/yr)',fontsize=60)
plt.savefig('lw.png')

### ------ ###

### --- SENSIBLE --- ###

fig=plt.figure(figsize=(40,10))
plin=plt.plot(np.linspace(0.5,83.5,84),sen_green-sen_blue,color=(51/255,160/255,44/255),linewidth=10)
fig=plt.fill_between(np.linspace(60,84,37),-320,320,color=(0.9,0.9,0.9))
fig=plt.bar(np.linspace(0.5,83.5,84),(sen_blue)*(poly_green-poly_blue)/poly_blue,color=(166/255,206/255,227/255),width=1.05)
plt.yticks(fontsize=40)
plt.ylim(-320,320)
plt.xlim(0,84)
plt.grid(axis='x')
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.ylabel('Heat anomaly (EJ/yr)',fontsize=60)
plt.savefig('sen.png')

### ------ ###

### --- LATENT --- ###

fig=plt.figure(figsize=(40,10))
plin=plt.plot(np.linspace(0.5,83.5,84),lat_green-lat_blue,color=(51/255,160/255,44/255),linewidth=10)
fig=plt.fill_between(np.linspace(60,84,37),-320,320,color=(0.9,0.9,0.9))
fig=plt.bar(np.linspace(0.5,83.5,84),(lat_blue)*(poly_green-poly_blue)/poly_blue,color=(166/255,206/255,227/255),width=1.05)
plt.yticks(fontsize=40)
plt.ylim(-320,320)
plt.xlim(0,84)
plt.grid(axis='x')
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.ylabel('Heat anomaly (EJ/yr)',fontsize=60)
plt.savefig('lat.png')

### ------ ###

### --- SHORTWAVE --- ###

fig=plt.figure(figsize=(40,10))
plin=plt.plot(np.linspace(0.5,83.5,84),sw_green-sw_blue,color=(51/255,160/255,44/255),linewidth=10)
fig=plt.fill_between(np.linspace(60,84,37),-320,320,color=(0.9,0.9,0.9))
fig=plt.bar(np.linspace(0.5,83.5,84),(sw_blue)*(poly_green-poly_blue)/poly_blue,color=(166/255,206/255,227/255),width=1.05)
plt.yticks(fontsize=40)
plt.ylim(-320,320)
plt.xlim(0,84)
plt.grid(axis='x')
plt.xticks(np.linspace(0,84,29),['                         08','','','',
                                 '              |          09','','','',
                                 '              |          10','','','',
                                 '              |          11','','','',
                                 '              |          12','','','',
                                 '              |          13','','','',
                                 '              |          14','','','',''],fontsize=40)
plt.ylabel('Heat anomaly (EJ/yr)',fontsize=60)
plt.savefig('sw.png')

### ------ ###

### --- Print out anomalies --- ###

ct=60
print('longwave:{}'.format(100*np.nanmean(lw_green[:ct]-lw_blue[:ct])/np.nanmean(lw_blue[:ct])))
print('shortwave:{}'.format(100*np.nanmean(sw_green[:ct]-sw_blue[:ct])/np.nanmean(sw_blue[:ct])))
print('sensible:{}'.format(100*np.nanmean(sen_green[:ct]-sen_blue[:ct])/np.nanmean(sen_blue[:ct])))
print('latent:{}'.format(100*np.nanmean(lat_green[:ct]-lat_blue[:ct])/np.nanmean(lat_blue[:ct])))
print('sum:{}'.format(np.nanmean(lw_green[:ct]+sen_green[:ct]+lat_green[:ct]-lw_blue[:ct]-sen_blue[:ct]-lat_blue[:ct])))
print('open water:{}'.format(np.nanmean(open_green[:ct]-open_blue[:ct])))
print('ice covered:{}'.format(np.nanmean(cov_green[:ct]-cov_blue[:ct])))
print('total:{}'.format(np.nanmean(open_green[:ct]+cov_green[:ct]-cov_blue[:ct]-open_blue[:ct])))

### ------ ###

### --- pie charts for supplement --- ###
net=-np.nanmean(sw_green[:ct]+lw_green[:ct]+sen_green[:ct]+lat_green[:ct]-sw_blue[:ct]-lw_blue[:ct]-sen_blue[:ct]-lat_blue[:ct])
srf=-np.nanmean(open_green[:ct]-open_blue[:ct])
res=srf-net
print('pie chart')
fig, ax = plt.subplots()
sizes = [100*net/srf, 100*res/srf]
labels = "{:.1f}%".format(100*net/srf), "{:.1f}%".format(100*res/srf)
ax.pie(sizes,labels=labels,radius=np.sqrt(srf)/4,colors=((90/255,180/255,172/255),(216/255,179/255,101/255)),textprops={'fontsize': 25})
plt.title('SIC < {:.0f}% | {:.1f} EJ/yr'.format(100*thresh,srf),fontsize=30)
plt.savefig('{}_threshold.png'.format(thresh))
#print('sw+lw+lat+sen:{}'.format(np.nanmean(sw_green[:ct]+lw_green[:ct]+sen_green[:ct]+lat_green[:ct]-sw_blue[:ct]-lw_blue[:ct]-sen_blue[:ct]-lat_blue[:ct])))
#print('open water:{}'.format(np.nanmean(open_green[:ct]-open_blue[:ct])))
