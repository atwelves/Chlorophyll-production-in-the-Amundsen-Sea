# This script plots upper ocean temperature profiles for GREEN and BLUE experiments, also mixed layers reported by Park et al. (2017)

# For Supplementary Figures 3 and 4

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import MITgcmutils as mitgcm
import MITgcmutils.mds as mds

# Amundsen Sea station data from Park et al. (2017)
mld_obs = [27,33,25,26,22,41,40,18,31,37,16,30,28,26,57,41,25,67,25,36,41,20,44]
sst_obs = [-1.81,-1.78,-1.79,-1.8,-1.68,-1.57,-1.12,-1.6,-1.38,-1.24,-1.67,-0.83,-1.78,-1.74,-1.29,-1.82,-1.17,-1.06,-1.5,-1.43,-1.15,-1.13,-1.32]
stat_id = [1,2,3,6,7,8,10,12,16,17,19,24,61,63,71,85,31,34,86,87,88,38,39]
lon_obs = 360-np.array([116.78,117.67,118.45,117.72,116.5,116.5,115,113.5,113,114,112.51,115.72,116.37,117.58,113.07,114.33,101.76,101.53,106.54,105,102.1,119,133.99])
lat_obs = [-71.66,-71.79,-71.95,-72.39,-72.85,-73.5,-73.25,-72.99,-73.5,-73.5,-74.2,-74.08,-72.46,-72.93,-73.82,-71.42,-75.09,-74.65,-73.81,-74.37,-74.86,-70.45,-71.58]

RC = np.squeeze(mds.rdmds('data/RC'))
rc = -RC

ds = xr.open_dataset('/media/twelves/My Passport/edinburgh_work/blue_temp.nc')
blue_temp = ds.physDiag.values
lon = ds.x.values
lat = ds.y.values
print(ds)
xr.Dataset.close(ds)

ds = xr.open_dataset('/media/twelves/My Passport/edinburgh_work/green_temp.nc')
green_temp = ds.physDiag.values
xr.Dataset.close(ds)

green_temp[green_temp==0] = np.nan
blue_temp[blue_temp==0] = np.nan

time_indx = 49

i = 0
match_lon = np.argmin(np.abs(lon_obs[i] - lon))
match_lat = np.argmin(np.abs(lat_obs[i] - lat))
print(green_temp[time_indx,:,match_lat,match_lon])
print(rc)

for i in range(0,len(lon_obs)):
    match_lon = np.argmin(np.abs(lon_obs[i] - lon))
    match_lat = np.argmin(np.abs(lat_obs[i] - lat))
    plt.figure(figsize=(40,9))
    plt.plot(np.linspace(np.nanmin([np.nanmin(green_temp),np.nanmin(blue_temp)]),np.nanmax([np.nanmax(green_temp),np.nanmax(blue_temp)]),10),np.linspace(mld_obs[i],mld_obs[i],10),'k',linestyle=':',linewidth=30,alpha=0.5)
    plt.plot(green_temp[time_indx,:,match_lat,match_lon],rc,color=(178/255,223/255,138/255),linewidth=30)
    plt.plot(blue_temp[time_indx,:,match_lat,match_lon],rc,color=(31/255,120/255,180/255),linewidth=30,linestyle='--')
    #plt.scatter(sst_obs[i],0)
    plt.yticks([0,50,100,150,200],['0','','100','','200'],fontsize=100)
    plt.xticks([-2,-1,0,1,2],fontsize=0)
    plt.grid()
    plt.xlim(-2.4,2.4)
    plt.ylim(0,200)
    plt.gca().invert_yaxis()
    plt.savefig('strat_park_station_{}'.format(stat_id[i]))

blue_diff = np.diff(blue_temp,1)
green_diff = np.diff(green_temp,1)
diff_rc = np.diff(rc)

print('shape of green_temp:{}'.format(np.shape(green_diff)))
print('shape of diff_rc:{}'.format(np.shape(diff_rc)))

for i in range(0,len(lon_obs)):
    match_lon = np.argmin(np.abs(lon_obs[i] - lon))
    match_lat = np.argmin(np.abs(lat_obs[i] - lat))
    plt.figure(figsize=(40,9))
    plt.plot(np.linspace(np.nanmin([np.nanmin(green_temp),np.nanmin(blue_temp)]),np.nanmax([np.nanmax(green_temp),np.nanmax(blue_temp)]),10),np.linspace(mld_obs[i],mld_obs[i],10),'k',linewidth=30,linestyle=':',alpha=0.5)
    plt.scatter(np.divide(np.diff(blue_temp[time_indx,:,match_lat,match_lon]),diff_rc),rc[:-1],s=3000,c=(31/255,120/255,180/255),alpha=0.7)
    plt.scatter(np.divide(np.diff(green_temp[time_indx,:,match_lat,match_lon]),diff_rc),rc[:-1],s=3000,c=(178/255,223/255,138/255),alpha=0.7)
    #plt.scatter(sst_obs[i],0)
    plt.yticks([0,50,100,150,200],['0','','100','','200'],fontsize=100)
    plt.xticks([-0.1,-0.05,0,0.05,0.1],fontsize=0)
    plt.xlim(-0.105,0.105)
    plt.ylim(0,200)
    plt.gca().invert_yaxis()
    plt.grid()
    plt.savefig('strat_delta_station_{}'.format(stat_id[i]))
