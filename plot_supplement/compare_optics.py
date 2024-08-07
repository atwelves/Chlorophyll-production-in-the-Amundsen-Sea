# This script compares euphotic depths measured with a PAR sensor and euphotic depths produced by applying light attenuation scheme
# from Manizza et al. (2005) to surface chlorophyll concentrations in the Amundsen Sea

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Amundsen Sea station observations from Park et al. 2017
chl_obs = [3.21,3.03,3.29,3.41,3.84,2.95,4.46,6.20,5.72,0.81,1.74,2.82,2.00,4.89,0.39,2.91,2.59,2.33,1.56,0.35]
zeu_obs = [27.9,25.6,24.7,22.8,20.9,14.8,18.1,18.1,17.6,41.8,28.2,30.7,15.3,14.8,26.3,33.1,26.3,37.7,40.0,15.5]
stat_id = [1,2,6,7,8,10,12,16,17,19,24,61,63,71,85,31,34,87,88,39]

# define some colormaps
rgb1    = [27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,217,217,217,217,117]
rgb2    = [158,158,158,158,158,158,158,158,158,158,158,158,158,158,158,95,95,95,95,112]
rgb3    = [119,119,119,119,119,119,119,119,119,119,119,119,119,119,119,2,2,2,2,179]

print(rgb1[0])

# calculate arrays of attenuation coefficients
k_red   = 0.225 + 0.037*np.power(chl_obs,0.629)
k_bg    = 0.0232 + 0.074*np.power(chl_obs,0.674)

k = k_red + k_bg

for i in range(0,len(chl_obs)):
    print('chl:{}'.format(chl_obs[i]))
    print('zeu:{}'.format(zeu_obs[i]))
    clr = (rgb1[i]/255,rgb2[i]/255,rgb3[i]/255)
    plt.figure(figsize=(40,35))
    plt.fill_between(np.linspace(-3.1,-2,10),0,72,color=(0.9,0.9,0.9))
    # calculate euphotic depth
    plt.plot(np.log10(0.5*np.exp(-k_red[i]*np.linspace(0,72,72))+0.5*np.exp(-k_bg[i]*np.linspace(0,72,72))),np.linspace(0,72,72),linewidth=20,color=clr)
    plt.scatter(-2,zeu_obs[i],s=2000,c=clr)
    plt.xlim(-3.1,0.1)
    plt.ylabel('Depth (m)',fontsize=150)
    plt.yticks([60,40,20,0],fontsize=150)
    #plt.xlabel('% of surface irradiance',fontsize=150)
    plt.xticks([-3,-2,-1,0],['0.1','1','10','100'],fontsize=0)
    plt.grid()
    plt.ylim(0,72)
    plt.gca().invert_yaxis()
    plt.title('Station {}'.format(stat_id[i]),fontsize=150)
    plt.savefig('optics_solver_station_{}.png'.format(stat_id[i]))

# Supplementary figure 1
plt.figure()
for i in range(0,len(chl_obs)):
    plt.plot(np.log10(0.5*np.exp(-k_red[i]*np.linspace(0,40,40))+0.5*np.exp(-k_bg[i]*np.linspace(0,40,40))),-np.linspace(0,40,40))
    plt.scatter(-2,-zeu_obs[i])
    plt.grid()
plt.savefig('optics_solver_combined.png'.format(i+1))
