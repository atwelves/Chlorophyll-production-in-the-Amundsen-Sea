# Plot pie charts to illustrate how choice of sea ice threshold impacts closure of the surface heat budget

import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# define range of thresholds
thresholds = np.linspace(0,25,6)

for th in thresholds:
  # find longwave, shortwave, sensible and latent heat fluxes
  net = lw+sw+sen+lat
  # total surface flux
  # find residual
  # Plot pie chart
  plt.figure(figsize=(20,20))
  pict = pie(net,surf,resid)
  plt.savefig('residual_pie_{}'.format(th))
