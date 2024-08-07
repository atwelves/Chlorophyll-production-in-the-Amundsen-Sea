# Chlorophyll-production-in-the-Amundsen-Sea

**Chlorophyll production in the Amundsen Sea boosts heat flux to atmopshere and weakens heat flux to ice shelves**

This repository contains the scripts needed to reproduce figures in a submission to Journal of Geophysical Research: Oceans.  The MITgcm model inputs are stored in netCDF format in the Zenodo repository **10.5281/zenodo.10830064**.  These plotting scripts mostly assume a structure where results from the GREEN experiment are stored in the directory hol_green and results from the BLUE experiment are stored in the directory hol_blue.

plot_figures:
  
  plot_bio.py concerns chlorophyll and euphotic depth and can be used to plot **Figure 1a, Figures 2a & 2b, and Supplementary Figures 2a and 2b**.  
  
  plot_sice.py concerns sea ice concentration and thickness and can be used to plot **Figures 4a, 4b, 4c & 4d**.
  
  plot_hov.py concerns ocean temperature and salinity and can be used to plot **Figure 1b, Figures 3a and 3b, Figure 8a, Figure 9a, and Figures 10a & 10b**.
  
  plot_pig.py concerns ice shelf melt rates and can be used to plot **Figure 8b, Figures 9b & 9c, and Figure 10c**.
  
  plot_poly.py concerns heat fluxes in units of W/m2 and can be used to plot **Figures 5a, 5b, 5c & 5d**.
  
  plot_poly_integrated.py conerns heat fluxes in units of EJ/yr and can be used to plot **Figures 6a, b, c, d, e & f**.

  plot_flux.py concerns heat fluxes and can be used to plot **Figure 8**

plot_supplement:
  compare_optics.py concerns chlorophyll and euphotic depth from *Park et al. (2017)* and can be used to plot **Supplementary Figure 1**.
  
  plot_supp5.py concerns ice concentration and heat fluxes and can be used to plot **Supplementary Figure 5**.
  
  plot_upper_temp.py concerns temperatures and can be used to plot **Supplementary Figure 3 and Supplementary Figure 4**.
