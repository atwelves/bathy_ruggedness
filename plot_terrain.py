# 16.02.14

### Load modules

import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

### --- Coordinates of region of interest --- ###

# model resolutions
resol = np.array([0.25,1])

# Gulf of Finland
jmin = np.array([0,620])
jmax = np.array([400,720])
imin = np.array([0,950])
imax = np.array([947,1230])

bathymetry = np.array([np.zeros((1+jmax[0]-jmin[0],1+imax[0]-imin[0])),np.zeros((1+jmax[1]-jmin[1],1+imax[1]-imin[1]))],dtype=object)

# filenames
filename1 = "../regional_model_comparisons/domain_files/gof_domain_cfg.nc"
filename2= "../regional_model_comparisons/domain_files/nordic_domain_cfg.nc"
files = [filename1,filename2]

### --- Read in bathymetry --- ###

for i in range(0,2):
    domain = xr.open_dataset(files[i])
    field  = domain.bathy_metry.values
    # region of interest
    field  = field[0,jmin[i]:jmax[i],imin[i]:imax[i]]
    # set land to nans
    field[field==0]=np.nan
    # assigna to object
    bathymetry[i] = field

### --- Hypsometry --- ### 

### Do this next week

### --- Terrain ruggedness index --- ###

# minimum neighbourhood size
min_buf = np.array([1,1])
# maximum neighbourhood size
max_buf = np.array([32,8]) 
# maximum neighbourhood sets mask?

tri_avg = np.array([np.zeros((2+max_buf[0]-min_buf[0])),np.zeros((2+max_buf[1]-min_buf[1]))],dtype=object)

for mod_iter in range(0,2):
    bathy   = bathymetry[mod_iter]
    tri     = np.zeros((np.shape(bathy)))
    mask    = np.ones((np.shape(bathy)))
    tri_temp = np.zeros((2+max_buf[mod_iter]-min_buf[mod_iter]))
    # loop over neighbourhood size
    for bf_iter in range(min_buf[mod_iter],max_buf[mod_iter]+1):
        #print(bf_iter)
        # reverse order
        bf = 1 + max_buf[mod_iter] - bf_iter
        #print(bf)
        bf_pts = np.square(1 + 2*bf) 
        #print(bf_pts)
        # loop over centre cells
        for jc in range(0,np.shape(bathy)[0]):
            for ic in range(0,np.shape(bathy)[1]):
                buffer = abs(bathy[jc-bf:jc+bf,ic-bf:ic+bf] - bathy[jc,ic])
                # exclude points already excluded at lower resolution
                # Terrain ruggedness index calculation
                tri[jc,ic]   = mask[jc,ic]*np.sum(buffer)/(bf_pts-1)
        # average over (sub) domain
        tri_temp[bf] = np.nanmean(tri)
        #print(tri_temp[bf])
        mask[np.isnan(tri)] = np.nan
    tri_avg[mod_iter] = tri_temp
# Plot outputs

plt.figure()
for mod_iter in range(0,2):
    endpoint = max_buf[mod_iter]
    xrange = np.linspace(min_buf[mod_iter],max_buf[mod_iter],1+max_buf[mod_iter]-min_buf[mod_iter])*resol[mod_iter]
    print(xrange)
    print(np.shape(xrange))
    tri_temp = tri_avg[mod_iter]
    print(np.shape(tri_temp))
    plt.plot(xrange,tri_temp[min_buf[mod_iter]:1+max_buf[mod_iter]])
plt.show()

plt.figure()
pcol = plt.pcolormesh(mask)
plt.show()

plt.figure()
pcol = plt.pcolormesh(tri)
cbar = plt.colorbar(pcol)
plt.show()
