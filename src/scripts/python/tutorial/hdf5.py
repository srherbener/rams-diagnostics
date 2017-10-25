#!/usr/bin/env python3

import h5py
import matplotlib.pyplot as plt
import numpy as np

# open file in read mode
Fname = '/Users/steve/Downloads/GMI_sample.HDF5'
f = h5py.File(Fname, mode='r')

# read in Tb data, index 6 is 37GHz H-Pol channel
lats = f['S1/Latitude'][:,:]
lons = f['S1/Longitude'][:,:]
t37h = f['S1/Tc'][:,:,6]

f.close()

# Mark out-of-range values
t37h[(t37h < 100.0)] = np.nan

# Mark pixels falling outside map range
t37h[(lats > 60.0) + (lats < 30.0) + (lons < 100.0) + (lons > 150.0)] = np.nan

# Plot
Fig = plt.figure()
plt.scatter(lons, lats, c=t37h, s=2, cmap=plt.cm.jet, edgecolor='', linewidth=0.)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.xlim((100.0, 150.0))
plt.ylim((30.0, 60.0))
plt.colorbar(label='K')
plt.title('37 GHz H-Pol Brightness Tempurature')

# show plot interactively
plt.show()

# save plot in file
Fig.savefig("/tmp/Tb_37GHz.png")
