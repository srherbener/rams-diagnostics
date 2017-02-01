#!/usr/bin/env python3

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# 1000 mb heights
data = np.load('/Users/steve/Downloads/jan1000mb.npz')
lon = data['lon']
lat = data['lat']
z = data['z']

# create mesh
lon2d, lat2d = np.meshgrid(lon, lat)

# create map, then place height data on map
Fig = plt.figure()

m = Basemap(projection = 'ortho', lat_0 = 20, lon_0 = 0)

m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color = 'black', linewidth = 0.5)
m.fillcontinents(color = '0.85',)
m.drawparallels(np.arange(-90,91, 30), color = '0.25', linewidth = 0.5)
m.drawmeridians(np.arange(-180,180, 30), color = '0.25', linewidth = 0.5)

# convert lon, lat to map coordinates
x, y = m(lon2d, lat2d)

# contours
cs = m.contour(x, y, z, levels = range(-180, 360, 30), colors = 'blue')
plt.clabel(cs, fmt = '%.0f', inline = True)

# show plot interactively
#plt.show()

# save plot in file
Fig.savefig("/tmp/heights_1000mb.png")
