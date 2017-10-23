#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0,2*np.pi,0.1*np.pi)
y = np.sin(x)

plt.plot(x,y)
Fig = plt.gcf()
plt.ylim(-1.1, 1.1)
plt.xlim(0,2*np.pi)
plt.xticks([ 0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi ],
     [ r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', 
       r'$\frac{3\pi}{2}$', r'$2\pi$' ], color = 'k',
       size = 'x-large')

# show plot interactively
plt.show()

# save plot in file
Fig.savefig("/tmp/plot_eg1.png")
