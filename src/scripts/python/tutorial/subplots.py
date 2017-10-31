#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

# 2 x 2 panels with shared axes
Fig, Ax = plt.subplots(2, 2, sharex = True, sharey = True)

# plots
x = np.linspace(0,2*np.pi,100)
Ax[0,0].plot(x, np.cos(x))
Ax[0,1].plot(x, np.sin(x))
Ax[1,0].plot(x, np.cos(2*x))
Ax[1,1].plot(x, np.sin(2*x))

Ax[0,0].set_xlim(0, 2*np.pi)

# show plot interactively
plt.show()

# save plot in file
Fig.savefig("/tmp/subplots.png")
