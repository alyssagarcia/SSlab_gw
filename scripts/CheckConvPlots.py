#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

h0=np.loadtxt('alyssaCHECK_sim_test.txt', usecols=(1), unpack=True)
plt.plot(h0)
plt.show()
