import numpy as np
from chainconsumer import ChainConsumer
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt


h0, post=np.loadtxt('gw170817-holz_test.txt', usecols=[0,1], unpack=True)
c = ChainConsumer()
c.add_chain([h0,post], parameters=["h0", "post"])
#fig = c.plotter.plot(display=True)
fig = c.plotter.plot_distributions()
plt.show()
