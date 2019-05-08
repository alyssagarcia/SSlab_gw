#!/usr/bin/env python

import numpy as np
from chainconsumer import ChainConsumer
from numpy.random import normal

"""
np.random.seed(0)
# Here we have some nice data, and then some bad data,
# where the last part of the chain has walked off, and the first part
# of the chain isn't agreeing with anything else!
data_good = normal(size=100000)
data_bad = data_good.copy()
data_bad += np.linspace(-0.5, 0.5, 100000)
data_bad[98000:] += 2
print('datagood',data_good)
print('databad',data_bad)

data=np.loadtxt('alyssaCHECK_sim_test.txt')
print('mydata', data[1])
# Lets load it into ChainConsumer, and pretend 10 walks went into making the chain
c = ChainConsumer()
c.add_chain(data_good, walkers=10, name="good")
c.add_chain(data_bad, walkers=10, name="bad")
print('c',c)


# Now, lets check our convergence using the Gelman-Rubin statistic
gelman_rubin_converged = c.diagnostic.gelman_rubin()
# And also using the Geweke metric
geweke_converged = c.diagnostic.geweke()

# Lets just output the results too
print(gelman_rubin_converged, geweke_converged)

"""
data=np.loadtxt('alyssaCHECK_sim_test.txt', usecols=6)

c = ChainConsumer()
c.add_chain(data, walkers=57)

print(data, len(data))
print(c)

gelman_rubin_converged = c.diagnostic.gelman_rubin()
# And also using the Geweke metric
geweke_converged = c.diagnostic.geweke()

# Lets just output the results too
print(gelman_rubin_converged, geweke_converged)

