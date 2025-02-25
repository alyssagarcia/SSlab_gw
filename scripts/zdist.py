import matplotlib.pyplot as plt
import numpy as np

dl, z=np.loadtxt('hsin-yu_data/testHYEvent2_HLVD.txt', usecols=[3,6], unpack=True)
dl2, z2=np.loadtxt('hsin-yu_data/testHYEvent2_O3.txt', usecols=[3,6], unpack=True)

"""
plt.plot(dl,z, 'ro')
plt.xlabel('dist')
plt.ylabel('z')
plt.title('37 sims')
#plt.show()
plt.savefig('HLVD_dist_vs_z.png')
"""

plt.plot(dl2,z2, 'ro')
plt.xlabel('dist')
plt.ylabel('z')
plt.title('8 sims')
#plt.show()
plt.savefig('O3_dist_vs_z.png')
