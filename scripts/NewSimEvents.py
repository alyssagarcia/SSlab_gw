#!/usr/bin/env python

import random
import math
#import scipy.constants 
import numpy as np
import argparse
import astropy.units as u
from astropy.cosmology import z_at_value , WMAP7
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM

#values from THE DISTANCE TO NGC 4993 paper
omegaM=0.24
h=0.73
omegaL=0.76
H0=73.24 
q=-0.53
table = []


parser = argparse.ArgumentParser()
parser.add_argument('run', choices=['O3', 'AligoDesign', 'AplusDesign'], help="Observing run")
parser.add_argument('--runTime', type=int, default=12,
                    help="How many months the detector will be observing. Default 12months") 

args=parser.parse_args()
#print(args.run)
#print(args.runTime)

time=args.runTime/12 #get time in units of years

if args.run == 'O3':
    #N = (1.54e3)*(0.007)*time*0.7 #N = rate/Gpc^3/yr * V Gpc^3 * Time yrs *epsilon(duty cycle)
    #updated (11/14/18)
    N = (1.54e3)*(0.0073)*time*0.7
    runs = int(round(N))
    #avgDist = 127 #Mpc
    #updated
    maxDist = 135 #Mpc
elif args.run == 'AligoDesign':
    #N = (1.54e3)*(0.03)*time*0.8
    N = (1.54e3)*(0.03)*time*0.75
    runs= int(round(N))
    #avgDist = 208
    maxDist = 190 #Mpc
elif args.run == 'AplusDesign':
    #N = (1.54e3)*(0.03)*time*1
    N = (1.54e3)*(0.14)*time*0.75
    runs= int(round(N))
    maxDist = 325
print("Number of events: "+str(runs))


outfile=open(str('test_1-16-19comove_')+args.run+str('.txt'), 'w')
outfile.write('#emid, gwid, prob, dl, derr, z, zerr, vp, vperr\n')

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
#cosmo = FlatLambdaCDM(H0=0.7, Om0=0.3)

#from marcelle event_sims.py
zmin=z_at_value(cosmo.luminosity_distance, 1 * u.Mpc) #10Mpc for BNS, 100 for BBH
zmax=z_at_value(cosmo.luminosity_distance, maxDist * u.Mpc)
comov_dist_min = cosmo.comoving_distance(zmin).value
comov_dist_max = cosmo.comoving_distance(zmax).value

comov_dist = np.linspace(comov_dist_min,comov_dist_max, num=5000) #array of evenly spaced distances
#comov_dist = np.linspace(1 ,maxDist, num=500000)

pdf=comov_dist**2 ###* mockmaps.selection_function(comov_dist,cosmo,mu=638.,sigma=50.)          
pdf=pdf/pdf.sum()
dist=np.random.choice(comov_dist, size=runs, replace=False,  p=pdf) #randomly choose N distances
dist_err = np.empty_like(dist)

for i in range(0,runs):
    #D=random.gauss(avgDist,100)
    #D = np.linspace(10, avgDist, num=runs) #used for runs labeled uniform
    D = dist[i]
    #z=z_at_value(cosmo.luminosity_distance, D*u.megaparsec, zmin=0, zmax=4) #used when not using comoving distanc
    z = z_at_value(cosmo.comoving_distance, D*u.Mpc)

    #emid, gwid, prob, dl, derr, z, zerr, vp, vperr
    outfile.write('%f %s %f %f %f %.3f %f %f %f\n' % (i, "gw"+str(i), 1, D, 0.25*D, z, 3e-5, 0, 300))


