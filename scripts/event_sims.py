import numpy as np
import random
import astropy.units as u
from astropy.cosmology import z_at_value 
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM, funcs
import scipy.interpolate as interpolate
from astropy.io import fits
from astropy.table import Table
import os
import argparse
import mockmaps
import glob, sys
import healpy as hp

parser = argparse.ArgumentParser()
parser.add_argument('--outname', type=str, default='mockmaps', help='Name of output file, no file extension')
parser.add_argument('--ligorun', choices=['o3', 'aligo', 'none'], default='none', help='What distribution to pick runs from')
parser.add_argument('--catalog', default='truth_hpix_Chinchilla-0Y3_v1.6_truth', help='What catalog to match RA/DEC/z to')
parser.add_argument('--rootdir', type=str, help='Path to top level work dir')
parser.add_argument('--nevents', type=int, default=3, help='Number of events to generate (default = 7)')
parser.add_argument('--skyarea', type=float, default=60., help='Approximate area enclosing 90% localization probability (in deg2, default is 60.)')
parser.add_argument('--disterr', type=float, default=0.15, help='Fractional error in distance.')
parser.add_argument('--H0', type=float, default=70., help='Hubble parameter (in km/s/Mpc, default is 70.)')
parser.add_argument('--Omegam', type=float, default=0.286, help='Omega_m value (default is 0.286)')
parser.add_argument('--cluster', type=bool, default=False, help='If True, put events in overdense regions only')
parser.add_argument('--central', type=bool, default=False, help='If True, put events on top of BCGs only')
parser.add_argument('--distmin', type=float, default=100., help='Max event distance, in Mpc. Will be overwritten if using ligorun files.')
parser.add_argument('--distmax', type=float, default=700., help='Min event distance, in Mpc. Will be overwritten if using ligorun files.')
parser.add_argument('--maplike', help='Filename of sky map that you want to the outputs to look like. If ommitted, will use NEST=True and NSIDE=1024')
parser.add_argument('--prob', type=float, default=0.90, help='Localization probability enclosed in skyarea')
parser.add_argument('--nomap', type=bool, default=False, help='Do not actually make any maps. Just create info files.')
parser.add_argument('--perturbd', type=bool, default=False, help='Perturb event distances.')
parser.add_argument('--perturbc', type=bool, default=False, help='Perturb event coords.')
parser.add_argument('--absmaglim', type=float, default=-17.23, help='Absmaglim')
parser.add_argument('--maglim', type=float, default=24., help='maglim')

args=parser.parse_args()

outname=args.outname

if args.rootdir == None:
    DIR_SOURCE = os.getcwd()
    DIR_MAIN = DIR_SOURCE.rsplit('/', 1)[0]
else:
    DIR_MAIN = args.rootdir

if args.ligorun == 'o3':
    ligo_z_distrib = np.loadtxt(DIR_MAIN+'/catalogs/stat_redshift_o3_120_60_hlv_bbh_10_10.txt', unpack=True)

if args.ligorun == 'aligo':
    ligo_z_distrib = np.loadtxt(DIR_MAIN+'/catalogs/stat_redshift_aLIGO_hlv_bbh_10_10.txt', unpack=True)

if args.ligorun == 'none':
    print "Input redshift distribution not provided." 
else: 
    np.random.shuffle(ligo_z_distrib)

print "Producing list of galaxy catalogs..."

galaxy_cats = glob.glob(DIR_MAIN+'/catalogs/'+args.catalog+'*.fits')

if len(galaxy_cats) == 0: 
    print "ERROR: No galaxy catalogs found!"
    sys.exit(1)


mock=False
nside=8
good_pixels = range(hp.nside2npix(nside))
nest=False
lon=0
lat=0
psi=0

# use only pixels in the DES footprint that are far away from the edges
if args.catalog == 'Chinchilla-0Y3_v1.6_truth.':
    mock=True
    good_pixels=[39,27,107,120,38,25,383,18,17,382,378,22,28,52,110,13,31,29,20,23,21,24,19,16,15,26,30,106,379,37,36,111,53,48,49,50,51,108,109]
    nside=8
    nest=True
    lon=-60
    lat=-180
    psi=0

# use only pixels far from the roi on the DES data
if args.catalog == 'GW_cat_hpx_':
    mock=False
    good_pixels =[10059, 10060, 10187, 10195, 10314, 10315, 10322, 10438, 10439, 10558, 10565, 10673, 10674, 10680, 10785, 10791, 10893, 10898, 10899, 10996, 10997, 11001, 11002, 11096, 11100, 11101, 11102, 11192, 11195, 11196, 11197]
    nside=32


good_cats=[]
bad_pixels=[]

for i in range(len(galaxy_cats)):
    pixel = 0 
    if mock:
        pixel = int(galaxy_cats[i].split('.')[-2].split('_')[0])
    else:
        pixel = int(galaxy_cats[i].split(args.catalog)[-1].split('.')[0])

    if pixel in good_pixels:
        good_cats.append(galaxy_cats[i])
    else:
        bad_pixels.append(pixel)

galaxy_cats=good_cats


np.random.shuffle(galaxy_cats)

hpix = np.zeros(len(galaxy_cats),dtype=int)

for i in range(len(galaxy_cats)):
    if mock: 
        hpix[i] = galaxy_cats[i].split('.')[-2].split('_')[0]
    else:
        hpix[i] = galaxy_cats[i].split(args.catalog)[-1].split('.')[0]
        

print "Producing list of events..."

H_0=args.H0
Omega_m=args.Omegam
cosmo = FlatLambdaCDM(H0=H_0, Om0=Omega_m)

outfile = open(DIR_MAIN+'/out/'+outname+'_events.txt','w')
outfile.write('# cosmo = FlatLambdaCDM\n')
outfile.write('# H0  = '+str(H_0)+'\n')
outfile.write('# Omegam  = '+str(Omega_m)+'\n')
outfile.write('# skyarea = '+str(args.skyarea)+'\n')
outfile.write('# prob = '+str(args.prob)+'\n')
outfile.write('# derrfrac = '+str(args.disterr)+'\n')
outfile.write('# cluster = '+str(args.cluster)+'\n')
outfile.write('# central = '+str(args.central)+'\n')
outfile.write('# glxycat = '+args.catalog+'\n')
outfile.write('# obsruntype = '+args.ligorun+'\n')
outfile.write('# mapslabel = '+outname+'\n')
outfile.write('# mapsdir = '+DIR_MAIN+'/skymaps'+'\n')
outfile.write('# EVENT_ID HOST_ID RA DEC Z DIST DIST_ERR ROI_AREA HPIX\n')

nevents=args.nevents
ra=np.zeros(nevents)
dec=np.zeros(nevents)
z=np.zeros(nevents)
host_id=np.zeros(nevents)
hpix=np.resize(hpix,nevents)
galaxy_cats=np.resize(galaxy_cats,nevents)
pix=np.zeros(nevents)
sky_area = np.full(nevents,args.skyarea) 
angular_size = np.full(nevents,mockmaps.apperture(args.skyarea))

zmin=z_at_value(cosmo.luminosity_distance, max(1.,args.distmin) * u.Mpc)
zmax=z_at_value(cosmo.luminosity_distance, args.distmax * u.Mpc)
comov_dist_min = cosmo.comoving_distance(zmin).value
comov_dist_max = cosmo.comoving_distance(zmax).value
comov_dist = np.linspace(comov_dist_min,comov_dist_max, num=5000)
pdf=comov_dist**2 ###* mockmaps.selection_function(comov_dist,cosmo,mu=638.,sigma=50.)
pdf=pdf/pdf.sum()
dist=np.random.choice(comov_dist, size=nevents, replace=False,  p=pdf)
dist_err = np.empty_like(dist)

i=0

while i in range(nevents):

    if args.ligorun == 'none':
        #event_z = z_at_value(cosmo.luminosity_distance, dist[i] * u.Mpc)
        event_z = z_at_value(cosmo.comoving_distance, dist[i] * u.Mpc)
    else:
        event_z = ligo_z_distrib[i]
        dist.fill(cosmo.luminosity_distance(event_z).value)
        dist_err.fill(dist[0]*args.disterr)

    print "Reading galaxy catalog hpix", hpix[i], "for event", i

    ucats = np.unique(galaxy_cats)

    gal_id, gal_ra, gal_dec, gal_z = mockmaps.read_galaxy_cat(galaxy_cats[i],args.central,args.cluster,mock,args.maglim,args.absmaglim)

    gal_pix = np.full_like(gal_id,hpix[i])

    nearest_gals = np.array([mockmaps.find_nearest(gal_z, event_z)])    

    for gal in nearest_gals:
        host_id[i] = gal_id[gal]
        ra[i] = gal_ra[gal]
        dec[i] = gal_dec[gal]
        z[i] = gal_z[gal]
        pix[i] = gal_pix[gal]
        dist[i]=cosmo.luminosity_distance(z[i]).value
        dist_err[i] = args.disterr * dist[i] 

        i=i+1


if args.nomap:
    roi90=-99.
    for i in range(nevents):
        outfile.write('%i %i %.10f %.10f %.10f %f %f %f %i\n' % (i, host_id[i], ra[i], dec[i], z[i], dist[i], dist_err[i], roi90, pix[i]))
else:
    print "Producing skymaps..."

    sample_map_file_name=args.maplike  

    maps=mockmaps.gaussian2d(ra=ra,dec=dec,sky_area=sky_area,dist=dist,dist_err=dist_err,
                             sample_map_file_name=sample_map_file_name,mock_map_name=outname,
                             rootdir=DIR_MAIN+'/skymaps/',plotdir=DIR_MAIN+'/plots/',
                             perturbd=args.perturbd,perturbc=args.perturbc)

    print "Computing ROI areas..."

    for i in range(nevents):
        m=maps[i]
        roi90=mockmaps.roi_area(m[0],args.prob)
        peak=np.argmax(m[0])
        rap,decp=hp.pix2ang(1024,peak,nest=True,lonlat=True)
        #print 'event', i+1 , 'of' , nevents
        #print 'true:', ra[i]    , dec[i]     , dist[i]
        #print 'peak:', rap      , decp       , m[1,peak]
        #print 'diff:', rap-ra[i], decp-dec[i], m[1,peak]-dist[i]
        #print '====='
        outfile.write('%i %i %.10f %.10f %.10f %f %f %f %i\n' % (i, host_id[i], ra[i], dec[i], z[i], dist[i], dist_err[i], roi90, pix[i]))

    print "Saving footprint plot..."

    gmap=np.zeros(hp.nside2npix(nside),dtype=int)
    gmap[bad_pixels]=1
    gmap[good_pixels]=2
    gmap[pix.astype(np.int)]=3
    
    mockmaps.plot(mapdata=gmap,outname=outname+'_footprint',nest=nest,
                  dirpath=DIR_MAIN+'/plots/',extension='.png',cmapname='Reds',
                  cbar=False,rot=[lon,lat,psi],grid=False)

outfile.close()

