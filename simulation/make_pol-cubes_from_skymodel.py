#!/usr/bin/env python
"""
make_pol-cubes_from_skymodel
based on: make_test_files

A simple script for creating a mock dataset with which to test pyrmsynth for
simple sky models.
Main contributors: Justin Bray, Henrik Junklewitz & Andreas Horneffer

"""

import datetime
import sys
import optparse

import pyfits
import numpy as np

# The model file is a text file containing model point parameters, with each
# parameter separated by spaces. One model definition per line!
# Each line contains:
#     Phi location (rad/m/m), Dec (pixel), RA (pixel), Stokes I, Q, U, V, specindex
# where the stokes Q and U parameters define the Faraday spectrum.
# The given Stokes I and V are simply assigned to all frequencies (alpha=0).

parser = optparse.OptionParser()
parser.add_option("--reffile", dest="REFFile",type="string",
                  help="The reference FITS file from which to grab model header")
parser.add_option("--model", dest="MODEL",type="string",
                  help="The model Faraday spectrum")
parser.add_option("--out-base", dest="OUTBASE",type="string",
                  help="The base name for output FITS files")
parser.add_option("--nsb", dest="NSB",type="int",default=10,
                  help="The number of subbands, one FITS file per SB (default: 10)")
parser.add_option("--nchan", dest="NCHAN",type="int",default=8,
                  help="The number of channels per subband (default: 8)")
parser.add_option("--freqstart", dest="FREQSTART",type="float",default=120e6,
                  help="The frequency of the first channel, in Hz (default: 120e6 Hz")
parser.add_option("--dfreq", dest="DFREQ",type="float",default=24.4e3,
                  help="The channel width, in Hz (default: 24.4e3 Hz)")
parser.add_option("--npix", dest="NPIX",type="int",default=100,
                  help="The number of pixles in RA and DEC (default: 100)")
parser.add_option("--noise", dest="NOISE",type="float",default=1.,
                  help="RMS of the noise to add to the data in Jy/beam (default: 1. Jy/beam)")
parser.add_option("--instq", dest="INSTQ",type="float",default=0.,
                  help="Percentage of Stokes-I to add to Stokes-Q for instrumental "
                  "polarization, can be negative. (default: 0.%)")
parser.add_option("--instu", dest="INSTU",type="float",default=0.,
                  help="Percentage of Stokes-I to add to Stokes-U for instrumental "
                  "polarization, can be negative. (default: 0.%)")
parser.add_option("--specmap", dest="SPECMAP",default=False,action="store_true",
                  help="Generate spectral index map (currently with strange header). "
                  "(default: don't write a spectral index map)")
parser.add_option("--specnoise", dest="SPECNOISE",type="float",default=0.,
                  help="Add (Gaussian-)noise to the spectral-index map (default: RMS=0.)")


(options, args) = parser.parse_args()
if not options.REFFile or not options.MODEL or not options.OUTBASE:
    parser.error("Need vaules for reffile, model and out-base to work!")

reffn = options.REFFile
mdl = np.loadtxt(options.MODEL)
outfn_base = options.OUTBASE
nsb = options.NSB
nchan = options.NCHAN 
freq0 = options.FREQSTART
dfreq =  options.DFREQ
nra = options.NPIX
ndec = options.NPIX
noise = options.NOISE
instq = options.INSTQ/100.
instu = options.INSTU/100.

cra = nra/2
cdec = ndec/2

c = 2.99792458e8

ref_head = pyfits.getheader(reffn)
ref_head.update('BUNIT','JY/BEAM')


for i in range(nsb):
    print "Working on subband",i,"out of",nsb
    sbcube = np.zeros((4, nchan, ndec, nra))
    for j in range(nchan):
        for k in range(len(mdl)):
            freq = freq0 + (i * nchan + j) * dfreq
            l2 = (c / freq) ** 2            
            #Stokes I
            stokesI = mdl[k][3] * (freq / freq0)**-mdl[k][7]
            sbcube[0, j, cdec + mdl[k][1], cra - mdl[k][2]] = stokesI
                
            #Stokes V
            sbcube[3, j, cdec + mdl[k][1], cra - mdl[k][2]] = \
                mdl[k][6] * (freq / freq0)**-mdl[k][7]

            # Stokes Q & U from Faraday spectrum
            val = complex(mdl[k][4], mdl[k][5])
            val = val * np.exp(2. * complex(0, 1) * mdl[k][0] * l2) 

            sbcube[1, j, cdec + mdl[k][1], cra - mdl[k][2]] = val.real \
                * (freq / freq0)**-mdl[k][7] + (instq*stokesI)
            sbcube[2, j, cdec + mdl[k][1], cra - mdl[k][2]] = val.imag \
                * (freq / freq0)**-mdl[k][7] + (instu*stokesI)

    if (noise > 0.):
        sbcube += np.random.normal(0.,noise,np.shape(sbcube))

    hdu = pyfits.PrimaryHDU(sbcube)
    hdu.header = ref_head.copy()

    hdu.header.update('ORIGIN', 'make_pol-cubes_from_skymodel.py', 'Origin of the data set')
    hdu.header.update('DATE', str(datetime.date.today()),
                      'Date when the file was created')
    hdu.header.update('NAXIS3', nchan, 'Length of the Frequency axis')
    hdu.header.update('CRPIX3', 1, 'Reference pixel')
    hdu.header.update('CRVAL3', freq0 + i * nchan * dfreq)
    hdu.header.update('CDELT3', dfreq)

    hdu.header.update('NAXIS1', nra)
    hdu.header.update('CRPIX1', nra / 2.)
    hdu.header.update('NAXIS2', ndec)
    hdu.header.update('CRPIX2', ndec / 2.)

    sbstr = ""
    if i / 100 == 0:
        sbstr += "0"
    if i / 10 == 0:
        sbstr += "0"
    sbstr += str(i)

    hdu_list = pyfits.HDUList([hdu])
    hdu_list.writeto(outfn_base + sbstr + '.fits', clobber=True)

    del hdu
    del hdu_list

if options.SPECMAP:
    print "Working on spectral index map"
    specmap = np.zeros((1, 1, ndec, nra))
    for k in range(len(mdl)):
        specmap[0, 0, cdec + mdl[k][1], cra - mdl[k][2]] = mdl[k][7]
    if options.SPECNOISE>0.:
        specmap += np.random.normal(0.,options.SPECNOISE,np.shape(specmap))

    hdu = pyfits.PrimaryHDU(specmap)
    hdu.header = ref_head.copy()
    
    hdu.header.update('ORIGIN', 'make_pol-cubes_from_skymodel.py', 'Origin of the data set')
    hdu.header.update('DATE', str(datetime.date.today()),
                      'Date when the file was created')
    hdu.header.update('NAXIS3', nchan, 'Length of the Frequency axis')
    hdu.header.update('CRPIX3', 1, 'Reference pixel')
    hdu.header.update('CRVAL3', freq0 + i * nchan * dfreq)
    hdu.header.update('CDELT3', dfreq)
    
    hdu.header.update('NAXIS1', nra)
    hdu.header.update('CRPIX1', nra / 2.)
    hdu.header.update('NAXIS2', ndec)
    hdu.header.update('CRPIX2', ndec / 2.)
    
    hdu_list = pyfits.HDUList([hdu])
    hdu_list.writeto(outfn_base + '.fits.specmap', clobber=True)


    
