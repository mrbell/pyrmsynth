#!/usr/bin/env python
"""
make_test_files

A simple script for creating a mock dataset with which to test the LOFAR
rmsynthesis.py script.
"""

import datetime
import sys

import pyfits
import numpy as np

if len(sys.argv) != 10:
    raise Exception("Incorrect number of inputs!")


# The model file is a text file containing model point parameters, with each
# parameter separated by spaces. One model definition per line!
# Each line contains:
#     Phi location (rad/m/m), Dec (pixel), RA (pixel), Stokes I, Q, U, V, specindex
# where the stokes Q and U parameters define the Faraday spectrum.
# The given Stokes I and V are simply assigned to all frequencies (alpha=0).

reffn = sys.argv[1]  # The reference FITS file from which to grab model header
mdl = np.loadtxt(sys.argv[2])  # The model Faraday spectrum
outfn_base = sys.argv[3]  # The base name for output FITS files
nchan = int(sys.argv[4])  # The number of channels per subband
nsb = int(sys.argv[5])  # The number of subbands, one FITS file per SB
freq0 = float(sys.argv[6])  # The frequency of the first channel, in Hz
dfreq = float(sys.argv[7])  # The channel width, in Hz
nra = int(sys.argv[8])  # The number of RA pixels
ndec = int(sys.argv[9])  # The number of Dec pixels

# TODO: convolve models with sky plane PSF

noise = True

# Pixel sizes and ref. location is determined from the header in the reffn
nra = 100  # Number of pixels in RA
cra = nra / 2
ndec = 100  # Number of pixels in Dec.
cdec = ndec / 2

freq0 = 120e6  # starting frequency
dfreq = 40e3  # channel width, applies to all channels

c = 2.99792458e8

ref_head = pyfits.getheader(reffn)

for i in range(nsb):
    sbcube = np.zeros((4, nchan, ndec, nra))
    for j in range(nchan):
        print 'j', j
        for k in range(len(mdl)):
            freq = freq0 + (i * nchan + j) * dfreq
            l2 = (c / freq) ** 2
            print 'k', k
            #Stokes I
            sbcube[0, j, cdec + mdl[k][1], cra - mdl[k][2]] = \
                mdl[k][3] * (freq / freq0)**-mdl[k][7]
                
            #Stokes V
            sbcube[3, j, cdec + mdl[k][1], cra - mdl[k][2]] = \
                mdl[k][6] * (freq / freq0)**-mdl[k][7]

            # Stokes Q & U from Faraday spectrum
            val = complex(mdl[k][4], mdl[k][5])
            val = val * np.exp(2. * complex(0, 1) * mdl[k][0] * l2) 

            sbcube[1, j, cdec + mdl[k][1], cra - mdl[k][2]] = val.real \
                * (freq / freq0)**-mdl[k][7]
            sbcube[2, j, cdec + mdl[k][1], cra - mdl[k][2]] = val.imag \
                * (freq / freq0)**-mdl[k][7]

    sbcube += np.random.normal(0,1,np.shape(sbcube))

    hdu = pyfits.PrimaryHDU(sbcube)
    hdu.header = ref_head.copy()

    hdu.header.update('ORIGIN', 'make_test_files.py', 'Origin of the data set')
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
