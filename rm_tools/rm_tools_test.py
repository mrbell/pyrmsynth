# -*- coding: utf-8 -*-
"""
Copyright 2012 Michael Bell, Henrik Junklewitz

This file is part of the pyrmsynth package.

pyrmsynth is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pyrmsynth is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pyrmsynth.  If not, see <http://www.gnu.org/licenses/>.
"""

#import datetime
import math
import sys

import pylab
import numpy

import rm_tools as R

# resolution in phi is about 40 rad/m/m with this frequency setup
NBands = 16
NFrequencies = 800
StartFreq = 1e9 + numpy.arange(NBands) * 2e8
StepFreq = 4e6


def progress(width, percent):
    """
    """
    marks = math.floor(width * (percent / 100.0))
    spaces = math.floor(width - marks)
    loader = '[' + ('=' * int(marks)) + (' ' * int(spaces)) + ']'
    sys.stdout.write("%s %d%%\r" % (loader, percent))
    if percent >= 100:
        sys.stdout.write("\n")
    sys.stdout.flush()


def rnd_test(mdl, ntrials=100, noise=1., niter=500, gain=0.1, cutoff=0.):
    """
    """

    nu = gennu(StartFreq, StepFreq, NFrequencies)

    # define image space
    nphi = 400
    dphi = 5.0
    phi = -0.5 * nphi * dphi + numpy.arange(nphi) * dphi

    cpks = numpy.zeros(ntrials)
    dpks = numpy.zeros(ntrials)

    rmsyn = R.RMSynth(nu, StepFreq, phi)
    rmcl = R.RMClean(rmsyn, niter, gain, cutoff)
    progress(20, 0.)
    for i in range(ntrials):
        ccloc, ccval, rmc, di = do_test(mdl, noise, niter, gain, cutoff, False,
                                        rms=rmsyn, rmc=rmcl)
        cpks[i] = max(abs(rmc.clean_map))
        dpks[i] = max(abs(di))
        progress(20, ((i + 1.) / ntrials) * 100.)

    plot_rnd_test(dpks, cpks, mdl)

    return dpks, cpks


def rnd_test2(mdl, ntrials=1000, noise=1., niter=500, gain=0.1, cutoff=0.):
    """
    """

    nu = gennu(StartFreq, StepFreq, NFrequencies)

    # define image space
    nphi = 400
    dphi = 5.0
    phi = -0.5 * nphi * dphi + numpy.arange(nphi) * dphi

    cpks1 = numpy.zeros(ntrials)
    dpks1 = numpy.zeros(ntrials)

    cpks2 = numpy.zeros(ntrials)
    dpks2 = numpy.zeros(ntrials)

    rmsyn = R.RMSynth(nu, StepFreq, phi)
    rmcl = R.RMClean(rmsyn, niter, gain, cutoff)
    progress(20, 0.)
    for i in range(ntrials):
        ccloc, ccval, rmc, di = do_test(mdl, noise, niter, gain, cutoff,
                                        False, rms=rmsyn, rmc=rmcl)
        cpks1[i] = abs(rmc.clean_map[170])
        dpks1[i] = abs(di[170])
        cpks2[i] = abs(rmc.clean_map[230])
        dpks2[i] = abs(di[230])
        progress(20, ((i + 1.) / ntrials) * 100.)

    plot_rnd_test(dpks1, cpks1, [mdl[1]])
    plot_rnd_test(dpks2, cpks2, [mdl[0]])

#    return dpks, cpks


def plot_rnd_test(dpks, cpks, mdl):
    """
    """

    pylab.figure()
    pylab.hist(dpks - abs(mdl[0][1]), bins=50, fc='r', label='Dirty image')
    pylab.hist(cpks - abs(mdl[0][1]), bins=50,
               hatch='/', fill=False, label='CLEAN image')
    pylab.legend()
    pylab.xlabel("Image peak - model")
    pylab.ylabel("Number")

    pylab.show()


def do_test(mdl, noise=0., niter=500, gain=0.1, cutoff=0.,
            doplot=True, rms=None, rmc=None):
    """
    A routine to test the rm_tools RMClean functionality

    input
    -----
    mdl:    a list of (phi, Q + iU) values
    noise:  the stdv of the Gaussian (white) noise to add to each channel of
            the data vector

    return
    ------
    ccloc:  a list of CLEAN component locations
    ccval:  a list of CLEAN component values
    map:    the restored CLEAN image
    phi:    coordinates for which the map are defined
    """

    # define nu range
    nu = gennu(StartFreq, StepFreq, NFrequencies)

    # define image space
    nphi = 400
    dphi = 5.0
    phi = -0.5 * nphi * dphi + numpy.arange(nphi) * dphi

    # compute l2
    l2 = convert_nu_to_l2(nu, StepFreq)
    l2 = numpy.flipud(l2)

    data = dft(mdl, noise, l2)

    # do rmsynth
    if rms is None:
        rms = R.RMSynth(nu, StepFreq, phi)
    di = rms.compute_dirty_image(data)

    if rmc is None:
        rmc = R.RMClean(rms, niter, gain, cutoff)
    else:
        rmc.reset()
    rmc.perform_clean(data)
    rmc.restore_clean_map()

    [ccloc, ccval] = sort_cc_list(rmc.cc_phi_list, rmc.cc_val_list)

    if doplot:
        plot_test(ccloc, ccval, rmc.clean_map, phi, mdl, di)

    return ccloc, ccval, rmc, di


def plot_test(ccloc, ccval, clean_map, phi, mdl, di):
    """
    function to plot test results
    """

    pylab.figure()
    pylab.plot(phi, abs(clean_map), 'k', linewidth=2, label="CLEAN spectrum")
    pylab.plot(phi, abs(di), ':k', linewidth=2, label="Dirty spectrum")
    for i in range(len(ccloc)):
        if i == 0:
            pylab.plot([ccloc[i], ccloc[i]], [0., abs(ccval[i])], 'k',
                       label="CLEAN component")
        else:
            pylab.plot([ccloc[i], ccloc[i]], [0., abs(ccval[i])], 'k')

    for i in range(len(mdl)):
        pylab.plot([mdl[i][0], mdl[i][0]], [0., abs(mdl[i][1])], ':r')
        if i == 0:
            pylab.plot(mdl[i][0], abs(mdl[i][1]), 'or', label="Model")
        else:
            pylab.plot(mdl[i][0], abs(mdl[i][1]), 'or')
    pylab.legend()
    pylab.xlim(phi[0] / 2., -1 * phi[0] / 2)
    pylab.xlabel("$\phi$")
    pylab.ylabel("FS")

    pylab.show()


def gennu(nu0, dnu, nnu):
    """
    generates a list of frequencies from the makems frequency definition format

    input
    -----
    nu0:    list of subband bins, given by the bottom edge of the subband
    dnu:    step size between bins
    nnu:    total number of bins (must be divisible by len(nu0))

    return
    ------
    nu:     list of frequencies
    """

    nchan = nnu / len(nu0)

    nu = numpy.zeros(nnu)

    for i in range(len(nu0)):
        for j in range(nchan):
            nu[i * nchan + j] = nu0[i] + dnu * (0.5 + j)

    return nu


def dft(mdl, noise, l2):
    """
    Performs the DFT from phi to lambda^2 space

    input
    -----
    mdl:     a list of (phi, Q+iU) values
    noise:   the stdv of the Gaussian (white) noise to add to each channel of
             the data vector
    l2:      the l2 vector on which to Fourier transform

    return
    ------
    data:    A complex vector containing noisy data
    """

    nl2 = len(l2)

    data = numpy.zeros(nl2, dtype=complex)
    #nphi = len(self.phi)
    #nl2 = len(self.l2)
    # trying this because I think it will be slightly faster than making
    # a bunch of self.whatever calls
    if noise > 0.:
        data.real = numpy.random.normal(0., noise, nl2)
        data.imag = numpy.random.normal(0., noise, nl2)

    for indx in range(nl2):
        temp = 2. * l2[indx]
        for jndx in range(len(mdl)):
            carg = mdl[jndx][0] * temp
            data[indx] = data[indx] + mdl[jndx][1] *\
                complex(numpy.cos(carg), numpy.sin(carg))

    return data


def convert_nu_to_l2(nu, dnu):
    """
    A routine to convert a frequency vector to a lambda^2 vector

    input
    -----
    nu:      numpy array containing a list of frequency points
    dnu:     channel width (single number applies to all nu)

    return
    ------
    l2:      numpy array with associated lambda^2 values
    """

    c2 = 8.9874e16  # speed of light squared

    # I do not use the approximate value listed in B&dB
    l2 = 0.5 * c2 * ((nu - 0.5 * dnu) ** -2 + (nu + 0.5 * dnu) ** -2)
    l2 = numpy.flipud(l2)

    return l2


def sort_cc_list(ccloc, ccval):
    """
    pass the cc_loc and cc_val lists from an RMClean object and arrange
    them so that all vals at the same location are added together.
    """

    rccloc = list()
    rccval = list()
    # Add up all of the components that are at the same location.
    while len(ccloc) > 0:
        lentry = ccloc.pop(0)[1]  # pop returns (indx, phi_val), select the val
        ventry = ccval.pop(0)

        # savetxt will not write the entire complex number, so i need to
        # store the real and imag parts separately.
#        entry = [entry[0], entry[1], entry[2], entry[3].real, entry[3].imag]
        # contains the indices of ccs that are at the same spot
        indices = list()
        for i in range(len(ccloc)):
            if ccloc[i][1] == lentry:
                indices.append(i)

        while len(indices) > 0:
            indx = indices.pop()
            #duploc = ccloc.pop(indx)
            dupval = ccval.pop(indx)
            ventry += dupval

        rccval.append(ventry)
        rccloc.append(lentry)

    return rccloc, rccval
