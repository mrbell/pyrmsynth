"""
rm_tools.py
Written by Michael Bell

Base classes for doing RM synthesis and RM CLEANing.  These classes operate
on single lines of sight.

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


# TODO: Check for consistency of weighting and l20 stuff
# TODO: Test that normalization is OK
# TODO: Check RMClean


import numpy
import math
from sys import stdout
import grid_tools as G
#import pylab

VERSION = '1.2.0'

toplot = True


def gaussian_taper(l2, l20, fwhm):
    """
    gaussian_taper(l2, l20, fwhm)

    Compute a gaussian taper to use as a weighting function.

    Inputs:
        l2-     The frequency axis over which the taper function should be
                computed.  Should be the same as is used to init the RMSynth
                class. Can be defined in lambda^2 or freq.
        l20-    Average of the gaussian.
        fwhm-   The FWHM of the gaussian defined in whatever units l2 is in.
    Outputs:
        1-      Gaussian taper weighting function with the same length as l2.
    """

    #w = numpy.zeros(len(l2))
    sdv = fwhm / 2. / math.sqrt(2. * math.log(2.))

    #for i in range(len(l2)):
        #w[i] = math.exp(-0.5 * (l2[i] - l20) ** 2. * sdv ** -2.)

    w = numpy.exp(-0.5 * (l2 - l20) ** 2. * sdv ** -2.)

    return w


class RMSynth:
    """
    An RM synthesis package whereby the RMSF and dirty map are computed
    explicitly (i.e. not using FFT). Once initialized, use compute_dirty_image
    to perform the inversion.
    """

    def __init__(self, nu, dnu, phi, weights=None):
        """
        RMSynth(nu, dnu, phi, weights=None, isl2=False)

        Initializes the RMSynth class.  For an image, this needs only be done
        once.  Each LOS may be inverted individually using the
        compute_dirty_image method.

        Inputs:
            nu-     A vector containing frequency(Hz) values. Either way they
                    must be ordered from lowest to highest value.
            dnu-    Width of the frequency channels (in Hz).
            phi-    The requested phi axis.  Use numpy.arange to construct.
            weights-A vector containing weights. Must be the same length as nu.
                    Values should be between 0 and 1.  If no weights vector is
                    given, it will be assumed that each value has weight 1.

        Outputs:
            None
        """
        # parameters used for gridding
        self.m = 6  # number of grid cells overwhich the GCF spans
        self.alpha = 1.5  # oversampling ratio

        self.dphi = phi[1] - phi[0]
        self.nphi = len(phi)

        self.dl2 = 1. / self.nphi / self.dphi / self.alpha
        self.nl2 = self.alpha * self.nphi

        stdout.write('RMSynth Initializing... ')
        stdout.flush()

        # equal weighting for all frequencies by default.
        if weights is None:
            weights = numpy.ones(len(nu))

        if len(nu) != len(weights):
            msg = 'The length of weight and frequency vectors must be equal!'
            raise IndexError(msg)

        # for now store the ungridded l2s (for diagnostics)
        self.l2_nonuni = self.convert_nu_to_l2(nu, dnu)

        self.weights_nonuni = numpy.flipud(weights)

        # another gridding parameter, derived from m and alpha
        self.beta = numpy.pi *\
            numpy.sqrt((self.m / self.alpha) ** 2. *
                       (self.alpha - 0.5) ** 2 - 0.8)

        self.l20 = self.compute_l20()
        self.l2i = self.l2_nonuni[0] - self.m * self.dl2 * numpy.pi

        # this axis is actually \lambda^2/\pi
        self.l2 = numpy.arange(0, self.nl2, 1) * self.dl2 + self.l2i / numpy.pi
        self.l2_beam = numpy.arange(0, self.nl2 * 2., 1) * self.dl2 / 2. \
            + self.l2i / numpy.pi

        self.phi = -self.dphi * self.nphi * 0.5 \
            + numpy.arange(0, self.nphi, 1) * self.dphi

        self.grid_corr = G.gridnorm(self.phi, self.dl2, self.m, self.beta)
        self.tndx_phi = int(0.5 * self.nphi * (self.alpha - 1))

        [self.rmsf, self.rmsf_phi] = self.compute_rmsf()

        stdout.write('Complete.\n')

    def compute_dirty_image(self, pol):
        """
        RMSynth.compute_dirty_image(pol)

        Inverts the given polarization vector to arrive at the dirty dispersion
        function.

        Inputs:
            pol-    An Nx1 numpy array containing the complex pol. emission
                    as a function of frequency (or lambda^2 if that is how the
                    freq. axis has been defined).  Must be in order from low to
                    high along freq. axis.
        Outputs:
            1-      A complex valued map of the dirty dispersion function.
        """

        pol = numpy.flipud(pol)
        pol = pol * self.weights_nonuni

        [l2_grid, w_grid] = G.grid_1d_complex(pol, self.l2_nonuni / numpy.pi,
                                              self.dl2, self.m, self.alpha)

        # Put the convolved points on a grid
        polgrid = G.sample_grid_complex(self.l2, l2_grid, w_grid)

        residual_map = self.K * numpy.fft.fftshift(numpy.fft.fft(polgrid))
        residual_map = residual_map[self. tndx_phi:self. tndx_phi + self.nphi]

        #for indx in range(self.nphi):

            #rot = 2. * (self.l20 - self.l2i) * self.phi[indx]
            ## shift results by the L20 and by the initial l2 value to
            ## account for the fact that we start at L2!=0 as FFT expects
            #residual_map[indx] = numpy.complex(math.cos(rot),
                #math.sin(rot)) * residual_map[indx]

        rot = 2. * (self.l20 - self.l2i) * self.phi
        # shift results by the L20 and by the initial l2 value to
        # account for the fact that we start at L2!=0 as FFT expects
        #residual_map = numpy.complex(numpy.cos(rot),
                                     #numpy.sin(rot)) * residual_map
        residual_map *= numpy.exp(complex(0, 1) * rot)

        return residual_map / self.grid_corr

    def compute_l20(self):
        """
        RMSynth.compute_l20()

        Calculate the weighted average of the lambda^2 axis. During the inverse
        process, the measured EVPA are derotated to this value rather than to
        lambda^2 = 0.  Doing so, the imag. part of the RMSF is better behaved
        throughout the central lobe (it will be flat at phi=0)

        Inputs:
            None
        Outputs:
            1-  The weighted average of the lambda^2 axis.
        """

        l20 = 0.

        for i in range(len(self.l2_nonuni)):
            l20 += self.weights_nonuni[i] * \
                self.l2_nonuni[i] / sum(self.weights_nonuni)

        return l20

    def compute_rmsf(self):
        """
        RMSynth.compute_rmsf()

        Inverts the given polarization vector to arrive at the dirty dispersion
        function.

        Inputs:
            None
        Outputs:
            1-  A map containing the RMSF.  A numpy array of complex values.
            2-  A map containing the RMSF phi axis.
        """

        # RMSF image size must be twice that of the df image for CLEANing
        phi_center = 0.5 * ((numpy.max(self.phi) + self.dphi) +
                            numpy.min(self.phi))
        phi_range = (numpy.max(self.phi) + self.dphi) - numpy.min(self.phi)

        rmsf_phi = phi_center - phi_range +\
            numpy.arange(2 * self.nphi) * self.dphi

        try:

            # Convolve weights with the GCF
            [l2_grid, w_grid] = G.grid_1d(self.weights_nonuni,
                                          self.l2_nonuni / numpy.pi,
                                          self.dl2 / 2., self.m, self.alpha)

            # Put the convolved points on a grid
            weights4rmsf = G.sample_grid(self.l2_beam, l2_grid, w_grid)

        except TypeError:
            print type(self.weights)
            print len(self.weights)
            print type(self.l2_nonuni)
            print len(self.l2_nonuni)
            print self.dphi
            print self.nphi * 2.
            raise

        rmsf = numpy.fft.fftshift(numpy.fft.fft(weights4rmsf))
        tndx_phi = int(self. nphi * (self.alpha - 1))
        rmsf = rmsf[tndx_phi:tndx_phi + 2 * self.nphi]

        #for indx in range(2 * self.nphi):
            #rot = 2. * (self.l20 - self.l2i) * rmsf_phi[indx]
            ## shift results by the L20 and by the initial l2 value to account
            ## for the fact that we start at L2!=0 as FFT expects
            #rmsf[indx] = numpy.complex(math.cos(rot), math.sin(rot)) *\
                #rmsf[indx]

        rot = 2. * (self.l20 - self.l2i) * rmsf_phi
        # shift results by the L20 and by the initial l2 value to account
        # for the fact that we start at L2!=0 as FFT expects
        rmsf *= numpy.exp(complex(0, 1) * rot)

        gridnorm = G.gridnorm(rmsf_phi, self.dl2 / 2., self.m, self.beta)

        rmsf = rmsf / gridnorm

        # The normalization should be as such that the rmsf has a peak
        # value of 1
        self.K = 1. / abs(rmsf[self.nphi])

        # normalize the rmsf
        rmsf = rmsf * self.K

        # adjust the normalization
        # because the rmsf has twice as many pixels as the image
        self.K = 2. * self.K

        return rmsf, rmsf_phi

    def convert_nu_to_l2(self, nu, dnu):
        """
        RMSynth.convert_nu_to_l2(nu, dnu)

        Given a vector of frequencies, with values that describe the central
        frequency of a channel each having width dnu, the associated lambda^2
        values are computed.

        Inputs:
            nu- An Nx1 numpy array containing frequency values.  These values
                are the central frequency of the channel. They must be given in
                Hz. Must be ordered from low to high values of frequency.
            dnu-A scalar containing the channel width in Hz.
        Outputs:
            1-  The associated lambda^2 axis in m^2.  The axis will be arranged
                from low to high values of lambda^2
        """

        c2 = 8.9874e16  # speed of light squared

        # I do not use the approximate value listed in B&dB
        l2 = 0.5 * c2 * ((nu - 0.5 * dnu) ** -2 + (nu + 0.5 * dnu) ** -2)
        l2 = numpy.flipud(l2)

        return l2


class RMClean:

    """
    Simple implementation of the RM Clean algorithm. Use perform_clean() to
    actually start the CLEAN algorithm.  Use restore_clean_map() to construct
    the final CLEAN map from the model components.
    """

    def __init__(self, rms, niter=500, gain=0.1, cutoff=0.):
        """
        RMClean(rms, niter=100, gain=0.1, cutoff=0)

        Initializes the RMClean class, parameters such as the number of
        iterations persist until the class is reinitialized.  So it is possible
        to, e.g., run perform_clean(), inspect the output, then continue
        cleaning by increasing the number of iterations. To restart on the same
        or another line of sight, the class should be reinitialized.

        Inputs:
            rms-    an instance of the rmsynth class
            niter-  number of CLEAN iterations
            gain-   loop gain
            cutoff- the algorithm stops when the absolute peak value of the
                    residual image is below this value
        Outputs:
            None
        """

        self.synth = rms
        self.niter = niter
        self.gain = gain
        self.cutoff = cutoff

        self.fwhm_restoring_sf = self.gauss_fit()

        self.current_iter = 0
        self.cc_phi_list = list()
        self.cc_val_list = list()

        self.residual_map = None
        self.clean_map = None

        # Compute the CLEAN beam and convert to data space
        # Store for use with all LOS
        clean_beam = self.compute_clean_beam()
        self.clean_beam_l2 = numpy.fft.ifft(numpy.fft.fftshift(clean_beam))

    def reset(self):
        """
        Re-initialize the class, simply getting rid of the residual and clean
        maps, as well as resetting the CC lists.  Global parameters (non-LOS
        dependent) do not get reset.
        """
        self.residual_map = None

        self.current_iter = 0
        self.cc_phi_list = list()
        self.cc_val_list = list()

        self.clean_map = None

    def compute_clean_beam(self):
        """
        RMClean.compute_clean_beam(len)

        Compute a gaussian CLEAN beam to be convolved with the CC list.

        Inputs:
            None
        Outputs:
            A map of gaussian clean beam of length len with the peak of the
                beam located at the central pixel.
        """

        sdev = self.fwhm_restoring_sf / 2. / math.sqrt(2. * math.log(2.))
        cphi = 0    # central phi value

        #nphi = len(self.synth.phi)

        #clean_beam = numpy.zeros(nphi, dtype=complex)

        #for i in range(nphi):
            #clean_beam[i] = numpy.complex(math.exp(-0.5 *
                #(self.synth.phi[i] - cphi) ** 2. * sdev ** -2.), 0)
        clean_beam = numpy.exp(-0.5 * complex(1,0) *
                (self.synth.phi - cphi) ** 2. * sdev ** -2.)

        return clean_beam

    def convolve_beam_w_cc_list(self):
        """
        RMClean.convolve_beam_w_cc_list()

        Convolve the clean component list with the gaussian CLEAN beam. Returns
        the restored image.

        Inputs:
            None
        Outputs:
            1-  Map of the CLEAN components convolved with the restoring beam

        """
        # map containing the CC point sources at the appropriate locations
        # to be convolved w/ CLEAN beam
        nphi = len(self.residual_map)

        ps_map = numpy.zeros(nphi, dtype=complex)

        for i in range(len(self.cc_phi_list)):
            ps_map[self.cc_phi_list[i][0]] = ps_map[self.cc_phi_list[i][0]] \
                + self.cc_val_list[i]

        ps_map_l2 = numpy.fft.ifft(numpy.fft.fftshift(ps_map))

        conv_image_l2 = ps_map_l2 * self.clean_beam_l2
        conv_image = numpy.fft.ifftshift(numpy.fft.fft(conv_image_l2)) * nphi

        return conv_image

    def find_peak_phi(self, cross_corr=False, img=None):
        """
        RMClean.find_peak_phi(cross_corr=False, img=None)

        Finds the index of the phi value at the peak of the residual image (or
        the image passed to the function). The peak can be found by a simple
        search algorithm, or by first cross correlating the image with the
        RMSF.

        Inputs:
            cross_corr- Perform a cross correlation between RMSF and res. map
                        prior to seeking for the peak ala Heald 09
            img-        Image to search through.  If none is given, the current
                        residual image will be used.
        Outputs:
            1-          The pixel index of the phi value at the map peak
        """
        phi_ndx = -1

        peak_res = 0.

        if img is None:
            img = self.residual_map

        if not cross_corr:
            for i in range(len(img)):
                if (abs(img[i]) > peak_res):
                    phi_ndx = i
                    peak_res = abs(img[i])
        else:
            temp_map_fft = numpy.fft(img)
            temp_map = numpy.ifft(temp_map_fft.conjugate() *
                                  numpy.fft(self.synth.rmsf))
            phi_ndx = self.find_peak_phi(img=temp_map)

        return phi_ndx

    def gauss_fit(self):
        """
        RMClean.gauss_fit()

        Calculates the FWHM of the CLEAN beam based on the size of the main
        peak of the RMSF.

        Inputs:
            None
        Outputs:
            1-  Returns the FWHM of the restoring spread function.
        """

        delta_l2 = self.synth.l2_nonuni[len(self.synth.l2_nonuni) - 1] -\
            self.synth.l2_nonuni[0]
        return 2 * math.sqrt(3) / delta_l2

    def perform_clean(self, pol=None):
        """
        RMClean.perform_clean(pol=None)

        Starts the CLEAN algorithm using parameters defined during class
        initialization.

        Inputs:
            pol-    A vector containing the complex polarized emission as a
                    function of freq. The vector must have the same length as
                    the number of frequencies that are provided and listed
                    from low to high frequency (or from low to high lambda^2 if
                    the freq. axis has been defined in these units). If CLEAN
                    has already been performed on a data set and one wishes to
                    continue cleaning, no pol data is needed, the algorithm
                    continues where it left off.
        Outputs:
            None
        """

        if pol is not None and (len(self.synth.l2_nonuni) != len(pol)):
            msg = 'Length of frequency and polarization vectors not equal!'
            raise IndexError(msg)

        if self.residual_map is None and pol is not None:
            self.residual_map = self.synth.compute_dirty_image(pol)
        elif self.residual_map is None and pol is None:
            msg = 'The class has not been initialized yet or no pol vector ' +\
                'has been given.  Nothing to clean.'
            raise ValueError(msg)

        peak_phi_indx = self.find_peak_phi()
        res_max = self.residual_map[peak_phi_indx]
        total_flux = 0.

        while (abs(res_max) >= self.cutoff) and (self.current_iter <
                                                 self.niter):

            cc_pol = res_max * self.gain
            self.cc_phi_list.append([peak_phi_indx,
                                     self.synth.phi[peak_phi_indx]])
            self.cc_val_list.append(cc_pol)
            total_flux = total_flux + cc_pol
            self.subtract_rmsf_from_residual(peak_phi_indx, cc_pol)

            self.current_iter = self.current_iter + 1
            peak_phi_indx = self.find_peak_phi()
            res_max = self.residual_map[peak_phi_indx]

    def restore_clean_map(self):
        """
        RMClean.restore_clean_map()

        Starts the CLEAN algorithm using parameters defined during class
        initialization.

        Inputs:
            None
        Outputs:
            1-  Final restored image including the residual map
        """

        self.clean_map = self.convolve_beam_w_cc_list()

        self.clean_map += self.residual_map

    def subtract_rmsf_from_residual(self, indx, val):
        """
        RMClean.restore_clean_map(indx, val)

        Subtracts the full RMSF, weighted by val, from the current residual
        image at the location defined by indx.

        Inputs:
            indx-   The pixel value where the peak of the RMSF should be
                    shifted prior to subtraction.
            val-    The value by which the RMSF should be weighted.
        Outputs:
            None
        """

        shft = self.synth.nphi - indx

        self.residual_map = self.residual_map -\
            val * self.synth.rmsf[shft:shft + self.synth.nphi]
