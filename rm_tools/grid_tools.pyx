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

import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
CTYPE = np.complex128
ctypedef np.float64_t DTYPE_t
ctypedef np.complex128_t CTYPE_t

cdef extern from "gsl/gsl_sf_bessel.h":
    double gsl_sf_bessel_I0(double x)

cdef extern from "math.h":
    double exp(double theta)
    double sqrt(double x)
    double ceil(double x)

@cython.boundscheck(False)
@cython.wraparound(False)

def sample_grid(np.ndarray[DTYPE_t, ndim=1] l, np.ndarray[DTYPE_t, ndim=1] l2,
                np.ndarray[DTYPE_t, ndim=1] d2):
    """
    """

##    d = np.zeros(len(l))
    cdef int Nl = l.shape[0]
    cdef int Nl2 = l2.shape[0]
    cdef int lndx = 0
    cdef np.ndarray[DTYPE_t, ndim=1] d = np.zeros(Nl)
    cdef float l0 = l[0]
    cdef float dl = l[1] - l0
    cdef Py_ssize_t i
    cdef float temp = 0.

    for i in range(Nl2):
        lndx = int((l2[i] - l0)/dl + 0.5)
        if lndx >= 0 and lndx < Nl:
            d[lndx] = d[lndx] + d2[i]
            
            
    return d


def sample_grid_complex(np.ndarray[DTYPE_t, ndim=1] l,
                        np.ndarray[DTYPE_t, ndim=1] l2,
                        np.ndarray[CTYPE_t, ndim=1] d2):
    """
    """

##    d = np.zeros(len(l))
    cdef int Nl = l.shape[0]
    cdef int Nl2 = l2.shape[0]
    cdef int lndx = 0
    cdef np.ndarray[CTYPE_t, ndim=1] d = np.zeros(Nl, dtype=CTYPE)
    cdef float l0 = l[0]
    cdef float dl = l[1] - l0
    cdef Py_ssize_t i
    cdef float temp = 0.

    for i in range(Nl2):
        lndx = int((l2[i] - l0)/dl + 0.5)
        if lndx >= 0 and lndx < Nl:
            d[lndx] = d[lndx] + d2[i]

    return d


def grid_1d(np.ndarray[DTYPE_t, ndim=1] d, np.ndarray[DTYPE_t, ndim=1] l, \
    double dx, int m, double alpha):
    """
    Given a data array d, sampled at locations l, this routine grids the data
    into pixels of size dx.  The gridded data array and the axis on which it is
    defined is returned.
    """

    cdef int N = d.shape[0]

    cdef np.ndarray[DTYPE_t, ndim=1] d2 = np.zeros(N*m)
    cdef np.ndarray[DTYPE_t, ndim=1] l2 = np.zeros(N*m)

    cdef Py_ssize_t indx, xndx, kndx

    cdef double val, x, xref, gcf_val, xg, temp, temp3

    # see Beatty et al. (2005)
    cdef double beta = np.pi*sqrt((m/alpha)**2.*(alpha - 0.5)**2 - 0.8)
    cdef double temp2 = (1./m/dx)
    for indx in range(N):

        val = d[indx]
        x = l[indx]
        xref = ceil((x - 0.5*m*dx)/dx)*dx

        for xndx in range(m):

            xg = xref + xndx*dx

            kndx = indx*m + xndx

            temp = xg-x
            temp3 = (2.*temp*temp2)**2.
            temp3 = sqrt(1 - temp3)
            temp3 = beta*temp3
            gcf_val = temp2*gsl_sf_bessel_I0(temp3)

            d2[kndx] = val*gcf_val
            l2[kndx] = xg

    return l2, d2


def grid_1d_complex(np.ndarray[CTYPE_t, ndim=1] d,
                    np.ndarray[DTYPE_t, ndim=1] l,
                    double dx, int m, double alpha):
    """
    Given a data array d, sampled at locations l, this routine grids the data
    into pixels of size dx.  The gridded data array and the axis on which it is
    defined is returned.
    """

    cdef int N = d.shape[0]

    cdef np.ndarray[CTYPE_t, ndim=1] d2 = np.zeros(N*m, dtype=CTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] l2 = np.zeros(N*m)

    cdef Py_ssize_t indx, xndx, kndx

    cdef DTYPE_t x, xref, gcf_val, xg, temp, temp3
    cdef CTYPE_t val

    # see Beatty et al. (2005)
    cdef double beta = np.pi*sqrt((m/alpha)**2.*(alpha - 0.5)**2 - 0.8)
    cdef double temp2 = (1./m/dx)
    for indx in range(N):

        val = d[indx]
        x = l[indx]
        xref = ceil((x - 0.5*m*dx)/dx)*dx

        for xndx in range(m):

            xg = xref + xndx*dx

            kndx = indx*m + xndx

            temp = xg-x
            temp3 = (2.*temp*temp2)**2.
            temp3 = sqrt(1 - temp3)
            temp3 = beta*temp3
            gcf_val = temp2*gsl_sf_bessel_I0(temp3)

            d2[kndx] = val*gcf_val
            l2[kndx] = xg

    return l2, d2



def gridnorm(np.ndarray[double, ndim=1] x, double dk, int W, double beta):

    cdef int xlen = x.shape[0]
    cdef np.ndarray[double, ndim=1] c
    cdef np.ndarray[double, ndim=1] temp
    cdef Py_ssize_t indx

    # avoid having to use the cmath sqrt function

    temp = np.sqrt(beta**2. - (np.pi*W*dk*x)**2.)
    c = (np.exp(temp) - np.exp(-1.*temp))/2./temp

    return c
