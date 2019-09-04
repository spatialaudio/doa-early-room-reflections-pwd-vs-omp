"""Identification of Early Reflections by Spherical Microphone Arrays"""

# Sascha Spors, University of Rostock

# source angles point towards source for plane wave/point source (= incidence angle of pwd)

import numpy as np
from numpy.core.umath_tests import inner1d
from scipy import signal
from scipy import special
import micarray
import sfs
from scipy.io import loadmat

# for modified blob_log function
from scipy.ndimage import gaussian_laplace
from skimage.util import img_as_float
from skimage.feature import peak_local_max
from skimage._shared.utils import assert_nD
from skimage.feature.blob import _prune_blobs


def mixing_time(L):
    """Predictor for physical mixing time of a rectangular room.

    Eq.(12) and (14) from
    Lindau, A.; Kosanke, L.; Weinzierl, S.: "Perceptual evaluation
    of model- and signal-based predictors of the mixing time in binaural
    room impulse responses", Journal of the AES, 60(11), 2012.
    """
    V = np.prod(L)
    A = 2*L[1]*L[2] + 2*L[0]*L[2] + 2*L[0]*L[1]
    tmp50  = 20.08 * V/A + 12
    tmp95  = 0.0117 * V + 50.1

    return tmp50, tmp95


def spherical_image_sources(N, k, azi, elev, R, x_src, x_mic, L, max_order, coeffs, setup='rigid'):
    """Surface pressure of a point source and its image sources on a sphere"""
    x_is, order = sfs.util.image_sources_for_box(x_src, L, max_order)
    strength = np.prod(coeffs**order, axis=1)
    xd = x_is - x_mic
    azi_is, elev_is, r_is = sfs.util.cart2sph(xd[:,0], xd[:,1], xd[:,2])

    p = np.zeros((len(k), (N+1)**2))
    Y_p = micarray.modal.angular.sht_matrix(N, azi, elev)
    for i in range(len(x_is)):
        bn = micarray.modal.radial.spherical_ps(N, k, R, r_is[i], setup)
        bn = micarray.modal.radial.repeat_n_m(bn)
        Y_pw = micarray.modal.angular.sht_matrix(N, azi_is[i], elev_is[i])

        if strength[i] != 0:
            p = p + strength[i] * np.squeeze(micarray.util.matdiagmul(np.conj(Y_pw), bn))

    return np.matmul(p, Y_p.T)


def upper_frequency_limit(N, R):
    """Approximate upper frequency limit of the spherical harmonics representation."""
    return N*343/(2*np.pi*R)


def schroeder_frequency(L, coeffs):
    """Calculate the Schroeder frequency."""
    V = np.prod(L)
    A = L[1]*L[2]*(coeffs[0] + coeffs[1]) + L[0]*L[2]*(coeffs[2] + coeffs[3]) + L[0]*L[1]*(coeffs[4] + coeffs[5])
    T = 13.8 * 4 * V / (343 * A)

    return 2000 * np.sqrt(T/V)


def pwd_matched_filter(N, k, azi, elev, R, weights, azi_pwd, elev_pwd, p, setup='rigid', window=None):
    """Plane wave decomposition by matched filter."""
    Y_p = micarray.modal.angular.sht_matrix(N, azi, elev, weights)
    Y_q = micarray.modal.angular.sht_matrix(N, azi_pwd, elev_pwd)
    bn = micarray.modal.radial.spherical_pw(N, k, R, setup)
    if window is not None:
        dn = bn * window
    dn = micarray.modal.radial.repeat_n_m(bn)
    A = 1/(4*np.pi) * np.matmul(micarray.util.matdiagmul(np.conj(Y_q), dn), Y_p.T)

    return np.squeeze(np.matmul(np.conj(A), np.expand_dims(p, 2)))
