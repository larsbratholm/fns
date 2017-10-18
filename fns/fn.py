# MIT License
# 
# Copyright (c) 2017 Lars Andersen Bratholm
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np

from .ffn import fiffn, fffn, fifafn, fobf

def iffn(features, seed = -1, npartitions = 1, nmax = 0, exact = True, brute = False):
    """ Iteratively perform the improved fast furthest neighbours algorithm
        described in Rahman & Rochan 2016 (10.1109/ICCITECHN.2016.7860222)
        to create a furthest neighbour ordering of a feature vector, given an
        initial seed index. If ``seed = -1`` a random index is chosen.
        The algorithm uses 2-norm as a distance metric.

        ``npartitions = 1`` corresponds to using the :math:`UB` upper bound from
        the paper, while ``npartitions = 2`` corresponds to using the :math:`UB2`
        upper bound etc..

        Specifying a positive integer value for :math:`nmax` results in the
        algorithm stopping after :math:`nmax` steps, resulting in the output
        only being the first :math:`nmax` furthest neighbour indices.

        The exact furthest neighbours, given the starting seed, is returned
        unless ``exact = False``. In that case the true distance between
        :math:`x` and :math:`y` (with dimensionality :math:`d`) are approximated as

        .. math::

            dist(x,y) \\approx d \\cdot \\left( (\\mu_{x} - \\mu_{y})^2 + \\sigma_x^2 + \\sigma_y^2)

        The algorithm is implemented using an OpenMP parallel Fortran routine.

        :param features: 2D array of samples - shape (N, d)
        :type features: numpy array
        :param seed: Index of the starting seed
        :type seed: integer
        :param npartitions: Algorithm parameter
        :type npartitions: integer
        :param nmax: Force the algorithm to stop after this number of iterations
        :type nmax: integer
        :param exact: To use exact distances
        :type exact: bool

        :return: Ordered indices of the furthest neighbour search
        :rtype: numpy array
    """
    try:
        if features.ndim != 2:
            raise ValueError('expected features of dimension=2')
    except AttributeError:
        raise AttributeError('expected features as numpy array')

    if int(seed) != seed:
        raise ValueError('expected integer seed')

    nsamples = features.shape[0]
    if seed > nsamples:
        raise ValueError('seed (%d) larger or equal to size of feature array (%d)' % (seed, nsamples))

    if seed <= 0:
        seed = np.random.randint(nsamples) + 1

    if int(nmax) != nmax:
        raise ValueError('expected integer nmax')

    if nmax <= 0:
        nmax = nsamples

    if nmax > nsamples:
        raise ValueError('nmax (%d) larger to size of feature array (%d)' % (nmax, nsamples))

    if int(npartitions) != npartitions:
        raise ValueError('expected integer npartitions')

    if npartitions <= 0:
        raise ValueError('expected positive number of partitions')

    nfeatures = features.shape[1]

    if brute:
        return fobf(features, seed, nmax)

    if npartitions >= nfeatures/2:
        raise ValueError('too many partitions')


    if exact not in [True, False]:
        raise ValueError('expected boolean value for parameter exact')

    if exact:
        if npartitions == 1:
            return fffn(features, seed, nmax)
        else:
            return fiffn(features, seed, npartitions, nmax)
    else:
        return fifafn(features, seed, npartitions, nmax)
