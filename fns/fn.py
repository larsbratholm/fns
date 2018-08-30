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

from .ffn import fiffn, fifafn, fobf, fifafn_l1, fobf_l1

def iffn(features, seed = -1, npartitions = 1, nmax = 0, nneighbours = -1):
    """ Iteratively perform the improved fast furthest nneighbours algorithm
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

        By default the furthest neighbour to all previously selected points are
        found in each iteration, however one can specify that only the last
        ``nneighbours`` points should be used, which results in the (now approximate) 
        algorithm having linear scaling.

        The algorithm is implemented using an OpenMP parallel Fortran routine.

        :param features: 2D array of samples - shape (N, d)
        :type features: numpy array
        :param seed: Index of the starting seed
        :type seed: integer
        :param npartitions: Algorithm parameter
        :type npartitions: integer
        :param nmax: Force the algorithm to stop after this number of iterations
        :type nmax: integer

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

    if npartitions >= nfeatures/2:
        raise ValueError('too many partitions')

    if nneighbours < 0:
        nneighbours = nsamples
    elif nneighbours == 0:
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamploes)
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamples)
    elif nneighbours > nsamples:
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamploes)
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamples)


    return fiffn(features, seed, npartitions, nmax, nneighbours)

def brutefn(features, seed = -1, metric = 'l2', nmax = 0, nneighbours = -1, store_memory = False):
    """ Brute force furthest neighbour ordering of a feature vector, given an
        initial seed index. If ``seed = -1`` a random index is chosen.
        The algorithm uses either l1-norm (``metric='l1'`` or ``metric='manhattan'``) or l2-norm
        (``metric='l2'`` or ``metric='euclidean'``) as a distance metric.

        Specifying a positive integer value for :math:`nmax` results in the
        algorithm stopping after :math:`nmax` steps, resulting in the output
        only being the first :math:`nmax` furthest neighbour indices.

        By default the furthest neighbour to all previously selected points are
        found in each iteration, however one can specify that only the last
        ``nneighbours`` points should be used, which results in the (now approximate) 
        algorithm having linear scaling.

        ``store_memory`` specifies if all distances should be precomputed or computed on the fly.
        Storing all distances in memory is much faster but might not be possible depending on the
        problem size.

        The algorithm is implemented using an OpenMP parallel Fortran routine.

        :param features: 2D array of samples - shape (N, d)
        :type features: numpy array
        :param seed: Index of the starting seed
        :type seed: integer
        :param metric: Distance metric
        :type metric: string
        :param nmax: Force the algorithm to stop after this number of iterations
        :type nmax: integer

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

    nfeatures = features.shape[1]

    if nneighbours < 0:
        nneighbours = nsamples
    elif nneighbours == 0:
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamples)
    elif nneighbours > nsamples:
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamples)

    if metric in ['l1', 'manhattan']:
        return fobf_l1(features, seed, nmax, nneighbours)
    elif metric in ['l2', 'euclidean']:
        return fobf(features, seed, nmax, nneighbours)
    else:
        raise ValueError('unknown metric %d (Allowed metrics are "l1", "manhattan", "l2", "euclidean")' % metric)

def fastfn(features, seed = -1, metric = 'l2', npartitions = 1, nmax = 0, approx_order = 4, nneighbours = -1):
    """ Turns out that while the improved fast furthest nneighbours algorithm
        described in Rahman & Rochan 2016 (10.1109/ICCITECHN.2016.7860222) skips a lot
        of distance calculations, it's actually not really any faster than the
        brute force algorithm.
        In many cases you just want the selected points to be far apart, and
        don't care as much if it's the exact furthest neighbour.
        This function uses approximations to the l1-norm (``metric='l1'`` or ``metric='manhattan'``) or l2-norm
        (``metric='l2'`` or ``metric='euclidean'``) to achieve this.

        For :math:`N` samples of :math:`d` dimensions, the mean :math:`\\mu` and
        and standard deviation :math:`\\sigma` is calculated for each sample.
        The squared l2 distance between samples :math:`X` and :math:`Y` is approximated as

        .. math::

            dist(X,Y)^2 = d \\cdot \\left((\\mu_{X} - \\mu_{Y})^2 + \\sigma_{X}^2 + \\sigma_{Y}^2 \\right)

        and the l1 distance as

        .. math::

            dist(X,Y) = d \\cdot \\left\\[ \\sqrt{\\frac{2}{\\pi}} \\cdot \\sqrt\\left(
                \\sigma_x^2 + \\sigma_y^2 \\right) \\cdot
                \\exp\\left(-\\frac{(\\mu_x - \\mu_y)^2}{2\\left(\\sigma_x^2+\\sigma_y^2\\right)}\\right) 
                + |\\mu_x - \\mu_y| \\text{erf}\\left( \\frac{|\\mu_x - \\mu_y||}{\\sqrt{2\\left(\\sigma_x^2 + \\sigma_y^2\\right)}} \\right) \\right\\]

        to create a furthest neighbour ordering of a feature vector, given an
        initial seed index. If ``seed = -1`` a random index is chosen.
        The algorithm uses 2-norm as a distance metric.

        ``npartitions`` splits up the features into smaller vectors of size :math:`d_{i}`, where

        .. math::

            d = \\sum_i^{n} d_{i}


        Specifying a positive integer value for :math:`nmax` results in the
        algorithm stopping after :math:`nmax` steps, resulting in the output
        only being the first :math:`nmax` furthest neighbour indices.

        ``approx_order`` is a parameter in how exponential functions are approximated and is
        set to a sensible default.

        By default the furthest neighbour to all previously selected points are
        found in each iteration, however one can specify that only the last
        ``nneighbours`` points should be used, which results in the
        algorithm having linear scaling.

        The algorithm is implemented using an OpenMP parallel Fortran routine.

        :param features: 2D array of samples - shape (N, d)
        :type features: numpy array
        :param seed: Index of the starting seed
        :type seed: integer
        :param metric: Distance metric
        :type metric: string
        :param npartitions: Algorithm parameter
        :type npartitions: integer
        :param nmax: Force the algorithm to stop after this number of iterations
        :type nmax: integer

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

    if npartitions >= nfeatures/2:
        raise ValueError('too many partitions')

    if nneighbours < 0:
        nneighbours = nsamples
    elif nneighbours == 0:
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamploes)
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamples)
    elif nneighbours > nsamples:
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamploes)
        raise ValueError('nneighbours must be between 0 and number of samples (%d)' % nsamples)

    if metric in ['l1', 'manhattan']:
        if approx_order < 1 or approx_order > 12:
            raise ValueError('expected approx_order to be in the range of 1-12')
        return fifafn_l1(features, seed, npartitions, nmax, approx_order, nneighbours)
    elif metric in ['l2', 'euclidean']:
        return fifafn(features, seed, npartitions, nmax, nneighbours)
    else:
        raise ValueError('unknown metric %d (Allowed metrics are "l1", "manhattan", "l2", "euclidean")' % metric)
