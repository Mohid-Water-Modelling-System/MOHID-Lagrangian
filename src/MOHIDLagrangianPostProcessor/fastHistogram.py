'''
This module contains a fast replacement for numpy's histogramdd and histogram2d.
Two changes were made. The first was replacing

np.digitize(a, b)

with 

np.searchsorted(b, a, "right")

This performance bug is explained on https://github.com/numpy/numpy/issues/2656
The speedup is around 2x for big number of bins (roughly >100).
It assumes that the bins are monotonic.

The other change is to allow lists of weight arrays. This is advantageous for
resampling as there is just one set of coordinates but several data arrays (=weights).
Therefore repeated computations are prevented.
'''


import numpy as np
from numpy import atleast_2d, asarray, zeros, ones, array, atleast_1d, arange,\
    isscalar, linspace, diff, empty, around, where, bincount, sort, log10,\
    searchsorted

def histogramdd(sample, bins=10, range=None, normed=False, weights=None):
    """
    Compute the multidimensional histogram of some data.

    Parameters
    ----------
    sample : array_like
        The data to be histogrammed. It must be an (N,D) array or data
        that can be converted to such. The rows of the resulting array
        are the coordinates of points in a D dimensional polytope.
    bins : sequence or int, optional
        The bin specification:

        * A sequence of arrays describing the bin edges along each dimension.
        * The number of bins for each dimension (nx, ny, ... =bins)
        * The number of bins for all dimensions (nx=ny=...=bins).

    range : sequence, optional
        A sequence of lower and upper bin edges to be used if the edges are
        not given explicitly in `bins`. Defaults to the minimum and maximum
        values along each dimension.
    normed : bool, optional
        If False, returns the number of samples in each bin. If True,
        returns the bin density ``bin_count / sample_count / bin_volume``.
    weights : array_like (N,), optional
        An array of values `w_i` weighing each sample `(x_i, y_i, z_i, ...)`.
        Weights are normalized to 1 if normed is True. If normed is False,
        the values of the returned histogram are equal to the sum of the
        weights belonging to the samples falling into each bin.
        Weights can also be a list of (weight arrays or None), in which case
        a list of histograms is returned as H.

    Returns
    -------
    H : ndarray
        The multidimensional histogram of sample x. See normed and weights
        for the different possible semantics.
    edges : list
        A list of D arrays describing the bin edges for each dimension.

    See Also
    --------
    histogram: 1-D histogram
    histogram2d: 2-D histogram

    Examples
    --------
    >>> r = np.random.randn(100,3)
    >>> H, edges = np.histogramdd(r, bins = (5, 8, 4))
    >>> H.shape, edges[0].size, edges[1].size, edges[2].size
    ((5, 8, 4), 6, 9, 5)

    """

    try:
        # Sample is an ND-array.
        N, D = sample.shape
    except (AttributeError, ValueError):
        # Sample is a sequence of 1D arrays.
        sample = atleast_2d(sample).T
        N, D = sample.shape
    
    if weights is None:
        W = None
    else:    
        try:
            # Weights is a 1D-array
            weights.shape
            W = -1
        except (AttributeError, ValueError):
            # Weights is a list of 1D-arrays or None's
            W = len(weights)

    if W == -1 and weights.ndim != 1:
        raise AttributeError('Weights must be a 1D-array, None, or a list of both')

    nbin = empty(D, int)
    edges = D*[None]
    dedges = D*[None]
    if weights is not None:
        if W == -1:
            weights = asarray(weights)
            assert weights.shape == (N,)
        else:
            for i in arange(W):
                if weights[i] is not None:
                    weights[i] = asarray(weights[i])
                    assert weights[i].shape == (N,)

    try:
        M = len(bins)
        if M != D:
            raise AttributeError(
                'The dimension of bins must be equal to the dimension of the '
                ' sample x.')
    except TypeError:
        # bins is an integer
        bins = D*[bins]

    # Select range for each dimension
    # Used only if number of bins is given.
    if range is None:
        # Handle empty input. Range can't be determined in that case, use 0-1.
        if N == 0:
            smin = zeros(D)
            smax = ones(D)
        else:
            smin = atleast_1d(array(sample.min(0), float))
            smax = atleast_1d(array(sample.max(0), float))
    else:
        smin = zeros(D)
        smax = zeros(D)
        for i in arange(D):
            smin[i], smax[i] = range[i]

    # Make sure the bins have a finite width.
    for i in arange(len(smin)):
        if smin[i] == smax[i]:
            smin[i] = smin[i] - .5
            smax[i] = smax[i] + .5

    # Create edge arrays
    for i in arange(D):
        if isscalar(bins[i]):
            if bins[i] < 1:
                raise ValueError(
                    "Element at index %s in `bins` should be a positive "
                    "integer." % i)
            nbin[i] = bins[i] + 2  # +2 for outlier bins
            edges[i] = linspace(smin[i], smax[i], nbin[i]-1)
        else:
            edges[i] = asarray(bins[i], float)
            nbin[i] = len(edges[i]) + 1  # +1 for outlier bins
        dedges[i] = diff(edges[i])
        if np.any(np.asarray(dedges[i]) <= 0):
            raise ValueError(
                "Found bin edge of size <= 0. Did you specify `bins` with"
                "non-monotonic sequence?")

    nbin = asarray(nbin)

    # Handle empty input.
    if N == 0:
        if W > 0:
            return [np.zeros(nbin-2) for _ in arange(W)], edges
        else:
            return np.zeros(nbin-2), edges

    # Compute the bin number each sample falls into.
    Ncount = {}
    for i in arange(D):
        # searchsorted is faster for many bins
        Ncount[i] = searchsorted(edges[i], sample[:, i], "right")
        #Ncount[i] = digitize(sample[:, i], edges[i])

    # Using digitize, values that fall on an edge are put in the right bin.
    # For the rightmost bin, we want values equal to the right
    # edge to be counted in the last bin, and not as an outlier.
    for i in arange(D):
        # Rounding precision
        mindiff = dedges[i].min()
        if not np.isinf(mindiff):
            decimal = int(-log10(mindiff)) + 6
            # Find which points are on the rightmost edge.
            not_smaller_than_edge = (sample[:, i] >= edges[i][-1])
            on_edge = (around(sample[:, i], decimal) == around(edges[i][-1], decimal))
            # Shift these points one bin to the left.
            Ncount[i][where(on_edge & not_smaller_than_edge)[0]] -= 1

    # Compute the sample indices in the flattened histogram matrix.
    ni = nbin.argsort()
    xy = zeros(N, int)
    for i in arange(0, D-1):
        xy += Ncount[ni[i]] * nbin[ni[i+1:]].prod()
    xy += Ncount[ni[-1]]

    # Compute the number of repetitions in xy and assign it to the
    # flattened histmat.
    if len(xy) == 0:
        if W > 0:
            return [np.zeros(nbin-2) for _ in arange(W)], edges
        else:
            return zeros(nbin-2, int), edges

    # Flattened histogram matrix (1D)
    # Reshape is used so that overlarge arrays
    # will raise an error.
    Wd = W if W > 0 else 1
    hists = [zeros(nbin, float).reshape(-1) for _ in arange(Wd)]
    for histidx, hist in enumerate(hists):
        weights_ = weights[histidx] if W > 0 else weights
        flatcount = bincount(xy, weights_)
        a = arange(len(flatcount))
        hist[a] = flatcount
    
        # Shape into a proper matrix
        hist = hist.reshape(sort(nbin))
        ni = nbin.argsort()
        for i in arange(nbin.size):
            j = ni.argsort()[i]
            hist = hist.swapaxes(i, j)
            ni[i], ni[j] = ni[j], ni[i]
    
        # Remove outliers (indices 0 and -1 for each dimension).
        core = D*[slice(1, -1)]
        hist = hist[core]
    
        # Normalize if normed is True
        if normed:
            s = hist.sum()
            for i in arange(D):
                shape = ones(D, int)
                shape[i] = nbin[i] - 2
                hist = hist / dedges[i].reshape(shape)
            hist /= s
    
        if (hist.shape != nbin - 2).any():
            raise RuntimeError(
                "Internal Shape Error: hist.shape != nbin-2 -> " + str(hist.shape) + " != " + str(nbin-2))
        
        hists[histidx] = hist
    
    if W in [None, -1]:
        return hists[0], edges
    else:
        return hists, edges

def histogram2d(x, y, bins=10, range=None, normed=False, weights=None):
    """
    Compute the bi-dimensional histogram of two data samples.

    Parameters
    ----------
    x : array_like, shape (N,)
        An array containing the x coordinates of the points to be
        histogrammed.
    y : array_like, shape (N,)
        An array containing the y coordinates of the points to be
        histogrammed.
    bins : int or [int, int] or array_like or [array, array], optional
        The bin specification:

          * If int, the number of bins for the two dimensions (nx=ny=bins).
          * If [int, int], the number of bins in each dimension
            (nx, ny = bins).
          * If array_like, the bin edges for the two dimensions
            (x_edges=y_edges=bins).
          * If [array, array], the bin edges in each dimension
            (x_edges, y_edges = bins).

    range : array_like, shape(2,2), optional
        The leftmost and rightmost edges of the bins along each dimension
        (if not specified explicitly in the `bins` parameters):
        ``[[xmin, xmax], [ymin, ymax]]``. All values outside of this range
        will be considered outliers and not tallied in the histogram.
    normed : bool, optional
        If False, returns the number of samples in each bin. If True,
        returns the bin density ``bin_count / sample_count / bin_area``.
    weights : array_like, shape(N,), optional
        An array of values ``w_i`` weighing each sample ``(x_i, y_i)``.
        Weights are normalized to 1 if `normed` is True. If `normed` is
        False, the values of the returned histogram are equal to the sum of
        the weights belonging to the samples falling into each bin.
        Weights can also be a list of (weight arrays or None), in which case
        a list of histograms is returned as H.

    Returns
    -------
    H : ndarray, shape(nx, ny)
        The bi-dimensional histogram of samples `x` and `y`. Values in `x`
        are histogrammed along the first dimension and values in `y` are
        histogrammed along the second dimension.
    xedges : ndarray, shape(nx,)
        The bin edges along the first dimension.
    yedges : ndarray, shape(ny,)
        The bin edges along the second dimension.

    See Also
    --------
    histogram : 1D histogram
    histogramdd : Multidimensional histogram

    Notes
    -----
    When `normed` is True, then the returned histogram is the sample
    density, defined such that the sum over bins of the product
    ``bin_value * bin_area`` is 1.

    Please note that the histogram does not follow the Cartesian convention
    where `x` values are on the abscissa and `y` values on the ordinate
    axis.  Rather, `x` is histogrammed along the first dimension of the
    array (vertical), and `y` along the second dimension of the array
    (horizontal).  This ensures compatibility with `histogramdd`.

    Examples
    --------
    >>> import matplotlib as mpl
    >>> import matplotlib.pyplot as plt

    Construct a 2D-histogram with variable bin width. First define the bin
    edges:

    >>> xedges = [0, 1, 1.5, 3, 5]
    >>> yedges = [0, 2, 3, 4, 6]

    Next we create a histogram H with random bin content:

    >>> x = np.random.normal(3, 1, 100)
    >>> y = np.random.normal(1, 1, 100)
    >>> H, xedges, yedges = np.histogram2d(y, x, bins=(xedges, yedges))

    Or we fill the histogram H with a determined bin content:

    >>> H = np.ones((4, 4)).cumsum().reshape(4, 4)
    >>> print H[::-1]  # This shows the bin content in the order as plotted
    [[ 13.  14.  15.  16.]
     [  9.  10.  11.  12.]
     [  5.   6.   7.   8.]
     [  1.   2.   3.   4.]]

    Imshow can only do an equidistant representation of bins:

    >>> fig = plt.figure(figsize=(7, 3))
    >>> ax = fig.add_subplot(131)
    >>> ax.set_title('imshow: equidistant')
    >>> im = plt.imshow(H, interpolation='nearest', origin='low',
                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

    pcolormesh can display exact bin edges:

    >>> ax = fig.add_subplot(132)
    >>> ax.set_title('pcolormesh: exact bin edges')
    >>> X, Y = np.meshgrid(xedges, yedges)
    >>> ax.pcolormesh(X, Y, H)
    >>> ax.set_aspect('equal')

    NonUniformImage displays exact bin edges with interpolation:

    >>> ax = fig.add_subplot(133)
    >>> ax.set_title('NonUniformImage: interpolated')
    >>> im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
    >>> xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
    >>> ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
    >>> im.set_data(xcenters, ycenters, H)
    >>> ax.images.append(im)
    >>> ax.set_xlim(xedges[0], xedges[-1])
    >>> ax.set_ylim(yedges[0], yedges[-1])
    >>> ax.set_aspect('equal')
    >>> plt.show()

    """
    try:
        N = len(bins)
    except TypeError:
        N = 1

    if N != 1 and N != 2:
        xedges = yedges = asarray(bins, float)
        bins = [xedges, yedges]
    hist, edges = histogramdd([x, y], bins, range, normed, weights)
    return hist, edges[0], edges[1]