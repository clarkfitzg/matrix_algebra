'''
trix.py

Utility functions for working with matrices in Numpy

Clark Fitzgerald
'''


from __future__ import division
import math
import itertools
import functools
import numpy as np


# TODO Change repr
def replicate(n, func, *args, **kwargs):
    '''
    Returns an iterator calling func n times.
    Similar to replicate in the R language.

    >>> list(replicate(3, pow, 2, 2))
    [4, 4, 4]

    '''
    pfunc = functools.partial(func, *args, **kwargs)
    return itertools.islice(iter(pfunc, 'NoSentinel'), n)


def matrixij(shape, ij):
    """Returns an array with 1 in position i, j, and zeros elsewhere

    >>> matrixij((2, 3), (0, 1))
    array([[0, 1, 0],
           [0, 0, 0]])

    """
    m = np.zeros(shape, dtype=int)
    m[ij] = 1
    return m


def matrixbasis(n, m):
    """Generates a basis for a vector space of n x m matrices

    >>> list(matrixbasis(1, 2))
    [matrix([[1, 0]]), matrix([[0, 1]])]

    """
    for i in range(n):
        for j in range(m):
            # Using matrix over array since matrix is hashable
            yield np.matrix(matrixij((n, m), (i, j)))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
