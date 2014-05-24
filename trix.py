'''
trix.py

Utility functions for working with matrices in Numpy

Clark Fitzgerald
'''


from __future__ import division
import math
import itertools
import numpy as np



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

    >>> basis = matrixbasis(1, 2)
    >>> basis.next()
    matrix([[1, 0]])
    >>> basis.next()
    matrix([[0, 1]])

    """
    for i in xrange(n):
        for j in xrange(m):
            # Converting to matrix since matrix is hashable
            yield np.matrix(matrixij((n, m), (i, j)))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
