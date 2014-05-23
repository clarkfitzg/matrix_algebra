'''
trix.py

Utility functions for working with matrices in Numpy

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
    """Returns a generator for a basis set for an vector space of n x m matrices

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


def mprint(a):
    '''
    Prints a matrix 

    Input string with variable name

    >>> m = np.matrix([[math.pi, 0, 23]])
    >>> mprint(m)

     m = matrix[[3.141, 0, 23]]
    
    '''
    print('\n {} = '.format(a))
    print(eval(a, globals()))



if __name__ == '__main__':
    import doctest
    doctest.testmod()
