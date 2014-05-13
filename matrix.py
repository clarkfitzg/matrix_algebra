'''
In this module we explore linear algebra using Python
following the text 'Matrix Algebra' by James Gentle.

Functions which have an equivalent in the Numpy stack are 
named Xfunc, where func is the name of the equivalent function. 
'''


import itertools
import operator
import numpy as np


def Xnorm(x):
    """L2 Norm for matrix or vector

    >>> a = np.array([3, 4])
    >>> Xnorm(a)
    5.0
    >>> np.linalg.norm(a)
    5.0

    """
    return np.sqrt(sum(x ** 2))


def Xdot(x, y):
    """Vector dot product

    >>> a = np.array([1, 2, 0])
    >>> b = np.array([3, 2, 9])
    >>> Xdot(a, b)
    7
    >>> np.dot(a, b)
    7

    """
    return sum(itertools.imap(operator.mul, x, y))


def Xnormalize(x):
    """Transform a vector so that it has L2 norm of 1

    >>> a = np.ones(4)
    >>> a_n = Xnormalize(a)
    >>> a_n
    array([ 0.5,  0.5,  0.5,  0.5])
    >>> np.linalg.norm(a_n)
    1.0
    >>> from sklearn.preprocessing import normalize
    >>> a.shape = (1, 4)  # Requires matrix
    >>> normalize(a)
    array([[ 0.5,  0.5,  0.5,  0.5]])

    """
    return x / np.linalg.norm(x)


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
    """Returns a basis set for n x m vector space of matrices

    >>> matrixbasis(1, 2)
    set(matrix([[ 1, 0]], matrix([[ 0, 1]]))

    """


if __name__ == '__main__':
    import doctest
    doctest.testmod()
