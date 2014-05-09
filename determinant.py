'''
Compute the determinant of a matrix
'''

import numpy as np

def drop_ij(m, i, j):
    '''
    Returns the submatrix formed by dropping row i and column j from m
    '''
    return np.delete(np.delete(m, i, axis=0), j, axis=1)


def det(m):
    '''
    Recursive definition -- Laplace expansion
    Returns determinant of matrix m
    '''
    if m.shape[0] != m.shape[1]:
        raise TypeError('Must input square matrix')

    # Stopping condition
    if m.shape == (1, 1):
        return m[0, 0]

    # Use Laplace expansion along first row
    ncol = m.shape[1]
    total = 0
    for j in range(ncol):
        total += m[0, j] * ((-1) ** (j)) * det(drop_ij(m, 0, j))
    
    return total


############################################################

m1 = np.matrix([[1, 2], [0, 3]])

if drop_ij(m1, 1, 0) == 2:
    print 'pass drop test'
else:
    print 'fail drop test'

if det(m1) == 3:
    print 'pass 2 x 2 test'
else:
    print 'fail 2 x 2 test'


m2 = np.matrix([1, 2, 3])

try:
    det(m2)
    print 'fail non-square test'
except TypeError:
    print 'pass non-square matrix input'

m10 = np.random.randn(10, 10)
