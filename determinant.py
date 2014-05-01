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
    # Stopping condition
    if m.ndim == 0:
        return m

    # Use Laplace expansion along first row
    ncol = m.shape[1]
    total = 0
    for j in range(ncol):
        total += m[0, j] * ((-1) ** (2 + j)) * det(drop_ij(m, 0, j))
    
    return total




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
