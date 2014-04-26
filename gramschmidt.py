'''
Perform gram-schmidt orthonormalization on a matrix
'''

import numpy as np
import numba
from numpy.linalg import norm


def gs_normal(m):
    '''
    Use gram-schmidt for orthonormalization
    Return a matrix whose columns are orthonormal
    '''
    n = m.shape[1]  # Number of column vectors
    for i in range(n):
        # Normalize the ith column
        m[:, i] = m[:, i] / norm(m[:, i])
        for j in range(i + 1, n):
            # Make the jth column orthogonal to ith
            j_on_i = np.dot(m[:, i], m[:, j]) * m[:, i]
            m[:, j] = m[:, j] - j_on_i 
    return m
    

# Can't seem to figure out how to get the performance out of Numba.
# It doesn't want to deal with np.sum. Does that mean I need to 
# write a for loop every time I take a sum? Yuck.

#@numba.jit(nopython=True)
@numba.jit('f8[:,:](f8[:,:])')  # Typed
def gs_fast(m):
    '''
    Use gram-schmidt for orthonormalization
    Return a matrix whose columns are orthonormal
    '''
    n = m.shape[1]  # Number of column vectors
    for i in range(n):
        # Normalize the ith column
        m_i_2 = 0
        for l in range(m.shape[0]):
            m_i_2 += m[l, i] ** 2
        m[:, i] = m[:, i] / np.sqrt(m_i_2)
        for j in range(i + 1, n):
            # Make the jth column orthogonal to ith
            inner_prod_ij = 0
            for k in range(m.shape[0]):
                inner_prod_ij += m[k, i] * m[k, j]
            j_on_i = inner_prod_ij * m[:, i]
            m[:, j] = m[:, j] - j_on_i 
    return m
    
############################################################

if __name__ == '__main__':
    
    m = 10 * np.random.randn(5, 3)
    
    print 'Original matrix:\n', m
    
    m = gs_normal(m)
    
    print '\n Normalized matrix:\n', m
    
    print '\n Dot products near 0 indicate orthogonality:'
    print m[:,0].dot(m[:,2]), m[:,1].dot(m[:,0]), m[:,1].dot(m[:,2])
    
    print '\n Lengths of the vectors should be near 1:'
    print norm(m, axis=0)
    
    # For speed testing
    # %timeit gs_normal(m10)
    m10 = np.random.randn(50, 10)
    # base- 2.1 ms      autojit- 2.4ms     typed jit- 2.2 ms
    m100 = np.random.randn(500, 100)
    # base - 237 ms     autojit- 263 ms     typed jit- 232ms
    
    # This line will compile the numba jit call
    #m2 = gs_fast(m10)
