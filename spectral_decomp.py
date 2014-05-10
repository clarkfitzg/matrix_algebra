'''
Exploring the spectral decomposition for matrices
'''

import numpy as np
import toolz


# Cache expensive eigen calculations
eigh = toolz.memoize(np.linalg.eigh)


class sym_matrix(np.ndarray):
    '''
    Class for real symmetric matrices
    '''

    def __new__(cls, input_array):
        m = np.asarray(input_array)
        if (m != m.T).any():
            raise ValueError('Matrix must be symmetric')
        return m.view(cls)


    @staticmethod
    def rand(n=3, maxint=10):
        '''
        Generate a random symmetric n x n matrix
        consisting of integers from 0 to maxint.
        '''
        m = np.random.randint(maxint // 2, size=(n, n))
        m = m + m.T     # Now m is symmetric
        return sym_matrix(m)


    # Why is this returning an instance of the parent class sym_matrix?
    @property
    def eigvec(self):
        '''
        Eigenvectors associated with matrix
        '''
        return eigh(self)[1]


    @property
    def eigval(self):
        '''
        Eigenvalues associated with matrix
        '''
        return eigh(self)[0]


a = sym_matrix.rand()

print a

#print('\n The matrices sum to the identity:')

projectors = [np.outer(vector, vector) for vector in a.eigvec]



############################################################
#
#t1 = np.random.randn(2, 3)
#
#try:
#    sym_matrix(t1)
#    print('fail non-square input')
#except ValueError:
#    print('pass non-square input')
#
#t2 = np.random.randn(3, 3, 3)
#
#
#try:
#    sym_matrix(t2)
#    print('fail high dimension input')
#except ValueError:
#    print('pass high dimension input')
