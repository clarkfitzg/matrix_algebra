'''
Exploring the spectral decomposition for matrices
'''

from itertools import imap
import numpy as np
import toolz
import utils


# Implement caching
eigh = toolz.memoize(np.linalg.eigh)


class sym_matrix(np.matrix):
    '''
    Class for real symmetric matrices
    '''

    # Rewrite this to __init__
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


    @property
    def eigvec(self):
        '''
        Eigenvectors associated with matrix
        '''
        return list(np.transpose(eigh(self)[1]))


    @property
    def eigval(self):
        '''
        Eigenvalues associated with matrix
        '''
        return list(eigh(self)[0])


    @property
    @toolz.memoize
    def projectors(self):
        '''
        Spectral projectors - 
        outer product of each eigenvector
        '''
        print '\n*** EVALUATING PROJECTOR FUNCTION ***\n'
        return [np.outer(v, v) for v in self.eigvec]


    @property
    @toolz.memoize
    def spec_decomp(self):
        '''
        Return matrix calculated as spectral decomposition
        '''
        mapper = imap(np.multiply, self.eigval, self.projectors)
        return sum(mapper)


a = sym_matrix.rand()

print a

print '\n The sum of the spectral projectors is:'
print sum(a.projectors)

print '\n Expressing a as the spectral decomposition:\n'
print a.spec_decomp

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
