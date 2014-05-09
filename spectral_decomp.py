'''
Exploring the spectral decomposition for matrices
'''

import numpy as np
from utils import lazy_property 


class sym_matrix(np.ndarray):
    '''
    Class for real symmetric matrices
    '''

    def __new__(cls, input_array):
        m = np.asarray(input_array)
        if m.shape[0] != m.shape[1] or len(m.shape) != 2:
            raise ValueError('Matrix must be square')
        return m.view(cls)


    @staticmethod
    def rand(n=3, maxint=10):
        '''
        Generate a random symmetric n x n matrix
        consisting of integers from 0 to maxint.
        '''
        m = np.random.randint(maxint // 2, size=(3, 3))
        m = m + m.T     # Now m is symmetric
        return sym_matrix(m)


    #@lazy_property
    def eigh(self):
        '''
        Compute eigenvalues and eigenvectors on an instance of sym_matrix
        '''
        print('Calculating Eigenvectors')
        #return np.linalg.eigh(self)
        self.eigval, self.eigvec = np.linalg.eigh(self)
        # Very strange- it's subclassing self.eigvec as sym_matrix


a = sym_matrix.rand()

print a

print('\n The matrices sum to the identity:')

############################################################

t1 = np.random.randn(2, 3)

try:
    sym_matrix(t1)
    print('fail non-square input')
except ValueError:
    print('pass non-square input')

t2 = np.random.randn(3, 3, 3)


try:
    sym_matrix(t2)
    print('fail high dimension input')
except ValueError:
    print('pass high dimension input')
