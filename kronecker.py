'''
Some general operations with vectors and matrices
From chapters 2 and 3 of Matrix Algebra
'''

import numpy as np

a = np.random.randint(6, size=(3, 4))
b = np.random.randint(6, size=(3, 4))

print '\n a = \n', a
print '\n b = \n', b

def kronecker_prod(A, B):
    '''
    Kronecker multiplication for matrices A, B
    '''
    bigB = np.tile(B, A.shape)
    bigA = np.repeat(np.repeat(A, B.shape[0], axis=0), B.shape[1], axis=1)
    # Funny - below we use Hadamard product as well
    return bigA * bigB

axb = kronecker_prod(a, b)
print '\n Kronecker Product axb = \n', axb
