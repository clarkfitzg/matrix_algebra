'''
Exploring properties of the Kronecker product
'''

import numpy as np

a = np.random.randint(6, size=(3, 4))
b = np.random.randint(6, size=(2, 3))


def mprint(a):
    '''
    Prints a matrix 

    Input string with variable name
    '''
    print('\n {} = '.format(a))
    print(eval(a, globals()))


mprint('a')
mprint('b')


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


# For inverses need square full rank matrices

a2 = np.array([[1, 0], [1, 1]])
b2 = np.array([[2, 3], [1, 7]])
