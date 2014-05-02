# P.75 of the text reads:
# From equations (3.78) and (3.79) we see that the matrices
# A and B are orthogonal to each other if and only if
# t(A)B and t(B)A are hollow (that is, they have 0s in all
# diagonal positions).

A <- matrix(c(1, 1, 1, 1), nrow=2)
B <- matrix(c(1, 0, -1, 0), nrow=2)
print('A = ')
print(A)
print('B = ')
print(B)

AdotB <- sum(A*B)
print(sprintf('A and B are orthogonal since A dot B = %.4f', AdotB))

tAB <- t(A) %*% B 
print('However, t(A)B is not hollow:')
print(tAB)
