import numpy as np
import sympy as sp

def W(n):
    z = np.zeros((n+1, n+1), dtype=object)
    for i in range(n+1):
        for j in range(0, i+1):
            if (i-j) % 2 == 0:
                c = 1 if j == 0 else 2
                z[i, j] = sp.Rational(c, 2**i)*sp.binomial(i, (i-j)//2)

    return z

def Wi(n):
    z = np.zeros((n+1, n+1), dtype=object)
    for i in range(n+1):
        if i == 0:
            z[0, 0] = 1
        else:
            for j in range(0, i+1):
                if (i-j) % 2 == 0:
                    z[i, j] = sp.Rational((-1)**(((i-j)//2)%2), i+j)*sp.binomial((i+j)//2, (i-j)//2)*2**j
            z[i] *= i
    return z

def A(n, a, b):
    z = np.zeros((n+1, n+1), dtype=object)
    for i in range(n+1):
        for j in range(0, i+1):
            z[i, j] = sp.binomial(i, i-j)*a**j*b**(i-j)
    return z

def BM(n):
    W0 = Wi(n-1)
    W1 = W(n-1)
    A0 = A(n-1, sp.S.Half, -sp.S.Half)
    A1 = A(n-1, sp.S.Half, sp.S.Half)
    return np.vstack((W0 @ A0 @ W1, W0 @ A1 @ W1)).reshape((2, n, n))

def CMtril(n):
    C = BM(n)
    D = np.zeros(n*n+n, dtype=object)
    k = 0
    for i in range(n):
        for j in range(i+1):
            D[k] = C[0, i, j]
            D[k+(n*n+n)//2] = C[1, i, j]
            k += 1
    return D

def CMtriu(n):
    C = BM(n)
    C0T = C[0].transpose()
    C1T = C[1].transpose()
    D = np.zeros(n*n+n, dtype=object)
    k = 0
    for i in range(n):
        for j in range(i, n):
            D[k] = C0T[i, j]
            D[k+(n*n+n)//2] = C1T[i, j]
            k += 1
    return D