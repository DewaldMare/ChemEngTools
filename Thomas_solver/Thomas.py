def thomas_solver(a, b, c, r):
    import numpy as np
    from numpy.linalg import norm
    '''Implementing th eThomas algorithm to solve a tridiagonal matrix with n unknows. a, b, c and d are the matrix entries:
    |b1 c1 0  0  0 | |x1|   |d1|
    |a2 b2 c2 0  0 | |x2|   |d2|
    |0  a3 b3 c3 0 | |x3| = |d3|
    |0  0  a4 b4 cn| |x4|   |d4|
    |0  0  0  an bn| |x5|   |d5| '''

    # Size of matrix
    n = len(b)
    #a has to have size n-1:
    if len(a) != n - 1:
        print('a should be of size n - 1. Input has length', len(a))
    #c has to have size n-1:
    if len(c) != n - 1:
        print('c should be of size n - 1. Input has length', len(c))
    #Stability check
    barr = np.array(b)
    aarr = np.array(a)
    carr = np.array(c)
    if norm(barr) <= norm(aarr) + norm(carr):
        print('The Thomas algorithm is not stable. For stability please ensure that the following condition is met: \n ||{b_{i}|| > ||a_{i}|| + ||c_{i}||')
   
    # Convert inputs to floats
    for i in range(len(a)):
        a[i] = float(a[i])
    for j in range(len(b)):
        b[j] = float(b[j])
    for k in range(len(c)):
        c[k] = float(c[k])
    for l in range(len(r)):
        r[l] = float(r[l])
    c.append(0.0)
    
    # Calculate gamma and rho
    gam = np.zeros(len(b))
    rho = np.zeros(len(b))
    gam[0] = c[0] / b[0]
    rho[0] = r[0] / b[0]
    for m in range(1, n):
        gam[m] = (c[m] / (b[m] - a[m - 1] * gam[m - 1]))
        rho[m] = (r[m] - a[m - 1] * rho[m - 1]) / (b[m] - a[m - 1] * gam[m - 1])

    #Solve through back substitution
    x = np.zeros(len(b))
    x[-1] = rho[-1]
    for k in range(len(b) - 1, -1, -1):
        x[k - 1] = rho[k - 1] - gam[k - 1] * x[k]
        
    return x