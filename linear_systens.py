from numpy.linalg import inv,norm
from numpy import array, zeros, diag, diagflat, dot, tril, triu, append, resize, identity

def breakcase(a) -> float:
    s = 0
    for i in range(len(a)):
        s += a[i][0] ** 2
    return s ** 0.5

def jacobi(A,b,N=25,x=None,epsi=0.01):
    F = []
    d = []

    for i in range(len(A)) :
        F.append(A[i]/A[i][i])
        d.append(b[i]/A[i][i])

    F = array(F)
    d = array(d).reshape(-1,1)

    x = zeros(d.shape)
    n = 0

    for _ in range(N):
        x = (identity(F.shape[0])-F).dot(x) + d
        if norm(A.dot(x)-b, 'fro') < epsi:
            break
        n+=1
    return x, n

def gauss_seidel(a:array,b:array,x=None,N=25,epsi=0.01):
    F = []
    d = []
    count = 0

    for i in range(len(A)) :
        F.append(A[i]/A[i][i])
        d.append(b[i]/A[i][i])

    F = array(F)
    d = array(d).reshape(-1,1)
    if x is None:
        x = zeros(d.shape) 
    n = 0

    U = triu(F,1) # Upper Triangle
    L = tril(F,-1) # Lower Triangle
    D = F-U-L # Main Diagonal
    LDI = inv(L + D)
    prev = x

    for _ in range(N):
        x = LDI.dot(d) - LDI.dot(U).dot(x)
        if breakcase(A.dot(x)-b) < epsi:
            break
        prev = x
        count+=1
    return x, count

A = array([[2.0,1.0],[5.0,7.0]])
b = array([11.0,13.0])
guess = array([1.0,1.0])
epsi = 0.001

sol,count = jacobi(A,b,x=guess,N=200,epsi=epsi)
print(f"Solução Jacobi: {sol} em {count} interações")

sol,count = gauss_seidel(A,b,N=200,epsi=epsi)
print(f"Solução Gauss-Seidel: {sol} em {count} interações")
