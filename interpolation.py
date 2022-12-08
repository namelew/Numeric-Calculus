from numpy import array,math,linalg
from functools import reduce

def finite_difference(x:float,xi:array,yi:array) -> float:
    h = xi[1]-xi[0]
    z = lambda x : (x - xi[0])/h
    fd = lambda n,i : yi[i] if n==0 else fd(n-1,i+1) - fd(n-1,i)
    return yi[0]+sum((fd(i,0)/math.factorial(i))*reduce(lambda a,b : a*b, [z(x)-j for j in range(i)], 1) for i in range(1,len(xi)))

def divided_difference(x:float,xi:array,yi:array) -> float:
    dd = lambda n,i : yi[i] if n==0 else (dd(n-1,i+1) - dd(n-1,i))/(xi[i+n]-xi[i])
    return yi[0]+sum(dd(i,0)*reduce(lambda a,b : a*b, [x-xi[j] for j in range(i)], 1) for i in range(1,len(xi)))

def lagrange(x:float,xi:array,yi:array) -> float:
    return sum(yi[i]*reduce(lambda a,b : a*b, [(x-xi[j])/(xi[i]-xi[j]) for j in range(len(xi)) if j!=i],1) for i in range(len(xi)))

def main():
    yi = array([1,2,3])
    xi = array([0,2,5])
    x = 7

    print(lagrange(x,xi,yi))

main()
