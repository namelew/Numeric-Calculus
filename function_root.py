
# Equações e Derivadas

from tkinter import E


def equation(x:float) -> float:
    return (x**3) - 2

def fderivate(x:float) -> float:
    return 3*(x**2)

def sderivate(x:float) -> float:
    return 6 * x

# Métodos númericos

def bissectionMethod(eq, gap:list, eps:float) -> tuple:
    x = (gap[0] + gap[1])/2
    y = eq(x)
    i = 0
    while(abs(y) > eps):
        if(y * eq(gap[0]) > 0):
            gap[0] = x
        else:
            gap[1] = x
        x = (gap[0] + gap[1])/2
        y = eq(x)
        i += 1
    return i,x,y

def newtonMethod(eq, fd, sd, gap:list, eps:float):
    recForm = True if fd(gap[0]) * sd(gap[0]) > 0 else False
    a = gap[0]
    b = gap[1]
    x = a - (eq(a)/fd(a)) if not recForm else b - (eq(b)/fd(b))
    y = eq(x)
    i = 0
    while(abs(y) > eps):
        print(i, x, y)
        x = x - (eq(x)/fd(x))
        y = eq(x)
        i += 1
    return i,x,y

def ropeMethod(eq, fd, sd, gap:list, eps:float):
    recForm = True if fd(gap[0]) * sd(gap[0]) > 0 else False
    a = gap[0]
    b = gap[1]
    x = ((a * eq(b)) - (b * eq(a)))/(eq(b) - eq(a))
    y = eq(x)
    i = 0
    while(abs(y) > eps and gap[0] >= a and gap[1] <= b):
        print(i, x, y)
        if recForm:
            x = ((x * eq(gap[1])) - (gap[1] * eq(x)))/(eq(gap[1]) - eq(x))
            gap[0] = x
        else:
            x = ((gap[0] * eq(x)) - (x * eq(gap[0])))/(eq(x) - eq(gap[0]))
            gap[1] = x
        y = eq(x)
        i += 1
    return i,x,y

def main():
    gap = [1,2]
    i,x,y = newtonMethod(equation, fderivate, sderivate, gap, 0.000001)
    print(i,x,y)

main()