import matplotlib.pyplot as plt
import numpy as np

# Equações e Derivadas

def equation(x:float) -> float:
    return (2*(x**3)) - (x**2) - (2.718281828459045 ** x) + 2

def fderivate(x:float) -> float:
    return 6*(x**2) - (2 * x) - (2.718281828459045 ** x)

def sderivate(x:float) -> float:
    return (12 * x) - (2.718281828459045 ** x) - 2

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
    a = gap[0] if gap[0] != 0 else 0.000000000001
    b = gap[1] if gap[1] != 0 else 0.000000000001
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
    a = gap[0] if gap[0] != 0 else 0.000000000001 
    b = gap[1] if gap[1] != 0 else 0.000000000001 
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

# Funções avulças
def haveRoot(eq, a:float, b:float) -> bool:
    return True if eq(a) * eq(b) < 0 else False

def haveOnlyOneRoot(eq, a:int, b:int) -> bool:
    prev = 0
    current = 0
    for i in range(a, b+1):
        current += eq(i)
        if(abs(current) < abs(prev)):
            return False
        prev = current
    return True

def plotFunctionGraph():
    # 100 linearly spaced numbers
    x = np.linspace(-5,7,100)

    # the function, which is y = x^2 here
    y = (2*(x**3)) - (x**2) - (2.718281828459045 ** x) + 2

    # setting the axes at the centre
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # plot the function
    plt.plot(x, y, 'r')

    # show the plot
    plt.show()

def main():
    plotFunctionGraph()

    gap = [5.9, 6]

    i,x,y = bissectionMethod(equation, gap, 0.0000001)
    print(f"Método da Bisseção: ({i}, {x}, {y:.6f})")
    i,x,y = ropeMethod(equation, fderivate, sderivate, gap, 0.0000001)
    print(f"Método das Cordas: ({i}, {x}, {y:.6f})")
    i,x,y = newtonMethod(equation, fderivate, sderivate, gap, 0.0000001)
    print(f"Método das Newtom: ({i}, {x}, {y:.6f})")

main()