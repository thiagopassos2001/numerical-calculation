def EulerMethod(df,x0,y0,a,b,n=1e6):

    '''Euler method is a first-order numerical procedure for solving ordinary differential equations with a given initial value
    Reference: https://en.wikipedia.org/wiki/Euler_method 
    Required inputs:
    'df' is first-order ordinary differential equation: df = g(x)
    'x0' and 'y0' are initial conditions
    'a' is left edge of interval
    'b' is right edge of interval
    'n=1e6' is a number of subdivisions in the interval
    '''

    h = (b - a) / n
    i_max = int(n) + 1

    x = []
    for i in range(i_max):
        xk = x0 + (i*h)
        x.append(xk)
    
    y = [y0]
    for i in range(1,i_max):
        yk = y[i-1] + (h*df(x[i-1],y[i-1]))
        y.append(yk)

    return x, y
