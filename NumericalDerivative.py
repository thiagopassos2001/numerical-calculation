# Read more in: https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods

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

def MidpointMethod(df,x0,y0,a,b,n=1e6):

    '''Midpoint method is a one-step method for numerically solving the differential equation
    Reference: https://en.wikipedia.org/wiki/Midpoint_method 
    Required inputs:
    'df' is first-order ordinary differential equation: df = g(x)
    'x0' and 'y0' are initial conditions
    'a' is left edge of interval
    'b' is right edge of interval
    'n=1e6' is a number of subdivisions in the interval
    '''

    h = (b - a) / n

    x = []
    for i in range(int(n) + 1):
        xk = x0 + (i*h)
        x.append(xk)

    y = [y0]
    for j in range(1,int(n) + 1):

        k1 = df(x[j-1],y[j-1])
        k2 = df(x[j-1] + (h/2),y[j-1] + (h*(k1/2)))
        yk = y[j-1] + (k2*h)

        y.append(yk)         
    
    return x, y

def HeunMethod(df,x0,y0,a,b,n=1e6,i_max=0,e_rel=1e-6):

    '''Heun method is a numerical procedure for solving ordinary differential equations (ODEs) with a given initial value
    Reference: https://en.wikipedia.org/wiki/Heun%27s_method
    Required inputs:
    'df' is first-order ordinary differential equation: df = g(x)
    'x0' and 'y0' are initial conditions
    'a' is left edge of interval
    'b' is right edge of interval
    'n=1e6' is a number of subdivisions in the interval
    '''

    h = (b - a) / n

    x = []
    for i in range(int(n) + 1):
        xk = x0 + (i*h)
        x.append(xk)

    y = [y0]
    for j in range(1,int(n) + 1):

        k1 = df(x[j-1],y[j-1])
        k2 = df(x[j-1] + h,y[j-1] + (h*k1))
        yk = y[j-1] + (((k1/2) + (k2/2))*h)

        i = i_max
        yk_1 = yk
        erk = 1

        while (i > 0) and (erk > e_rel):
            k2 = df(x[j-1] + h,yk_1)
            yk = y[j-1] + (((k1/2) + (k2/2))*h)
            
            erk = abs((yk - yk_1)/yk)
            yk_1 = yk
            i = i - 1

        y.append(yk)

    return x, y

def RalstonMethod(df,x0,y0,a,b,n=1e6):
    
    '''
    Reference: https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
    Required inputs:
    'df' is first-order ordinary differential equation: df = g(x)
    'x0' and 'y0' are initial conditions
    'a' is left edge of interval
    'b' is right edge of interval
    'n=1e6' is a number of subdivisions in the interval
    '''

    h = (b - a) / n

    x = []
    for i in range(int(n) + 1):
        xk = x0 + (i*h)
        x.append(xk)

    y = [y0]
    for j in range(1,int(n) + 1):

        k1 = df(x[j-1],y[j-1])
        k2 = df(x[j-1] + ((3/4)*h),y[j-1] + ((3*h*k1)/4))
        yk = y[j-1] + ((k1 + (2*k2))*(h/3))

        y.append(yk)         
    
    return x, y

def RungeKutta4Method(df,x0,y0,a,b,n=1e6):

    '''Runge–Kutta fourth-order method using the general formula with s=4 evaluated, as explained above, at the starting point, the midpoint and the end point of any interval
    Reference: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Derivation_of_the_Runge%E2%80%93Kutta_fourth-order_method
    Required inputs:
    'df' is first-order ordinary differential equation: df = g(x)
    'x0' and 'y0' are initial conditions
    'a' is left edge of interval
    'b' is right edge of interval
    'n=1e6' is a number of subdivisions in the interval
    '''
    
    h = (b - a) / n

    x = []
    for i in range(int(n) + 1):
        xk = x0 + (i*h)
        x.append(xk)

    y = [y0]
    for j in range(1,int(n) + 1):

        k1 = df(x[j-1],y[j-1])
        k2 = df(x[j-1] + ((1/2)*h),y[j-1] + (k1*((1/2)*h)))
        k3 = df(x[j-1] + ((1/2)*h),y[j-1] + (k2*((1/2)*h)))
        k4 = df(x[j-1] + h, y[j-1] + (k3*h))
        yk = y[j-1] + ((k1 + (2*k2) + (2*k3) + k4)*(h/6))

        y.append(yk)         
    
    return x, y

def FiniteDifference1Method(p,q,g,x0,y0,x1,y1,a,b,n):
    '''
    p, q ang g functions are the following format y" + p(x)y' + q(x)y = g(x)
    'x0' and 'y0' are the initial solution for y(x) - Boundary Value Problem
    'x1' and 'y1' are the initial solution for y(x) - Boundary Value Problem
    'a' and 'b' are the limits of method (a=x0, b=x1)
    'n' is a number of subdivisions in the interval
    '''
    h = (b - a)/n
    i_max = int(n) + 1
    
    x = []
    for i in range(i_max):
        xk = x0 + (i*h)
        x.append(xk)

    # Boundary starting point
    y = [y0]

    # Inside point
    A, b = TridiagonalMatrix(p,q,g,x0,y0,x1,y1,x,h)
    A, b = GaussElimination(A,b)
    y = y + SolveUpperTriangular(A,b)

    # Boundary end point
    y.append(y1)

    return x, y
