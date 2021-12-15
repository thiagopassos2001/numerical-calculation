# Read more in: https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods

#Auxiliary Methods
def TridiagonalMatrix(p,q,g,x0,y0,x1,y1,x_list,h):
    '''
    p, q ang g functions are the following format y" + p(x)y' + q(x)y = g(x)
    'x0' and 'y0' are the initial solution for y(x) - Boundary Value Problem
    'x1' and 'y1' are the initial solution for y(x) - Boundary Value Problem
    'a' and 'b' are the limits of method (a=x0, b=x1)
    'list_x' is the list of discretized x values
    'h' is step close to 0
    '''

    n = len(x_list[1:-1])
    A = []
    b = []

    Awk = [(1 - (h*p(i)*0.5)) for i in x_list[1:-1]]
    Aek = [(1 + (h*p(i)*0.5)) for i in x_list[1:-1]]
    Apk = [(((h**2)*q(i)) - 2) for i in x_list[1:-1]]
    bk = [((h**2)*g(i)) for i in x_list[1:-1]]

    for i in range(0,n):
        A_line = []
        for j in range(0,n):
            if (i == j):
                A_line.append(Apk[i])
            elif ((j-1) == i):
                A_line.append(Aek[i])
            elif ((i-1) == j):
                A_line.append(Awk[i])
            else:
                A_line.append(0)
        
        A.append(A_line)
        b.append(bk[i])
    
    b[0] = b[0] - (1 - (h*p(x0)*0.5))*y0
    b[-1] = b[-1] - (1 + (h*p(x1)*0.5))*y1

    return A, b

def GaussElimination(A,b):
    '''
    A is coefficients
    b is result of system
    '''
    n = len(b)
    x = [0]*n

    for k in range(1,n,1):
        
        for i in range(k+1,n+1,1):
            Mult = A[i-1][k-1]/A[k-1][k-1]
            A[i-1][k-1] = 0
            b[i-1] = b[i-1] - Mult*b[k-1]
            
            for j in range(k+1,n+1,1):
                A[i-1][j-1] = A[i-1][j-1] - Mult*A[k-1][j-1] 

    return A,b

def SolveUpperTriangular(U,b):
    '''
    U is upper triangular matrix
    b is result of system
    '''
    n = len(b)
    x = [0]*n
    x[-1] = b[-1] / U[-1][n-1]

    for i in range(n-1, 0,-1):
        sum_temporary = 0
        for j in range(i + 1,n + 1):
            sum_temporary = sum_temporary + (U[i-1][j-1] * x[j-1])
        
        x[i-1] = (b[i-1] - sum_temporary) / U[i-1][i-1]

    return x

# Main functions
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
