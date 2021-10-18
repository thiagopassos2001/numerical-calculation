#Auxiliary Methods

def VerificationBolzanoWeierstrassTheorem(function, a, b):

    '''Check if exist a root in interval ([a,b])
    Return True or False
    Reference: https://en.wikipedia.org/wiki/Bolzano%E2%80%93Weierstrass_theorem

    Required inputs:
    'function' is analyzed function
    'a' is left edge of interval
    'b' is right edge of interval
    'xk' if a value value within interval'''

    if function(a) * function(b) <= 0:
        return True
    else:
        return False

def BinaryChoice(function, a, b, xk):

    '''Selects the interval ([a,xk] or [xk,b]) that contains the root
    Return a, xk or xk, b
    Reference: https://en.wikipedia.org/wiki/Binary_search_algorithm
    
    Required inputs:
    'function' is analyzed function
    'a' is left edge of interval
    'b' is right edge of interval
    'xk' if a value value within interval'''

    if VerificationBolzanoWeierstrassTheorem(function, a, xk):
        return a, xk
    else:
        return xk, b

#Interval Methods

def Bisection(function,a,b,i_max=1e6,e_abs=0.001,e_rel=0.001,logs=False):

    '''Calculates the numerical root by the bisection method
    Reference: https://en.wikipedia.org/wiki/Bisection_method

    Required inputs:
    'function' is analyzed function
    'a' is left edge of interval
    'b' is right edge of interval

    Optional keyword arguments:
    Halting problem. Set 'False' to don't use the stop criterion
    'i_max' is maximum number of iterations of the method
    'e_abs' is a absolute error between approximate root and 0
    'e_rel' is a relative error between current root approximation and previous root approximation
    
    Show data. Set True for for show data for each iteration
    'logs'
    '''

    i = i_max
    erk = 1
    xk_1 = a

    if VerificationBolzanoWeierstrassTheorem(function, a, b):
        while ((i > 0) and (erk > e_rel) and (abs(function(xk_1)) > e_abs)):
            
            #Bisection update
            xk = (a + b)/2
            a, b = BinaryChoice(function, a, b, xk)
            #

            i = i - 1
            erk = abs((xk_1 - xk) / xk)
            
            if logs:
                print({'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))})
            
            xk_1 = xk

        return {'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))}

    else:
        print('No root found in interval')

def FalsePosition(function,a,b,i_max=1e6,e_abs=0.001,e_rel=0.001,logs=False):

    '''Calculates the numerical root by the false position method
    Reference: https://en.wikipedia.org/wiki/Regula_falsi

    Required inputs:
    'function' is analyzed function
    'a' is left edge of interval
    'b' is right edge of interval

    Optional keyword arguments:
    Halting problem. Set 'False' to don't use the stop criterion
    'i_max' is maximum number of iterations of the method
    'e_abs' is a absolute error between approximate root and 0
    'e_rel' is a relative error between current root approximation and previous root approximation
    
    Show data. Set True for for show data for each iteration
    'logs'
    '''

    i = i_max
    erk = 1
    xk_1 = a

    if VerificationBolzanoWeierstrassTheorem(function, a, b):
        while ((i > 0) and (erk > e_rel) and (abs(function(xk_1)) > e_abs)):

            #FalsePosition update
            m_inclination = (function(b) - function(a)) / (b - a)
            y_intersection = function(b) - (m_inclination * b)
            xk = (-y_intersection)/m_inclination
            a, b = BinaryChoice(function, a, b, xk)
            #

            i = i - 1
            erk = abs((xk_1 - xk) / xk)

            if logs:
                print({'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))})

            xk_1 = xk

        return {'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))}

    else:
        print('No root found in interval')

#Open Methods

def FixedPoint(function,phi,x0,i_max=1e6,e_abs=0.001,e_rel=0.001,logs=False):

    '''Calculates the numerical root by the fixed point method
    Reference: https://en.wikipedia.org/wiki/Fixed-point_iteration

    Required inputs:
    'function' is analyzed function
    'phi' is phi function (check reference)
    'x0' is initial value

    Optional keyword arguments:
    Halting problem. Set 'False' to don't use the stop criterion
    'i_max' is maximum number of iterations of the method
    'e_abs' is a absolute error between approximate root and 0
    'e_rel' is a relative error between current root approximation and previous root approximation
    
    Show data. Set True for for show data for each iteration
    'logs'
    '''

    i = i_max
    erk = 1
    xk_1 = x0

    while ((i > 0) and (erk > e_rel) and (abs(function(xk_1)) > e_abs)):
        
        #FixedPoint update
        xk = phi(xk_1)
        #

        i = i - 1
        erk = abs((xk_1 - xk) / xk)

        if logs:
            print({'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))})

        xk_1 = xk

    return {'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))}

def NewtonRaphson(function,derivative_function,x0,i_max=1e6,e_abs=0.001,e_rel=0.001,logs=False):

    '''Calculates the numerical root by the Newton-Raphson method
    Reference: https://en.wikipedia.org/wiki/Newton's_method

    Required inputs:
    'function' is analyzed function
    'derivative_function' is derivative function (check reference)
    'x0' is initial value

    Optional keyword arguments:
    Halting problem. Set 'False' to don't use the stop criterion
    'i_max' is maximum number of iterations of the method
    'e_abs' is a absolute error between approximate root and 0
    'e_rel' is a relative error between current root approximation and previous root approximation
    
    Show data. Set True for for show data for each iteration
    'logs'
    '''
    
    i = i_max
    erk = 1
    xk_1 = x0

    while ((i > 0) and (erk > e_rel) and (abs(function(xk_1)) > e_abs)):

        #NewtonRaphson update
        xk = xk_1 - (function(xk_1)/derivative_function(xk_1))
        #

        i = i - 1
        erk = abs((xk_1 - xk) / xk)

        if logs:
            print({'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))})

        xk_1 = xk

    return {'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))}

def Secant(function,x0,d=0.00001,x_1=False,i_max=1e6,e_abs=0.001,e_rel=0.01,logs=False):

    '''Calculates the numerical root by the secant method
    Reference: https://en.wikipedia.org/wiki/Secant_method

    Required inputs:
    'function' is analyzed function
    'x0' is initial value
    'd' is infinitesimal value 
    'x_1' is second value. If False, x_1 = x0 + d

    Optional keyword arguments:
    Halting problem. Set 'False' to don't use the stop criterion
    'i_max' is maximum number of iterations of the method
    'e_abs' is a absolute error between approximate root and 0
    'e_rel' is a relative error between current root approximation and previous root approximation
    
    Show data. Set True for for show data for each iteration
    'logs'
    '''
    
    i = i_max
    erk = 1
    xk_1 = x0
 
    if x_1 == False:
        xk_2 = xk_1 + d
    else:
        xk_2 = x_1

    while ((i > 0) and (erk > e_rel) and (abs(function(xk_1)) > e_abs)):

        # Secant update
        xk = xk_1 - (function(xk_1)*(xk_2 - xk_1)) / (function(xk_2) - function(xk_1))
        #

        i = i - 1
        erk = abs((xk_1 - xk) / xk)

        if logs:
            print({'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))})

        xk_2 = xk_1
        xk_1 = xk

    return {'root':xk, 'iterations':i_max - i, 'e_rel':erk, 'e_abs':abs(function(xk))}
