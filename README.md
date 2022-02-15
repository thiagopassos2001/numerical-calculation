# Cálculo Numérico (Numerical Calculation)

**Descrição**: Contém funções para a resolução numérica de problemas de álgebra e cálculo, aplicando métodos tradicionais em matemática. Desenvolvido durante a disciplina de  Métodos Numérico para Engenharia Civil (Contains functions for numerically solving of algebra and calculus problems, applying traditional methods in mathematics. Developed during the course of Numerical Methods for Civil Engineering).

## Conteúdo:

### Raízes Numéricas - NumericalRoot.py

* BinaryChoice

Example:
```python
function = lambda x: x + 2
a = 0
b = -5
xk = (a + b)/2

a, b = BinaryChoice(function, a, b, xk)

print(a,b)

OUTPUT: -2.5, 0
```

* VerificationBolzanoWeierstrassTheorem

Example:
```python
function = lambda x: x + 2
a = 0
b = 5

if VerificationBolzanoWeierstrassTheorem(function, a, b):
    print('yes, the root is in the range')
else:
    print('no, the root is not in the range')
    
OUTPUT: 'no, the root is not in the range'
```

* Bisection

Example:
```python
f = lambda x:(x**3) + (2*x) - 4
Bisection(f,1,5,i_max=False,e_abs=False,e_rel=0.01)

OUTPUT: 
{'e_abs': 0.0011019706726074219,
 'e_rel': 0.006622516556291391,
 'iterations': 9,
 'root': 1.1796875}
```

* False Position

Example:
```python
f = lambda x:(x**3) + (2*x) - 4
FalsePosition(f,1,5,i_max=100,e_abs=False,e_rel=0.01)

OUTPUT:
{'e_abs': 0.28503459703327927,
 'e_rel': 0.00897367452059495,
 'iterations': 7,
 'root': 1.1320673171230085}
```

* FixedPoint
* NewtonRaphson
* Secant

### Solução EDOs de 1ª Ordem - NumericalDerivative.py

* EulerMethod
* MidpointMethod
* HeunMethod
* RalstonMethod
* RungeKutta4Method
* FiniteDifference1Method
