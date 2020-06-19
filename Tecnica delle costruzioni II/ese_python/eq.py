###ESERCIZIO TDSII ###
import numpy as np
import sympy as sy
#from sympy import symbols ,solve ,Eq, diff,pprint
import matplotlib.pyplot as plt
import pylab

## SIMBOLICO ##
e, p, A, B,t, E= sy.symbols ('e p A B t E') # costanti d'integrazione


## Soluzione dell'equazione differenziale della trave con carico distribuito ##
Equa=E*(e-p)-A-B*p**6

## condizioni al contorno ##
cc1=sy.Eq(Equa,0)

## RISOLVO IL SISTEMA
sol=sy.solve(cc1,p)
print(sol)
#v=v.subs(c1,sol[c1]).subs(c2,sol[c2]).subs(c3,sol[c3]).subs(c4,sol[c4])