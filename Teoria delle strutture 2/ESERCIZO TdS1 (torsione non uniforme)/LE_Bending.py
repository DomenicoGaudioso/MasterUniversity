import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

## simbolico

c1, c2, c3, c4 = sy.symbols ('c1 c2 c3 c4') # costanti d'integrazione
p1=sy.symbols('p1') # carico distribuito in direzione 1
E,J2,H =sy.symbols ('E J2 H') # Modulo di Young & Momento d'inerzia asse 2
z=sy.symbols ('z') # ascissa curvilinea & altezza

## Soluzione dell'equazione differenziale della trave con carico distribuito

v1=(p1*z**4/24+c1*z**3/6+c2*z**2/2+c3*z+c4)*1/(E*J2) #spostamento tratto 1
r2=sy.diff(v1,z) #rotazione
M2=E*J2*sy.diff(v1,z,z) #Momento
T2=E*J2*sy.diff(v1,z,z,z) #Taglio

## condizioni al contorno

cc1=sy.Eq(v1.subs(z,0),0) #spostamento nullo in s=0
cc2=sy.Eq(r2.subs(z,0),0) #rotazione nulla in s=0
cc3=sy.Eq(M2.subs(z,H),0) #Momento in s=h
cc4=sy.Eq(T2.subs(z,H),0) #Taglio in s=h

## RISOLVO IL SISTEMA

sol=sy.solve((cc1,cc2,cc3,cc4),(c1,c2,c3,c4))
v1=v1.subs(c1,sol[c1]).subs(c2,sol[c2]).subs(c3,sol[c3]).subs(c4,sol[c4])
r2=r2.subs(c1,sol[c1]).subs(c2,sol[c2]).subs(c3,sol[c3]).subs(c4,sol[c4])
M2=M2.subs(c1,sol[c1]).subs(c2,sol[c2]).subs(c3,sol[c3]).subs(c4,sol[c4])
T1=T2.subs(c1,sol[c1]).subs(c2,sol[c2]).subs(c3,sol[c3]).subs(c4,sol[c4])

#print(sol)
print('v1=', sy.simplify(v1) )
print('r2=', sy.simplify(r2) )
print('M2=', sy.simplify(M2) )
print('T2=', sy.simplify(T1) )