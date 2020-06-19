import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

## simbolico

A, B, C, D = sy.symbols ('A B C D') # costanti d'integrazione

w,m=sy.symbols ('w m') # rotazione e momento distribuito

# Modulo di Young, Momento d'inerzia settoriale ,constante di torsione
E,Jw=sy.symbols ('E Jw')
z,H=sy.symbols ('z H ') # ascissa curvilinea & altezza
a = sy.Symbol('a', real=True) # rigidezza torsionale

## Soluzione dell'equazione differenziale della trave con momento torcente distribuito

w=A*sy.cosh(a*z)+B*sy.sinh(a*z)+C*z+D-m*z**2/(2*E*Jw*a**2) #rotazione
dw=sy.diff(w,z) #ingobbamento
Mt=a**2*E*Jw*dw-E*Jw*sy.diff(w,z,z,z) #Momento torcente
Bw=E*Jw*sy.diff(w,z,z,z) #Bimomento

## condizioni al contorno

cc1=sy.Eq(w.subs(z,0),0) #rotazione nulla in z=0
cc2=sy.Eq(dw.subs(z,0),0) #ingobbamento nullo in z=0
cc3=sy.Eq(Mt.subs(z,H),0) #Momento torcente nullo in z=H
cc4=sy.Eq(Bw.subs(z,H),0) #Bimomento nullo in z=H

## risolvo il sistema
sol=sy.solve((cc1,cc2,cc3,cc4),(A,B,C,D))
w=w.subs(A,sol[A]).subs(B,sol[B]).subs(C,sol[C]).subs(D,sol[D])
dw=dw.subs(A,sol[A]).subs(B,sol[B]).subs(C,sol[C]).subs(D,sol[D])
Mt=Mt.subs(A,sol[A]).subs(B,sol[B]).subs(C,sol[C]).subs(D,sol[D])
Bw=Bw.subs(A,sol[A]).subs(B,sol[B]).subs(C,sol[C]).subs(D,sol[D])
Mts=-sy.diff(Bw,z)


#print(sol)
print('w=', sy.simplify(w))
print('dw=', sy.simplify(dw))
print('Mt=', sy.simplify(Mt))
print('Mts=', sy.simplify(Mts))
print('Bw=', sy.simplify(Bw))