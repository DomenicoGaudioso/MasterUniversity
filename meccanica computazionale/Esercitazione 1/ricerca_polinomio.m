clear all
close all

syms s B C D E F L
cc= sym('c', [1 4])

A=1
B=5

j=0
phi(s)=A*s^(4+(j))+B*s^(3+(j))+C*s^(2+(j))+D*s^(1+(j))+E*s^(j)

phi(s)=phi
dphi(s)=diff(phi)


cc(1)=dphi(0)==0
cc(2)=dphi(L)==0
cc(3)=phi(0)==0
cc(4)=phi(L)==0

[MatA, b] =equationsToMatrix([cc],[B C D E])

Costanti=linsolve(MatA,b)

phi=subs(phi,C,Costanti(1))
phi=subs(phi,D,Costanti(2))
phi=subs(phi,E,Costanti(3))
phi=subs(phi,F,Costanti(4))

pretty(phi)

