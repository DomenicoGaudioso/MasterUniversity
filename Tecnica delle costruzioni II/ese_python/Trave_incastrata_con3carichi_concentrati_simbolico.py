###TRAVE INCASTRO-INCASTRO CON TRE CARICHI CONCENTRATI ###
import numpy as np
import sympy
from sympy import symbols ,solve ,Eq, diff,pprint
import matplotlib.pyplot as plt
import pylab

## SIMBOLICO ##
c1, c2, c3, c4, c5, c6, c7, c8 = symbols ('c1 c2 c3 c4 c5 c6 c7 c8')
v1,v2=symbols ('v1 v2')
E, J=symbols ('E J')
s, a, P=symbols ('s a P')

## LINEA ELASTICA TRATTO 1 ##
v1=(c1*s**3/6+c2*s**2/2+c3*s+c4)*1/(E*J) #spostamento tratto 1
dv1=diff(v1,s) #rotazione tratto 1
M1=-E*J*diff(v1,s,s) #Momento 1
T1=-E*J*diff(v1,s,s,s) #Taglio tratto 1

## LINEA ELASTICA TRATTO 2 ##
v2=(c5*s**3/6+c6*s**2/2+c7*s+c8)*1/(E*J) #spostamento tratto 2
dv2=diff(v2,s) #rotazione tratto 2
M2=-E*J*diff(v2,s,s) #Momento tratto 2
T2=-E*J*diff(v2,s,s,s) #Taglio tratto 2

## CONDIZIONI A CONTORNO TRATTO 1 ##
cc1=Eq(v1.subs(s,0),0) #spostamento nullo in A
cc2=Eq(dv1.subs(s,0),0) #rotazione nulla in A
cc3=Eq(v1.subs(s,a),v2.subs(s,a)) #congruenza sullo spostamento in C
cc4=Eq(T1.subs(s,a),3*P/2) #Taglio in C
cc5=Eq(dv1.subs(s,a),dv2.subs(s,a)) #congruenza sulla rotazione in C

## CONDIZIONI A CONTORNO TRATTO 2 ##
cc6=Eq(dv2.subs(s,2*a),0) #rotazione nulla in B
cc7=Eq(T2.subs(s,2*a),P/2) #Taglio in B
cc8=Eq(M1.subs(s,a),M2.subs(s,a)) #congruenza sul momento in C

## RISOLVO IL SISTEMA
sol=solve((cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8),(c1,c2,c3,c4,c5,c6,c7,c8))
# print(sol.keys ())
vf1=v1.subs(c1,sol[c1]).subs(c2,sol[c2]).subs(c3,sol[c3]).subs(c4,sol[c4])
vf2=v2.subs(c5,sol[c5]).subs(c6,sol[c6]).subs(c7,sol[c7]).subs(c8,sol[c8])
Mf1=M1.subs(c1,sol[c1]).subs(c2,sol[c2])
Mf2=M2.subs(c5,sol[c5]).subs(c6,sol[c6])
Tf1=T1.subs(c1,sol[c1])
Tf2=T2.subs(c5,sol[c5])

print( Mf1.subs(s,0))

## PROBLEMA NUMERICO
En=(210000)*1000
Jn=float(179768785.01)/(1000**4)
an=8.5/4.
Pn=133.61 #Allo SLU
#Pn=94.51 #RARA
#Pn=78.57 #Quasi Perm
L=8.5/2.
delta=0.96
vn1=vf1.subs(E,En).subs(J,Jn).subs(a,an).subs(P,Pn)
vn2=vf2.subs(E,En).subs(J,Jn).subs(a,an).subs(P,Pn)
Mn1=Mf1.subs(E,En).subs(J,Jn).subs(a,an).subs(P,Pn)
Mn2=Mf2.subs(E,En).subs(J,Jn).subs(a,an).subs(P,Pn)
Tn1=Tf1.subs(E,En).subs(J,Jn).subs(a,an).subs(P,Pn)
Tn2=Tf2.subs(E,En).subs(J,Jn).subs(a,an).subs(P,Pn)
print(Mn1.subs(s,0)) # momento all'incastro
print(Mn2.subs(s,2*an)) # momento in mezzeria

## PLOT
s1=np.arange( 0.00, an , 0.005) #ascissa tratto 1
s2=np.arange( an, L , 0.005) #ascissa tratto 2

y1=np.zeros(len(s1))
y2=np.zeros(len(s2))
M1y=np.zeros(len(s1))
M2y=np.zeros(len(s2))
T1y=np.zeros(len(s1))
T2y=np.zeros(len(s2))

for i in range(0,(len(s1))):
      y1[i]=vn1.subs(s,s1[i])
      M1y[i]=Mn1.subs(s,s1[i])
      T1y[i]=Tn1.subs(s,s1[i])
for i in range(0,(len(s2))):
      y2[i]=vn2.subs(s,s2[i])
      M2y[i]=Mn2.subs(s,s2[i])
      T2y[i]=Tn2.subs(s,s2[i])

# plot deformata
pylab.subplot(212)
pylab.plot(s1,y1*1000,'-g', linewidth=1.5)
pylab.plot(s2,y2*1000,'-g', linewidth=1.5)
pylab.xlabel('s [m]')
pylab.ylabel('v(s) [mm]')
pylab.xlim(0, L)
pylab.ylim(0, max(y2)*1000+10)
pylab.grid(True)
pylab.gca().invert_yaxis()
pylab.title('Deformata')
# plot momento
pylab.subplot(221)
pylab.plot(s1,M1y,'-b',linewidth=1.5)
pylab.fill_between(s1,M1y, color="skyblue", alpha=0.4)
pylab.plot(s2,M2y,'-b',linewidth=1.5)
pylab.fill_between(s2,M2y, color="skyblue", alpha=0.4)
pylab.xlabel('s [m]')
pylab.ylabel('M(s) [KNm]')
pylab.xlim(0, L)
pylab.grid(True)
pylab.gca().invert_yaxis()
pylab.title('Momento')
# plot taglio
pylab.subplot(222)
pylab.plot(s1,T1y,'-r',linewidth=1.5)
pylab.fill_between(s1,T1y, color="red", alpha=0.4)
pylab.plot(s2,T2y,'-r',linewidth=1.5)
pylab.fill_between(s2,T2y, color="red", alpha=0.4)
pylab.xlabel('s [m]')
pylab.ylabel('T(s) [KN]')
pylab.xlim(0, L)
pylab.grid(True)
pylab.gca().invert_yaxis()
pylab.title('Taglio')
pylab.show ()