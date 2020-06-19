import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
import pylab as py
from mpl_toolkits.mplot3d import Axes3D

## SEZIONE AD H
# Dati sezione
B=300.
H=B
ta=10.
tp=20.
# coordinate
x_neg=list(np.linspace(-B/2,0, num=11,endpoint=True))
x_pos=list(np.linspace(0,B/2, num=11,endpoint=True))
y=np.linspace(-H/2,H/2, num=11,endpoint=True)

## DATI TRAVE
L=8000.
z=np.linspace(0,L, num=11,endpoint=True)

## MESHGRID
X_pos, Z = np.meshgrid(x_pos,z)
X_neg,Z=np.meshgrid(x_neg,z)
Y,Z=np.meshgrid(y,z)

##SOLLECITAZIONE
T1=50*(L-Z)
T2=50*(L-Z)
M2=50*(L*Z-Z**2/2)
M1=50*(L*Z-Z**2/2)
Mts=500000
##Ascisse
#Ascisse curviline sulla sezione per S1 e S2
n1sx=B/2-np.float32(X_pos)
n1dx=B/2+np.float32(X_neg)
n2=H/2+Y
#Ascisse curviline sulla sezione per Sphi (Momento statico settoriale)
s1=X_pos
s2=Y
s3=X_neg
## CALCOLO FUNZIONE D'INGOBBAMENTO
psi1=-s1*H/2
psi2=0
psi3=-s1*H/2

## CALCOLO MOMENTI STATICI
#Momento statico S1 ali superiori
S1_1sx=tp*n1sx*(-H/2)
S1_1dx=tp*n1dx*(-H/2)
# Momento statico S1 anima
S1_2=tp*B*(-H/2)-n2*ta*(H/2-n2/2)
#Momento statico S1 ali inferiori
S1_3sx=tp*(n1sx)*(H/2)
S1_3dx=tp*(n1dx)*(H/2)

#Momento statico S2 ali superiori
S2_1sx=tp*n1sx*(B/2-n1sx/2)
S2_1dx=tp*n1dx*(B/2-n1dx/2)
# Momento statico S2 anima
S2_2=0
#Momento statico S2 ali inferiori
S2_3sx=tp*n1sx*(B/2-n1sx/2)
S2_3dx=tp*n1dx*(B/2-n1dx/2)

#Momento statico Settoriale ali sup
Spsi_s1sx=-H*tp*(B**2/4-s1**2)/4
Spsi_s1dx=+H*tp*(-B**2/4+s3**2)/4
#Momento statico Settoriale ali inf
Spsi_s3sx=+H*tp*(B**2/4-s1**2)/4
Spsi_s3dx=-H*tp*(-B**2/4+s3**2)/4
#Momento statico Settoriale anima
Spsi_s2=0

##CALCOLO AREA
A=B*H-2*(B-ta)*(H-2*tp)

## CALCOLO MOMENTI D'INERZIA
# Momento d'inerzia rispetto l'asse 1
I1=(B*H**3-2*(B-ta)*(H-2*tp)**3)/12
# Momento d'inerzia rispetto l'asse 2
I2=(H*B**3-2*(H-2*tp)*(B-ta)**3)/12
# Momento d'inerzia polare
Io=I1+I2
# Momento d'inerzia settoriale
Ipsi=(H**2*B**3*tp)/24


## CALCOLO TENSIONI TANGENZIALI
# tensioni (T2) tangenziali ali superiori
tauT2_1sx=T2*S1_1sx/(I1*tp)
tauT2_1dx=T2*S1_1dx/(I1*tp)
# tensioni (T2) tangenziali ali inferiori
tauT2_3sx=T2*S1_3sx/(I1*tp)
tauT2_3dx=T2*S1_3dx/(I1*tp)
# tensioni (T2) tangenziali anima
tauT2_2=T2*S1_2/(I1*ta)

# tensioni (T1) tangenziali ali superiori
tauT1_1sx=T1*S2_1sx/(I1*tp)
tauT1_1dx=T1*S2_1dx/(I1*tp)
# tensioni (T1) tangenziali ali inferiori
tauT1_3sx=T1*S2_3sx/(I1*tp)
tauT1_3dx=T1*S2_3dx/(I1*tp)
# tensioni (T1) tangenziali anima
tauT1_2=T1*S2_2/(I1*ta)

# tensioni (Mts) tangenziali ali superiori
tauMts_s1sx=Mts*Spsi_s1sx/(Ipsi*tp)
tauMts_s1dx=Mts*Spsi_s1dx/(Ipsi*tp)
# tensioni (Mts) tangenziali ali inferiori
tauMts_s3sx=Mts*Spsi_s3sx/(Ipsi*tp)
tauMts_s3dx=Mts*Spsi_s3dx/(Ipsi*tp)
# tensioni (Mts) tangenziali anima
tauMts_s2=Mts*Spsi_s2/(Ipsi*ta)

## CALCOLO TENSIONI NORMALI
# tensioni normali dovuti a M1 sulla piattabanda
sigmazM1_1=(M1/I1)*Y[:,1]
# tensioni normali dovuti a M1 sull'anima
sigmazM1_2=(M1/I1)*Y

# tensioni normali dovuti a M2 sulla piattabanda
sigmazM2_sx=-(M2/I1)*X_pos
sigmazM2_dx=-(M2/I1)*X_neg

### PLOT
plt.ion()
plt.show ()
## Tensioni per T2
plt.figure().gca(projection='3d')
delta=20
#ali sup
plt.title('$\Gamma_{T_2}$ sezione ad H [$N/mm^2$]')
plt.contour(X_pos, Z, tauT2_1sx, delta, zdir='z', offset=H/2)
plt.contour(X_neg, Z, tauT2_1dx, delta, zdir='z', offset=H/2)
A1=plt.contourf(X_pos, Z, tauT2_1sx, delta, zdir='z', offset=H/2,cmap='jet')
A2=plt.contourf(X_neg, Z, tauT2_1dx, delta, zdir='z', offset=H/2 ,cmap='jet')
# plt.quiver(X_pos,Z,tauT2_1sx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
# plt.quiver(X_neg,Z,-tauT2_1dx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(A1)
#ali inf
plt.subplot(1,3,2)
plt.contour(X_pos, Z, -tauT2_3sx, delta, zdir='z', offset=-H/2)
plt.contour(X_neg, Z, -tauT2_3dx, delta, zdir='z', offset=-H/2)
B1=plt.contourf(X_pos, Z, -tauT2_3sx, delta, zdir='z', offset=-H/2,cmap='jet')
B2=plt.contourf(X_neg, Z, -tauT2_3dx, delta, zdir='z', offset=-H/2,cmap='jet')
#plt.quiver(X_pos,Z,tauT2_3sx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
#plt.quiver(X_neg,Z,-tauT2_3dx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(B1)
#anima
plt.subplot(1,3,3)
plt.contour(Y, Z, tauT2_2, delta)
C=plt.contourf(Y, Z, tauT2_2, delta,cmap='jet')
#plt.quiver(Y,Z,tauT2_2,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' y  (mm)')
plt.ylabel(' z (mm)')
plt.gca().invert_yaxis()
plt.colorbar(C)

##  Tensioni per T1
plt.figure()
delta=20
#ali sup
plt.subplot(1,3,(1))
plt.suptitle('$\Gamma_{T_1}$ sezione ad H [$N/mm^2$]')
plt.contour(X_pos, Z, -tauT1_1sx, delta)
plt.contour(X_neg, Z, -tauT1_1dx, delta)
A1=plt.contourf(X_pos, Z, -tauT1_1sx, delta,cmap='jet')
A2=plt.contourf(X_neg, Z, -tauT1_1dx, delta,cmap='jet')
plt.quiver(X_pos,Z,tauT1_1sx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.quiver(X_neg,Z,tauT1_1dx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.title('$\Gamma_{xz}$ ali superiori [$N/mm^2$]')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(A1)
#ali inf
plt.subplot(1,3,2)
plt.contour(X_pos, Z, -tauT1_3sx, delta)
plt.contour(X_neg, Z, -tauT1_3dx, delta)
B1=plt.contourf(X_pos, Z, -tauT1_3sx, delta,cmap='jet')
B2=plt.contourf(X_neg, Z, -tauT1_3dx, delta,cmap='jet')
plt.quiver(X_pos,Z,tauT1_3sx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.quiver(X_neg,Z,tauT1_3dx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.title('$\Gamma_{xz}$ ali inferiori [$N/mm^2$]')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(B1)
#anima
# plt.subplot(1,3,3)
# plt.title('$\Gamma_{yz}$ anima [$N/mm^2$]')
# plt.xlabel(' y  (mm)')
# plt.ylabel(' z (mm)')
# plt.gca().invert_yaxis()

## Tensioni tangenziali per Mts
plt.figure()
delta=20
#ali sup
plt.subplot(1,3,(1))
plt.suptitle('$\Gamma_{Mt^s}$ sezione ad H [$N/mm^2$]')
plt.contour(X_pos, Z, -tauMts_s1sx, delta)
plt.contour(X_neg, Z, -tauMts_s1dx, delta)
A1=plt.contourf(X_pos, Z, -tauMts_s1sx, delta,cmap='jet')
A2=plt.contourf(X_neg, Z, -tauMts_s1dx, delta,cmap='jet')
plt.quiver(X_pos,Z,tauMts_s1sx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.quiver(X_neg,Z,tauMts_s1dx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.title('$\Gamma_{xz}$ ali superiori [$N/mm^2$]')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(A1)
#ali inf
plt.subplot(1,3,2)
plt.contour(X_pos, Z, tauMts_s3sx, delta)
plt.contour(X_neg, Z, tauMts_s3dx, delta)
B1=plt.contourf(X_pos, Z, tauMts_s3sx, delta,cmap='jet')
B2=plt.contourf(X_neg, Z, tauMts_s3dx, delta,cmap='jet')
plt.quiver(X_pos,Z,tauMts_s3sx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.quiver(X_neg,Z,tauMts_s3dx,0*Z,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.title('$\Gamma_{xz}$ ali inferiori [$N/mm^2$]')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(B1)
#anima
# plt.subplot(1,3,3)
# plt.title('$\Gamma_{yz}$ anima [$N/mm^2$]')
# plt.xlabel(' y  (mm)')
# plt.ylabel(' z (mm)')
# plt.gca().invert_yaxis()

# plt.plot(x_pos,np.transpose(tauMts_s3sx))
# plt.plot(x_neg,-np.transpose(tauMts_s3dx))
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()

## Tensioni normali per M1
plt.figure()
delta=20
#ali sup
plt.subplot(1,3,(1))
plt.suptitle('$\sigma_{M1}$ sezione ad H [$N/mm^2$]')
plt.contour(X_pos, Z, sigmazM1_1, delta)
plt.contour(X_neg, Z, sigmazM1_1, delta)
A1=plt.contourf(X_pos, Z, sigmazM1_1, delta,cmap='jet')
A2=plt.contourf(X_neg, Z, sigmazM1_1, delta,cmap='jet')
plt.quiver(X_pos,Z,0*Z,sigmazM1_1,color='k', units='inches',width=0.01)
plt.quiver(X_neg,Z,0*Z,sigmazM1_1,color='k', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.title('$\Gamma_{xz}$ ali superiori [$N/mm^2$]')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(A1)
#ali inf
plt.subplot(1,3,2)
plt.contour(X_pos, Z, sigmazM1_1, delta)
plt.contour(X_neg, Z, sigmazM1_1, delta)
A1=plt.contourf(X_pos, Z, sigmazM1_1, delta,cmap='jet')
A2=plt.contourf(X_neg, Z, sigmazM1_1, delta,cmap='jet')
plt.quiver(X_pos,Z,0*Z,-sigmazM1_1,color='k',pivot='tip', units='inches',width=0.01)
plt.quiver(X_neg,Z,0*Z,-sigmazM1_1,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.title('$\Gamma_{xz}$ ali superiori [$N/mm^2$]')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(A1)
#anima
plt.subplot(1,3,3)
plt.contour(Y, Z, sigmazM1_2, delta)
C=plt.contourf(Y, Z, sigmazM1_2, delta,cmap='jet')
plt.quiver(Y,Z,0*Z,sigmazM1_2,color='k',units='inches',width=0.01)
plt.title('$\sigma_{yz}$ anima [$N/mm^2$]')
plt.xlabel(' y  (mm)')
plt.ylabel(' z (mm)')
plt.gca().invert_yaxis()
plt.colorbar(C)

## Tensioni normali per M2
plt.figure()
delta=20
#ali sup
plt.subplot(1,3,(1))
plt.suptitle('$\sigma_{M2}$ sezione ad H [$N/mm^2$]')
plt.contour(X_pos, Z, sigmazM2_sx, delta)
plt.contour(X_neg, Z, -sigmazM2_dx, delta)
A1=plt.contourf(X_pos, Z, sigmazM2_sx, delta,cmap='jet')
A2=plt.contourf(X_neg, Z, -sigmazM2_dx, delta,cmap='jet')
plt.quiver(X_pos,Z,0*Z,sigmazM2_sx,color='k', units='inches',width=0.01)
plt.quiver(X_neg,Z,0*Z,sigmazM2_dx,color='k', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.title('$\Gamma_{xz}$ ali superiori [$N/mm^2$]')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(A1)
#ali inf
plt.subplot(1,3,2)
plt.contour(X_pos, Z, sigmazM2_sx, delta)
plt.contour(X_neg, Z, -sigmazM2_dx, delta)
A1=plt.contourf(X_pos, Z, sigmazM2_sx, delta,cmap='jet')
A2=plt.contourf(X_neg, Z, -sigmazM2_dx, delta,cmap='jet')
plt.quiver(X_pos,Z,0*Z,sigmazM2_sx,color='k', units='inches',width=0.01)
plt.quiver(X_neg,Z,0*Z,sigmazM2_dx,color='k', units='inches',width=0.01)
plt.xlabel(' x  (mm)')
plt.ylabel(' z (mm)')
plt.title('$\Gamma_{xz}$ ali superiori [$N/mm^2$]')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.colorbar(A1)








