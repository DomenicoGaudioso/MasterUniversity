import numpy as np
import sympy
from sympy import symbols ,solve ,Eq, diff,pprint
import matplotlib.pyplot as plt
import math
import pylab

##INPUT##
De=406.4 # mm
t=16 # mm
fyd=224 # S235 [N/mm2]
fck=25 # C25/30
fcd=0.85*fck/1.5 # C25/30 [N/mm2]
##Parametri calcolati
#Parametri geometrici
Di=De-2*t # Diametro interno
Re=De/2
Ri=Di/2 # Raggio interno
Rm=(De-t)/2 # Raggio medio acciao
Aa=np.pi*(De**2-Di**2)/4 #Area acciaio
Ac=np.pi*(Di)**2/4 # Area cls

#Parametri sui materiali
fcd2=0.85*fcd #Resistenza a compressione per stress-block (cls)
## FORZA NORMALE (PURA) No confinamento [PUNTO A]
Nprd=(Ac*fcd+Aa*fyd)/1000 # Forza normale resistente [KN]

## FLESSIONE RETTA [PUNTO B]
# alpha compreso tra 0 e pi/2
i=np.linspace(0,De/2,1000)
#Per trovare l'asse neutro
plt.figure(1)
plt.ion()
plt.show()
for n in range(0,len(i)):
    yn=i[n]
    alphan=(np.arccos(1-(yn)/Re))  # angolo in corrispondenza dell'asse neutro
    alphanc=(np.arccos(1-0.8*yn/Ri))  # angolo in corrispondenza di 0.8*yn
    Fc=Ri**2*(alphanc-np.sin(alphanc)*np.cos(alphanc))*fcd2  # Forza di compressione nel calcestruzzo
    Fa1=2*alphan*Rm*t*fyd # Forza di compressione nell'acciaio
    Fa2=(2*np.pi*Rm-(2*alphan)*Rm)*t*fyd # Forza di trazione nell'acciaio
    F=Fc+Fa1-Fa2 # Equilibrio
    plt.plot(F,yn,marker='o',markersize=2, color='red')
    if 0.00<=F<=0.6:
        print(yn)
        break
# Calcolo Mrd per FLESSIONE SEMPLICE
yga1=Re*np.sin(alphan)/alphan # baricentro acciaio compresso
yga2=yn+Re*np.sin(np.pi-alphan)/(np.pi-alphan) # baricentro acciaio teso
ygc=t+2*Ri*np.sin(alphanc)**3/(3*(alphanc-np.sin(alphanc)*np.cos(alphanc))) #baricentro calcestruzzo compresso
Mprd=(Fa1*(yn-yga1)+Fc*(yn-ygc)+Fa2*(yga2-yn))/(1000**2) # Momento resistente [KNm]

##[PUNTO C] Momento resistente plastico pari a quello ottenuto per flessione semplice e sforzo normale resistente di progetto della sola porzione di calcestruzzo della sezione composta
Npm_rd=fcd*Ac/1000
MprdC=Mprd
##[PUNTO D] Momento resistente plastico Massimio e sforzo normale pari a metà di quello ottenuto al punto C
alphaD=(np.arccos(1-0.8*Ri/Ri))
NrdD=0.5*Npm_rd
Fc_D=Ri**2*(alphaD-np.sin(alphaD)*np.cos(alphaD))*fcd2
Fa=Aa/2*fyd
yga=2*Re/np.pi # baricentro acciaio teso e compresso
ygcD=t+2*Ri*np.sin(alphaD)**3/(3*(alphaD-np.sin(alphaD)*np.cos(alphaD))) #baricentro calcestruzzo
Mrd_D=(Fa*(Re-yga)+Fa*(Re-yga)+Fc*(Re-ygcD))/(1000**2)
print(Mrd_D,Mprd)

##VERIFICA A PRESSOFLESSIONE
mu=0.9
plt.figure()
plt.ion()
plt.show()
MRD=[0,mu*MprdC,mu*Mrd_D,mu*Mprd]
NRD=[Nprd,Npm_rd,NrdD,0]
plt.plot(MRD,NRD,marker='o',markersize=4, markerfacecolor='k', color='blue',lw=1)


























