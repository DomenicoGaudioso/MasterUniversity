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
fcd=14.17  # C25/30 [N/mm2]
##
fcd2=0.85*fcd
Re=De/2
Di=De-2*t # Diametro interno
Ri=Di/2 # Raggio interno
Rm=(De-t)/2 # Raggio medio acciao

## FLESSIONE RETTA
# alpha compreso tra 0 e pi/2
i=np.linspace(0,De/2,100000000)
#Per trovare l'asse neutro
for n in range(0,len(i)):
    yn=i[n]
    alphan=(np.arccos(1-(yn)/Re))  # angolo in corrispondenza dell'asse neutro
    alphanc=(np.arccos(1-yn/Re))  # angolo in corrispondenza di 0.8*yn
    Fc=Ri**2*(alphanc-np.sin(alphanc)*np.cos(alphanc))*fcd2  # Forza di compressione nel calcestruzzo
    Fa1=2*alphan*Rm*t*fyd # Forza di compressione nell'acciaio
    Fa2=(2*np.pi*Rm-(2*alphan)*Rm)*t*fyd # Forza di trazione nell'acciaio
    F=Fc+Fa1-Fa2 # Equilibrio
    if 0.00<=F<=0.4:
        print(yn)






















