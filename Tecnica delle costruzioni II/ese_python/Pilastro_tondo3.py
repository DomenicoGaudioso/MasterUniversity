import numpy as np
import sympy
from sympy import symbols ,solve ,Eq, diff,pprint
import matplotlib.pyplot as plt
import pylab

##INPUT##
De=406.60 # mm
t=16.00 # mm
fyd=224 # S235 [N/mm2]
fcd=14.17  # C25/30 [N/mm2]
##
fcd2=0.85*fcd
Di=De-2*t # Diametro interno
Ri=Di/2 # Raggio interno
Rm=(De-t)/2 # Raggio medio acciao

## FLESSIONE RETTA
i=np.linspace(t,De,De/100)
for n in range(0,len(i)):

alpha=np.arccos(1-0.85*(yn-t)/Ri)  # angolo in corrispondenza dell'asse neutro
Fc=2*Ri**2*(alpha-np.sin(alpha)*np.cos(alpha))*fcd2  # Forza di compressione nel calcestruzzo
Fa1=2*alpha*Rm*t*fyd # Forza di compressione nell'acciaio
Fa2=(2*np.pi*Rm-(2*np.pi-2*alpha)*Rm)*t*fyd # Forza di trazione nell'acciaio
F=Fc+Fa1-Fa2 # Equilibrio














