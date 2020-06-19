###TRAVE INCASTRO-INCASTRO CON TRE CARICHI CONCENTRATI ###
import numpy as np
import sympy
from sympy import symbols ,solve ,Eq, diff,pprint
import matplotlib.pyplot as plt
import pylab

##INPUT##
De=406.6 # mm
t=16 # mm
fyd=224 # S235 [N/mm2]
fcd=14.17  # C25/30 [N/mm2]
##
fcd2=0.85*fcd
Di=De-2*t # Diametro interno
Ri=Di/2 # Raggio interno
Rm=(De-t)/2 # Raggio medio acciao

## FLESSIONE RETTA
#def alpha(yn):  # angolo in corrispondenza dell'asse neutro
#	return np.arccos(1-0.85*(yn-t)/Ri)

def F(alpha):
	return 2*Ri**2*(alpha-np.sin(alpha)*np.cos(alpha))*fcd2+2*alpha*Rm*t*fyd-(2*np.pi*Rm-(2*np.pi-2*alpha)*Rm)*t*fyd


def d_F(alpha):
	#return 2*Ri**2*(1-np.sin(alpha)**2+np.cos(alpha)**2)*fcd2+2*Rm*t*fyd+2*Rm*t*fyd
	return sympy.diff(F,x)



def tangenti(alpha_k, scarto, f, df):
    #controllo sulla vicinanza del risultato allo zero
	while f(alpha_k) > scarto:
		alpha_k = alpha_k - ( f(alpha_k) / df(alpha_k) )
		tangenti(alpha_k, scarto, f, df)
	return alpha_k

init_val = 0
scarto = 0.01

alpha_k = tangenti(init_val, scarto, F, d_F)
print(alpha_k)
print(F(alpha_k))
print("E' la soluzione corretta?"),
ans = F(alpha_k) <= scarto
print(ans)

