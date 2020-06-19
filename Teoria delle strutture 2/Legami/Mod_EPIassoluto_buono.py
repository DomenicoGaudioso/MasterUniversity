import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

E=1 # Modulo elastico [N/mm2]
A=E*np.sin(np.pi/4) # tensione di snervamento [N/mm2]
t=np.arange( 0.,2.5*np.pi,0.0001)  # range di tempo [s]
eyk=A/E # epslon snervamento
B=E

e=np.sin(t) # funzione epslon armonica
# e=np.sin(t)*np.exp(t/20) # funzione epslon armonica

sigma=np.zeros(len(t))
p=np.zeros(len(t))
f=np.zeros(len(t))
f[0]=-A

def sigmaf(e,p):
    return E*(e-p)


for i in range(1,(len(t))):
    p[i]=p[i-1]
    sigma[i]=sigmaf(e[i],p[i]) # legame costitutivo molla
    f[i]=np.abs(sigmaf(e[i],p[i]))-(A+B*np.abs(p[i])) #funzione di snervamento

    if f[i] <= 0: # FASE ELASTICA
         p[i]=p[i-1]
         sigma[i]=sigmaf(e[i],p[i])

    elif f[i]>0 and sigma[i]>=0  : # SCORRIMENTO PLASTICO 1 #and  e[i]>eyk

        p[i]=(e[i]*E-A)/(E+B)

        sigma[i]=sigmaf(e[i],p[i]) # legame costitutivo molla

        f[i]=sigma[i]-(A+B*(p[i])) #funzione di snervamento

    elif f[i]>0 and sigma[i]<=0 : # SCORRIMENTO PLASTICO 2 #e[i] < -eyk

        p[i]=(e[i]*E+A)/(E+B)

        sigma[i]=sigmaf(e[i],p[i]) # legame costitutivo molla

        f[i]=-(sigma[i])-(A+B*(-p[i])) #funzione di snervamento



exec (open ("Legamiplot.py") .read ())
plt.suptitle('Modello elasto-plastico incrudente $\sigma^*=(A+B\cdot|p|)$',fontsize=20,color='darkred',fontstyle='italic')

file = open('datiscritti_prova4.txt', 'w')  # Apertura file in scrittura ('w')

for i in range(len(t)):
                                     # Preparazione stringa da scrivere
    st = str(round(sigma[i], 2)) + '\t' + str(round(e[i], 2)) + '\n' # ho arrotondato con round
    file.write(st)                   # Scrittura su file
file.close() # Chiusura file