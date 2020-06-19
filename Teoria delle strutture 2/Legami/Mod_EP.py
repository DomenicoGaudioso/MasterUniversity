import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

E=1 # Modulo elastico [N/mm2]
A=E*np.sin(np.pi/4) # tensione di snervamento [N/mm2]
t=np.arange( 0.,2.5*np.pi,0.01)  # range di tempo [s]
ek=A/E # epslon snervamento
pk=ek-A/E # p snervamento

e=np.sin(t) # funzione epslon armonica
#e=0.01*np.sin(t)*np.exp(-t/10) # funzione epslon armonica

sigma=np.zeros(len(t))
p=np.zeros(len(t))
f=np.zeros(len(t))
f[0]=-A

for i in range(1,(len(t))):
    p[i]=p[i-1]
    sigma[i]=E*(e[i]-p[i]) # legame costitutivo molla
    f[i]=np.abs(sigma[i])-A #funzione di snervamento

    if f[i] < 0: # FASE ELASTICA
        p[i]=p[i-1]
        sigma[i]=E*(e[i]-p[i]) # legame costitutivo molla

    elif f[i]>0: # SCORRIMENTO PLASTICO
        if sigma[i]>0 :
            p[i]=e[i]-A/E
            sigma[i]=E*(e[i]-p[i]) # legame costitutivo molla
            f[i]=(sigma[i])-A #funzione di snervamento

        if sigma[i]<0 :
            p[i]=e[i]+A/E
            sigma[i]=E*(e[i]-p[i]) # legame costitutivo molla
            f[i]=-(sigma[i])-A #funzione di snervamento

exec (open ("Legamiplot.py") .read ())
plt.suptitle('Modello elasto-plastico perfetto $\sigma^*= A = sin(\pi/4)$',fontsize=20,color='darkred',fontstyle='italic')

file = open('datimodelloEP.txt', 'w')  # Apertura file in scrittura ('w')

for i in range(len(t)):
# Preparazione stringa da scrivere
# ho arrotondato con round
    st = str(round(sigma[i], 2)) + '\t' + str(round(e[i], 2)) + '\n'
    file.write(st)   # Scrittura su file
file.close() # Chiusura file