import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

E=1. # Modulo elastico [N/mm2]
fyk=1. # tensione di snervamento [N/mm2]
t=np.arange( 0.,np.pi/2,0.01)  # renge di tempo [s]
eyk=fyk/E # epslon snervamento
pyk=eyk-fyk/E # p snervamento

e=np.sin(t) # funzione epslon armonica
#e=0.01*np.sin(t)*np.exp(-t/20) # funzione epslon armonica

sigma=np.zeros(len(t))
p=np.zeros(len(t))
f=np.zeros(len(t))
f[0]=-fyk

for i in range(1,(len(t))):
    p[i]=p[i-1]
    sigma[i]=E*(e[i]-p[i]) # legame costitutivo molla
    f[i]=np.abs(sigma[i])-fyk #funzione di snervamento

    if f[i] < 0: # FASE ELASTICA
        p[i]=p[i-1]
        sigma[i]=E*(e[i]-p[i]) # legame costitutivo molla

    elif f[i]>0: # SCORRIMENTO PLASTICO
        if sigma[i]>0 :
            p[i]=e[i]-fyk/E
            sigma[i]=E*(e[i]-p[i]) # legame costitutivo molla
            f[i]=(sigma[i])-fyk #funzione di snervamento

        if sigma[i]<0 :
            p[i]=e[i]+fyk/E
            sigma[i]=E*(e[i]-p[i]) # legame costitutivo molla
            f[i]=-(sigma[i])-fyk #funzione di snervamento



## PLOT ##
plt.ion()
plt.show ()
fig, axs = plt.subplots(2, 2)
plt.suptitle('Modello elasto-plastico perfetto',fontsize=20,color='darkred',fontstyle='italic')
plt.subplots_adjust(top=0.89,
bottom=0.11,
left=0.08,
right=0.965,
hspace=0.44,
wspace=0.265)
#Plot epslon in funzione del tempo
axs[0, 0].plot(t,e,'-g',linewidth=1)
axs[0, 0].grid(True)
axs[0, 0].set_title('Diagramma $\epsilon_{(t)}$ - t')
axs[0, 0].set_xlim(min(t), max(t))
axs[0, 0].set_ylim(min(e)-0.002, max(e)+0.002)
axs[0, 0].set_xlabel('t [s]')
axs[0, 0].set_ylabel('$\epsilon_{(t)}$')
#Plot sigma in funzione del tempo
axs[1, 0].plot(t,sigma,'-r',linewidth=1)
axs[1, 0].grid(True)
axs[1, 0].set_title('Diagramma $\sigma_{(t)}$ - t')
axs[1, 0].set_xlim(min(t), max(t))
axs[1, 0].set_ylim(min(sigma)-20, max(sigma)+20)
axs[1, 0].set_xlabel('t [s]')
axs[1, 0].set_ylabel('$\sigma_{(t)}$ $[N/mm^2]$')
#Plot sigma in funzione di epslon
axs[0, 1].plot(e,sigma,'-b',linewidth=1)
axs[0, 1].grid(True)
axs[0, 1].set_title('Diagramma $\sigma_{(t)}$- $\epsilon_{(t)}$')
axs[0, 1].set_xlim(min(e)-0.005, max(e)+0.005)
axs[0, 1].set_ylim(min(sigma)-50, max(sigma)+50)
axs[0, 1].set_xlabel('$\epsilon_{(t)}$ [-]')
axs[0, 1].set_ylabel('$\sigma_{(t)}$ $[N/mm^2]$')
#Plot f in funzione di t
axs[1, 1].plot(t,f,'-k',linewidth=1)
axs[1, 1].grid(True)
axs[1, 1].set_title('Funzione di snervamento')
axs[1, 1].set_xlim(min(t), max(t))
axs[1, 1].set_ylim(max(f)+20, min(f)-20)
axs[1, 1].set_xlabel('$t [s]$')
axs[1, 1].set_ylabel('$f_{(t)} [-]$')
axs[1, 1].invert_yaxis()
