## PLOT ##
plt.ion()
plt.show ()
fig, axs = plt.subplots(4, 2)
plt.subplots_adjust(top=0.88,
bottom=0.07,
left=0.07,
right=0.99,
hspace=0.96,
wspace=0.18)

#Plot epslon in funzione del tempo
axs[0, 0].plot(t,e,'-g',linewidth=1)
axs[0, 0].grid(True)
axs[0, 0].set_title('Diagramma $\epsilon_{(t)}$ - t')
axs[0, 0].set_xlim(min(t), max(t))
# axs[0, 0].set_ylim(min(e)-0.002, max(e)+0.002)
axs[0, 0].set_xlabel('t [s]')
axs[0, 0].set_ylabel('$\epsilon_{(t)}$')

#Plot sigma in funzione del tempo
axs[1, 0].plot(t,sigma,'-r',linewidth=1)
axs[1, 0].grid(True)
axs[1, 0].set_title('Diagramma $\sigma_{(t)}$ - t')
axs[1, 0].set_xlim(min(t), max(t))
# axs[1, 0].set_ylim(min(sigma)-20, max(sigma)+20)
axs[1, 0].set_xlabel('t [s]')
axs[1, 0].set_ylabel('$\sigma_{(t)}$ $[N/mm^2]$')

#Plot sigma in funzione di epslon

col=["k", "m", "c", "g" , "y", "b" ]

ax3 = plt.subplot(122)
ax3.plot(e,sigma,'-',color=col[5],linewidth=1.5)
ax3.grid(True)
ax3.set_title('Diagramma $\sigma_{(t)}$- $\epsilon_{(t)}$')
# ax3.set_xlim(min(e)-0.001, max(e)+max(e)/10)
# ax3.set_ylim(min(sigma)-20, max(sigma)+20)
ax3.set_xlabel('$\epsilon_{(t)}$ [-]')
ax3.set_ylabel('$\sigma_{(t)}$ $[N/mm^2]$')

#Plot f in funzione di t
axs[3, 0].plot(t,f,color='fuchsia',linewidth=1)
axs[3, 0].grid(True)
axs[3, 0].set_title('Funzione di snervamento')
axs[3, 0].set_xlim(min(t), max(t))
# axs[3, 0].set_ylim(min(f)-20, max(f)+20)
axs[3, 0].set_xlabel('$t [s]$')
axs[3, 0].set_ylabel('$f_{(t)} [-]$')

#Plot p in funzione di t
axs[2, 0].plot(t,p,'-k',linewidth=1)
axs[2, 0].grid(True)
axs[2, 0].set_title('deformazione plastica [p]')
axs[2, 0].set_xlim(min(t), max(t))
# axs[2, 0].set_ylim(min(p)-0.005, max(p)+0.005)
axs[2, 0].set_xlabel('$t [s]$')
axs[2, 0].set_ylabel('$p_{(t)} [-]$')
