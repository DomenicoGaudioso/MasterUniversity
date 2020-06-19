## plot Linea elastica Torsion
plt.ion()
plt.show ()
fig, axs = plt.subplots(2, 2)
plt.subplots_adjust(top=0.885,
bottom=0.075,
left=0.095,
right=0.99,
hspace=0.515,
wspace=0.26)
plt.suptitle('TORSIONE MISTA $\omega^{IV}-a^2\cdot\omega^{II}={m}/({E\cdot J_{\psi}})$',fontsize=20,color='darkred',fontstyle='italic')

# Plot Momento torcente
axs[0, 0].plot(zn,Mtn, color='red', linewidth=1,label='Mt totale')
axs[0, 0].plot(zn,Mtsn, color='m', linewidth=1,label='Mt secondario')
axs[0, 0].plot(zn,Mtpn, color='k', linewidth=1,label='Mt primario')
axs[0, 0].grid(True)
axs[0, 0].set_title('Momento torcente a z=H')
axs[0, 0].set_xlim(min(zn), max(zn))
# axs[1, 0].set_ylim(min(Mtn), max(Mtn))
axs[0, 0].set_xlabel('z [m]')
axs[0, 0].set_ylabel('Mt (z) [KN*m]')
axs[0, 0].legend(loc='upper right')


# Plot rotazione

axs[0, 1].plot(zn,np.degrees(wn),'-b',linewidth=1)
axs[0, 1].grid(True)
axs[0, 1].set_title('rotazione $\omega$')
axs[0, 1].set_xlim(min(zn), max(zn))
# axs[0, 1].set_ylim(min(dwn), max(dwn))
axs[0, 1].set_xlabel('$z [m]$')
axs[0, 1].set_ylabel('$\omega$ [°] ')
axs[0, 1].invert_yaxis()


# Plot zoom momento torcente

axs[1, 0].plot(zn,Mtn, color='red', linewidth=1,label='Mt totale')
axs[1, 0].plot(zn,Mtsn, color='m', linewidth=1,label='Mt secondario')
axs[1, 0].plot(zn,Mtpn, color='k', linewidth=1,label='Mt primario')
axs[1, 0].grid(True)
axs[1, 0].set_title('Momento torcente')
axs[1, 0].set_xlim(min(zn), max(zn))
# axs[1, 0].set_ylim(min(Mtn), max(Mtn))
axs[1, 0].set_xlabel('z [m]')
axs[1, 0].set_ylabel('Mt (z) [KN*m]')
axs[1, 0].legend(loc='upper right')




# Plot Bimomento

axs[1, 1].fill_between(zn,Bwn,facecolor='greenyellow', edgecolor='k', linewidth=1)
axs[1, 1].grid(True)
axs[1, 1].set_title('Bimomento $B_{\omega}$')
axs[1, 1].set_xlim(min(zn), max(zn))
# axs[1, 1].set_ylim(min(Bwn), max(Bwn))
axs[1, 1].set_xlabel('$z [m]$')
axs[1, 1].set_ylabel('$B_{\omega} [KN]$')