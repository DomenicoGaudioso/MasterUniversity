## plot Linea elastica Bending
plt.ion()
plt.show ()
fig, axs = plt.subplots(2, 2)
plt.subplots_adjust(top=0.885,
bottom=0.075,
left=0.07,
right=0.99,
hspace=0.515,
wspace=0.26)
plt.suptitle('FLESSIONE SEMPLICE    $(EJ_{2}\cdot a_1^{IV}-p_1=0)$',fontsize=20,color='darkred',fontstyle='italic')

# Plot spostamento

axs[0, 0].plot(zn,v1n*1000,'-g',linewidth=1)
axs[0, 0].grid(True)
axs[0, 0].set_title('Spostamento v1')
axs[0, 0].set_xlim(min(zn), max(zn))
# axs[0, 0].set_ylim(min(v1n)*1000, max(v1n)*1000)
axs[0, 0].set_xlabel('z [m]')
axs[0, 0].set_ylabel('$v1_{(z)}$ [mm]')
axs[0, 0].invert_yaxis()

#Plot rotazione

axs[0, 1].plot(zn,np.degrees(r2n),'-b',linewidth=1)
axs[0, 1].grid(True)
axs[0, 1].set_title('rotazione $\phi_{2} $')
axs[0, 1].set_xlim(min(zn), max(zn))
# axs[0, 1].set_ylim(min(r2n), max(r2n))
axs[0, 1].set_xlabel('$z [m]$')
axs[0, 1].set_ylabel('$\phi_{2}$ [°]')
axs[0, 1].invert_yaxis()

# axs[1, 0].plot(zn,T1n,linewidth=1)
axs[1, 0].fill_between(zn,T1n,facecolor='lightsalmon', edgecolor='k', linewidth=1)
axs[1, 0].grid(True)
axs[1, 0].set_title('Taglio T1')
axs[1, 0].set_xlim(min(zn), max(zn))
# axs[1, 0].set_ylim(min(T1n), max(T1n))
axs[1, 0].set_xlabel('z [m]')
axs[1, 0].set_ylabel('T2(z) [KN]')

# Plot Momento

axs[1, 1].fill_between(zn,M2n,facecolor='skyblue', edgecolor='k', linewidth=1)
axs[1, 1].grid(True)
axs[1, 1].set_title('Momento M2')
axs[1, 1].set_xlim(min(zn), max(zn))
# axs[1, 1].set_ylim(min(M2n), max(M2n))
axs[1, 1].set_xlabel('$z [m]$')
axs[1, 1].set_ylabel('$M2_{(z)} [KN*m]$')