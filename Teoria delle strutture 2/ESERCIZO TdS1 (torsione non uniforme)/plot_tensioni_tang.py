## Tensioni per T1
plt.figure()
delta=20

# plt.plot(e1[0,:],tau_T1_1[0,:])
# plt.plot(e2[0,:],tau_T1_2[0,:])

plt.subplots_adjust(top=0.88,
bottom=0.11,
left=0.05,
right=0.9,
hspace=0.45,
wspace=0.23)

plt.subplot(2,1,2)
plt.suptitle('Tensioni tangenziali $\Gamma_{T_1}$ [$KN/m^2$]',fontsize=15,color='darkred',fontstyle='italic')
plt.contour(Z, e1, tau_T1_1, delta)
A1=plt.contourf(Z, e1, tau_T1_1 , delta,cmap='jet')
plt.quiver(Z,e1,0,tau_T1_1*0.001,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' z  (m)')
plt.ylabel(' $\eta_1 \,\,(m)$')
plt.title('$\Gamma_{\eta z}\,\,\,\, Tratto 1\,\,\,\, [KN/m^2]$')
plt.colorbar(A1)

plt.subplot(2,1,1)
plt.contour(Z, e2, tau_T1_2, delta)
A2=plt.contourf(Z, e2, tau_T1_2 , delta,cmap='jet')
plt.quiver(Z,e2,0,tau_T1_2*0.001,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' z  (m)')
plt.ylabel('$\eta_2 \,\,(m)$')
plt.title('$\Gamma_{\eta z}\,\,\,\, Tratto 2\,\,\,\, [KN/m^2]$')
plt.colorbar(A2)

## Tensioni per Mts
plt.figure()


# plt.fill_between(e1[0,:],tau_Mts_2[0,:])
# plt.gca().invert_yaxis()
# plt.fill_between(e2[0,:],tau_Mts_1[0,:])
# plt.gca().invert_yaxis()

plt.subplot(2,1,2)
plt.suptitle('Tensioni tangenziali $\Gamma_{Mt^s}$ [$KN/m^2$]',fontsize=15,color='darkred',fontstyle='italic')
plt.contour(Z, e1, tau_Mts_2, delta)
A1=plt.contourf(Z, e1, tau_Mts_2, delta, cmap= ('jet'))
plt.quiver(Z,e1,0,-tau_Mts_2*0.001,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' z  (m)')
plt.ylabel(' $\eta_1 \,\,(m)$')
plt.title('$\Gamma_{sz}\,\,\,\, Tratto 1 \,\,\,\, [KN/m^2]$')
plt.colorbar(A1)

plt.subplot(2,1,1)
plt.contour(Z, e2, tau_Mts_1, delta)
A2=plt.contourf(Z, e2, tau_Mts_1, delta , cmap= ('jet'))
plt.quiver(Z,e2,0,-tau_Mts_1*0.001,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' z  (m)')
plt.ylabel('$\eta_2 \,\,(m)$')
plt.title('$\Gamma_{sz}\,\,\,\, Tratto 2\,\,\,\, [KN/m^2]$')
plt.colorbar(A2)

## Tensioni tangenziali totali
plt.figure()


# plt.fill_between(e1[0,:],tau_1[0,:])
# plt.fill_between(e2[0,:],tau_2[0,:])

plt.subplots_adjust(top=0.88,
bottom=0.11,
left=0.05,
right=0.9,
hspace=0.45,
wspace=0.23)

plt.subplot(2,1,2)
plt.suptitle('Tensioni tangenziali totali $\Gamma_{sz}$ [$KN/m^2$]',fontsize=15,color='darkred',fontstyle='italic')
plt.contour(Z, e1, tau_1, delta)
A1=plt.contourf(Z, e1, tau_1, delta, cmap= ('jet'))
plt.quiver(Z,e1,0,tau_1*0.001,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' z  (m)')
plt.ylabel(' $\eta_1 \,\,(m)$')
plt.title('$\Gamma_{sz}\,\,\,\, Tratto 1 \,\,\,\, [KN/m^2]$')
plt.colorbar(A1)

plt.subplot(2,1,1)
plt.contour(Z, e2, tau_2, delta)
A2=plt.contourf(Z, e2, tau_2, delta , cmap= ('jet'))
plt.quiver(Z,s2,0,tau_2*0.001,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' z  (m)')
plt.ylabel('$\eta_2 \,\,(m)$')
plt.title('$\Gamma_{sz}\,\,\,\, Tratto 2\,\,\,\, [KN/m^2]$')
plt.colorbar(A2)

## PLOT tensioni a Z=0
fig, ax = plt.subplots(2, 6)
plt.subplots_adjust(top=0.83,
bottom=0.1,
left=0.045,
right=0.98,
hspace=0.0,
wspace=0.0)

plt.suptitle(' Tensioni Tangenziali $[KN/m^2]$',fontsize=20,color='k',fontstyle='italic')

# tensioni dovuti a T1
ax[0,1].set_title('Tensioni tangenziali dovuti a T1',fontsize=12,color='darkred',fontstyle='italic')
ax[0,0].set_axis_off()
ax[1,1].set_axis_off()

ax[1,0].fill_betweenx(e1[1,:],tau_T1_1[1,:],facecolor='red', interpolate=True,alpha=0.5)# tratto 1
ax[1,0].plot(tau_T1_1[1], e1[1,:] ,'-k', linewidth=1)
ax[1,0].plot(0*e1[1,:],e1[1,:] ,'-k', linewidth=2)
ax[1,0].invert_xaxis()

ax[0,1].plot(e2[1,:], tau_T1_2[1],'-k', linewidth=1)
ax[0,1].plot(e2[1,:], 0*e2[1,:] ,'-k', linewidth=2)
ax[0,1].fill_between(e2[1], tau_T1_2[1],facecolor='red', interpolate=True,alpha=0.5) # tratto 2

# Tensioni per Mts
ax[0,3].set_title('Tensioni tangenziali dovuti a $Mt^s$',fontsize=12,color='y',fontstyle='italic')
ax[0,2].set_axis_off()
ax[1,3].set_axis_off()


ax[1,2].fill_betweenx(e1[1,:],tau_Mts_2[1,:],facecolor='y', interpolate=True,alpha=0.5)# tratto 1
ax[1,2].plot(tau_Mts_2[1], e1[1,:] ,'-k', linewidth=1)
ax[1,2].plot(0*e1[1,:],e1[1,:] ,'-k', linewidth=2)
ax[1,2].invert_xaxis()

ax[0,3].plot(e2[1,:], tau_Mts_1[1],'-k', linewidth=1)
ax[0,3].plot(e2[1,:], 0*e2[1,:] ,'-k', linewidth=2)
ax[0,3].fill_between(e2[1], tau_Mts_1[1],facecolor='y', interpolate=True,alpha=0.5) # tratto 2

# Tensioni totali
ax[0,5].set_title('Tensioni tangenziali totali',fontsize=12,color='g',fontstyle='italic')
ax[0,4].set_axis_off()
ax[1,5].set_axis_off()


ax[1,4].fill_betweenx(e1[1,:],tau_1[1,:],facecolor='darkred', interpolate=True,alpha=0.5)# tratto 1
ax[1,4].plot(tau_1[1], e1[1,:] ,'-k', linewidth=1)
ax[1,4].plot(0*e1[1,:],e1[1,:] ,'-k', linewidth=2)
ax[1,4].invert_xaxis()

ax[0,5].plot(e2[1,:], tau_2[1],'-k', linewidth=1)
ax[0,5].plot(e2[1,:], 0*e2[1,:] ,'-k', linewidth=2)
ax[0,5].fill_between(e2[1], tau_2[1],facecolor='darkred', interpolate=True,alpha=0.5) # tratto 2
