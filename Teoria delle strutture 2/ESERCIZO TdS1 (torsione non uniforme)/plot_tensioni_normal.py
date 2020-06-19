## Tensioni normali per M2
plt.figure()
delta=20

# plt.fill_between(e1[0,:],sigma_M2_1[0,:])
# plt.fill_between(s1[0,:],sigma_M2_2[0,:])
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()

plt.subplots_adjust(top=0.88,
bottom=0.11,
left=0.05,
right=0.9,
hspace=0.45,
wspace=0.23)

plt.subplot(2,1,2)
plt.suptitle('Tensioni normali $\sigma_{z}^{M_2}$ [$KN/m^2$]',fontsize=15,color='darkred',fontstyle='italic')
plt.contour(Z, e1, sigma_M2_1, delta)
A1=plt.contourf(Z, e1, sigma_M2_1 , delta,cmap='jet_r')
plt.quiver(Z,e1,-sigma_M2_1*0.001,0,color='k',pivot='tip', units='inches',width=0.01)
plt.xlabel(' z  (m)')
plt.ylabel(' $x1=l/2 \,\,(m)$')
plt.title('$\sigma_{z}^{M_2}\,\,\,\, Tratto 1\,\,\,\, [KN/m^2]$')
plt.colorbar(A1)

plt.subplot(2,1,1)
plt.suptitle('Tensioni normali $\sigma_{z}^{M_2}$ [$KN/m^2$]',fontsize=15,color='darkred',fontstyle='italic')
plt.contour(Z, s1, sigma_M2_2, delta)
A2=plt.contourf(Z, s1, sigma_M2_2 , delta,cmap='jet_r')
plt.quiver(Z,s1,-sigma_M2_2*0.001,0,color='k',pivot='tip', units='inches',width=0.01)
plt.gca().invert_yaxis()
plt.xlabel(' z  (m)')
plt.ylabel(' $x1 \,\,(m)$')
plt.title('$\sigma_{z}^{M_2}\,\,\,\, Tratto 2\,\,\,\, [KN/m^2]$')
plt.colorbar(A2)

## Tensioni normali per Bimomento
plt.figure()

# plt.fill_between(s1[0,:],sigma_B_1[0,:])
# plt.fill_between(s2[0,:],sigma_B_2[0,:])
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()

plt.subplots_adjust(top=0.88,
bottom=0.11,
left=0.05,
right=0.9,
hspace=0.45,
wspace=0.23)

plt.subplot(2,1,2)
plt.suptitle('Tensioni normali $\sigma_{z}^{B}$ [$KN/m^2$]',fontsize=15,color='darkred',fontstyle='italic')
plt.contour(Z, s2, sigma_B_2, delta)
A1=plt.contourf(Z, s2, sigma_B_2 , delta,cmap='jet')
plt.quiver(Z,s2,-sigma_B_2*0.001,0,color='k',pivot='tip', units='inches',width=0.01)
plt.gca().invert_yaxis()
plt.xlabel(' z  (m)')
plt.ylabel(' $s2\,\,(m)$')
plt.title('$\sigma_{z}^{B}\,\,\,\, Tratto 1\,\,\,\, [KN/m^2]$')
plt.colorbar(A1)

plt.subplot(2,1,1)
plt.suptitle('Tensioni normali $\sigma_{z}^{B}$ [$KN/m^2$]',fontsize=15,color='darkred',fontstyle='italic')
plt.contour(Z, s1, sigma_B_1, delta)
A2=plt.contourf(Z, s1, sigma_B_1 , delta,cmap='jet')
plt.quiver(Z,s1,-sigma_B_1*0.001,0,color='k',pivot='tip', units='inches',width=0.01)
plt.gca().invert_yaxis()
plt.xlabel(' z  (m)')
plt.ylabel(' $s1\,\,(m)$')
plt.title('$\sigma_{z}^{B}\,\,\,\, Tratto 2\,\,\,\, [KN/m^2]$')
plt.colorbar(A2)

## Tensioni normali totali
plt.figure()

# plt.fill_between(s1[0,:],sigma_1[0,:])
# plt.fill_between(s2[0,:],sigma_2[0,:])
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()

plt.subplots_adjust(top=0.88,
bottom=0.11,
left=0.05,
right=0.9,
hspace=0.45,
wspace=0.23)

plt.subplot(2,1,2)
plt.suptitle('Tensioni normali totali $\sigma_{z}$ [$KN/m^2$]',fontsize=15,color='darkred',fontstyle='italic')
plt.contour(Z, s1, sigma_1, delta)
A1=plt.contourf(Z, s1, sigma_1 , delta,cmap='jet_r')
plt.quiver(Z,s1,-sigma_1*0.001,0,color='k',pivot='tip', units='inches',width=0.01)
plt.gca().invert_yaxis()
plt.xlabel(' z  (m)')
plt.ylabel(' $s1\,\,(m)$')
plt.title('$\sigma_{z}\,\,\,\, Tratto 1\,\,\,\, [KN/m^2]$')
plt.colorbar(A1)

plt.subplot(2,1,1)
plt.suptitle('Tensioni normali $\sigma_{z}$ [$KN/m^2$]',fontsize=15,color='darkred',fontstyle='italic')
plt.contour(Z, s2, sigma_2, delta)
A2=plt.contourf(Z, s2, sigma_2 , delta,cmap='jet_r')
plt.quiver(Z,s2,-sigma_2*0.001,0,color='k',pivot='tip', units='inches',width=0.01)
plt.gca().invert_yaxis()
plt.xlabel(' z  (m)')
plt.ylabel(' $s2\,\,(m)$')
plt.title('$\sigma_{z}\,\,\,\, Tratto 2\,\,\,\, [KN/m^2]$')
plt.colorbar(A2)

## PLOT tensioni a Z=0
fig, ax = plt.subplots(2, 6)
plt.subplots_adjust(top=0.83,
bottom=0.1,
left=0.045,
right=0.98,
hspace=0.0,
wspace=0.0)

plt.suptitle(' Tensioni Normali $[KN/m^2]$',fontsize=20,color='k',fontstyle='italic')

# tensioni dovuti a M2
ax[0,1].set_title('Tensioni normali dovuti a M2',fontsize=12,color='b',fontstyle='italic')
ax[0,0].set_axis_off()
ax[1,1].set_axis_off()

ax[1,0].fill_betweenx(e1[1,:],sigma_M2_1[1,:],facecolor='b', interpolate=True,alpha=0.5)# tratto 1
ax[1,0].plot(sigma_M2_1[1], e1[1,:] ,'-k', linewidth=1)
ax[1,0].plot(0*e1[1,:],e1[1,:] ,'-k', linewidth=2)

ax[0,1].plot(e2[1,:], sigma_M2_2[1],'-k', linewidth=1)
ax[0,1].plot(e2[1,:], 0*e2[1,:] ,'-k', linewidth=2)
ax[0,1].fill_between(e2[1], sigma_M2_2[1],facecolor='b', interpolate=True,alpha=0.5) # tratto 2
ax[0,1].invert_yaxis()

# Tensioni per B
ax[0,3].set_title('Tensioni normali dovuti a $Bw$',fontsize=12,color='y',fontstyle='italic')
ax[0,2].set_axis_off()
ax[1,3].set_axis_off()


ax[1,2].fill_betweenx(e1[1,:],sigma_B_2[1,:],facecolor='gold', interpolate=True,alpha=1)# tratto 1
ax[1,2].plot(sigma_B_2[1], e1[1,:] ,'-k', linewidth=1)
ax[1,2].plot(0*e1[1,:],e1[1,:] ,'-k', linewidth=2)
ax[1,2].invert_xaxis()



ax[0,3].plot(e2[1,:], sigma_B_1[1],'-k', linewidth=1)
ax[0,3].plot(e2[1,:], 0*e2[1,:] ,'-k', linewidth=2)
ax[0,3].fill_between(e2[1], sigma_B_1[1],facecolor='gold', interpolate=True,alpha=1) # tratto 2

# Tensioni totali
ax[0,5].set_title('Tensioni tangenziali totali',fontsize=12,color='b',fontstyle='italic')
ax[0,4].set_axis_off()
ax[1,5].set_axis_off()


ax[1,4].fill_betweenx(e1[1,:],sigma_1[1,:],facecolor='darkblue', interpolate=True,alpha=0.8)# tratto 1
ax[1,4].plot(sigma_1[1], e1[1,:] ,'-k', linewidth=1)
ax[1,4].plot(0*e1[1,:],e1[1,:] ,'-k', linewidth=2)

ax[0,5].plot(e2[1,:], sigma_2[1],'-k', linewidth=1)
ax[0,5].plot(e2[1,:], 0*e2[1,:] ,'-k', linewidth=2)
ax[0,5].fill_between(e2[1], sigma_2[1],facecolor='darkblue', interpolate=True,alpha=0.8) # tratto 2
ax[0,5].invert_yaxis()

for i in range(0,19) :
    kk=sigma_2[1,i]
    k=e2[1,i]
    ax[0,5].plot(k, kk,'-k', linewidth=1)


