import numpy as np
import matplotlib.pyplot as plt

## INPUT ##
L = 10
s1 =  L*0.0
P = 1

def Bending( s, xp, L, P ):
    Va = P*(L - xp)/L
    Vb = P*xp/L
    if xp >= s:
        return Va*s
    elif xp < s:
        return Vb*( L - s )

# LINEA D'INFLUENZA DI s1 ##
xp = np.linspace( 0, L, int(L/0.1) )

def LdI( s, xp, L, p ):
    Msi = []
    for xpi in xp:
        Msi.append( Bending( s, xpi, L, P ) )
    return Msi

Ms1 = LdI(s1, xp, L, P )

# DIAGRAMMA DEI MASSIMI E DEI MINIMI ##
xs = np.linspace( 0, L,int(L/0.1))

def diagrammaMaxMin( s, xp, L, P ):
    maxMsi = []
    minMsi = []
    for xsi in s:  
        Msi = LdI(xsi, xp, L, P )
        maxMsi.append( max(Msi) )
        minMsi.append( min(Msi) )
        
    return  maxMsi, minMsi 

Max, Min = diagrammaMaxMin( xs, xp, L, P )   

## PLOT ##
plt.ion()
plt.show()

## Plot linea d'influenza di s1 ##    
plt.figure()
plt.plot(xp, Ms1)
s_str = str( s1 )
plt.title('Linea di influenza (LdI) Momento flettente in s1 = '+ s_str + ' m', fontsize=15,color='darkred',fontstyle='italic' )
plt.xlabel('xp  (m)')
plt.ylabel(' M (KNm)')
plt.xlim(0, L)
plt.ylim(min(Ms1)- 0.1, max(Ms1)+0.1)
plt.gca().invert_yaxis()
plt.grid()

## Plot diagramma di max e minimi ##    
plt.figure()
plt.plot(xs, Max, Min)

plt.title('Diagramma dei Massimi e dei Minimi del Momento flettente ', fontsize=15,color='darkred',fontstyle='italic' )
plt.xlabel('xs  (m)')
plt.ylabel(' M (KNm)')
plt.xlim(0, L)
plt.ylim(min(Max)- 0.1, max(Max)+0.1)
plt.gca().invert_yaxis()
plt.grid()
