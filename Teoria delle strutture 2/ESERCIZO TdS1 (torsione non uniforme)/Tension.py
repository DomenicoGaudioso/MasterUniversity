exec (open ("Numerico.py") .read ())
from mpl_toolkits import mplot3d

## MESHGRID
ee1 = (np.linspace (0,b, num=len(zn), endpoint=True))
ee2 = (np.linspace (0,l/2, num=len(zn), endpoint=True))


Z , T1nn = np.meshgrid(zn,T1n)
Z , Mtsnn = np.meshgrid(zn,Mtsn)
Z , Bwnn = np.meshgrid(zn,Bwn)
Z , M2nn = np.meshgrid(zn,M2n)
e1, Z = np.meshgrid(ee1,zn)
e2, Z = np.meshgrid(ee2,zn)

s2=b-e1
s1=l/2-e2
#x1 = (np.linspace (-l/2,0, num=92, endpoint=True))

## funzione ingobamento

psi_1 =-s1*(b+d/2)
psi_2 =-(b+d/2)*l/2-s2*l/2

## Calcolo momenti statici
# Momento statico S1
S1_1 = t*e1*(-(e1+d)*0.5)
S1_2 = t*b*(-(b+d)/2)+e2*t*(-(b+d))

# Momento statico S2
S2_1 = t*e1*l/2
S2_2 = t*b*l/2+t*e2*(l/2-e2/2)

# Momento statico settoriale
Spsi_1 = -t*l/2*((b+d/2)*b+b**2/2)-t*(b+d/2)*(l**2/8-s1**2/2)
Spsi_2 = -t*l/2*((b+d/2)*(b-s2)+(b**2/2-s2**2/2))

## Calcolo tensioni tangenziali
# tensioni tangenziali dovuti a (T1)
tau_T1_1 = - T1nn*S2_1/(I2*t)
tau_T1_2 = - T1nn*S2_2/(I2*t)

# tensioni tangenziali dovuti a (Mts)
tau_Mts_1 =  Mtsnn*Spsi_1/(Iw*t)
tau_Mts_2 =  Mtsnn*Spsi_2/(Iw*t)

# tensioni tangenziali totali
tau_1=tau_T1_1+tau_Mts_2
tau_2=tau_T1_2+tau_Mts_1

## Calcolo tensioni normali
# tensioni normali dovuti a (M2)
sigma_M2_1 = -M2nn/I2*l/2
sigma_M2_2 = -M2nn/I2*s1

# tensioni normali dovuti a (B)
sigma_B_1 = Bwnn*psi_1/Iw
sigma_B_2 = Bwnn*psi_2/Iw

# tensioni normali totali
sigma_1 = sigma_M2_1 + sigma_B_2
sigma_2 = sigma_M2_2 + sigma_B_1

### PLOT
exec (open ("plot_tensioni_tang.py") .read ())
exec (open ("plot_tensioni_normal.py") .read ())

