## Dati sezione

B=2.45          # [m]
L=5.8           # [m]
d=3.10          # [m]
t=0.30          # [m]
Hn=91.5         # [m]

b=B-t/2         # [m]
l=L-t           # [m]

## Prop Materiale

E1 = 25000 # MPa
G1 = 10800 # Mpa
En=E1*1000 # KN/m2
Gn=G1*1000 # KN/m2

## Area

An=6*(B-t)+2*L

## Momenti d'inerzia

# [m**4] Momento d'inerzia asse 1
I1=L*(2*B+d)**3/12-2*(l/2-t)*(2*(B-t)+d)**3/12-3*t*d**3/12

#I1=(L*t**3/12+L*t*(B-t/2+d/2)**2)*2+6*((B-t)**3*t/12+(B-t)*t*((B-t)/2+d/2)**2)

# [m**4] Momento d'inerzia asse 2
I2=2*t*L**3/12+2*t**3*(B-t)/12+4*t**3*(B-t)/12+4*t*(B-t)*((L-t)/2)**2

# [m**4] Momento d'inerzia polare
I0=I1+I2

# [m**4] Momento d'inerzia settoriale
Iw=t*(b+d/2)**2*l**3/6+t*l**2*((b+d/2)**2*b+b**3/3+b**2*(b+d/2))
print(Iw)

# [m**4] momento d'inerzia torsionale
It=6*(b*t**3)/3+2*(l*t**3)/3

print('A=', round(An, 3), 'm^2')
print('I1=', round(I1, 3), 'm^4')
print('I2=', round(I2, 3), 'm^4')
print('I0=', round(I0, 3), 'm^4')
print('Iw=', round(Iw, 3), 'm^4')
print('It=', round(It, 3), 'm^4')

