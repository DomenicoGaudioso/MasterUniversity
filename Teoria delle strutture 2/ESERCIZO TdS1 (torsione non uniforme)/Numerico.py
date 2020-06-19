exec (open ("LE_Bending.py") .read ())
exec (open ("LE_Torsion.py") .read ())
exec (open ("Prop_sez.py") .read ())

# Problema Numerico

p1n=1.2*16.5 # KN/m
q=I0/It
print( 'cazzo =' , q)
an=np.sqrt(Gn*I0/(En*Iw*q))
mn=p1n*4.25

zn=np.linspace( 0.00, Hn , 20) # ascissa

v1n=np.zeros(len(zn)) # spostamento in direzione 1
r2n=np.zeros(len(zn)) # rotazione in direzione 2
T1n=np.zeros(len(zn)) # Taglio in direzione 1
M2n=np.zeros(len(zn)) # Momento in direzione 2

wn=np.zeros(len(zn)) # rotazione torzionale
dwn=np.zeros(len(zn)) # ingobbamento
Bwn=np.zeros(len(zn)) # Bimomento
Mtsn=np.zeros(len(zn)) # Momento torcente secondario
Mtpn=np.zeros(len(zn)) # Momento torcente primario
Mtn=np.zeros(len(zn)) # Momento torcente totale

for i in range(0,(len(zn))):
    v1n[i]=np.float(v1.subs(z,zn[i]).subs(E,En).subs(J2,I2).subs(p1,p1n).subs(H,Hn))
    r2n[i]=np.float(r2.subs(z,zn[i]).subs(E,En).subs(J2,I2).subs(p1,p1n).subs(H,Hn))
    M2n[i]=np.float(M2.subs(z,zn[i]).subs(E,En).subs(J2,I2).subs(p1,p1n).subs(H,Hn))
    T1n[i]=np.float(T1.subs(z,zn[i]).subs(E,En).subs(J2,I2).subs(p1,p1n).subs(H,Hn))
    wn[i]=np.float(w.subs(z,zn[i]).subs(E,En).subs(Jw,Iw).subs(a,an).subs(H,Hn).subs(m,mn))
    dwn[i]=np.float(dw.subs(z,zn[i]).subs(E,En).subs(Jw,Iw).subs(a,an).subs(H,Hn).subs(m,mn))
    Mtsn[i]=np.float(Mts.subs(z,zn[i]).subs(E,En).subs(Jw,Iw).subs(a,an).subs(H,Hn).subs(m,mn))
    Mtn[i]=np.float(Mt.subs(z,zn[i]).subs(E,En).subs(Jw,Iw).subs(a,an).subs(H,Hn).subs(m,mn))
    Bwn[i]=np.float(Bw.subs(z,zn[i]).subs(E,En).subs(Jw,Iw).subs(a,an).subs(H,Hn).subs(m,mn))
    # Mtpn[i]=Gn*I0*dwn[i]/q
    Mtpn[i]=Mtn[i]-Mtsn[i]


print('v1{max}=', round(max(v1n)*1000,3), 'mm' )
print('spostamento massimo consentito=', 'H/400=',Hn*1000/400, 'mm' )
print('M2{max}=', round(max(M2n),3),'KN*m' )
print('T1{max}=', round(min(T1n),3),'KN')
print('w{max}=', round(max(wn),3),'1/m')
print('ingobbamento{max}=', round(min(dwn),3),'1/m')
print('Mtp{max}=', round(max(Mtpn),3),'KN*m' )
print('Mts{max}=', round(max(Mtsn),3),'KN*m')
print('Mtn{max}=', round(max(Mtn),3),'KN*m')

## plot
exec (open ("Bending_plot.py") .read ())
exec (open ("Torsion_plot.py") .read ())