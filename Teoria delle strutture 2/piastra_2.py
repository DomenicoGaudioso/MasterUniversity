import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

b=1
p2=5 #[KN/m]
D=1
a=1

yy = np.linspace(0, b)


yy = np.linspace(0, b)

# F=(m*sy.pi/l)**2+(n*sy.pi/b)**2-sy.sqrt(p1*(m*sy.pi/l)**2/D+p2*(n*sy.pi/b)**2/D)

def w(x,y):
    return np.sin(np.pi*m*x/l)*np.sin(np.pi*n*y/b)

def p(n,m) :
    # return D*np.pi**2*(((m/l)**2+(n/b)**2)**2)/((m/l)**2)-p2*((n/b)**2)/((m/l)**2)
        return D*np.pi**2*(((m/l)**2+(n/b)**2)**2)/((m/l)**2+a*(n/b)**2)

plt.ion()
plt.show()


plt.figure()
plt.ylim(0, 80)
plt.grid(True)
plt.title('Curve instabilità',fontsize=20,color='darkred',fontstyle='italic')
plt.xlabel('l/b ')
plt.ylabel('p [KN/m]')

c=1000

for i in range(1,2):
    n=i
    for j in range(1,6):
        m=j
        l=np.linspace(0.1*b,6*b,c)
        pn=p(n,m)

        plt.plot(l/b,pn,'-',label= (n,m)  ,linewidth=1)
        plt.xlim(0, max(l/b))
        plt.legend(title='n,m:',loc='lower right')

