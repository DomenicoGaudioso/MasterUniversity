import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

n=1
b=1
yy = np.linspace(0, b)

def w(x, y):
    return np.sin(np.pi*m*x/l)*np.sin(np.pi*n*y/b)

plt.ion()
plt.show()

fig = plt.figure(figsize=plt.figaspect(0.5))

d=3
c=3
A=c*d


for x in range(1,A+1):
    l=1*b
    m=x*l/b
    ax = fig.add_subplot(d, c, x, projection='3d')
    xx = np.linspace(0, l)

    X, Y = np.meshgrid(xx, yy)
    W=w(X,Y)


    ax.plot_surface(X, Y, W, rstride=2, cstride=2,
                    cmap=plt.cm.jet, edgecolor='k')
    ax.set_axis_off()

    ax.set_title(m,fontsize=10,color='darkred',fontstyle='italic')

plt.suptitle('Deformata instabile di una piastra (al variare di m) con carico lungo y in direzione x ',fontsize=20,color='darkred',fontstyle='italic')

c=1000


def k(m,l):
    return (m*b/l+l/(m*b))**2

plt.figure()
plt.ylim(0, 8)
plt.grid(True)
plt.title('Curve instabilità',fontsize=20,color='darkred',fontstyle='italic')
plt.xlabel('l/b ')
plt.ylabel('k')

for i in range(1,A+1):
    m=i
    l=np.linspace(0.1*b,6*b,c)
    kn=k(m,l)

    plt.plot(l/b,kn,'-',label= m,linewidth=1)
    plt.xlim(0, max(l/b))
    plt.legend(title='m:',loc='lower right')




