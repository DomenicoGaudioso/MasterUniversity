import numpy as np
import math
geom=np.array([[0.0,0.0],[4.0,4.0],[8.0,0.0]]) #nodes coordinates x and y
conn=np.array([0,1],[1,2]]) # element connectivity

nnd=cord.shape[0] or size(a, axis=1) #number of nodes
nel=conn.shape[0] or size(a, axis=1) #number of elements
nne=2 #number of nodes per element
nodof=3 #number of degrees of freedom per node (2 traslazioni e una rotazione)

E=np.zeros(nel) # Young's modulus [KN/m^2]
E[0]=21E7 #steel
E[1]=21E7 #steel

A=np.zeros(nel) #Area section [m^2]
A[0]=0.3*0.6
A[1]=0.3*0.6

I=np.zeros(nel) # inertia [m^3]
I[0]=0.3*(0.6**3)/12 #rectangle
I[1]=0.3*(0.6**3)/12 #..

#BOUNDARY CONDITIONS 
nf=np.zeros((nnd,nodof)) #tutti gli elementi sono incastrati impostando ''zero sia sugli spost. che sulle rotazioni''
# se voglio un nodo a cerniera basta fare nf[0,]=1 nella riga dove voglio la cerniera

L=np.zeros(nel)#element length
Theta=np.zeros(nel)# angle formed by the element with the horizontal
KK=np.zeros((nne*nnd,nne*nnd)) #dimensions of the global stiffness matrix
for e in range(0,nel):
    node_1=conn[e,0] #node start
    node_2=conn[e,1] #node end

    x1=geom[node_1,0]; x2=geom[node_2,0]
    y1=geom[node_1,1]; y2=geom[node_1,1]
   
    
    L=math.sqrt((x2-x1)**2+(y2-y1)**2) #element length
    Theta[e]=math.atan2((j[1]-i[1]),(j[0]-i[0])) # angle formed by the element with the horizontal
    
    Q=np.matrix([[math.cos(Theta[e]),-math.sin(Theta[e]), 0],
                 [math.sin(Theta[e]),math.cos(Theta[e]),0],
                 [0,0,1]])# transfer matrix
    
    T=np.matrix([[math.cos(Theta[e]),-math.sin(Theta[e]), 0,0,0,0],
                 [math.sin(Theta[e]), math.cos(Theta[e]), 0,0,0,0],
                 [0,0,1,0,0,0],
                 [0,0,0,math.cos(Theta[e]),-math.sin(Theta[e]),0],
                 [0,0,0,math.sin(Theta[e]),math.cos(Theta[e]),0],
                 [0,0,0,0,0,1]])
    

  # STIFFNESS MATRIX IN LOCAL
    Kl=np.matrix([[(A[e]*E[e]/L),0,0,-(A[e]*E[e]/L),0,0],
                        [0,(12*E[e]*I[e]/(L**3)),(6*E[e]*I[e]/(L**2)),0,-(12*E[e]*I[e]/(L**3)),(6*E[e]*I[e]/(L**2))],
                        [0,(6*E[e]*I[e]/(L**2)),(4*E[e]*I[e]/(L)),0,-(6*E[e]*I[e]/(L**2)),(2*E[e]*I[e]/(L))],
                        [-(A[e]*E[e]/L),0,0,(A[e]*E[e]/L),0,0],
                        [0,-(12*E[e]*I[e]/(L**3)),-(6*E[e]*I[e]/(L**2)),0,(12*E[e]*I[e]/(L**3)),-(6*E[e]*I[e]/(L**2))],
                        [0,(6*E[e]*I[e]/(L**2)),(2*E[e]*I[e]/(L)),0,-(6*E[e]*I[e]/(L**2)),(4*E[e]*I[e]/(L))]])

    # STIFFNESS MATRIX IN GLOBAL
    Kg=T*Kl*T.transpose()

    #STEERING VECTOR
# è il vettore guida contenente il numero di gradi di libertà dell'elemento
#g=np.array([nf[i,0],nf[i,1],nf[i,2],nf[j,0],nf[j,1],nf[j,2]])
#KK[range((1*i),(1*i+1)),range((1*i),(1*i+1))]=KK[range((1*i),(1*i+1)),range((1*i),(1*i+1))]+Kg[range((0),(1)),range((0),(1))]

print((conn[0,0]).ravel().nonzero())
