%%  Computational Mechanics - TERZA ESERCITAZIONE
% risoluzione di elementi triangolari a 6 nodi [Funzioni di forma
% quadratiche]
%% Clear workspace and close any open windows
clear all
close all
%% Set Command Window output display format
format bank;
%% data of the problem                      
Matricola=506682;             % my university number
H=(Matricola /100000);       % [m] length
b=1;                               % [m] length
B=3*b;                            % [m] length
vu=0.2;                           % Poisson's ratio
rho= 2400/1000;               % density [kg/m3]  diviso 1000  perchè poi voglio le forze di volume in KN
E1=3E7;                          % Young's modulus concrete [KN/m^2]
t=1;                               % [m] Thickness
%% data structure
nnd =25;                 % Number of nodes in system
nel = 8;                   % Number of elements
nne = 6;                  % Number of nodes for element
ndof =2;                  % Number of degrees of freedom for node
sdof = nnd*ndof;      % Number of degrees of freedom for element
%% Nodes cFoordinates X and Y
NX=3;
NY=2;

dx=B/(2*NX); % discretizzazione lungo X
dy=H/(2*NY); % discretizzazione lungo Y

i=0:dy:H; % discretizzazione lungo Y
j=0:dx:B; % discretizzazione lungo X
%% nodi
 e=1:5;
 y(e)=i;    x(e)=j(1);
 
 e=6:10;   
 y(e)=i;    x(e)=j(2);
 
 e=11:15;   
 y(e)=i;    x(e)=j(3);
 
 e=16:19;
 i=0:dy:(H-dy);
 y(e)=i;    x(e)=j(4);
 
 e=20:22;
  i=0:dy:(H-2*dy);
 y(e)=i;    x(e)=j(5);
 
 e=23:24;
 i=0:dy:(H-3*dy);
 y(e)=i;    x(e)=j(6);
 
 e=25;
 i=0:dy:(H-4*dy);
 y(e)=i;   x(e)=j(7);
%% Element connectivity
e=1:2; % elementi 
I=1:2:3;
J=11:2:13;
K=3:2:5;
L=2:2:4;
M=6:2:8;
N=7:2:9;
conn(e,1)=I;  conn(e,2)=M;  conn(e,3)=J; conn(e,4)=N;  conn(e,5)=K;  conn(e,6)=L; 

e=3:4; % elementi
I=13:2:15;
J=3:2:5;
K=11:2:13;
L=12:2:14;
M=8:2:10;
N=7:2:9;
conn(e,1)=I;  conn(e,2)=M;  conn(e,3)=J; conn(e,4)=N;  conn(e,5)=K;  conn(e,6)=L; 

e=5:6; % elementi 
I=11:2:13;
J=20:2:22;
K=13:2:15;
L=12:2:14;
M=16:2:18;
N=17:2:19;
conn(e,1)=I;  conn(e,2)=M;  conn(e,3)=J; conn(e,4)=N;  conn(e,5)=K;  conn(e,6)=L; 

e=7; % elementi
I=22;
J=13;
K=20;
L=21;
M=18;
N=17;
conn(e,1)=I;  conn(e,2)=M;  conn(e,3)=J; conn(e,4)=N;  conn(e,5)=K;  conn(e,6)=L; 

e=8; % elementi 
I=20;
J=25;
K=22;
L=21;
M=23;
N=24;
conn(e,1)=I;  conn(e,2)=M;  conn(e,3)=J; conn(e,4)=N;  conn(e,5)=K;  conn(e,6)=L;  

for i=1:nel
I=conn(i,1); 
M=conn(i,2);
J=conn(i,3);
N=conn(i,4); 
K=conn(i,5);
L=conn(i,6);
xx1=[x(L);x(I)];
yy1=[y(L);y(I)];
S1(i)=plot(xx1,yy1,'K-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','k');
hold on
xx=[x(I);x(M);x(J);x(N);x(K);x(L)];
yy=[y(I);y(M);y(J);y(N);y(K);y(L)];
S(i)=plot(xx,yy,'K-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','k');
axis([-0.1 1 -0.06  5.8]);
axis equal
axis off
hold on
end
%% Elasticity matrix
E=E1/((1+vu)*(1-2*vu))*[1-vu  -vu   0.  ;...  %Matrix elastic %stato piano di deformazione
                                     -vu  1-vu   0.  ;...
                                      0      0 (1-2*vu)/2];                               
%% matrix stiffnes composition              
 syms   X Y Xi Etha xI yI xJ yJ xK yK xL yL xM yM xN yN  real
%% Shape functions
lamda=1-Xi-Etha;

NI(Xi,Etha)=-lamda*(1-2*lamda);
NM(Xi,Etha)=4*Xi*lamda;
NJ(Xi,Etha)=-Xi*(1-2*Xi);
NN(Xi,Etha)=4*Xi*Etha;
NK(Xi,Etha)=-Etha*(1-2*Etha);
NL(Xi,Etha)=4*Etha*lamda;
% 
X=NI(Xi,Etha)*xI+NM(Xi,Etha)*xM+NJ(Xi,Etha)*xJ+NN(Xi,Etha)*xN+NK(Xi,Etha)*xK+NL(Xi,Etha)*xL;
Y=NI(Xi,Etha)*yI+NM(Xi,Etha)*yM+NJ(Xi,Etha)*yJ+NN(Xi,Etha)*yN+NK(Xi,Etha)*yK+NL(Xi,Etha)*yL;
%
Jacob(Xi,Etha)=[diff(X,Xi) diff(X,Etha);diff(Y,Xi) diff(Y,Etha)]; %Jacobiano
%
Nglobal(Xi,Etha)=[NI(Xi,Etha) 0 NM(Xi,Etha) 0 NJ(Xi,Etha) 0 NN(Xi,Etha) 0 NK(Xi,Etha) 0 NL(Xi,Etha) 0;...
                   0 NI(Xi,Etha) 0 NM(Xi,Etha) 0 NJ(Xi,Etha) 0 NN(Xi,Etha) 0 NK(Xi,Etha) 0 NL(Xi,Etha)];
%% Strain-displacement matrices
Diff_I=[diff(NI(Xi,Etha),Xi) diff(NI(Xi,Etha),Etha)]';
Diff_M=[diff(NM(Xi,Etha),Xi) diff(NM(Xi,Etha),Etha)]';
Diff_J=[diff(NJ(Xi,Etha),Xi) diff(NJ(Xi,Etha),Etha)]';
Diff_N=[diff(NN(Xi,Etha),Xi) diff(NN(Xi,Etha),Etha)]';
Diff_K=[diff(NK(Xi,Etha),Xi) diff(NK(Xi,Etha),Etha)]';
Diff_L=[diff(NL(Xi,Etha),Xi) diff(NL(Xi,Etha),Etha)]';

diffN1=inv(Jacob(Xi,Etha))*Diff_I;
diffN2=inv(Jacob(Xi,Etha))*Diff_M;
diffN3=inv(Jacob(Xi,Etha))*Diff_J;
diffN4=inv(Jacob(Xi,Etha))*Diff_N;
diffN5=inv(Jacob(Xi,Etha))*Diff_K;
diffN6=inv(Jacob(Xi,Etha))*Diff_L;

BI=[diffN1(1) 0;... 
       0 diffN1(2);... 
         diffN1(2) diffN1(1)]; 

 BM=[diffN2(1) 0;... 
         0 diffN2(2);... 
         diffN2(2) diffN2(1)]; 
     
  BJ=[diffN3(1) 0;... 
         0 diffN3(2);... 
     diffN3(2) diffN2(1)]; 
 
   BN=[diffN4(1) 0;... 
         0 diffN4(2);... 
     diffN4(2) diffN4(1)]; 
 
   BK=[diffN5(1) 0;... 
         0 diffN5(2);... 
     diffN5(2) diffN5(1)]; 
 
   BL=[diffN6(1) 0;... 
         0 diffN6(2);... 
     diffN6(2) diffN6(1)]; 
 
 Bsym=[BI BM BJ BN BK BL];
 %% Elasticity matrix
E=E1/((1+vu)*(1-2*vu))*[1-vu  -vu   0.  ;...  %Matrix elastic %stato piano di deformazione
                                     -vu  1-vu   0.  ;...
                                      0     0 (1-2*vu)/2];   
 %% passaggio da sym a numerico
  for e=1:nel
        I=conn(e,1); 
        M=conn(e,2);
        J=conn(e,3);
        N=conn(e,4); 
        K=conn(e,5);
        L=conn(e,6);

 Jac(:,:,e)=subs(Jacob(Xi,Etha),[xI yI xM yM xJ yJ xN yN xK yK xL yL],[x(I) y(I) x(M) y(M) x(J) y(J) x(N) y(N) x(K) y(K) x(L) y(L)]);      
 BB(:,:,e)=subs(Bsym,[xI yI xM yM xJ yJ xN yN xK yK xL yL],[x(I) y(I) x(M) y(M) x(J) y(J) x(N) y(N) x(K) y(K) x(L) y(L)]);
 
 Ke(:,:,e)=int(int(BB(:,:,e)'*E*BB(:,:,e)*det(Jac(:,:,e)),Etha,[0 1-Xi]),Xi,[0 1]);
 KK(:,:,e)=double(Ke(:,:,e));
  end

  %% Global stiffness matrix
Kg=zeros(2*nnd, 2*nnd);
for e = 1:nel
    I=conn(e,1); 
    M=conn(e,2);
    J=conn(e,3);
    N=conn(e,4); 
    K=conn(e,5);
    L=conn(e,6);
    
 % transfer of the submatrices to their position in the stiffness matrix global
    Kg(2*I-1:2*I,2*I-1:2*I) = Kg(2*I-1:2*I,2*I-1:2*I) + KK(1:2,1:2,e);
    Kg(2*I-1:2*I,2*M-1:2*M) = Kg(2*I-1:2*I,2*M-1:2*M) + KK(1:2,3:4,e);
    Kg(2*I-1:2*I,2*J-1:2*J) = Kg(2*I-1:2*I,2*J-1:2*J)+ KK(1:2,5:6,e);
    Kg(2*I-1:2*I,2*N-1:2*N) = Kg(2*I-1:2*I,2*N-1:2*N)+ KK(1:2,7:8,e);    
    Kg(2*I-1:2*I,2*K-1:2*K) = Kg(2*I-1:2*I,2*K-1:2*K) + KK(1:2,9:10,e);        
    Kg(2*I-1:2*I,2*L-1:2*L) =  Kg(2*I-1:2*I,2*L-1:2*L)+ KK(1:2,11:12,e); 
    
    Kg(2*M-1:2*M,2*I-1:2*I) = Kg(2*M-1:2*M,2*I-1:2*I) + KK(3:4,1:2,e);
    Kg(2*M-1:2*M,2*M-1:2*M) = Kg(2*M-1:2*M,2*M-1:2*M) + KK(3:4,3:4,e);
    Kg(2*M-1:2*M,2*J-1:2*J) = Kg(2*M-1:2*M,2*J-1:2*J)+ KK(3:4,5:6,e);
    Kg(2*M-1:2*M,2*N-1:2*N) = Kg(2*M-1:2*M,2*N-1:2*N)+ KK(3:4,7:8,e);    
    Kg(2*M-1:2*M,2*K-1:2*K) = Kg(2*M-1:2*M,2*K-1:2*K) + KK(3:4,9:10,e);        
    Kg(2*M-1:2*M,2*L-1:2*L) =  Kg(2*M-1:2*M,2*L-1:2*L)+ KK(3:4,11:12,e);   
    
    Kg(2*J-1:2*J,2*I-1:2*I) = Kg(2*J-1:2*J,2*I-1:2*I) + KK(5:6,1:2,e);
    Kg(2*J-1:2*J,2*M-1:2*M) = Kg(2*J-1:2*J,2*M-1:2*M) + KK(5:6,3:4,e);
    Kg(2*J-1:2*J,2*J-1:2*J) = Kg(2*J-1:2*J,2*J-1:2*J)+ KK(5:6,5:6,e);
    Kg(2*J-1:2*J,2*N-1:2*N) = Kg(2*J-1:2*J,2*N-1:2*N)+ KK(5:6,7:8,e);    
    Kg(2*J-1:2*J,2*K-1:2*K) = Kg(2*J-1:2*J,2*K-1:2*K) + KK(5:6,9:10,e);        
    Kg(2*J-1:2*J,2*L-1:2*L) =  Kg(2*J-1:2*J,2*L-1:2*L)+ KK(5:6,11:12,e);   
    
    Kg(2*N-1:2*N,2*I-1:2*I) = Kg(2*N-1:2*N,2*I-1:2*I) + KK(7:8,1:2,e);
    Kg(2*N-1:2*N,2*M-1:2*M) = Kg(2*N-1:2*N,2*M-1:2*M) + KK(7:8,3:4,e);
    Kg(2*N-1:2*N,2*J-1:2*J) = Kg(2*N-1:2*N,2*J-1:2*J)+ KK(7:8,5:6,e);
    Kg(2*N-1:2*N,2*N-1:2*N) = Kg(2*N-1:2*N,2*N-1:2*N)+ KK(7:8,7:8,e);    
    Kg(2*N-1:2*N,2*K-1:2*K) = Kg(2*N-1:2*N,2*K-1:2*K) + KK(7:8,9:10,e);        
    Kg(2*N-1:2*N,2*L-1:2*L) =  Kg(2*N-1:2*N,2*L-1:2*L)+ KK(7:8,11:12,e);   
 
    Kg(2*K-1:2*K,2*I-1:2*I) = Kg(2*K-1:2*K,2*I-1:2*I) + KK(9:10,1:2,e);
    Kg(2*K-1:2*K,2*M-1:2*M) = Kg(2*K-1:2*K,2*M-1:2*M) + KK(9:10,3:4,e);
    Kg(2*K-1:2*K,2*J-1:2*J) = Kg(2*K-1:2*K,2*J-1:2*J)+ KK(9:10,5:6,e);
    Kg(2*K-1:2*K,2*N-1:2*N) = Kg(2*K-1:2*K,2*N-1:2*N)+ KK(9:10,7:8,e);    
    Kg(2*K-1:2*K,2*K-1:2*K) = Kg(2*K-1:2*K,2*K-1:2*K) + KK(9:10,9:10,e);        
    Kg(2*K-1:2*K,2*L-1:2*L) =  Kg(2*K-1:2*K,2*L-1:2*L)+ KK(9:10,11:12,e);   
     
    Kg(2*L-1:2*L,2*I-1:2*I) = Kg(2*L-1:2*L,2*I-1:2*I) + KK(11:12,1:2,e);
    Kg(2*L-1:2*L,2*M-1:2*M) = Kg(2*L-1:2*L,2*M-1:2*M) + KK(11:12,3:4,e);
    Kg(2*L-1:2*L,2*J-1:2*J) = Kg(2*L-1:2*L,2*J-1:2*J)+ KK(11:12,5:6,e);
    Kg(2*L-1:2*L,2*N-1:2*N) = Kg(2*L-1:2*L,2*N-1:2*N)+ KK(11:12,7:8,e);    
    Kg(2*L-1:2*L,2*K-1:2*K) = Kg(2*L-1:2*L,2*K-1:2*K) + KK(11:12,9:10,e);        
    Kg(2*L-1:2*L,2*L-1:2*L) =  Kg(2*L-1:2*L,2*L-1:2*L)+ KK(11:12,11:12,e);      
end

%% forze di superficie
% carico distribuito uniforme
s=H:-H/NY:0;
qtr=45*s/H; % carico triangolare

for e=1:nel
   f1_load(:,:,e)=zeros(nne*ndof,1);
end   

% Nload(X,Y)=[NI(X,Y) 0 0 0 0 0 0 0 NK(X,Y) 0 NL(X,Y) 0;...
%                  0 NI(X,Y) 0 0 0 0 0 0 0 NK(X,Y) 0 NL(X,Y)];            
for e=1:2  % poichè il carico agisce solo sugli elementi da 1 a 2
    
    I=conn(e,1); 
    M=conn(e,2);
    J=conn(e,3);
    N=conn(e,4); 
    K=conn(e,5);
    L=conn(e,6);
     
    L(e)=y(K)-y(I);
% Rqtr(e)=((qtr(e)-qtr(e+1))*L(e))/2; % risultante carico triangolare 
   
   q(:,:,e)=[5+qtr(e+1); 0];  

% N_load(:,:,e) = int( int(Nload(X,Y),Y,[y(J) (y(J)+(X-x(J))*(y(K)-y(J))/(x(K)-x(J)))]),X,[x(I) x(J)]); 
% N_load(X,Y)=subs(Nload(X,Y),[ xI yI xJ yJ xK yK],[ x(I) y(I) x(J) y(J) x(K) y(K)]); 
%  f_sup(:,:,e)=int(N_load(0,Y)'*q(:,:,e),Y,[0 dy])*t
    f1_load(:,:,e)=[q(1,e)/6; 0; 0; 0; 0; 0; 0; 0; q(1,e)/6; 0; q(1,e)*2/3; 0]*L(e);
end

F_load=zeros(ndof*nnd,1);

for e=1:nel
    
    I=conn(e,1); 
    M=conn(e,2);
    J=conn(e,3);
    N=conn(e,4); 
    K=conn(e,5);
    L=conn(e,6);
     
    % transfer of the submatrices to their position in the stiffness matrix global
    F_load(2*I-1:2*I,1) = F_load(2*I-1:2*I,1) + f1_load(1:2,1,e);

    F_load(2*M-1:2*M,1) = F_load(2*M-1:2*M,1) + f1_load(3:4,1,e);

    F_load(2*J-1:2*J,1) = F_load(2*J-1:2*J,1) + f1_load(5:6,1,e);
    
    F_load(2*N-1:2*N,1) = F_load(2*N-1:2*N,1) + f1_load(7:8,1,e);
      
    F_load(2*K-1:2*K,1) = F_load(2*K-1:2*K,1) + f1_load(9:10,1,e);
    
    F_load(2*L-1:2*L,1) = F_load(2*L-1:2*L,1) + f1_load(11:12,1,e);
end   

%% forze di volume
g=9.80665; %accelerazione di gravità [m/s^2]

% Nglobal(X,Y)=[NI(X,Y) 0 NM(X,Y) 0 NJ(X,Y) 0 NN(X,Y) 0 NK(X,Y) 0 NL(X,Y) 0;...
%                    0 NI(X,Y) 0 NM(X,Y) 0 NJ(X,Y) 0 NN(X,Y) 0 NK(X,Y) 0 NL(X,Y)];
            
for e=1:nel
    
    I=conn(e,1); 
    M=conn(e,2);
    J=conn(e,3);
    N=conn(e,4); 
    K=conn(e,5);
    L=conn(e,6);
    
  % Nnum(:,:,e)=subs(Nglobal(Xi,Etha),[xI yI xM yM xJ yJ xN yN xK yK xL yL],[x(I) y(I) x(M) y(M) x(J) y(J) x(N) y(N) x(K) y(K) x(L) y(L)]);
    Nnum(:,:,e)=Nglobal(Xi,Etha)

   f_volume(:,:,e)=double(int(int(Nnum(:,:,e)'*[0;-rho*g]*det(Jac(:,:,e)),Etha,[0 1-Xi]),Xi,[0 1]));
 end   

F_volume=zeros(ndof*nnd,1);

for e=1:nel;
    
    I=conn(e,1); 
    M=conn(e,2);
    J=conn(e,3);
    N=conn(e,4); 
    K=conn(e,5);
    L=conn(e,6);
     
% transfer of the submatrices to their position in the stiffness matrix global

F_volume(2*I-1:2*I,1) =F_volume(2*I-1:2*I,1) + f_volume(1:2,1,e);

F_volume(2*M-1:2*M,1) = F_volume(2*M-1:2*M,1) + f_volume(3:4,1,e);

 F_volume(2*J-1:2*J,1) = F_volume(2*J-1:2*J,1) + f_volume(5:6,1,e);
    
  F_volume(2*N-1:2*N,1) = F_volume(2*N-1:2*N,1) + f_volume(7:8,1,e);
      
 F_volume(2*K-1:2*K,1) = F_volume(2*K-1:2*K,1) + f_volume(9:10,1,e);
    
 F_volume(2*L-1:2*L,1) =F_volume(2*L-1:2*L,1) + f_volume(11:12,1,e);
end   

F=F_load+F_volume;
%fare una copia F  in Fstar e azzerare le righe in corrispondenza dei
%gradi di libertà vincolati.
%% Boundary conditions
Kstar = Kg;
Fstar=F;

for e=[1 6 11 16 20 23 25] % nodi vincolati
    Kstar(e*ndof-1,:)=0;
    Kstar(:,e*ndof-1)=0;
    
    Kstar(e*ndof,:)=0;
    Kstar(:,e*ndof)=0;

    Kstar(e*ndof-1,e*ndof-1)=1;
    Kstar(e*ndof,e*ndof)=1;

    Fstar(e*ndof-1)=0;
    Fstar(e*ndof)=0;
end
%% Solution
v= inv(Kstar)*Fstar;   % nodal diplacements 
f= Kg*v;                   % reaction forces in the nodes
r =(f-F);                    % reaction force in external constraints
%% diplacements 
i=1:2:ndof*nnd;
vx=v(i);  % nodal diplacements in X directions
fx=f(i);  % reaction forces in the nodes in X directions 
rx=r(i);  %reazioni vincolari in X directions 
j=2:2:ndof*nnd;
vy=v(j);  % nodal diplacements in Y directions 
fy=f(j);  % reaction forces in the nodes in Y directions
ry=r(j);  %reazioni vincolari in in Y directions

%verifica 1 % massa deve essere uguale a Ry
massa=rho*g*((B+b)*H)/2;
Ry=sum(ry);
% verifica 2 % la risultante del carico trapezoidale deve essere uguale
% alla reazione orizzontale
R_trapezio=(45+b)*H/2;
Rx=sum(rx);

for e=1:nel
    
I=conn(e,1); 
M=conn(e,2);
J=conn(e,3);
N=conn(e,4); 
K=conn(e,5);
L=conn(e,6);

        u(:,:,e)=[vx(I);...
                    vy(I);...
                    vx(M);...
                    vy(M);...
                    vx(J);...
                    vy(J);...
                    vx(N);...
                    vy(N);...
                    vx(K);...
                    vy(K);...
                    vx(L);...
                    vy(L)];
        
    epslon(:,:,e)=BB(:,:,e)*u(:,:,e) ;    
    
    sigma(:,:,e)=E*epslon(:,:,e);
    
sigmaZ(:,:,e)=vu*(sigma(1,:,e)+sigma(2,:,e));
    
    sigma_VM(:,:,e)=sqrt((sigma(1,:,e)^2+sigma(2,:,e)^2+sigmaZ(:,:,e)^2)-(sigma(1,:,e)*sigma(2,:,e)+sigma(2,:,e)*sigmaZ(:,:,e)+sigmaZ(:,:,e)*sigma(1,:,e))+3*sigma(3,:,e)^2);
end

%% Plot 
%deformata
for i=1:nel
I=conn(i,1); 
M=conn(i,2);
J=conn(i,3);
N=conn(i,4); 
K=conn(i,5);
L=conn(i,6);

c=50000 %costante di amplificazione
xx2=[x(L)+c*vx(L);x(I)+c*vx(I)];
yy2=[y(L)+c*vy(L);y(I)+c*vy(I)];
S2(i)=plot(xx2,yy2,'g-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','g');
hold on
xx3=[x(I)+c*vx(I);x(M)+c*vx(M);x(J)+c*vx(J);x(N)+c*vx(N);x(K)+c*vx(K);x(L)+c*vx(L)];
yy3=[y(I)+c*vy(I);y(M)+c*vy(M);y(J)+c*vy(J);y(N)+c*vy(N);y(K)+c*vy(K);y(L)+c*vy(L)];
S3(i)=plot(xx3,yy3,'g-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','g');
title(['\fontname{Courier}\fontsize{17}Plot Structure & Deformed '],'color','K');
axis([-0.1 1 -0.06  5.8]);
axis equal
axis off
hold on
end
%Tensioni
for e=1:nel
I=conn(e,1); 
M=conn(e,2);
J=conn(e,3);
N=conn(e,4); 
K=conn(e,5);
L=conn(e,6);

XX(:,:,e) = linspace(x(I),x(J));
YY(:,:,e) = y(J)+(XX(:,:,e)-x(J))*(y(K)-y(J))/(x(K)-x(J));
end
