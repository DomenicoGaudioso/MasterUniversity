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

I=[ 1 3 13 15 11 13 22 20]';
J=[ 11  13 3 5 20 22 13 25]';
K=[ 3 5 11 13 13 15 20 22]';
L=[ 6 8 8 10 16 18 18 23]';
M=[ 7 9 7 9 17 19 17 24]';
N=[ 2 4 12 14 12 14 21 21]';

for e=1:nel
conn(e,1)=I(e);  conn(e,2)=L(e);  conn(e,3)=J(e); conn(e,4)=M(e);  conn(e,5)=K(e);  conn(e,6)=N(e); 
end


for i=1:nel
 figure(1)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
I=conn(i,1); 
L=conn(i,2);
J=conn(i,3);
M=conn(i,4); 
K=conn(i,5);
N=conn(i,6);
xx1=[x(N);x(I)];
yy1=[y(N);y(I)];
S1(i)=plot(xx1,yy1,'K-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','k');
hold on
xx=[x(I);x(L);x(J);x(M);x(K);x(N)];
yy=[y(I);y(L);y(J);y(M);y(K);y(N)];
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
 syms   X Y  xI yI xJ yJ xK yK xL yL xM yM xN yN  real
%% Shape functions
for e=1:nel
    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);
    
P = [ 1 X Y X.^2 X.*Y Y.^2 0 0 0 0 0 0;...
      0 0 0 0 0 0 1 X Y X.^2 X.*Y Y.^2 ];    

C= [ 1 x(I) y(I) x(I).^2 x(I).*y(I) y(I)^2 0 0 0 0 0 0;              % Shape function calulation
                  0 0 0 0 0 0 1 x(I) y(I) x(I).^2 x(I).*y(I) y(I).^2;
                  1 x(J) y(J) x(J).^2 x(J).*y(J) y(J).^2 0 0 0 0 0 0;
                  0 0 0 0 0 0 1 x(J) y(J) x(J).^2 x(J).*y(J) y(J).^2;
                  1 x(K) y(K) x(K).^2 x(K).*y(K) y(K).^2 0 0 0 0 0 0;
                  0 0 0 0 0 0 1 x(K) y(K) x(K).^2 x(K).*y(K) y(K).^2;
                  1 x(L) y(L) x(L).^2 x(L).*y(L) y(L).^2 0 0 0 0 0 0;
                  0 0 0 0 0 0 1 x(L) y(L) x(L).^2 x(L).*y(L) y(L).^2;
                  1 x(M) y(M) x(M).^2 x(M).*y(M) y(M).^2 0 0 0 0 0 0;
                  0 0 0 0 0 0 1 x(M) y(M) x(M).^2 x(M).*y(M) y(M).^2;
                  1 x(N) y(N) x(N).^2 x(N).*y(N) y(N).^2 0 0 0 0 0 0;
                  0 0 0 0 0 0 1 x(N) y(N) x(N).^2 x(N).*y(N) y(N).^2 ];
              
  N=P*inv(C);
 
  NI(:,:,e)=N(1,1); NJ(:,:,e)=N(1,3); NK(:,:,e)=N(1,5); NL(:,:,e)=N(1,7); NM(:,:,e)=N(1,9); NN(:,:,e)=N(1,11);
  
  Nglobal(:,:,e)=[NI(:,:,e) 0 NJ(:,:,e) 0 NK(:,:,e) 0 NL(:,:,e) 0 NM(:,:,e) 0 NN(:,:,e) 0;...
                       0 NI(:,:,e) 0 NJ(:,:,e) 0 NK(:,:,e) 0 NL(:,:,e) 0 NM(:,:,e) 0 NN(:,:,e)];
  
  BI(:,:,e)=[diff(NI(:,:,e),X) 0; 0 diff(NI(:,:,e),Y); diff(NI(:,:,e),Y) diff(NI(:,:,e),X)];
  BJ(:,:,e)=[diff(NJ(:,:,e),X) 0; 0 diff(NJ(:,:,e),Y); diff(NJ(:,:,e),Y) diff(NJ(:,:,e),X)];
  BK(:,:,e)=[diff(NK(:,:,e),X) 0; 0 diff(NK(:,:,e),Y); diff(NK(:,:,e),Y) diff(NK(:,:,e),X)];
  BL(:,:,e)=[diff(NL(:,:,e),X) 0; 0 diff(NL(:,:,e),Y); diff(NL(:,:,e),Y) diff(NL(:,:,e),X)];
  BM(:,:,e)=[diff(NM(:,:,e),X) 0; 0 diff(NM(:,:,e),Y); diff(NM(:,:,e),Y) diff(NM(:,:,e),X)];
  BN(:,:,e)=[diff(NN(:,:,e),X) 0; 0 diff(NN(:,:,e),Y); diff(NN(:,:,e),Y) diff(NN(:,:,e),X)];
  
  BB(:,:,e)=[BI(:,:,e) BJ(:,:,e) BK(:,:,e) BL(:,:,e) BM(:,:,e) BN(:,:,e)];
  
  K0(:,:,e)=BB(:,:,e)'*E*BB(:,:,e);
  
  K1(:,:,e)=int(K0(:,:,e),Y,[y(J) y(J)+(X-x(J))*(y(K)-y(J))/(x(K)-x(J))]);
  Ke(:,:,e)=int(K1(:,:,e),X,[x(I) x(J)]);
  KK(:,:,e)=double(Ke(:,:,e));
end

%% Global stiffness matrix
Kg=zeros(2*nnd, 2*nnd);
for e = 1:nel
    
    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);
%     
%  sctr = [ 2*I-1 2*I 2*J-1 2*J 2*K-1 2*K 2*L-1 2*L 2*M-1 2*M 2*N-1 2*N ];
%  K(sctr,sctr) = K(sctr,sctr) + KK(:,:,e);
    
 % transfer of the submatrices to their position in the stiffness matrix global
    Kg(2*I-1:2*I,2*I-1:2*I) = Kg(2*I-1:2*I,2*I-1:2*I) + KK(1:2,1:2,e);
    Kg(2*I-1:2*I,2*J-1:2*J) = Kg(2*I-1:2*I,2*J-1:2*J) + KK(1:2,3:4,e);
    Kg(2*I-1:2*I,2*K-1:2*K) = Kg(2*I-1:2*I,2*K-1:2*K)+ KK(1:2,5:6,e);
    Kg(2*I-1:2*I,2*L-1:2*L) = Kg(2*I-1:2*I,2*L-1:2*L)+ KK(1:2,7:8,e);    
    Kg(2*I-1:2*I,2*M-1:2*M) = Kg(2*I-1:2*I,2*M-1:2*M)+ KK(1:2,9:10,e);        
    Kg(2*I-1:2*I,2*N-1:2*N) =  Kg(2*I-1:2*I,2*N-1:2*N)+ KK(1:2,11:12,e); 
    
    Kg(2*J-1:2*J,2*I-1:2*I) = Kg(2*J-1:2*J,2*I-1:2*I) + KK(3:4,1:2,e);
    Kg(2*J-1:2*J,2*J-1:2*J) = Kg(2*J-1:2*J,2*J-1:2*J) + KK(3:4,3:4,e);
    Kg(2*J-1:2*J,2*K-1:2*K) = Kg(2*J-1:2*J,2*K-1:2*K)+ KK(3:4,5:6,e);
    Kg(2*J-1:2*J,2*L-1:2*L) = Kg(2*J-1:2*J,2*L-1:2*L)+ KK(3:4,7:8,e);    
    Kg(2*J-1:2*J,2*M-1:2*M) = Kg(2*J-1:2*J,2*M-1:2*M) + KK(3:4,9:10,e);        
    Kg(2*J-1:2*J,2*N-1:2*N) =  Kg(2*J-1:2*J,2*N-1:2*N) + KK(3:4,11:12,e);   
    
    Kg(2*K-1:2*K,2*I-1:2*I) = Kg(2*K-1:2*K,2*I-1:2*I)+ KK(5:6,1:2,e);
    Kg(2*K-1:2*K,2*J-1:2*J) = Kg(2*K-1:2*K,2*J-1:2*J) + KK(5:6,3:4,e);
    Kg(2*K-1:2*K,2*K-1:2*K) = Kg(2*K-1:2*K,2*K-1:2*K)+ KK(5:6,5:6,e);
    Kg(2*K-1:2*K,2*L-1:2*L) = Kg(2*K-1:2*K,2*L-1:2*L)+ KK(5:6,7:8,e);    
    Kg(2*K-1:2*K,2*M-1:2*M) = Kg(2*K-1:2*K,2*M-1:2*M) + KK(5:6,9:10,e);        
    Kg(2*K-1:2*K,2*N-1:2*N) =  Kg(2*K-1:2*K,2*N-1:2*N) + KK(5:6,11:12,e);   
    
    Kg(2*L-1:2*L,2*I-1:2*I) = Kg(2*L-1:2*L,2*I-1:2*I) + KK(7:8,1:2,e);
    Kg(2*L-1:2*L,2*J-1:2*J) = Kg(2*L-1:2*L,2*J-1:2*J) + KK(7:8,3:4,e);
    Kg(2*L-1:2*L,2*K-1:2*K) = Kg(2*L-1:2*L,2*K-1:2*K)+ KK(7:8,5:6,e);
    Kg(2*L-1:2*L,2*L-1:2*L) = Kg(2*L-1:2*L,2*L-1:2*L)+ KK(7:8,7:8,e);    
    Kg(2*L-1:2*L,2*M-1:2*M) = Kg(2*L-1:2*L,2*M-1:2*M) + KK(7:8,9:10,e);        
    Kg(2*L-1:2*L,2*N-1:2*N) =  Kg(2*L-1:2*L,2*N-1:2*N) + KK(7:8,11:12,e);   
 
    Kg(2*M-1:2*M,2*I-1:2*I) = Kg(2*M-1:2*M,2*I-1:2*I)  + KK(9:10,1:2,e);
    Kg(2*M-1:2*M,2*J-1:2*J) = Kg(2*M-1:2*M,2*J-1:2*J) + KK(9:10,3:4,e);
    Kg(2*M-1:2*M,2*K-1:2*K) = Kg(2*M-1:2*M,2*K-1:2*K)+ KK(9:10,5:6,e);
    Kg(2*M-1:2*M,2*L-1:2*L) = Kg(2*M-1:2*M,2*L-1:2*L)+ KK(9:10,7:8,e);    
    Kg(2*M-1:2*M,2*M-1:2*M) = Kg(2*M-1:2*M,2*M-1:2*M)+ KK(9:10,9:10,e);        
    Kg(2*M-1:2*M,2*N-1:2*N) =  Kg(2*M-1:2*M,2*N-1:2*N)+ KK(9:10,11:12,e);   
     
    Kg(2*N-1:2*N,2*I-1:2*I) = Kg(2*N-1:2*N,2*I-1:2*I)  + KK(11:12,1:2,e);
    Kg(2*N-1:2*N,2*J-1:2*J) = Kg(2*N-1:2*N,2*J-1:2*J) + KK(11:12,3:4,e);
    Kg(2*N-1:2*N,2*K-1:2*K) = Kg(2*N-1:2*N,2*K-1:2*K)+ KK(11:12,5:6,e);
    Kg(2*N-1:2*N,2*L-1:2*L) = Kg(2*N-1:2*N,2*L-1:2*L)+ KK(11:12,7:8,e);    
    Kg(2*N-1:2*N,2*M-1:2*M) = Kg(2*N-1:2*N,2*M-1:2*M) + KK(11:12,9:10,e);        
    Kg(2*N-1:2*N,2*N-1:2*N) =  Kg(2*N-1:2*N,2*N-1:2*N)+ KK(11:12,11:12,e);      
end
%% Forze di superficie
s=H:-H/NY:0; %ascissa per il carico triangolare agente sulla parete
qtr=45*s/H; % carico triangolare

for e=1:nel
   f_load(:,:,e)=zeros(12,1);
end

for i=1:NY
for e=1:NY  % poichè il carico agisce solo sugli elementi di questa serie

    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);

    L=(y(K)-y(I));
    
    ql(i)=(qtr(i)-qtr(i+1)) ;  % max carico triangolare a quella quota x i
    
      Nq(:,:,e)=[NI(:,:,e) 0 0 0 NK(:,:,e) 0 0 0 0 0 NN(:,:,e) 0;...
                       0 NI(:,:,e) 0 0 0 NK(:,:,e) 0 0 0 0 0 NN(:,:,e)];

    Nq1(:,:,e)=subs(Nq(:,:,e),X,x(I)); %poichè tramite gauss trasformo l'integrale d'area in uno perimetrale  
       
       qtriangle=ql(1)*(L-Y)/L;
    %forza distribuita +  carico triangolare che varia con Y   
    q1(:,:,e)=[5+ql(i)+qtriangle; 0];
    f_load0(:,:,e)=Nq1(:,:,e)'*q1(:,:,e);
    f_load(:,:,e)=int(f_load0(:,:,e),Y,[y(I) y(K)])
 end
end

F_load=zeros(ndof*nnd,1);

for e=1:nel
   
    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);
     
    % transfer of the submatrices to their position in the stiffness matrix global
    F_load(2*I-1:2*I,1) = F_load(2*I-1:2*I,1) + f_load(1:2,1,e);

    F_load(2*J-1:2*J,1) = F_load(2*J-1:2*J,1) + f_load(3:4,1,e);

    F_load(2*K-1:2*K,1) = F_load(2*K-1:2*K,1) + f_load(5:6,1,e);
    
    F_load(2*L-1:2*L,1) = F_load(2*L-1:2*L,1) + f_load(7:8,1,e);
      
    F_load(2*M-1:2*M,1) = F_load(2*M-1:2*M,1) + f_load(9:10,1,e);
    
    F_load(2*N-1:2*N,1) = F_load(2*N-1:2*N,1) + f_load(11:12,1,e);
end  

%% forze di volume
g=9.80665; %accelerazione di gravità [m/s^2]
  
F_volume=zeros(ndof*nnd,1);
for e=1:nel
    
    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);
    
f_volume0(:,:,e)=Nglobal(:,:,e)'*[0;-rho*g];
f_volume1(:,:,e)=int(f_volume0(:,:,e),Y,[y(J) (y(J)+(X-x(J))*(y(K)-y(J))/(x(K)-x(J)))]);
f_volume(:,:,e)=int(f_volume1(:,:,e),X,[x(I) x(J)]);

% transfer of the submatrices to their position in the stiffness matrix global
F_volume(2*I-1:2*I,1) =F_volume(2*I-1:2*I,1) + f_volume(1:2,1,e);

F_volume(2*J-1:2*J,1) = F_volume(2*J-1:2*J,1) + f_volume(3:4,1,e);

 F_volume(2*K-1:2*K,1) = F_volume(2*K-1:2*K,1) + f_volume(5:6,1,e);
    
  F_volume(2*L-1:2*L,1) = F_volume(2*L-1:2*L,1) + f_volume(7:8,1,e);
      
 F_volume(2*M-1:2*M,1) = F_volume(2*M-1:2*M,1) + f_volume(9:10,1,e);
    
 F_volume(2*N-1:2*N,1) =F_volume(2*N-1:2*N,1) + f_volume(11:12,1,e);
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
R_trap=(50+b)*H/2;
Rx=sum(rx);

for e=1:nel
    
    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);

        u(:,:,e)=[vx(I);...
                    vy(I);...
                    vx(J);...
                    vy(J);...
                    vx(K);...
                    vy(K);...
                    vx(L);...
                    vy(L);...
                    vx(M);...
                    vy(M);...
                    vx(N);...
                    vy(N)];
    digits(2)
    
    epslon(:,:,e)=BB(:,:,e)*u(:,:,e);    
    
    sigma(:,:,e)=vpa(E*epslon(:,:,e));
    
    sigmaZ(:,:,e)=vu*(sigma(1,:,e)+sigma(2,:,e));

    sigma_VM(:,:,e)=vpa(simplify(sqrt((sigma(1,:,e)^2+sigma(2,:,e)^2+sigmaZ(:,:,e)^2)-(sigma(1,:,e)*sigma(2,:,e)+sigma(2,:,e)*sigmaZ(:,:,e)+sigmaZ(:,:,e)*sigma(1,:,e))+3*sigma(3,:,e)^2)));
end

%% Plot 
%deformata
for i=1:nel
    I=conn(i,1); 
    L=conn(i,2);
    J=conn(i,3);
    M=conn(i,4); 
    K=conn(i,5);
    N=conn(i,6);

c=20000 %costante di amplificazione
xx2=[x(N)+c*vx(N);x(I)+c*vx(I)];
yy2=[y(N)+c*vy(N);y(I)+c*vy(I)];
S2(i)=plot(xx2,yy2,'g-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','g');
hold on
xx3=[x(I)+c*vx(I);x(L)+c*vx(L);x(J)+c*vx(J);x(M)+c*vx(M);x(K)+c*vx(K);x(N)+c*vx(N)];
yy3=[y(I)+c*vy(I);y(L)+c*vy(L);y(J)+c*vy(J);y(M)+c*vy(M);y(K)+c*vy(K);y(N)+c*vy(N)];
S3(i)=plot(xx3,yy3,'g-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','g');
title(['\fontname{Courier}\fontsize{17}Plot Structure & Deformed '],'color','K');
axis([-0.1 1 -0.06  5.8]);
axis equal
axis off
hold on
end

%% Plot the stresses
cc=0.02; % parametro che scala il meshgrid
cc1=4
% sigma XX

for e=1:nel    
figure(2)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';

    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);
    
 if y(K)>y(I)
[X1,Y1] = meshgrid(x(I):cc:x(J),y(I):cc:y(K));
 else 
[X1,Y1] = meshgrid(x(I):-cc:x(J),y(I):-cc:y(K));
 end
sigma_X0(X,Y)=sigma(1,:,e);
sigma_X=sigma_X0(X1,Y1);
sigma_X((X1-x(J))/(x(K)-x(J))<(Y1-y(J))/(y(K)-y(J))) = nan;
[C,h] = contourf(X1,Y1,sigma_X,abs(x(J)-x(I))*cc1,'EdgeColor','k','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','K')
title(['\fontname{Courier}\fontsize{15}Stresses in the x direction (\sigma_x_x) '],'color','K');
zlabel('\sigma_x_x(x,y)','LineWidth',2)
xlabel('x','LineWidth',2)
ylabel('y','LineWidth',2)
axis([-0.1 1 -0.06  5.8]);
axis equal
c=colorbar
colormap(jet)
c.Label.String = '\sigma_x_x [KN/m^2]'
axis off
hold on
end

% sigma YY
for e=1:nel
    figure(3)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);
    
 if y(K)>y(I)
[X1,Y1] = meshgrid(x(I):cc:x(J),y(I):cc:y(K));
 else 
[X1,Y1] = meshgrid(x(I):-cc:x(J),y(I):-cc:y(K));
 end
sigma_Y0(X,Y)=sigma(2,:,e);
sigma_Y=sigma_Y0(X1,Y1);
sigma_Y((X1-x(J))/(x(K)-x(J))<(Y1-y(J))/(y(K)-y(J))) = nan;
[C,h] = contourf(X1,Y1,sigma_Y,abs(x(J)-x(I))*cc1,'EdgeColor','k','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','K')
title(['\fontname{Courier}\fontsize{15}Stresses in the y direction (\sigma_y_y) '],'color','K');
zlabel('\sigma_y_y(x,y)','LineWidth',2)
xlabel('x','LineWidth',2)
ylabel('y','LineWidth',2)
axis([-0.1 1 -0.06  5.8]);
axis equal
c=colorbar
colormap(jet)
c.Label.String = '\sigma_y_y [KN/m^2]'
axis off
hold on
end

% tau XY
for e=1:nel
    figure(4)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);
    
 if y(K)>y(I)
[X1,Y1] = meshgrid(x(I):cc:x(J),y(I):cc:y(K));
 else 
[X1,Y1] = meshgrid(x(I):-cc:x(J),y(I):-cc:y(K));
 end
tau_XY0(X,Y)=sigma(3,:,e);
tau_XY=tau_XY0(X1,Y1);
tau_XY((X1-x(J))/(x(K)-x(J))<(Y1-y(J))/(y(K)-y(J))) = nan;
[C,h] = contourf(X1,Y1,tau_XY,abs(x(J)-x(I))*cc1,'EdgeColor','k','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','K')
title(['\fontname{Courier}\fontsize{15}Stresses in the xy direction (\tau_x_y) '],'color','K');
zlabel('\tau_x_y','LineWidth',2)
xlabel('x','LineWidth',2)
ylabel('y','LineWidth',2)
axis([-0.1 1 -0.06  5.8]);
axis equal
c=colorbar
colormap(jet)
c.Label.String = '\tau_x_y [KN/m^2]'
axis off
hold on
end

% sigma Z
for e=1:nel
    figure(5)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6);    
 if y(K)>y(I)
[X1,Y1] = meshgrid(x(I):cc:x(J),y(I):cc:y(K));
 else 
[X1,Y1] = meshgrid(x(I):-cc:x(J),y(I):-cc:y(K));
 end
sigma_Z0(X,Y)=sigmaZ(:,:,e);
sigma_Z=sigma_Z0(X1,Y1);
sigma_Z((X1-x(J))/(x(K)-x(J))<(Y1-y(J))/(y(K)-y(J))) = nan;
[C,h] = contourf(X1,Y1,sigma_Z,abs(x(J)-x(I))*cc1,'EdgeColor','k','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','K')
title(['\fontname{Courier}\fontsize{15}Stresses in the z direction (\sigma_z_z) '],'color','K');
zlabel('\sigma_z_z','LineWidth',2)
xlabel('x','LineWidth',2)
ylabel('y','LineWidth',2)
axis([-0.1 1 -0.06  5.8]);
axis equal
c=colorbar
colormap(jet)
c.Label.String = '\sigma_z_z [KN/m^2]'
axis off
hold on
end

% sigma Von Mises
for e=1:nel
    figure(6)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
    I=conn(e,1); 
    L=conn(e,2);
    J=conn(e,3);
    M=conn(e,4); 
    K=conn(e,5);
    N=conn(e,6); 
 if y(K)>y(I)
[X1,Y1] = meshgrid(x(I):cc:x(J),y(I):cc:y(K));
 else 
[X1,Y1] = meshgrid(x(I):-cc:x(J),y(I):-cc:y(K));
 end
sigma_VM0(X,Y)= sigma_VM(:,:,e);
sigma_VM1=sigma_VM0(X1,Y1);
sigma_VM1((X1-x(J))/(x(K)-x(J))<(Y1-y(J))/(y(K)-y(J))) = nan;
[C,h] = contourf(X1,Y1,sigma_VM1,abs(x(J)-x(I))*cc1,'EdgeColor','k','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','K')
title(['\fontname{Courier}\fontsize{15}Von Mises stresses (\sigma_V_M) '],'color','K');
zlabel('\sigma_V_M','LineWidth',2)
xlabel('x','LineWidth',2)
ylabel('y','LineWidth',2)
axis([-0.1 1 -0.06  5.8]);
axis equal
c=colorbar
colormap(jet)
c.Label.String = '\sigma_V_M [KN/m^2]'
axis off
hold on
end