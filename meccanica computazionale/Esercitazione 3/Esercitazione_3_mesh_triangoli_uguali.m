%%  Computational Mechanics - TERZA ESERCITAZIONE
% stress triangle (CST) element
%% Clear workspace and close any open windows
clear all
close all
%% Set Command Window output display format
format bank;
%% ALTER THE NEXT LINES TO CHOOSE AN OUTPUT FILE FOR THE RESULTS
% Open file for output of results
% fid = fopen('Ese2_results.txt','w'); 
%% ALTER THE NEXT LINE TO CHOOSE AN INPUT FILE
% disp('Results printed in file : truss_1_results.txt');
%% data of the problem                      
Matricola=506682;             % my university number
H=(Matricola /100000)        % [m] length
b=1;                                % [m] length
B=3*b;                             % [m] length
vu=0.2;                            % Poisson's ratio
rho= 2400;                        % density [kg/m3]
E1=3E7;                            % Young's modulus concrete [KN/m^2]
t=1                                  % [m] Thickness
%% data structure
nnd = 25;                 % Number of nodes in system
nel = 32;     % Number of elements
nne = 3;                   % Number of nodes for element
ndof =2;                   % Number of degrees of freedom for node
sdof = nnd*ndof;       % Number of degrees of freedom for element
%% Nodes cFoordinates X and      Y
NX=6
NY=4

i=0:H/NY:H; % discretizzazione lungo Y
j=0:B/NX:B; % discretizzazione lungo X

dx=B/NX; % discretizzazione lungo X
dy=H/NY; % discretizzazione lungo Y

 e=1:5;
 y(e)=i;    x(e)=j(1);
 
 e=6:10;   
 y(e)=i;    x(e)=j(2);
 
 e=11:15   
 y(e)=i;    x(e)=j(3);
 
 e=16:19;
 i=0:H/4:(H-(H/4));
 y(e)=i;    x(e)=j(4);
 
 e=20:22;
 i=0:H/4:(H-(2*H/4));
 y(e)=i;   x(e)=j(5);
 
 e=23:24;
 i=0:H/4:(H-(3*H/4));
 y(e)=i;   x(e)=j(6);
 
 e=25;
 i=0:H/4:(H-(4*H/4));
 y(e)=i;   x(e)=j(7);
 
%% Element connectivity
e=1:4; % elementi 
I=1:4;
J=2:5;
K=6:9;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K; 

e=5:8; % elementi
I=6:9;
J=2:5;
K=7:10;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;

e=9:12; % elementi
I=6:9;
J=7:10;
K=11:14;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;

e=13:16; % elementi 
I=11:14;
J=7:10;
K=12:15;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;

e=17:20; % elementi 
I=11:14;
J=12:15;
K=16:19;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;

e=21:23; % elementi 
I=16:18;
J=12:14;
K=17:19;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;

e=24:26; % elementi 
I=16:18;
J=17:19;
K=20:22;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;

e=27:28; % elementi 
I=20:21;
J=17:18;
K=21:22;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;

e=29:30; % elementi 
I=20:21;
J=21:22;
K=23:24;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;

e=31; % elementi 
I=23;
J=21;
K=24;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;

e=32; % elementi 
I=23;
J=24;
K=25;
conn(e,1)=I;  conn(e,2)=J;  conn(e,3)=K;
%% Elasticity matrix
E=E1/((1+vu)*(1-2*vu))*[1-vu -vu 0.  ;...  %Matrix elastic
                                     -vu 1-vu 0.  ;...
                                     0  0 (1-2*vu)/2];
                                 
%% Load vector in [KN]
%% forze di superficie

s=H:-H/4:0;
qtr=45*s/H; % carico triangolare

for e=1:nel
   f_load(:,:,e)=zeros(6,1);
   
for e=1:4  % poichè il carico agisce solo sugli elementi da 1 a 4
    
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
    
    L(e)=y(J)-y(I);
    
    Rqtr(e)=((qtr(e)-qtr(e+1))*L(e))/2; % risultante carico triangolare  
     
    f_load(:,:,e)=[5/2+((qtr(e+1)*L(e)/2)+Rqtr(e)/3); 0;5/2+((qtr(e+1)*L(e)/2)+Rqtr(e)/6);0;0;0]*t*L(e);  
end
end

F_load=zeros(ndof*nnd,1);

for e=1:nel
    
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
     
    % transfer of the submatrices to their position in the stiffness matrix global
    F_load(2*I-1:2*I,1) = F_load(2*I-1:2*I,1) + f_load(1:2,1,e);

    F_load(2*J-1:2*J,1) = F_load(2*J-1:2*J,1) + f_load(3:4,1,e);

    F_load(2*K-1:2*K,1) = F_load(2*K-1:2*K,1) + f_load(5:6,1,e);
end   

%% forze di volume
g=9.80665 %accelerazione di gravità [m/s^2]

for e=1:nel
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
    
    % Area elements
    A=1/2*abs(det([1 x(I) y(I); 1  x(J)  y(J); 1  x(K)  y(K)]));
    
    % Forze di volume
    f_volume(:,:,e)=-(t/3)*[0; rho*g*A;0;rho*g*A;0;rho*g*A];
end

F_volume=zeros(ndof*nnd,1);

for e=1:nel;
    
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
     
% transfer of the submatrices to their position in the stiffness matrix global
F_volume(2*I-1:2*I,1) = F_volume(2*I-1:2*I,1) + f_volume(1:2,1,e);
   
F_volume(2*J-1:2*J,1) = F_volume(2*J-1:2*J,1) + f_volume(3:4,1,e);
  
F_volume(2*K-1:2*K,1) = F_volume(2*K-1:2*K,1) + f_volume(5:6,1,e);
end   

F=F_load+F_volume;
%fare una copia F  in Fstar e azzerare le righe in corrispondenza dei
%gradi di libertà vincolati.
%% matrix stiffnes composition              
for e=1:nel
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    
    % Area elements
    A=1/2*abs(det([1 x(I) y(I); 1  x(J)  y(J); 1  x(K)  y(K)]));
    
    % costant in Shape functions
    aI=(x(J)*y(K)-x(K)*y(J))/(2*A);
    bI=(y(J)-y(K))/(2*A);
    cI=(x(K)-x(J))/(2*A);
    
    aJ=(x(K)*y(I)-x(I)*y(K))/(2*A);
    bJ=(y(K)-y(I))/(2*A);
    cJ=(x(I)-x(K))/(2*A);
    
    aK=(x(I)*y(J)-x(J)*y(I))/(2*A);
    bK=(y(I)-y(J))/(2*A);
    cK=(x(J)-x(I))/(2*A);
    
    % Strain-displacement matrices
    BI = [ bI 0; 0 cI; cI bI ];
    BJ = [ bJ 0; 0 cJ; cJ bJ ];
    BK = [ bK 0; 0 cK; cK bK ];
    BB = [ BI BJ BK ];
    
    % Stiffness matrix
    Ke(:,:,e) = BB'*E*BB*A*t;
end

%% Global stiffness matrix
Kg = zeros(2*nnd,2*nnd); % dimensions of the global stiffness matrix

for e = 1:nel;
    I = conn(e,1); % node I
    J = conn(e,2); % node J
    K= conn(e,3); % node K
    
 % transfer of the submatrices to their position in the stiffness matrix global
    Kg(2*I-1:2*I,2*I-1:2*I) = Kg(2*I-1:2*I,2*I-1:2*I) + Ke(1:2,1:2,e);
    Kg(2*I-1:2*I,2*J-1:2*J) = Kg(2*I-1:2*I,2*J-1:2*J) + Ke(1:2,3:4,e);
    Kg(2*I-1:2*I,2*K-1:2*K) = Kg(2*I-1:2*I,2*K-1:2*K)+ Ke(1:2,5:6,e);
    
    Kg(2*J-1:2*J,2*I-1:2*I) = Kg(2*J-1:2*J,2*I-1:2*I) + Ke(3:4,1:2,e);
    Kg(2*J-1:2*J,2*J-1:2*J) = Kg(2*J-1:2*J,2*J-1:2*J) + Ke(3:4,3:4,e);
    Kg(2*J-1:2*J,2*K-1:2*K) = Kg(2*J-1:2*J,2*K-1:2*K)+ Ke(3:4,5:6,e);
    
    Kg(2*K-1:2*K,2*I-1:2*I) = Kg(2*K-1:2*K,2*I-1:2*I) + Ke(5:6,1:2,e);
    Kg(2*K-1:2*K,2*J-1:2*J) = Kg(2*K-1:2*K,2*J-1:2*J) + Ke(5:6,3:4,e);
    Kg(2*K-1:2*K,2*K-1:2*K) = Kg(2*K-1:2*K,2*K-1:2*K) + Ke(5:6,5:6,e);
end

%% Boundary conditions
Kstar = Kg;
Fstar=F

for e=[1 6 11 16 20 23 25] % nodi vincolati
    Kstar(e*ndof-1,:)=0;
    Kstar(:,e*ndof)=0;

    Kstar(e*ndof-1,e*ndof-1)=1;
    Kstar(e*ndof,e*ndof)=1;

    Fstar(e*ndof-1)=0
    Fstar(e*ndof)=0
end

%verifica
massa=rho*g*((B+b)*H)/2
%% Solution
v= inv(Kstar)*Fstar;   % nodal diplacements 
f= Kg*v;               % reaction forces in the nodes
r =(f-Fstar);              % reaction force in external constraints
%% diplacements 

i=1:2:ndof*nnd;
vx=v(i);  % nodal diplacements in X directions in global 
fx=f(i);  % reaction forces in the nodes in X directions in global
j=2:2:ndof*nnd;
vy=v(j);  % nodal diplacements in Y directions in global
fy=f(j);  % reaction forces in the nodes in Y directions in global

 i=1:nnd
rr=sum(fy(i))

for e=1:nel
       I = conn(e,1); % node I
       J = conn(e,2); % node J
       K= conn(e,3); % node K
    
        u(:,:,e)=[vx(I);...
                    vy(I);...
                    vx(J);...
                    vy(J);...
                    vx(K);...
                    vy(K)]
                
    epslon(:,:,e)=BB*u(:,:,e)     
    
    sigma(:,:,e)=E*epslon(:,:,e)
end

%% Plot the stresses in the x_direction
figure(2)

x_stress=sigma(1,:);
cmin=min(x_stress);
cmax=max(x_stress);
caxis([ cmin cmax])
geom=[x' y']
patch('Faces', conn, 'Vertices',geom, 'FaceVertexCData', x_stress',...
    'Facecolor','flat','Marker','o')
axis equal
colorbar
%


%% TRIANGULATION
figure(1)
W=[x',y']
T=[conn]
TR = triangulation(T,W)
TR.Points
TR.ConnectivityList

triplot(TR,'LineWidth',0.1,'Color',[0.1 0.2 0])
axis equal

IC = incenter(TR) % centro del triangolo
hold on

numtri = size(TR,1);
trilabels = arrayfun(@(A) {sprintf('T%d', A)}, (1:numtri)');
Htl = text(IC(:,1),IC(:,2),trilabels,'FontWeight','bold', ...
'HorizontalAlignment','center','Color',[0 0 0],'FontSize',3.5);
hold on

%% Print node
for i=1:nnd
text(x(i)+0.1,y(i)+0.06,{(i)},...
    'HorizontalAlignment','right',...
    'FontSize',3.5)
axis off
hold on
end
%% Print force
a=0.050
X=[-0.3 -0.3 -0.3-a*5 -0.3-a*50]
Y=[0 H H 0]
hpatch=patch(X,Y,'w','Facealpha',0.2,'EdgeColor', 'r');

for s=0:H/12:H
j=a*((45/H)*(H-s)+5)
plot(-0.33,s,'>r','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
p0=[-0.3-j s];
p1=[-0.3 s];
vectarrow(p0,p1)
axis off
hold on
end

figure(3)

%%print deformata
WW=[x'+vx,y'+0.00001*vy]
TT=[conn]
TR=triangulation(TT,WW)
TR.Points
TR.ConnectivityList

triplot(TR,'LineWidth',1.5)
axis equal
hold on






