%% Elementi mappati
% Computational Mechanics
clear all
close all
%% Set Command Window output display format
format bank;
%% data of the problem                      
Matricola=506682;             % my university number
H=(Matricola /100000);       % [m] length
b=1;                                % [m] length
B=3*b;                            % [m] length
vu=0.2;                            % Poisson's ratio
rho= 2400/1000;                % density [kg/m3]  ''diviso 1000  perchè poi voglio che le forze di volume moltiplicate per l'area siano in KN''
E1=3E7;                           % Young's modulus concrete [KN/m^2]
t=1;                                 % [m] Thickness
%% D
NX=3; % numero di discretizzazioni in x
NY=NX; % numero di discretizzazioni in y

dx=B/NX;  % dimensione lungo x dei triangolini della mesh
dy=H/NY;  % dimensione lungo y dei triangolini della mesh
%% Data structure
nnd = (NX+1)*(NY+1); % Number of nodes in system
nel = NX*NY;      % Number of elements
nne = 4;                    % Number of nodes for element
ndof =2;                    % Number of degrees of freedom for node
sdof = nnd*ndof;        % Number of degrees of freedom for element  
%% cordinate 
for s=0:NY
b_h(s+1)=(B-b)*(s*dy)/H+b; %lunghezza della base del trapezio alle varie quote
end
b_h=fliplr(b_h);

for i=1:NX+1
dxx(i)= b_h(i)/NY;  % larghezza triangolini alle varie quote
end

%coordinate dei punti nodali
  n=0;
 i=1:NX+1;
 y(i)=n*dx ; x(i)=(i-min(i))*dxx(n+1);
     
    while length(y)< (NX+1)*(NY+1) 
i=i+NX+1;
n=n+1;
    y(i)=n*dy ; x(i)=(i-min(i))*dxx(n+1);
    end 
 %% Element connectivity
for j=1:NX
  for  i=1:NY
    e=((j-1)*NY)+i;
    conn(e,1)=(j-1)*(NY+1)+i;
    conn(e,2)=(j-1)*(NY+1)+i+1;
    conn(e,3)=(j)*(NY+1)+i+1;
    conn(e,4)=(j)*(NY+1)+i;
  end  
end
%% Plot figure
figure(1)
 for e=1:nel
         I = conn(e,1);
         J = conn(e,2);
         K = conn(e,3);
         L=conn(e,4);
        
        xx1=[x(L);x(I)];
        yy1=[y(L);y(I)];
        S1=plot(xx1,yy1,'K-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','k');
        hold on
        xx=[x(I);x(J);x(K);x(L)];
        yy=[y(I);y(J);y(K);y(L)];
        S=plot(xx,yy,'K-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','k');
        axis([-0.1 1 -0.06  5.8]);
        axis equal
        axis off
        hold on        
 end 
%% Symbol definition
syms Xi Etha X1 Y1 xI yI xJ yJ xK yK xL yL  real
%% Shape functions
NI(Xi,Etha)=((1-Xi)*(1-Etha))/4;
NJ(Xi,Etha)=((1+Xi)*(1-Etha))/4;
NK(Xi,Etha)=((1+Xi)*(1+Etha))/4;
NL(Xi,Etha)=((1-Xi)*(1+Etha))/4;

X=NI(Xi,Etha)*xI+NJ(Xi,Etha)*xJ+NK(Xi,Etha)*xK+NL(Xi,Etha)*xL;
Y=NI(Xi,Etha)*yI+NJ(Xi,Etha)*yJ+NK(Xi,Etha)*yK+NL(Xi,Etha)*yL;

Jacob(Xi,Etha)=[diff(X,Xi) diff(Y,Xi);diff(X,Etha) diff(Y,Etha)]; %Jacobiano

N(Xi,Etha)=[NI(Xi,Etha) 0 NJ(Xi,Etha) 0 NK(Xi,Etha) 0 NL(Xi,Etha) 0;0 NI(Xi,Etha) 0 NJ(Xi,Etha) 0 NK(Xi,Etha) 0 NL(Xi,Etha)];
%% Strain-displacement matrices
Diff_I=[diff(NI(Xi,Etha),Xi) diff(NI(Xi,Etha),Etha)]';
Diff_J=[diff(NJ(Xi,Etha),Xi) diff(NJ(Xi,Etha),Etha)]';
Diff_K=[diff(NK(Xi,Etha),Xi) diff(NK(Xi,Etha),Etha)]';
Diff_L=[diff(NL(Xi,Etha),Xi) diff(NL(Xi,Etha),Etha)]';

% B1=[diff(NI(Xi,Etha),Xi) 0; 0 diff(NI(Xi,Etha),Etha); diff(NI(Xi,Etha),Etha) diff(NI(Xi,Etha),Xi)]
% B2=[diff(NJ(Xi,Etha),Xi) 0; 0 diff(NJ(Xi,Etha),Etha); diff(NJ(Xi,Etha),Etha) diff(NJ(Xi,Etha),Xi)]
% B3=[diff(NK(Xi,Etha),Xi) 0; 0 diff(NK(Xi,Etha),Etha); diff(NK(Xi,Etha),Etha) diff(NK(Xi,Etha),Xi)]
% B4=[diff(NL(Xi,Etha),Xi) 0; 0 diff(NL(Xi,Etha),Etha); diff(NL(Xi,Etha),Etha) diff(NL(Xi,Etha),Xi)]
% 
% BB2=[B1 B2 B3 B4]

diffN1=inv(Jacob(Xi,Etha))*Diff_I;
diffN2=inv(Jacob(Xi,Etha))*Diff_J;
diffN3=inv(Jacob(Xi,Etha))*Diff_K;
diffN4=inv(Jacob(Xi,Etha))*Diff_L;

Bsym=[diffN1(1) 0 diffN2(1) 0 diffN3(1) 0 diffN4(1) 0;...
     0 diffN1(2) 0 diffN2(2) 0 diffN3(2) 0 diffN4(2);...
     diffN1(2) diffN1(1)  diffN2(2)  diffN2(1)  diffN3(2)  diffN3(1)  diffN4(2) diffN4(1)];
 %% Elasticity matrix
E=E1/((1+vu)*(1-2*vu))*[1-vu  -vu   0.  ;...  %Matrix elastic %stato piano di deformazione
                                     -vu  1-vu   0.  ;...
                                      0     0 (1-2*vu)/2];   
 %% passaggio da sym a numerico
  for e=1:nel
         I = conn(e,1);
         J = conn(e,2);
         K = conn(e,3);
         L=conn(e,4);
         
 Jac(:,:,e)=subs(Jacob(Xi,Etha),[xI yI xJ yJ xK yK xL yL],[x(I) y(I) x(J) y(J) x(K) y(K) x(L) y(L)]);      
 BB(:,:,e)=subs(Bsym,[xI yI xJ yJ xK yK xL yL],[x(I) y(I) x(J) y(J) x(K) y(K) x(L) y(L)]);
  end
%   for e=1:nel
%   Ke(:,:,e)=int(int(BB(:,:,e)'*E*BB(:,:,e)*det(Jac(:,:,e)),Xi,[-1 1]),Etha,[-1 1]);
%   end
  %% Numerical integral (Gauss)    
   DeltaXi=[0.555 0.888 0.555];
   DeltaEtha=[0.555 0.888 0.555];
    
   i=[-0.775  0  +0.775];
   j=[-0.775  0  +0.775];
for e=1:nel                                        
Function(:,:,e)=(BB(:,:,e)'*E*BB(:,:,e)*det(Jac(:,:,e)));

K1(:,:,e)=subs(Function(:,:,e),[Xi Etha],[i(1) j(1)])*DeltaXi(1)*DeltaEtha(1);
K2(:,:,e)=subs(Function(:,:,e),[Xi Etha],[i(2) j(2)])*DeltaXi(2)*DeltaEtha(2);
K3(:,:,e)=subs(Function(:,:,e),[Xi Etha],[i(3) j(3)])*DeltaXi(3)*DeltaEtha(3);

Ke(:,:,e)=double(K1(:,:,e)+K2(:,:,e)+K3(:,:,e));
end
 %% Global stiffness matrix
Kg = zeros(2*nnd,2*nnd); % dimensions of the global stiffness matrix

for e = 1:nel;
         I = conn(e,1);
         J = conn(e,2);
         K = conn(e,3);
         L=conn(e,4);
         
 % transfer of the submatrices to their position in the stiffness matrix global
    Kg(2*I-1:2*I,2*I-1:2*I) = Kg(2*I-1:2*I,2*I-1:2*I) + Ke(1:2,1:2,e);
    Kg(2*I-1:2*I,2*J-1:2*J) = Kg(2*I-1:2*I,2*J-1:2*J) + Ke(1:2,3:4,e);
    Kg(2*I-1:2*I,2*K-1:2*K) = Kg(2*I-1:2*I,2*K-1:2*K)+ Ke(1:2,5:6,e);
    Kg(2*I-1:2*I,2*L-1:2*L)=Kg(2*I-1:2*I,2*L-1:2*L)+Ke(1:2,7:8,e);
    
    Kg(2*J-1:2*J,2*I-1:2*I) = Kg(2*J-1:2*J,2*I-1:2*I) + Ke(3:4,1:2,e);
    Kg(2*J-1:2*J,2*J-1:2*J) = Kg(2*J-1:2*J,2*J-1:2*J) + Ke(3:4,3:4,e);
    Kg(2*J-1:2*J,2*K-1:2*K) = Kg(2*J-1:2*J,2*K-1:2*K)+ Ke(3:4,5:6,e);
    Kg(2*J-1:2*J,2*L-1:2*L)=Kg(2*J-1:2*J,2*L-1:2*L)+Ke(3:4,7:8,e);
    
    Kg(2*K-1:2*K,2*I-1:2*I) = Kg(2*K-1:2*K,2*I-1:2*I) + Ke(5:6,1:2,e);
    Kg(2*K-1:2*K,2*J-1:2*J) = Kg(2*K-1:2*K,2*J-1:2*J) + Ke(5:6,3:4,e);
    Kg(2*K-1:2*K,2*K-1:2*K) = Kg(2*K-1:2*K,2*K-1:2*K) + Ke(5:6,5:6,e);
    Kg(2*K-1:2*K,2*L-1:2*L)=Kg(2*K-1:2*K,2*L-1:2*L)+Ke(5:6,7:8,e);

    Kg(2*L-1:2*L,2*I-1:2*I) = Kg(2*L-1:2*L,2*I-1:2*I) + Ke(7:8,1:2,e);
    Kg(2*L-1:2*L,2*J-1:2*J) = Kg(2*L-1:2*L,2*J-1:2*J) + Ke(7:8,3:4,e);
    Kg(2*L-1:2*L,2*K-1:2*K) = Kg(2*L-1:2*L,2*K-1:2*K) + Ke(7:8,5:6,e);
    Kg(2*L-1:2*L,2*L-1:2*L)=Kg(2*L-1:2*L,2*L-1:2*L)+Ke(7:8,7:8,e);   
end
%% Forze di superficie
s=H:-H/NY:0; %ascissa per il carico triangolare agente sulla parete
qtr=45*s/H; % carico triangolare

for e=1:nel
   f_load(:,:,e)=zeros(ndof*nne,1);
end

for i=1:NY
for e=(i-1)*NY+1  % poichè il carico agisce solo sugli elementi di questa serie
    
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
     L=conn(e,4);  % node L 
     
    L=y(L)-y(I);
    
    ql(i)=(qtr(i)-qtr(i+1)) ;  % max carico triangolare a quella quota x i
    
    %forza trapezoidale 
    f_load(:,:,e)=[5/2+((qtr(i+1)/2)+ql(i)/6);0;0; 0;0;0;5/2+((qtr(i+1)/2)+ql(i)/3);0]*t*L;       
end
end

F_load=zeros(ndof*nnd,1);

for e=1:nel
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
     L=conn(e,4);  % node L 
     
    % transfer of the submatrices to their position in the stiffness matrix global
    F_load(2*I-1:2*I,1) = F_load(2*I-1:2*I,1) + f_load(1:2,1,e);
    F_load(2*J-1:2*J,1) = F_load(2*J-1:2*J,1) + f_load(3:4,1,e);
    F_load(2*K-1:2*K,1) = F_load(2*K-1:2*K,1) + f_load(5:6,1,e);
    F_load(2*L-1:2*L,1) = F_load(2*L-1:2*L,1) + f_load(7:8,1,e);  
end   

%% Forze di volume
g=9.80665; %accelerazione di gravità [m/s^2]

   DeltaXi=[0.555 0.888 0.555];
   DeltaEtha=[0.555 0.888 0.555];
    
   i=[-0.775  0  +0.775];
   j=[-0.775  0  +0.775];

for e=1:nel
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
     L=conn(e,4);  % node L 
     
    NN(:,:,e)=subs(N(Xi,Etha),[xI yI xJ yJ xK yK xL yL],[x(I) y(I) x(J) y(J) x(K) y(K) x(L) y(L)]);

%     NN1(:,:,e)=subs(NN(:,:,e)'*[0;-rho*g]*det(Jac(:,:,e)),[Xi Etha],[i(1) j(1)])*DeltaXi(1)*DeltaEtha(1);
%     NN2(:,:,e)=subs(NN(:,:,e)'*[0;-rho*g]*det(Jac(:,:,e)),[Xi Etha],[i(2) j(2)])*DeltaXi(2)*DeltaEtha(2);    
%     NN3(:,:,e)=subs(NN(:,:,e)'*[0;-rho*g]*det(Jac(:,:,e)),[Xi Etha],[i(3) j(3)])*DeltaXi(3)*DeltaEtha(3);
    % Forze di volume
%     f_volume(:,:,e)=t*double( NN1(:,:,e)+ NN2(:,:,e)+ NN3(:,:,e));

% chiedere perchè questo integrale numerico non va bene

   f_volume(:,:,e)=double(int(int(NN(:,:,e)'*[0;-rho*g]*det(Jac(:,:,e)),Xi,[-1 1]),Etha,[-1 1]));
end
F_volume=zeros(ndof*nnd,1);
for e=1:nel;
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
     L=conn(e,4);  % node L     
     
% transfer of the submatrices to their position in the stiffness matrix global
F_volume(2*I-1:2*I,1) = F_volume(2*I-1:2*I,1) + f_volume(1:2,1,e); 
F_volume(2*J-1:2*J,1) = F_volume(2*J-1:2*J,1) + f_volume(3:4,1,e);
F_volume(2*K-1:2*K,1) = F_volume(2*K-1:2*K,1) + f_volume(5:6,1,e);
F_volume(2*L-1:2*L,1) = F_volume(2*L-1:2*L,1) + f_volume(7:8,1,e);
end   
F=F_load+F_volume;
%% Boundary conditions
Kstar = Kg;
Fstar=F;
for e=1:NX+1 % nodi vincolati
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
massa=double(rho*g*((B+b)*H)/2);
Ry=sum(ry);
% verifica 2 % la risultante del carico trapezoidale deve essere uguale
% alla reazione orizzontale
R_trapezio=(50+5)*H/2;
Rx=sum(rx);

for e=1:nel
    
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    L=conn(e,4);    

        u(:,:,e)=[vx(I);...
                    vy(I);...
                    vx(J);...
                    vy(J);...
                    vx(K);...
                    vy(K);...
                    vx(L);...
                    vy(L)];
        
    epslon(:,:,e)=BB(:,:,e)*u(:,:,e);    
    
    sigma(:,:,e)=E*epslon(:,:,e);
    
    sigmaZ(:,:,e)=vu*(sigma(1,:,e)+sigma(2,:,e));
    
    sigma_VM(:,:,e)=sqrt((sigma(1,:,e)^2+sigma(2,:,e)^2+sigmaZ(:,:,e)^2)-(sigma(1,:,e)*sigma(2,:,e)+sigma(2,:,e)*sigmaZ(:,:,e)+sigmaZ(:,:,e)*sigma(1,:,e))+3*sigma(3,:,e)^2);
end

%% cordinate del baricentro dell'elemento
for e=1:nel
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    L=conn(e,4);   
 % coordinate del baricentro   
    xc=(x(I)+x(J)+x(K)+x(L))/4; 
    yc=(y(I)+y(J)+y(K)+y(L))/4;
    
    a=(x(K)-x(I))/2;   b=(y(K)-y(J))/2;
 % Ascisse
 Xi1=(X1-xc)/a;%/a; %compresa tra [-1,1]
 Etha1=(Y1-yc)/b;%/b; %compresa tra [-1,1]

sigmaXY(:,:,e)=subs(sigma(:,:,e),[Xi Etha],[Xi1 Etha1]);
sigmaXY_Z(:,:,e)=subs(sigmaZ(:,:,e),[Xi Etha],[Xi1 Etha1]);
sigmaXY_VM(:,:,e)=subs(sigma_VM(:,:,e),[Xi Etha],[Xi1 Etha1]);
end
%% Plot 
%deformata

for e=1:nel
    fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    L=conn(e,4);    

c=5000%costante di amplificazione
xx2=[x(L)+c*vx(L);x(I)+c*vx(I)];
yy2=[y(L)+c*vy(L);y(I)+c*vy(I)];
S2(i)=plot(xx2,yy2,'g-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','g');
hold on
xx3=[x(I)+c*vx(I);x(J)+c*vx(J);x(K)+c*vx(K);x(L)+c*vx(L)];
yy3=[y(I)+c*vy(I);y(J)+c*vy(J);y(K)+c*vy(K);y(L)+c*vy(L)];
S3(i)=plot(xx3,yy3,'g-o','LineWidth',1.5,'MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','g');
title(['\fontname{Courier}\fontsize{17}Plot Structure & Deformed '],'color','K');
axis([-0.1 1 -0.06  5.8]);
axis equal
axis off
hold on
end

%% Plot the stresses
cc=0.01; % parametro che scala il meshgrid
ccc=[-100:10:100]
% sigma XX

for e=1:nel
    figure(2)
    fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    L=conn(e,4);  
    
[X2,Y2] = meshgrid(x(L):cc:x(J),y(I):cc:y(L));

sigma_X0(X1,Y1)=sigmaXY(1,:,e);
sigma_X=sigma_X0(X2,Y2);
if x(I)==x(L)
sigma_X((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_X((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
sigma_X((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_X((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,sigma_X,'EdgeColor','none','LineWidth',0.1);
clabel(C,h,'FontSize',5,'Color','k');
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
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    L=conn(e,4);  
    
[X2,Y2] = meshgrid(x(L):cc:x(J),y(I):cc:y(L));
sigma_Y0(X1,Y1)=sigmaXY(2,:,e);
sigma_Y=sigma_Y0(X2,Y2);
if x(I)==x(L)
sigma_Y((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_Y((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
sigma_Y((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_Y((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,sigma_Y,'EdgeColor','k','LineWidth',0.1);
clabel(C,h,'FontSize',5,'Color','k');
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
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    L=conn(e,4);  
    
[X2,Y2] = meshgrid(x(L):cc:x(J),y(I):cc:y(L));
tau_XY0(X1,Y1)=sigmaXY(3,:,e);
tau_XY=tau_XY0(X2,Y2);
if x(I)==x(L)
tau_XY((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
tau_XY((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
tau_XY((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
tau_XY((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,tau_XY,'EdgeColor','k','LineWidth',0.1);
clabel(C,h,'FontSize',5,'Color','k');
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
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    L=conn(e,4);  
    
[X2,Y2] = meshgrid(x(L):cc:x(J),y(I):cc:y(L));
sigmaXY_Z0(X1,Y1)=sigmaXY_Z(:,:,e);
sigma_Z=sigmaXY_Z0(X2,Y2);
if x(I)==x(L)
sigma_Z((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_Z((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
sigma_Z((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_Z((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,sigma_Z,'EdgeColor','k','LineWidth',0.1);
clabel(C,h,'FontSize',5,'Color','k');
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
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    L=conn(e,4);  
    
[X2,Y2] = meshgrid(x(L):cc:x(J),y(I):cc:y(L));
sigmaXY_VM0(X1,Y1)=sigmaXY_VM(:,:,e);
sigma_VM1=sigmaXY_VM0(X2,Y2);
if x(I)==x(L)
sigma_VM1((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_VM1((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
sigma_VM1((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_VM1((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,sigma_VM1,'EdgeColor','k','LineWidth',0.1);
clabel(C,h,'FontSize',5,'Color','k');
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