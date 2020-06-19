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
NX=2; % numero di discretizzazioni in x
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
 for e=1:nel
figure(1)
    fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none'; 
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
 %% simbolico
 syms   X Y  xI yI xJ yJ xK yK xL yL   real
%% Elasticity matrix
E=E1/((1+vu)*(1-2*vu))*[1-vu  -vu   0.  ;...  %Matrix elastic %stato piano di deformazione
                                 -vu  1-vu   0.  ;...
                                  0    0 (1-2*vu)/2];  
%% Shape functions                                 
for e=1:nel
 I = conn(e,1);
 J = conn(e,2);
 K = conn(e,3);
 L=conn(e,4);
    
P = [ 1 X Y  X.*Y 0 0 0 0 ;...
      0 0 0 0 1 X Y  X.*Y];    

C= [ 1 x(I) y(I) x(I).*y(I)  0 0 0 0 ;              % Shape function calulation
      0 0 0 0 1 x(I) y(I) x(I).*y(I);
      1 x(J) y(J) x(J).*y(J)  0 0 0 0 ;              
      0 0 0 0 1 x(J) y(J) x(J).*y(J) ;
      1 x(K) y(K) x(K).*y(K)  0 0 0 0 ;              
      0 0 0 0 1 x(K) y(K) x(K).*y(K)   ;
      1 x(L) y(L) x(L).*y(L)  0 0 0 0 ;              
      0 0 0 0 1 x(L) y(L) x(L).*y(L) ];

  N=P*inv(C);
 
  NI(:,:,e)=N(1,1); NJ(:,:,e)=N(1,3); NK(:,:,e)=N(1,5); NL(:,:,e)=N(1,7); 
  
  Nglobal(:,:,e)=[NI(:,:,e) 0 NJ(:,:,e) 0 NK(:,:,e) 0 NL(:,:,e) 0;...
                       0 NI(:,:,e) 0 NJ(:,:,e) 0 NK(:,:,e) 0 NL(:,:,e)];
  
  BI(:,:,e)=[diff(NI(:,:,e),X) 0; 0 diff(NI(:,:,e),Y); diff(NI(:,:,e),Y) diff(NI(:,:,e),X)];
  BJ(:,:,e)=[diff(NJ(:,:,e),X) 0; 0 diff(NJ(:,:,e),Y); diff(NJ(:,:,e),Y) diff(NJ(:,:,e),X)];
  BK(:,:,e)=[diff(NK(:,:,e),X) 0; 0 diff(NK(:,:,e),Y); diff(NK(:,:,e),Y) diff(NK(:,:,e),X)];
  BL(:,:,e)=[diff(NL(:,:,e),X) 0; 0 diff(NL(:,:,e),Y); diff(NL(:,:,e),Y) diff(NL(:,:,e),X)];
  
  BB(:,:,e)=[BI(:,:,e) BJ(:,:,e) BK(:,:,e) BL(:,:,e) ];
  
  K0(:,:,e)=BB(:,:,e)'*E*BB(:,:,e);
  
  % integrale sull'elemento genitore
x0=x(L)+(x(I)-x(L))/2; x1=x(K)+(x(J)-x(K))/2;

Ke0(:,:,e)=int(K0(:,:,e),Y,[y(I) y(K)]);
Ke(:,:,e)=int(Ke0(:,:,e),X,[x0 x1]);
end
Ke(:,:,e)=double(Ke(:,:,e));
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
     
  YY2=simplify(y(J)+(X-x(J))*(y(K)-y(J))/(x(K)-x(J)));
  YY1=simplify(y(I)+(X-x(I))*(y(L)-y(I))/(x(L)-x(I)));
  
  f_volume1(:,:,e)=Nglobal(:,:,e)'*[0;-rho*g];
% integrale sull'elemento genitore
x0=x(L)+(x(I)-x(L))/2; x1=x(K)+(x(J)-x(K))/2;

f_volume0(:,:,e)=int(f_volume1(:,:,e),Y,[y(I) y(K)]);
f_volume(:,:,e)=int(f_volume0(:,:,e),X,[x0 x1]);
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
figure(1)    
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
%% Plot 
%deformata

for e=1:nel
figure(1)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';  

    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
    L=conn(e,4);    

c=5000; %costante di amplificazione
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
cc1=8; %parametro che scala le linee delle tensioni
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

sigma_X0(X,Y)=sigma(1,:,e);
sigma_X=sigma_X0(X2,Y2);
if x(I)==x(L)
sigma_X((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_X((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
sigma_X((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_X((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,sigma_X,(x(J)-x(I))*cc1,'EdgeColor','K','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','k');
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
sigma_Y0(X,Y)=sigma(2,:,e);
sigma_Y=sigma_Y0(X2,Y2);
if x(I)==x(L)
sigma_Y((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_Y((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
sigma_Y((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_Y((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,sigma_Y,(x(J)-x(I))*cc1,'EdgeColor','k','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','k');
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
tau_XY0(X,Y)=sigma(3,:,e);
tau_XY=tau_XY0(X2,Y2);
if x(I)==x(L)
tau_XY((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
tau_XY((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
tau_XY((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
tau_XY((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,tau_XY,(x(J)-x(I))*cc1,'EdgeColor','k','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','k');
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
sigma_Z0(X,Y)=sigmaZ(:,:,e);
sigma_Z=sigma_Z0(X2,Y2);
if x(I)==x(L)
sigma_Z((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_Z((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
sigma_Z((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_Z((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,sigma_Z,(x(J)-x(I))*cc1,'EdgeColor','k','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','k');
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
sigmaXY_VM0(X,Y)=sigma_VM(:,:,e);
sigma_VM1=sigmaXY_VM0(X2,Y2);
if x(I)==x(L)
sigma_VM1((X2-x(I))/(x(L)-x(I))<(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_VM1((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan; 
else
sigma_VM1((X2-x(I))/(x(L)-x(I))>(Y2-y(I))/(y(L)-y(I))) = nan;
sigma_VM1((X2-x(J))/(x(K)-x(J))<(Y2-y(J))/(y(K)-y(J))) = nan;
end
[C,h] = contourf(X2,Y2,sigma_VM1,(x(J)-x(I))*cc1,'EdgeColor','k','LineWidth',0.1);
%clabel(C,h,'FontSize',5,'Color','k');
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