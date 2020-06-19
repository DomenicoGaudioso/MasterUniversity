%%  Computational Mechanics - TERZA ESERCITAZIONE
% stress triangle (CST) element
%% Clear workspace and close any open windows
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
%% Nodes cFoordinates X and      Y
NX=10; % numero di discretizzazioni in x
NY=NX; % numero di discretizzazioni in y

dx=B/NX;  % dimensione lungo x dei triangolini della mesh
dy=H/NY;  % dimensione lungo y dei triangolini della mesh

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
for j=1:NY
    i=1:NX;
    e=(2*(j-1)*NY)+i;
    conn(e,1)=(j-1)*(NY+1)+i;
    conn(e,2)=(j-1)*(NY+1)+i+1;
    conn(e,3)=(j)*(NY+1)+i+1;
    
    e=2*(j-1)*NY+NY+i;
    conn(e,1)=(j)*(NY+1)+i+1;
    conn(e,2)=(j)*(NY+1)+i;
    conn(e,3)=(j-1)*(NY+1)+i;
end  
%% Data structure
nnd = (NX+1)*(NY+1); % Number of nodes in system
nel = length(conn);      % Number of elements
nne = 3;                    % Number of nodes for element
ndof =2;                    % Number of degrees of freedom for node
sdof = nnd*ndof;        % Number of degrees of freedom for element   
%% Elasticity matrix
E=E1/((1+vu)*(1-2*vu))*[1-vu  -vu   0.  ;...  %Matrix elastic %stato piano di deformazione
                                     -vu  1-vu   0.  ;...
                                      0    0 (1-2*vu)/2];                            
%% Forze di superficie
s=H:-H/NY:0; %ascissa per il carico triangolare agente sulla parete
qtr=45*s/H; % carico triangolare

for e=1:nel
   f_load(:,:,e)=zeros(6,1);
end

for i=1:NY
for e=2*(i-1)*NY+NY+1  % poichè il carico agisce solo sugli elementi di questa serie
    
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
     
    L=y(J)-y(K);
    
    ql(i)=(qtr(i)-qtr(i+1)) ;  % max carico triangolare a quella quota x i
    
    %forza  distribuiton alle varie quote senza la componente triangolare 
    f_load(:,:,e)=[0;0;5/2+((qtr(i+1)/2)+ql(i)/6); 0;5/2+((qtr(i+1)/2)+ql(i)/3);0]*t*L;      
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

%% Forze di volume
g=9.80665; %accelerazione di gravità [m/s^2]

for e=1:nel
     I = conn(e,1); % node I
     J = conn(e,2); % node J
     K= conn(e,3); % node K
    
    % Area elements
    A(e)=1/2*abs(det([1 x(I) y(I); 1  x(J)  y(J); 1  x(K)  y(K)]));
    
    % Forze di volume
    f_volume(:,:,e)=-(t/3)*[0; rho*g*A(e);0;rho*g*A(e);0;rho*g*A(e)];
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
%% matrix stiffnes composition              
for e=1:nel
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);
     
%  costant in Shape functions
    aI=x(J)*y(K) - x(K)*y(J);
    bI=y(J)-y(K);
    cI=x(K)-x(J);

    aJ =x(K)*y(I) - x(I)*y(K);
    bJ=y(K)-y(I);
    cJ=x(I)-x(K);

    aK=x(I)*y(J) - x(J)*y(I);
    bK=y(I)-y(J);
    cK=x(J)-x(I);

    % Strain-displacement matrices
    BI(:,:,e)= [ bI 0; 0 cI; cI bI ]/(2*A(e));
    BJ(:,:,e) = [ bJ 0; 0 cJ; cJ bJ ]/(2*A(e));
    BK(:,:,e) = [ bK 0; 0 cK; cK bK ]/(2*A(e));
    BB(:,:,e)= [ BI(:,:,e) BJ(:,:,e) BK(:,:,e) ];
    
    % Stiffness matrix
    Ke(:,:,e) = BB(:,:,e)'*E*BB(:,:,e)*A(e)*t;
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
massa=rho*g*((B+b)*H)/2;
Ry=sum(ry);
% verifica 2 % la risultante del carico trapezoidale deve essere uguale
% alla reazione orizzontale
R_trapezio=(50+5)*H/2
Rx=sum(rx)

for e=1:nel
    
    I = conn(e,1);
    J = conn(e,2);
    K = conn(e,3);

        u(:,:,e)=[vx(I);...
                    vy(I);...
                    vx(J);...
                    vy(J);...
                    vx(K);...
                    vy(K)];
        
    epslon(:,:,e)=BB(:,:,e)*u(:,:,e) ;    
    
    sigma(:,:,e)=E*epslon(:,:,e);
    
    sigmaZ(:,:,e)=vu*(sigma(1,:,e)+sigma(2,:,e));
    
    sigma_VM(:,:,e)=sqrt((sigma(1,:,e)^2+sigma(2,:,e)^2+sigmaZ(:,:,e)^2)-(sigma(1,:,e)*sigma(2,:,e)+sigma(2,:,e)*sigmaZ(:,:,e)+sigmaZ(:,:,e)*sigma(1,:,e))+3*sigma(3,:,e)^2);
end
%% TRIANGULATION
figure(6)
W=[x',y'];
T=[conn];
TR = triangulation(T,W);
TR.Points;
TR.ConnectivityList;
IC = incenter(TR); % centro del triangolo

%% Mappaggio
Mappax=containers.Map();
Mappay=containers.Map();

for i=1:nel
    Centro=IC(i,:);
    stress_x=sigma(1,i);
    stress_x=int32(stress_x/10)*10;
    Mappax(num2str(stress_x))=[];
    Mappay(num2str(stress_x))=[];
end
    
for i=1:nel;
    Centro=IC(i,:);
    stress_x=sigma(1,i);
    stress_x=int32(stress_x/10)*10;
    Mappax(num2str(stress_x)) = [Mappax(num2str(stress_x)), Centro(1)];
    Mappay(num2str(stress_x)) = [Mappay(num2str(stress_x)), Centro(2)];
end
figure(8)
chiavi = keys(Mappax); ;% sono anche le chiavi di Mappay

 hb=bar(1:10,[1:10;2:11]');
 get(hb,'DisplayName');

for i=1:length(chiavi);
    chiave= cell2mat(chiavi(i));
    sz=10;
    scatter(Mappax(chiave), Mappay(chiave),sz,'filled');
    title(['\fontname{Courier}\fontsize{15} Mapping \sigma_x_x '],'color','K');
    axis([-0.1 1 -0.06  5.8]);
    axis equal
    axis off
    hold on
end
% una mappa è un array ma ansichè accederlo solo con indici, si può accedere anche con stringhe
%come nel nostro caso
%% POST PROCESSING
% disegno struttura
figure(1)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
triplot(TR,'LineWidth',0.1,'Color',[0.1 0.2 0])
title(['\fontname{Courier}\fontsize{17}Plot Structure & Deformed '],'color','K');
axis([-0.1 1 -0.06  5.8]);
axis equal
hold on
% Print force
a=0.050;
X=[-0.3 -0.3 -0.3-a*5 -0.3-a*50];
Y=[0 H H 0];
hpatch=patch(X,Y,'w','Facealpha',0.1,'EdgeColor', [0.8 0.2 0]);

for s=0:H/NY:H
j=a*((45/H)*(H-s)+5);
plot(-0.33,s,'>','MarkerSize',2,'MarkerEdgeColor', [1 0.2 0],'MarkerFaceColor', [1 0.2 0])
p0=[-0.3-j s];
p1=[-0.3 s];
vectarrow(p0,p1);
axis off
hold on
end

c=20000; % costante di amplificazione della deformata
%%print deformata
xd=x'+c*vx;
yd=y'+c*vy;

WW=[xd,yd];
TT=[conn];
TR1=triangulation(TT,WW);
triplot(TR1,'LineWidth',0.1,'Color','g')
Bordo=freeBoundary(TR1)
axis equal
hold on
plot(xd(Bordo),yd(Bordo),'Color',[0 0.5 0.1],'LineWidth',2)
%% Plot the stresses in the x_direction
c1=15 % controllo le curve di livello che voglio
geom=[x' y'];
figure(2)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
x_stress=sigma(1,:);
cmin=min(x_stress);
cmax=max(x_stress);
caxis([ cmin cmax]);
patch('Faces', conn, 'Vertices',geom, 'FaceVertexCData', x_stress',...
    'Facecolor','flat','LineWidth',0.1,'EdgeColor','none');
hold on
c=colorbar
colormap(jet(c1)) 
c.Label.String = '\sigma_x_x  [KN/m^2]'
Bordo2=freeBoundary(TR)
plot(x(Bordo2),y(Bordo2),'Color','k','LineWidth',1)
title(['\fontname{Courier}\fontsize{15} Stresses in the x direction (\sigma_x_x)'],'color','K');
axis([-0.1 1 -0.06  5.8]);
axis equal
axis off
hold on
%% Plot the stresses in the y_direction
figure(3)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
y_stress=sigma(2,:);
cmin=min(y_stress);
cmax=max(y_stress);
caxis([ cmin cmax]);
patch('Faces', conn, 'Vertices',geom, 'FaceVertexCData', y_stress',...
    'Facecolor','flat','LineWidth',0.1,'EdgeColor','none');
hold on
c=colorbar
colormap(jet(c1))  %fa tipo delle curve di livello
c.Label.String = '\sigma_y_y  [KN/m^2]'
plot(x(Bordo2),y(Bordo2),'Color','k','LineWidth',1)
title(['\fontname{Courier}\fontsize{15}Stresses in the y direction (\sigma_y_y) '],'color','K');
axis([-0.1 1 -0.06  5.8]);
axis equal
axis off
hold on
%% Plot the stresses in the xy_direction
figure(4)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
xy_stress=sigma(3,:);
cmin=min(xy_stress);
cmax=max(xy_stress);
caxis([cmin cmax]);
patch('Faces', conn, 'Vertices',geom, 'FaceVertexCData', xy_stress',...
    'Facecolor','flat','LineWidth',0.1,'EdgeColor','none');
hold on
c=colorbar
colormap(jet(c1)) 
c.Label.String = '\tau_x_y  [KN/m^2]'
plot(x(Bordo2),y(Bordo2),'Color','k','LineWidth',1)
title(['\fontname{Courier}\fontsize{15}Stresses in the xy direction (\tau_x_y) '],'color','K');
axis([-0.1 1 -0.06  5.8]);
axis equal
axis off
hold on
%% Plot the stresses in the z_direction
figure(5)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
z_stress=sigmaZ(1,:)
cmin=min(z_stress);
cmax=max(z_stress);
caxis([ cmin cmax]);
patch('Faces', conn, 'Vertices',geom, 'FaceVertexCData', z_stress',...
    'Facecolor','flat','LineWidth',0.1,'EdgeColor','none');
hold on
c=colorbar
colormap(jet(c1)) 
c.Label.String = '\sigma_z [KN/m^2]'
plot(x(Bordo2),y(Bordo2),'Color','k','LineWidth',1)
title(['\fontname{Courier}\fontsize{15}Stresses in the z direction (\sigma_z) '],'color','K');
axis([-0.1 1 -0.06  5.8]);
axis equal
axis off
hold on
%% Plot the Von Mises stresses
figure(6)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
vm_stress=sigma_VM(1,:)
cmin=min(vm_stress);
cmax=max(vm_stress);
caxis([ cmin cmax]);
patch('Faces', conn, 'Vertices',geom, 'FaceVertexCData', vm_stress',...
    'Facecolor','flat','LineWidth',0.1,'EdgeColor','none');
hold on
c=colorbar
colormap(jet(c1)) 
c.Label.String = '\sigma_V_M [KN/m^2]'
plot(x(Bordo2),y(Bordo2),'Color','k','LineWidth',1)
title(['\fontname{Courier}\fontsize{15}Von Mises stresses (\sigma_V_M) '],'color','K');
axis([-0.1 1 -0.06  5.8]);
axis equal
axis off
hold on