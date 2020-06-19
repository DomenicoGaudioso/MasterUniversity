%%  Computational Mechanics - SECONDA ESERCITAZIONE
% 3D TRUSS STRUCTURE
%% Clear workspace and close any open windows
clear all
close all
%% Set Command Window output display format
format bank;
%% data of the problem                      
a=3;        % [m] length 
E=21E7;                               % Young's modulus steel [KN/m^2]
%% data structure
np=2  % Number of level;
nnd = 4*(np+1);    % Number of nodes:
c=8 % numero di elementi per piano
nel = c*np;              % Number of elements:
nne =2 ;              % Number of nodes for element:
nodof =3 ;             % Number of degrees of freedom for node
eldof = nne*nodof;     % Number of degrees of freedom for element
%% Nodes cFoordinates X and      Y
for i=0:1:np
    x(4*i+1)=0; z(4*i+1)=0; y(4*i+1)=a*i;
    x(4*i+2)=0; z(4*i+2)=a; y(4*i+2)=a*i;
    x(4*i+3)=a; z(4*i+3)=a; y(4*i+3)=a*i;
    x(4*i+4)=a; z(4*i+4)=0; y(4*i+4)=a*i;
end
%% Element connectivity
for n=0:np
    % elementi verticali
conn(c*n+1,1)=4*n+1; conn(c*n+1,2)=4*n+5; 
conn(c*n+2,1)=4*n+2; conn(c*n+2,2)=4*n+6; 
conn(c*n+3,1)=4*n+3; conn(c*n+3,2)=4*n+7; 
conn(c*n+4,1)=4*n+4; conn(c*n+4,2)=4*n+8; 
    % elementi orizzontali
conn(c*n+5,1)=4*n+5; conn(c*n+5,2)=4*n+6; 
conn(c*n+6,1)=4*n+6; conn(c*n+6,2)=4*n+7; 
conn(c*n+7,1)=4*n+7; conn(c*n+7,2)=4*n+8; 
conn(c*n+8,1)=4*n+8; conn(c*n+8,2)=4*n+5; 
    % elementi diagonali
% conn(12*n+9,1)=4*n+1; conn(12*n+9,2)=4*n+6; 
% conn(12*n+10,1)=4*n+2; conn(12*n+10,2)=4*n+7; 
% conn(12*n+11,1)=4*n+3; conn(12*n+11,2)=4*n+8; 
% conn(12*n+12,1)=4*n+4; conn(12*n+12,2)=4*n+5; 
end
%% Plot structure 
figure(1)
pos0 = [0.06 0.03 0.94 0.94];
subplot('Position',pos0)
for i=1:nel
node_1=conn(i,1); % node start
node_2=conn(i,2); % node end
xx=[x(node_1),x(node_2)];
yy=[y(node_1),y(node_2)];
zz=[z(node_1),z(node_2)];
S(i)=plot3(xx,yy,zz,'K-o','LineWidth',1.5,'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5],'DisplayName','structure');
axis equal
hold on
end
%% Area section [m^2]
A_c=0.005381;    % Area correnti orizzontali_IPE300 [m^2]

A=ones(nel,1)*A_c;
%% Element properties
for e = 1:nel
    I = conn(e,1); % node start
    J = conn(e,2); % node end
    % Element length
    L(e) = sqrt((x(J)-x(I))^2 + (y(J)-y(I))^2 +(z(J)-z(I))^2);
   %versori 
   cosxX=(x(J)-x(I))/L(e); cosxY=(y(J)-y(I))/L(e); cosxZ=(z(J)-z(I))/L(e);
   
   xp=(x(J)-x(I))/2+2; yp=(y(J)-y(I))/2+2; zp=(z(J)-z(I))/2; % punto p nel piano xz
   
   zx=(y(J)-y(I))*(zp-z(I))-(z(J)-z(I))*(yp-y(I));
   zy=(z(J)-z(I))*(xp-x(I))-(x(J)-x(I))*(zp-z(I));
   zz=(x(J)-x(I))*(yp-y(I))-(y(J)-y(I))*(xp-x(I));
   Z=sqrt(zx^2 +zy^2+zz^2);
   coszX=zx/Z; coszY=zy/Z; coszZ=zz/Z;
   
   yx=cosxY*coszZ-cosxZ*coszY;
   yy=cosxZ*coszX-cosxX*coszZ;
   yz=cosxX*coszY-cosxY*coszX;
   Y=sqrt(yx^2 +yy^2+yz^2);
   cosyX=yx/Y; cosyY=yy/Y; cosyZ=yz/Y;
   
   % rotation matrix
   Q(:,:,e) = [ cosxX  cosxY  cosxZ;...
                   cosyX   cosyY  cosyZ;...
                   coszX   coszY   coszZ];
    % stiffness matrix in local
    KLoc(:,:,e)=((E*A(e))/L(e))*[1  0  0  -1  0  0;...
                                              0  0  0   0  0  0;...
                                              0  0  0   0  0  0;...
                                             -1  0  0   1  0  0;...
                                               0  0  0   0  0  0;...
                                               0  0  0   0  0  0] ;
    % transfer matrix
    T(:,:,e)  = [ Q(:,:,e)  zeros(3,3); zeros(3,3) Q(:,:,e)];
    % stiffness matrix in global
    KGlob(:,:,e) = T(:,:,e)*KLoc(:,:,e)*T(:,:,e)';
end
%% Global stiffness matrix
K = zeros(nodof*nnd,nodof*nnd); % dimensions of the global stiffness matrix
for e = 1:nel;
    I = conn(e,1); % node start
    J = conn(e,2); % node end
    
 % transfer of the submatrices to their position in the stiffness matrix
 % global
    K(3*I-2:3*I,3*I-2:3*I) = K(3*I-2:3*I,3*I-2:3*I)+KGlob(1:3,1:3,e);
    K(3*I-2:3*I,3*J-2:3*J) = K(3*I-2:3*I,3*J-2:3*J)+KGlob(1:3,4:6,e);
    K(3*J-2:3*J,3*I-2:3*I) = K(3*J-2:3*J,3*I-2:3*I)+KGlob(4:6,1:3,e);
    K(3*J-2:3*J,3*J-2:3*J) = K(3*J-2:3*J,3*J-2:3*J)+KGlob(4:6,4:6,e);
end
%% Load vector in [KN]
for e=1:nel
   p(:,:,e)=zeros(6,1);
end
    %forza  distribuiton alle varie quote senza la componente triangolare 
    p(:,:,15)=[50;60;0;50;50;0];
    p(:,:,np*5)=[0;0;-70;0;0;-70];
    p(:,:,np*7)=[0;0;-70;0;0;-70];
    
P=p; % Load vector

F_load=zeros(nodof*nnd,1);
for e=1:nel
    I = conn(e,1); % node start
    J = conn(e,2); % node end
     
    % transfer of the submatrices to their position in the stiffness matrix global
    F_load(3*I-2:3*I,1) = F_load(3*I-2:3*I,1) + P(1:3,1,e);

    F_load(3*J-2:3*J,1) = F_load(3*J-2:3*J,1) + P(4:6,1,e);
end 
%% Boundary conditions
Kstar = K;
Fstar=F_load;

for e=1:4 % nodi vincolati
    %annullo ux
    Kstar(e*nodof-2,:)=0;
    Kstar(:,e*nodof-2)=0;
    %annullo uy
    Kstar(e*nodof-1,:)=0;
    Kstar(:,e*nodof-1)=0;
    %annullo uz
    Kstar(e*nodof,:)=0;
    Kstar(:,e*nodof)=0;
    
    Kstar(e*nodof-2,e*nodof-2)=1; % inserisco 1 ad ux
    Kstar(e*nodof-1,e*nodof-1)=1; % inserisco 1 ad uy
    Kstar(e*nodof,e*nodof)=1; % inserisco uno 1ad uz

%     Fstar(e*nodof-1)=0;
%     Fstar(e*nodof)=0;
end
% Kstar(1:12,:) = zeros(4,2*nnd);
% Kstar(:,1:12) = zeros(2*nnd,4);
% primo incastro
% Kstar(1,1) = 1; %ux bloccata
% Kstar(2,2) = 1; %uy bloccata
% Kstar(3,3) = 1; % rotazione bloccata
% % secondo incastro
% Kstar(4,4) = 1;
% Kstar(5,5) = 1;
% Kstar(6,6) = 1;
% % terzo incastro
% Kstar(7,7) = 1;
% Kstar(8,8) = 1;
% Kstar(9,9) = 1;
% % quarto incastro
% Kstar(10,10) = 1;
% Kstar(11,11) = 1;
% Kstar(12,12) = 1;
%% Solution
v= inv(Kstar)*Fstar;  % nodal diplacements 
f= K*v;               % reaction forces in the nodes
r =(f-Fstar);         % reaction force in external constraints
% %% Normal forces


