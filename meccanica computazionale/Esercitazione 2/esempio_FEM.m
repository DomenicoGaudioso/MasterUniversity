%% TRUSS STRUCTURE
clear all
close all
%% Problem data
n=4            % numero di nodi
m=4            % numero di aste
h=3.5          % width [m]
l=5            % height
Q=50           % load [KN]
E = 21E7;      % Young's modulus steel [KN/m^2] 
A=20e-4        % Cross-section area [m^2]
%% Nodal coords
x(1)=0.0;  y(1)=0.0;
x(2)=l;    y(2)=0.0;
x(3)=0.0;  y(3)=h;
x(4)=l;    y(4)=h;
%% Load vector
p=zeros(2*n,1)
p(5)=Q
%% Connectivity matrix
conn(1,1)=1; conn(1,2)=3;  % Asta 1
conn(2,1)=1; conn(2,2)=4;  % Asta 2
conn(3,1)=2; conn(3,2)=4;  % Asta 3
conn(4,1)=3; conn(4,2)=4;  % Asta 4
%% Element properties
for e=1:m
    i=conn(e,1);
    j=conn(e,2);
    L(e)=sqrt((x(j)-x(i)).^2+(y(j)-y(i)).^2);
    alpha(e)=atan2(y(j)-y(i),x(j)-x(i));                                 % rotazione dell x segnato rispetto alla x del riferimento globale
    Kloc(:,:,e)=((E*A)/L(e))*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0]     % matrice rigidezza elastica in locale
    Q=[cos(alpha(e)) -sin(alpha(e));sin(alpha(e)) cos(alpha(e))]        % matrice di rotazione
    T=[ Q zeros(2,2); zeros(2,2) Q];                                    % matrice trasferimento
    Kglob(:,:,e)=T*Kloc(:,:,e)*T'                                           % matrice di rigidezza globale
end
%% Global stiffness matrix
K=zeros(2*n*2*n) % final dimension of the matrix
for e=1:m
    i=conn(e,1);
    j=conn(e,2);
    
   K(2*i-1:2*i,2*i-1:2*i)=K(2*i-1:2*i,2*i-1:2*i)+Kglob(1:2,1:2,e); 
   K(2*i-1:2*i,2*j-1:2*j)=K(2*i-1:2*i,2*j-1:2*j)+Kglob(1:2,3:4,e); 
   K(2*j-1:2*j,2*i-1:2*i)=K(2*j-1:2*j,2*i-1:2*i)+Kglob(3:4,1:2,e); 
   K(2*j-1:2*j,2*j-1:2*j)=K(2*j-1:2*j,2*j-1:2*j)+Kglob(3:4,3:4,e); 
end
%% Restraints
Kstar=K;

Kstar(1:4,:)=zeros(4,2*n);
Kstar(:,1:4)=zeros(2*n,4);

Kstar(1,1)=1;
Kstar(2,2)=1;
Kstar(3,3)=1;
Kstar(4,4)=1;

pstar=p;
pstar(1:4)=0;

%% solution
u=inv(Kstar)*pstar;% spostamenti effettivi di tutti i nodi del sistema

%% Reactions
f=K*u;
r=f-p