%% Truss structure
% Computational Mechanics - Prof. Paolo S. Valvo - Lessons 10-17.04.2018

clear all
close all

%% Problem data

n = 4;           % Number of nodes
m = 4;           % Number of bars
a = 5.0;         % (m) Length
h=3.5            % altezza
% E=21E7;          % Young's modulus steel [KN/m^2]
syms A  E real
%% Nodal coords

x(1) = 0.0;     y(1) = 0.0;
x(2) = a;       y(2) = 0.0;
x(3) = 0.0;     y(3) = h;
x(4) = a;       y(4) = h;

%% Load vector

p = zeros(2*n,1)
p(5) = 300;

%% Connectivity matrix

conn(1,1) = 1;  conn(1,2) = 3;
conn(2,1) = 1;  conn(2,2) = 4;
conn(3,1) = 2;  conn(3,2) = 4;
conn(4,1) = 3;  conn(4,2) = 4;

%% Element properties

for e = 1:m
    I = conn(e,1);
    J = conn(e,2);
    L(e) = sqrt((x(J)-x(I))^2 + (y(J)-y(I))^2);
    alpha(e) = atan2(y(J)-y(I), x(J)-x(I));
    KLoc(:,:,e) = E*A/L(e)*[ 1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0 ];
    Q = [ cos(alpha(e)) -sin(alpha(e)); sin(alpha(e)) cos(alpha(e)) ];
    T = [ Q zeros(2,2); zeros(2,2) Q ];
    KGlob(:,:,e) = T*KLoc(:,:,e)*T';
end
digits(2)
vpa(KGlob(:,:,3))
%% Global stiffness matrix

K = zeros(2*n,2*n);

for e = 1:m
    I = conn(e,1);
    J = conn(e,2);

    K(2*I-1:2*I,2*I-1:2*I) = K(2*I-1:2*I,2*I-1:2*I) + KGlob(1:2,1:2,e);
    K(2*I-1:2*I,2*J-1:2*J) = K(2*I-1:2*I,2*J-1:2*J) + KGlob(1:2,3:4,e);
    K(2*J-1:2*J,2*I-1:2*I) = K(2*J-1:2*J,2*I-1:2*I) + KGlob(3:4,1:2,e);
    K(2*J-1:2*J,2*J-1:2*J) = K(2*J-1:2*J,2*J-1:2*J) + KGlob(3:4,3:4,e);

end
    
%% Restraints

Kstar = K;

Kstar(1:4,:) = zeros(4,2*n);
Kstar(:,1:4) = zeros(2*n,4);

Kstar(1,1) = 1;
Kstar(2,2) = 1;
Kstar(3,3) = 1;
Kstar(4,4) = 1;

pstar = p;
%% Solution
u = inv(Kstar)*pstar

f = K*u

r = f - p