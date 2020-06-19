%% Simple two-bar truss structure
% Computational Mechanics - Prof. Paolo S. Valvo - Lesson 20.03.2017

%% Clear workspace and close any open windows
clear all
close all

%% Symbolic calculations
syms E A L alpha real

KeLoc = E*A/L*[ 1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0 ]
Q = [ cos(alpha) -sin(alpha); sin(alpha) cos(alpha) ]
Te = [ Q zeros(2,2); zeros(2,2) Q ]
Ke = Te*KeLoc*Te'

%% Numerical data
% Geometry

alpha1 = 0.0; % (rad)
alpha2 = pi/6; % (rad)

Ln = 2.00; % (m)
L1 = Ln;
L2 = Ln/cos(alpha2)

D1 = 101.6/1000; % (m)
t1 = 4/1000; % (m)
A1 = pi*D1*t1

D2 = 60.3/1000; % (m)
t2 = 4/1000; % (m)
A2 = pi*D2*t2

% Materials

Es = 210E9; % (Pa)
E1 = Es;
E2 = Es;

alphat = 12e-6; % (1/C)

% Actions
Q = 5000.0; % (N)
t = 30.0; % (C)

% Stiffness matrices
Ke1 = double(subs(Ke, [E A L alpha], [E1 A1 L1 alpha1]))
Ke2 = double(subs(Ke, [E A L alpha], [E2 A2 L2 alpha2]))

K = Ke1(3:4,3:4) + Ke2(3:4,3:4)

pB = [ 0 -Q ]'

uB = inv(K)*pB
