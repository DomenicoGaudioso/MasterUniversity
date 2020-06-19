%% Formulation of 3D beam element
% Computational Mechanics - Prof. Paolo S. Valvo - Lesson 24.04.2018

%% Beam element

clear all
close all

syms L EA EIy EIz GJt x real

%% Shape functions

NIx = 1 - x/L;
NIy = 1 - 3*x^2/L^2 + 2*x^3/L^3;
NIz = x - 2*x^2/L + x^3/L^2;

NJx = x/L;
NJy = 3*x^2/L^2 - 2*x^3/L^3;
NJz = -x^2/L + x^3/L^2;

%% Element matrices

% u(x) = [ ux uy uz phix phiy phiz ]'
%
% d = [ epsilon chiy chiz thetax ]'
%
% S = [ d/dx 0      0      0    0 0; ...
%       0    0     -d2/dx2 0    0 0; ...
%       0    d2/dx2 0      0    0 0; ...
%       0    0      0      d/dx 0 0 ]

N(x) = [ NIx 0   0   0   0   0   NJx 0   0   0   0   0; ...
         0   NIy 0   0   0   NIz 0   NJy 0   0   0   NJz; ...
         0   0   NIy 0  -NIz 0   0   0   NJy 0  -NJz 0; ...
	     0   0   0   NIx 0   0   0   0   0   NJx 0   0; ...
         0  -diff(NIy,x) 0 0 0 -diff(NIz,x) 0 -diff(NJy,x) 0 0 0 -diff(NJz,x); ...
         0   0 diff(NIy,x) 0 diff(NIz,x) 0 0 0 diff(NJy,x) 0 diff(NJz,x) 0 ]

B(x) = [ diff(NIx,x) 0 0 0 0 0 diff(NJx,x) 0 0 0 0 0; ...
         0 0 -diff(NIy,x,2) 0 diff(NIz,x,2) 0 0 0 -diff(NJy,x,2) 0 diff(NJz,x,2) 0; ...
         0 diff(NIy,x,2) 0 0 0 diff(NIz,x,2) 0 diff(NJy,x,2) 0 0 0 diff(NJz,x,2); ...
         0 0 0 diff(NIx,x) 0 0 0 0 0 diff(NJx,x) 0 0 ]

D = [ EA 0   0   0; ...
      0  EIy 0   0; ...
      0  0   EIz 0; ...
      0  0   0   GJt ]

%% Stiffness matrix

KLoc = int(B(x)'*D*B(x),x,0,L)

