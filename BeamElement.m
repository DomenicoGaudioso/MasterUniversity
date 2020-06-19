%% Formulation of plane beam element
% Computational Mechanics - Prof. Paolo S. Valvo - Lesson 24.04.2018

%% Beam element

clear all
close all

syms H L EA EI x px py mz alphat tminus tplus tm Deltat real

%% Shape functions

NIx = 1 - x/L;
NIy = 1 - 3*x^2/L^2 + 2*x^3/L^3;
NIz = x - 2*x^2/L + x^3/L^2;

NJx = x/L;
NJy = 3*x^2/L^2 - 2*x^3/L^3;
NJz = -x^2/L + x^3/L^2;

%% Element matrices

N(x) = [ NIx 0 0 NJx 0 0; ...
         0 NIy NIz 0 NJy NJz; ...
         0 diff(NIy,x) diff(NIz,x) 0 diff(NJy,x) diff(NJz,x) ]

B(x) = [ diff(NIx,x) 0 0 diff(NJx,x) 0 0; ...
         0 diff(NIy,x,2) diff(NIz,x,2) 0 diff(NJy,x,2) diff(NJz,x,2) ]

D = [ EA 0; 0 EI ]

%% Stiffness matrix

KLoc = int(B(x)'*D*B(x),x,0,L)

%% Uniform load equivalent nodal forces

p = [ px py mz ]'
fload = -int(N(x)'*p,x,0,L)

%% Linear temperature nodal forces

% tm = (tminus + tplus)/2;
% Deltat = tplus - tminus;
dt = [ alphat*tm -alphat*Deltat/H ]';
ftemp = -int(B(x)'*D*dt,x,0,L)

