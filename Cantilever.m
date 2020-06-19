%% Uniformly loaded cantilever beam solved by different approaches
% Computational Mechanics - Prof. Paolo S. Valvo - Lesson 06.03.2017

%% Clear workspace and close any open windows
clear all
close all

%% Problem data
E = 210e9       % (N/m^2) Young's modulus
L = 2;          % (m) length
I = 317.8e-8    % (m^4) moment of inertia (IPE 120)
q = 3000.0      % (N/m) distributed load
s = 0:L/100:L   % (m) abscissa (for plots)

%% Exact solution
vex = (1/24*q*s.^4-1/6*q*L*s.^3+1/4*q*L^2*s.^2)/(E*I);

%% Rigid-segment discretisation
a = L/4
Q = q*a
k0 = E*I/a
A = [3*k0 -k0 0 0; -k0 2*k0 -k0 0; 0 -k0 2*k0 -k0; 0 0 -k0 k0]
b = [7/2; 5/2; 3/2; 1/2]*Q*a
theta = A\b
vrig = zeros(1,101)
for i = 1:26
    vrig(i) = theta(1)*s(i)
end
for i = 27:101
    vrig(i) = vrig(26) + theta(2)*(s(i)-a)
end
for i = 52:101
    vrig(i) = vrig(51) + theta(3)*(s(i)-2*a)
end
for i = 77:101
    vrig(i) = vrig(76) + theta(4)*(s(i)-3*a)
end

%% 1st-order Rayleigh-Ritz solution
c = 32*(1-2/pi)/pi^4*q*L^4/(E*I)
vRR1 = c*(1-cos(pi*s/(2*L)))

%% 2nd-order Rayleigh-Ritz solution
c1 = 96/pi^3*(3*pi-7)/(9*pi^2-16)*q*L^4/(E*I)
c2 = 6/pi^4*(3*pi^2-16*pi+32)/(9*pi^2-16)*q*L^4/(E*I)
vRR2 = c1*(1-cos(pi*s/(2*L)))+c2*(1-cos(pi*s/L))

%% Plot deformed shapes
plot(s,vex,'r-','LineWidth',2)
title('Cantilever beam deflection')
xlabel('s (m)')
ylabel('v(s) m')
axis ij
hold on

plot(s,vRR1,'b.','LineWidth',2)
plot(s,vRR2,'-.','LineWidth',2,'color',[0 0.5 0])
plot(s,vrig,'--','LineWidth',2,'color',[1 0.5 0])

legend('exact','1st-order R-R','2nd-order R-R','rigid segments')
