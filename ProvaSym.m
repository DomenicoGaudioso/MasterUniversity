%% Uniformly loaded simply supported beam - Analytical solution
% Computational Mechanics - Dr. Luca Taglialegne - Lesson 13.03.2017

clear all
close all

syms s EI q A B C D L

v4(s)=q/EI;
v3(s)=int(v4) + A;
v2(s)=int(v3) + B;
v1(s)=int(v2) + C;
v(s)=int(v1) + D

cc1 = v(0)==0
cc2 = v(L)==0
cc3 = v2(0)==0
cc4 = v2(L)==0

[MatA, b] = equationsToMatrix([cc1 cc2 cc3 cc4],[A B C D])

Costanti = linsolve(MatA,b)

v(s) = subs(v, A, Costanti(1));
v(s) = subs(v, B, Costanti(2));
v(s) = subs(v, C, Costanti(3));
v(s) = subs(v, D, Costanti(4))

L_plot = 10;
s_plot = linspace(0,L_plot,50)

v(s) = subs(v,q,3000)
v(s) = subs(v,L,L_plot)
v(s) = subs(v,EI,1*1)

v_plot = double(v(s_plot))

plot(s_plot,v_plot)