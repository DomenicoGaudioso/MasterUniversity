
clear all
close all

syms s A B C D E F L
cc= sym('c', [1 4])
s_plot=0:0.2:12
phi1(s)=1-cos(2*(pi*(s)/(12)))
phi2(s)=((sin((pi*(s)/(12)))).^2)*3
phi3(s)=(1-cos(4*(pi*(s)/(12))).^2)*3
G3=subs(phi3,s,s_plot)
G=subs(phi1,s,s_plot)
G2=subs(phi2,s,s_plot)
plot(s_plot,G3)
hold on
plot(s_plot,G)
hold on
plot(s_plot,G2)
hold on
phi(s)=phi
dphi(s)=diff(phi)


cc(1)=dphi(0)==0
cc(2)=dphi(L)==0
cc(3)=phi(0)==0
cc(4)=phi(L)==0

%[MatA, b] =equationsToMatrix([cc],[A B C D ])

%Costanti=linsolve(MatA,b)

%phi=subs(phi,C,Costanti(1))
%phi=subs(phi,D,Costanti(2))
%phi=subs(phi,E,Costanti(3))
%phi=subs(phi,F,Costanti(4))

%pretty(phi)