clear all
close all
syms s  A  B  C D real
%% TRAVE INCASTRATA CON TRE CARICHI CONCENTRATI
P=133.61; %KN
delta=0.97
Mo=251; %KNm
L=8.5; %m
a=L/4; %m
E=(210000/1000)*1000^2 %modulo elastico acciaio KN/m2
J=179768785.01/1000^3 %inerzia fessurata m3
M1(s)=3*P/2*s-Mo;
M2(s)=3*P/2*(a+s)-P*s-Mo;
M3(s)=3*P/2*(2*a+s)-P*(a+s)-P*s-Mo;
M4(s)=3*P/2*(3*a+s)-P*(2*a+s)-P*(a+s)-P*s-Mo;
% plot per i valori
figure(1)
a1=0:0.1:a;
a2=max(a1):0.1:2*a;
a3=max(a2):0.1:3*a;
a4=max(a3):(L-max(a3))/21:4*a;
area(a1,M1(a1),'FaceColor', 'b','LineWidth',1.5)
hold on
area(a2,M2(a1),'FaceColor', 'b','LineWidth',1.5)
hold on
area(a3,M3(a1),'FaceColor', 'b','LineWidth',1.5)
hold on
area(a4,M4(a1),'FaceColor', 'b','LineWidth',1.5)
hold on
alpha(.4)
title({'diagramma del momento ridistribuito'})
xlabel('s  [m]')
ylabel('M(z)   [KN]')
axis ij
hold on
%% PRIMO TRATTO
v1(s)=M1/(E*J) %curvatura
v2(s)=int(v1,s)+A %rotazione 
v3(s)=int(v2,s)+B %spostamento
%condizioni al contorno
cc1=v3(0)==0;
cc2=v2(0)==0;
%% SECONDO TRATTO
vv1(s)=M2/(E*J); %curvatura
vv2(s)=int(vv1,s)+C; %rotazione 
vv3(s)=int(vv2,s)+D; %spostamento
%condizioni al contorno
cc3=vv2(a)==0;
cc4=v3(a)==vv3(0);
[MatA, b] = equationsToMatrix([cc1 cc2 cc3 cc4 ],[A B C D]);
Costanti = linsolve(MatA,b);
figure(2)
f1(s) = subs(v3, [A B ],[ Costanti(1) Costanti(2)]);
f2(s)=  subs(vv3, [C D],[ Costanti(3) Costanti(4)]);
plot(a1,f1(a1)*1000)
hold on
plot(a2,f2(a1)*1000)
title({'freccia'})
xlabel('s  [m]')
ylabel('f(z)   [mm]')
hold on
ftot=double(f2(a))*1000

