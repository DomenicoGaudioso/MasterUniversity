clear all
close all
syms x real
%% Propietà della sezione
As=1520;
n=15;
c=56;
d=1144;
b=300;
B=1200;
H=1200;
s=300;
fck=25 %resistenza cilindrica cls C25/30
fctm=0.3*fck^(2/3) %resistenza media a trazione semplice per classi <c50/60
%% Stato limite delle tensioni
E1=0.5*b*x^2+n*As*(x-c)-n*As*(d-x)==0; %anima della sezione compressa
S=solve([E1],[x]);
x1=[double(S(1)); double(S(2))]
E2=0.5*B*x^2+n*As*(x-c)-n*As*(d-x)==0; %anima della piattabanda compressa
S2=solve([E2],[x]);
x2=[double(S2(1));double(S2(2))]
%% Stato limite fessurazione
% STADIO I
%Asse neutro nell'anima
n=6;
n1=0.5;
SxI=b*x^2*0.5-((H-s-x)^2*b*0.5+B*s*(H-0.5*s-x))*n1+n*As*(x-c)-n*As*(d-x)==0
S3=solve([SxI],[x]);
x3=[double(S3(1));double(S3(2))]
IxI=(b*x3(2)^3)/3+(b*(H-s-x3(2))^3*n1)/3+n1*(B*s^3)/12+B*s*n1*(H-s*0.5-x3(2))^2+n*As*(x3(2)-c)^2+n*As*(d-x3(2))^2
MfI=fctm*IxI /(1.2*n1*(H-x3(2))) % Momento di prima fessurazione
%Asse neutro nella piattabanda
Sx1I=B*x^2*0.5-((s-x)^2*B*0.5+b*(H-s)*((H-s)*0.5+s-x))*n1+n*As*(x-c)-n*As*(d-x)==0
S4=solve([Sx1I],[x]);
x4=[double(S4(1));double(S4(2))]
Ix1I=(B*x4(2)^3)/3+(B*(s-x4(2))^3*n1)/3+n1*(b*(H-s)^3)/12+b*(H-s)*n1*((H-s)*0.5+s-x4(2))^2+n*As*(x4(2)-c)^2+n*As*(d-x4(2))^2
Mf1I=fctm*Ix1I /(1.2*n1*(H-x3(2))) % Momento di prima fessurazione
