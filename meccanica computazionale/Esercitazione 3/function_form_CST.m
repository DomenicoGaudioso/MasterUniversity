%% Stiffness matrix for constant stress triangle (CST) element
% Computational Mechanics - Prof. Paolo S. Valvo - Lesson 27.03.2017
clear all
close all
%% Symbol definition
syms aI bI cI aJ bJ cJ aK bK cK X Y xI yI xJ yJ xK yK E1 nu t real
%% cordinate elemento 
x(1)=2; y(1)=2;
x(2)=6; y(2)=4;
x(3)=4; y(3)=8;
%% Connessioni
conn(1,1)=1;
conn(1,2)=2;
conn(1,3)=3;
%nodi elementi
I=conn(1,1);
J=conn(1,2);
K=conn(1,3);
%% Shape functions
A=1/2*(det([1 xI yI; 1  xJ  yJ; 1  xK  yK]));

NI(X,Y) = (aI + bI*X +cI*Y)/(2*A);
NJ(X,Y) = (aJ + bJ*X +cJ*Y)/(2*A);
NK(X,Y) = (aK + bK*X +cK*Y)/(2*A);

equI = [ NI(xI,yI) == 1 NI(xJ,yJ) == 0 NI(xK,yK) == 0 ];
equJ = [ NJ(xI,yI) == 0 NJ(xJ,yJ) == 1 NJ(xK,yK) == 0 ];
equK = [ NK(xI,yI) == 0 NK(xJ,yJ) == 0 NK(xK,yK) == 1 ];

[aI bI cI] = solve(equI, [aI bI cI]);
[aJ bJ cJ] = solve(equJ, [aJ bJ cJ]);
[aK bK cK] = solve(equK, [aK bK cK]);

%% Strain-displacement matrices
BI = [ bI 0; 0 cI; cI bI ];
BJ = [ bJ 0; 0 cJ; cJ bJ ];
BK = [ bK 0; 0 cK; cK bK ];
B = [ BI BJ BK ];

%% Shape functions no syms
aI=subs(aI,[xJ xK yJ yK],[x(J) x(K) y(J) y(K)]);
bI=subs(bI,[yJ yK],[y(J) y(K)]);
cI=subs(cI,[xJ xK],[x(J) x(K)]);

aJ=subs(aJ,[xI xK yI yK],[x(I) x(K) y(I) y(K)]);
bJ=subs(bJ,[yI yK],[y(I) y(K)]);
cJ=subs(cJ,[xI xK],[x(I) x(K)]);

aK=subs(aK,[xI xJ yI yJ],[x(I) x(J) y(I) y(J)]);
bK=subs(bK,[yI yJ],[y(I) y(J)]);
cK=subs(cK,[xI xJ],[x(I) x(J)]);
%% Strain-displacement matrices no syms
BI = [ bI 0; 0 cI; cI bI ];
BJ = [ bJ 0; 0 cJ; cJ bJ ];
BK = [ bK 0; 0 cK; cK bK ];
B = [ BI BJ BK ];

A=1/2*(det([1 x(I) y(I); 1  x(J)  y(J); 1  x(K)  y(K)]));

NI(X,Y) = (aI + bI*X +cI*Y)/(2*A);
NJ(X,Y) = (aJ + bJ*X +cJ*Y)/(2*A);
NK(X,Y) = (aK + bK*X +cK*Y)/(2*A);
%% plottaggi
figure(1)
subplot(2,2,2)
W=[x' y'];
T=[conn];
TR=triangulation(T,W);
triplot(TR,'--','LineWidth',1)
axis equal
hold on
view(3)

I=[x(1) x(2) x(3)];
J=[y(1) y(2) y(3)];
K=[1 0 0];
fill3(I,J,K,'r','Facealpha',0.4)
title('NI(x,y)')
axis equal
str1 = 'I';
str2 = 'J';
str3='K';
text(I+0.5,J+0.5,{str1 str2 str3},'HorizontalAlignment','center','Color','k','FontSize',10);
hold on

subplot(2,2,3)
triplot(TR,'--','LineWidth',1)
hold on
view(3)

I=[x(1) x(2) x(3)];
J=[y(1) y(2) y(3)];
K=[0 1 0];
fill3(I,J,K,'y','Facealpha',0.4)
title('NJ(x,y)')
axis equal
text(I+0.5,J+0.5,{str1 str2 str3},'HorizontalAlignment','center','Color','k','FontSize',10);
hold on

subplot(2,2,4)
triplot(TR,'--','LineWidth',1)
hold on
view(3)

I=[x(1) x(2) x(3)];
J=[y(1) y(2) y(3)];
K=[0 0 1];
fill3(I,J,K,'b','Facealpha',0.4)
title('NK(x,y)')
axis equal
text(I+0.5,J+0.5,{str1 str2 str3},'HorizontalAlignment','center','Color','k','FontSize',10);
hold on

subplot(2,2,1)
str1 = 'NI(x,y)=(aI + bI*x +cI*y)/(2*A)';
str2 = 'NJ(x,y)=(aJ + bJ*x +cJ*y)/(2*A)';
str3='NK(x,y)=(aK + bK*x +cK*y)/(2*A)';
text([0.5 0.5 0.5],[0.2 0.4 0.6],{str1 str2 str3},'HorizontalAlignment','center','Color','k','FontSize',10);
title('Shape functions, (CST) element')
axis equal
axis ij
axis off