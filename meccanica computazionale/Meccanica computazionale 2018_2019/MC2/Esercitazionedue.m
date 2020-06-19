%% ESERCITAZIONE 2 MC 2018_2019

%Readme

%Nelle righe 33,34,35,160 ci sono i valori per confronto con SAP

%Run-->viene visualizzata la figura con i nodi assegnati e vengono...
%...visualizzati i grafici di Taglio e Momento delle "aste"

%Visualizzare Tabella spostamenti nodali-->  digitare TU--> (m)
%Visualizzare Tabella sforzo normale aste--> digitare TN--> (KN)
%Visualizzare Tabella reazioni vincolari-->  digitare TR--> (KN)

clear all
close all

%% Dati del problema

n=25;                % Numero nodi
m=72;                % Numero aste
rho=78.5;           % (KN/m^3) Peso specifico acciaio
E=210e6;             % (KN/m^2) Modulo di Young
a=3;                 % (m) Lunghezza aste orizzontali
%b=a/2*(2)^(1/2);    % (m) Altezza nodo aste
do=114.3e-3;                                %(m)   Diametro aste orizzontali
to=10e-3;                                   %(m)   Diametro aste orizzontali
%Ao=(pi*do^2/4-pi*(do-2*to)^2/4);           %(m^2) Area aste orizzontali
di=139.7e-3;                                %(m)   Diametro aste inclinate
ti=10e-3;                                   %(m)   Diametro aste inclinate
%Ai=(pi*di^2/4-pi*(di-2*ti)^2/4);           %(m^2) Area asteinclinate
t=25;                                       %(°C)  Variazione termica
alpha = 12e-6;                              %(1/C) Coefficiente di dilatazione termica
Ai=4.075e-3;                                %Area aste inclinate SAP
Ao=3.277e-3;                                %Area aste orizzontali SAP
b=2.12;                                     %Altezza nodo aste SAP             

P1=-10;    %(KN) Carico concentrato nodo 17 
P2=-20;    %(KN) Carico concentrato nodi 18,20
P3=-30;    %(KN) Carico concentrato nodi 19,21,23
P4=-40;    %(KN) Carico concentrato nodi 22,24
P5=-50;    %(KN) Carico concentrato nodo 25

for e=1:24
    A(e)=Ao;
end
for e=61:72
    A(e)=Ao;
end
for e=25:60
    A(e)=Ai;
end
    
%% Coordinate dei nodi

x(1)=0;       y(1)=0;      z(1)=0;
x(2)=a;       y(2)=0;      z(2)=0;
x(3)=2*a;     y(3)=0;      z(3)=0;
x(4)=3*a;     y(4)=0;      z(4)=0;
x(5)=0;       y(5)=a;      z(5)=0;
x(6)=a;       y(6)=a;      z(6)=0;
x(7)=2*a;     y(7)=a;      z(7)=0;
x(8)=3*a;     y(8)=a;      z(8)=0;
x(9)=0;       y(9)=2*a;    z(9)=0;
x(10)=a;      y(10)=2*a;   z(10)=0;
x(11)=2*a;    y(11)=2*a;   z(11)=0;
x(12)=3*a;    y(12)=2*a;   z(12)=0;
x(13)=0;      y(13)=3*a;   z(13)=0;
x(14)=a;      y(14)=3*a;   z(14)=0;
x(15)=2*a;    y(15)=3*a;   z(15)=0;
x(16)=3*a;    y(16)=3*a;   z(16)=0;
x(17)=a/2;    y(17)=a/2;   z(17)=b;
x(18)=3/2*a;  y(18)=a/2;   z(18)=b;
x(19)=5/2*a;  y(19)=a/2;   z(19)=b;
x(20)=a/2;    y(20)=3/2*a; z(20)=b;
x(21)=3/2*a;  y(21)=3/2*a; z(21)=b;
x(22)=5/2*a;  y(22)=3/2*a; z(22)=b;
x(23)=a/2;    y(23)=5/2*a; z(23)=b;
x(24)=3/2*a;  y(24)=5/2*a; z(24)=b;
x(25)=5/2*a;  y(25)=5/2*a; z(25)=b;

%% Matrice di connettività

conn(1,1)=1;     conn(1,2)=2;
conn(2,1)=2;     conn(2,2)=3;
conn(3,1)=3;     conn(3,2)=4;
conn(4,1)=5;     conn(4,2)=6;
conn(5,1)=6;     conn(5,2)=7;
conn(6,1)=7;     conn(6,2)=8;
conn(7,1)=9;     conn(7,2)=10;
conn(8,1)=10;    conn(8,2)=11;
conn(9,1)=11;    conn(9,2)=12;
conn(10,1)=13;   conn(10,2)=14;
conn(11,1)=14;   conn(11,2)=15;
conn(12,1)=15;   conn(12,2)=16;
conn(13,1)=1;    conn(13,2)=5;
conn(14,1)=5;    conn(14,2)=9;
conn(15,1)=9;    conn(15,2)=13;
conn(16,1)=2;    conn(16,2)=6;
conn(17,1)=6;    conn(17,2)=10;
conn(18,1)=10;   conn(18,2)=14;
conn(19,1)=3;    conn(19,2)=7;
conn(20,1)=7;    conn(20,2)=11;
conn(21,1)=11;   conn(21,2)=15;
conn(22,1)=4;    conn(22,2)=8;
conn(23,1)=8;    conn(23,2)=12;
conn(24,1)=12;   conn(24,2)=16;
conn(25,1)=1;    conn(25,2)=17;
conn(26,1)=2;    conn(26,2)=17;
conn(27,1)=2;    conn(27,2)=18;
conn(28,1)=3;    conn(28,2)=18;
conn(29,1)=3;    conn(29,2)=19;
conn(30,1)=4;    conn(30,2)=19;
conn(31,1)=5;    conn(31,2)=17;
conn(32,1)=6;    conn(32,2)=17;
conn(33,1)=6;    conn(33,2)=18;
conn(34,1)=7;    conn(34,2)=18;
conn(35,1)=7;    conn(35,2)=19;
conn(36,1)=8;    conn(36,2)=19;
conn(37,1)=5;    conn(37,2)=20;
conn(38,1)=6;    conn(38,2)=20;
conn(39,1)=6;    conn(39,2)=21;
conn(40,1)=7;    conn(40,2)=21;
conn(41,1)=7;    conn(41,2)=22;
conn(42,1)=8;    conn(42,2)=22;
conn(43,1)=9;    conn(43,2)=20;
conn(44,1)=10;   conn(44,2)=20;
conn(45,1)=10;   conn(45,2)=21;
conn(46,1)=11;   conn(46,2)=21;
conn(47,1)=11;   conn(47,2)=22;
conn(48,1)=12;   conn(48,2)=22;
conn(49,1)=9;    conn(49,2)=23;
conn(50,1)=10;   conn(50,2)=23;
conn(51,1)=10;   conn(51,2)=24;
conn(52,1)=11;   conn(52,2)=24;
conn(53,1)=11;   conn(53,2)=25;
conn(54,1)=12;   conn(54,2)=25;
conn(55,1)=13;   conn(55,2)=23;
conn(56,1)=14;   conn(56,2)=23;
conn(57,1)=14;   conn(57,2)=24;
conn(58,1)=15;   conn(58,2)=24;
conn(59,1)=15;   conn(59,2)=25;
conn(60,1)=16;   conn(60,2)=25;
conn(61,1)=17;   conn(61,2)=18;
conn(62,1)=18;   conn(62,2)=19;
conn(63,1)=20;   conn(63,2)=21;
conn(64,1)=21;   conn(64,2)=22;
conn(65,1)=23;   conn(65,2)=24;
conn(66,1)=24;   conn(66,2)=25;
conn(67,1)=17;   conn(67,2)=20;
conn(68,1)=20;   conn(68,2)=23;
conn(69,1)=18;   conn(69,2)=21;
conn(70,1)=21;   conn(70,2)=24;
conn(71,1)=19;   conn(71,2)=22;
conn(72,1)=22;   conn(72,2)=25;

%% Matrice di rigidezza

for e = 1 : m
    I = conn (e,1);
    J = conn (e,2);
    L(e) = sqrt((x(J)-x(I))^2+(y(J)-y(I))^2+(z(J)-z(I))^2);     %(m) Lunghezza elemento
    %L(e)=2.99907;                                               %(m) Lunghezza elemento SAP
    exloc(:,e)=1/L(e)*[x(J)-x(I)   y(J)-y(I)  z(J)-z(I)];       % Versore ex locale
    eyloc=zeros(3,m);
    if exloc(1,e)^2+exloc(2,e)~=0
       eyloc(1,e)=-exloc(2,e)/sqrt(exloc(2,e)^2+exloc(1,e)^2);  %Versore ey locale
       eyloc(2,e)=exloc(1,e)/sqrt(exloc(2,e)^2+exloc(1,e)^2);
     else
       eyloc(1,e)=exloc(3,e); 
    end
    ezloc(:,e)=cross(exloc(1:3,e),eyloc(1:3,e));                      %Versore ez locale
    Q(:,:,e)=[exloc(1:3,e) eyloc(1:3,e) ezloc(1:3,e)];                %Matrice di trasferimento
    T(:,:,e) = [Q(:,:,e)  zeros(3,3); zeros(3,3)  Q(:,:,e)];          %Matrice di trasferimento dell'elemento
    Kloc(:,:,e) = E*A(e)/L(e)*[1 0 0 -1 0 0;0 0 0 0 0 0;0 0 0 0 0 0;  %Matrice di rigidezza elemento in locale
                               -1 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 0 0]; 
    Kglob(:,:,e)=T(:,:,e)*Kloc(:,:,e)*T(:,:,e)';                      %Matrice di rigidezza dell'elemento in globale   
end

K=zeros(3*n,3*n);
for e = 1:m                                               
    I=conn(e,1);
    J=conn(e,2);
    K(3*I-2:3*I,3*I-2:3*I) = K(3*I-2:3*I,3*I-2:3*I) + Kglob(1:3,1:3,e);
    K(3*I-2:3*I,3*J-2:3*J) = K(3*I-2:3*I,3*J-2:3*J) + Kglob(1:3,4:6,e);
    K(3*J-2:3*J,3*I-2:3*I) = K(3*J-2:3*J,3*I-2:3*I) + Kglob(4:6,1:3,e);
    K(3*J-2:3*J,3*J-2:3*J) = K(3*J-2:3*J,3*J-2:3*J) + Kglob(4:6,4:6,e);
end
 
%% Vettore dei carichi

for e=1:m
Pp(e)=A(e)*L(e)*rho;      %(KN) Peso proprio elemento
end

ppn=zeros(3*n,1);         %(KN) Peso proprio nodale 
for e =1:m
    I =conn(e,1);
    J=conn(e,2);
    ppn(3*I)=ppn(3*I)+Pp(e)/2;
    ppn(3*J)=ppn(3*J)+Pp(e)/2;
end

pc=zeros(3*n,1);                 %Vettore carichi concentrati
pc(3*17,1)=P1;                   %Carico concentrato nodo 17
pc(3*18,1)=P2;                   %Carico concentrato nodo 18
pc(3*20,1)=P2;                   %Carico concentrato nodo 20
pc(3*19,1)=P3;                   %Carico concentrato nodo 19
pc(3*21,1)=P3;                   %Carico concentrato nodo 21
pc(3*23,1)=P3;                   %Carico concentrato nodo 23
pc(3*22,1)=P4;                   %Carico concentrato nodo 22
pc(3*24,1)=P4;                   %Carico concentrato nodo 24
pc(3*25,1)=P5;                   %Carico concentrato nodo 25

te=zeros(3*n,1);                  %Vettore carichi termici
te(3*17-2)=-E*Ao*alpha*t;         %Carico termico nodo 17
te(3*19-2)=E*Ao*alpha*t;          %Carico termico nodo 19
te(3*20-2)=-E*Ao*alpha*t;         %Carico termico nodo 20
te(3*22-2)=E*Ao*alpha*t;          %Carico termico nodo 22
te(3*23-2)=-E*Ao*alpha*t;         %Carico termico nodo 23
te(3*25-2)=E*Ao*alpha*t;          %Carico termico nodo 25

p=pc+te-ppn;                      %Vettore totale carichi
%% Vincoli

Kv=K;
Kv(1:3,1:3*n)=zeros;                     %Vincolo nodo 1
Kv((4*3)-2:4*3,1:3*n)=zeros;             %Vincolo nodo 4
Kv((13*3)-2:13*3,1:3*n)=zeros;           %Vincolo nodo 13
Kv((16*3)-2:16*3,1:3*n)=zeros;           %Vincolo nodo 16

pv=p;
pv(1:3,1)=0;                              %Vincolo nodo 1
pv((4*3)-2:4*3,1)=0;                      %Vincolo nodo 4
pv((13*3)-2:13*3,1)=0;                    %Vincolo nodo 13
pv((16*3)-2:16*3,1)=0;                    %Vincolo nodo 16

for i=1:3
    Kv(i,i)=1;
end
for i=(4*3)-2:4*3
    Kv(i,i)=1;
end
for i=(13*3)-2:13*3
    Kv(i,i)=1;
end
for i=(16*3)-2:16*3
    Kv(i,i)=1;
end
%% Soluzione

u=inv(Kv)*pv;           %(m)  Spostamento nodi
r=K*u-p;                %(KN) Reazioni vincolari

for e=1:n
    unodi(:,:,e)=u(e*3-2:e*3,1);                %Vettore degli spostamenti in globale dei singoli nodi
end
for e = 1 : m
    I = conn (e,1);
    J = conn (e,2);
    uglob(:,:,e)=[u(I*3-2:I*3);u(J*3-2:J*3)];   %Spostamento estremi elemento in globale
end
for e=1:m
    uloc(:,:,e)=T(:,:,e)'*uglob(:,:,e);         %Spostamenti estremi elemento in locale
end
for e=1:m
N(:,:,e)=E*A(e)/L(e)*(uloc(4,e)-uloc(1,e));
end

N(:,:,61)=N(:,:,61)-E*Ao*alpha*t;
N(:,:,62)=N(:,:,62)-E*Ao*alpha*t;
N(:,:,63)=N(:,:,63)-E*Ao*alpha*t;
N(:,:,64)=N(:,:,64)-E*Ao*alpha*t;
N(:,:,65)=N(:,:,65)-E*Ao*alpha*t;
N(:,:,66)=N(:,:,66)-E*Ao*alpha*t;

%% Tabella Spostamenti

ux=zeros(n,1);
uy=zeros(n,1);
uz=zeros(n,1);
for i=1:n
ux(i) = [unodi(1,1,i)];
uy(i) = [unodi(2,1,i)];
uz(i) = [unodi(3,1,i)];
end
Nodo= [1:n]';

TU=table(Nodo,ux,uy,uz);       %Tabella Spostamenti

%% Tabella reazioni vincolari

Nodo=[1 4 13 16]';
Rvx=[r(1*3-2,1) r(4*3-2,1) r(13*3-2,1) r(16*3-2,1)]';
Rvy=[r(1*3-1,1) r(4*3-1,1) r(13*3-1,1) r(16*3-1,1)]';
Rvz=[r(1*3,1) r(4*3,1) r(13*3,1) r(16*3,1)]';
TR=table(Nodo,Rvx,Rvy,Rvz);   %Tabella Reazioni vincolari

%% Tabella sforzi normali

for e=1:m
Win(:,:,e)=[0 0 Pp(e)/2]';                   % Forza nodale in globale
Winloc(:,:,e)=Q(:,:,e)'*Win(:,:,e);          % Forza nodale in locale
end

for e=1:m
    Nt(:,:,e)=[N(e)-Winloc(1,1,e); N(e); N(e)+Winloc(1,1,e)];
end
Ni=zeros(m,1);
Nmedio=zeros(m,1);
Nj=zeros(m,1);
for i=1:m
    Ni(i)=Nt(1,:,i);
    Nmedio(i)=Nt(2,:,i);
    Nj(i)=Nt(3,:,i);
end

Nodoi=zeros(m,1);
Nodoj=zeros(m,1);
for i=1:m
Nodoi(i)= conn(i,1);
Nodoj(i)=conn(i,2);
end

TN=table(Nodoi,Nodoj,Ni,Nmedio,Nj);         %Tabella Sforzi normali

%% Grafico Taglio e Momento

for e=1:m
qv(:,:,e)=[0 0 Pp(e)/L(e)]';        % (KN/m) Carico distribuito elemento sistema globale (z)
qvloc(:,:,e)=Q(:,:,e)'*qv(:,:,e);   % (KN/m) Carico distribuito elemento sistema locale
end

syms sloc

figure;
subplot(2,2,1)
slocn=0:L(1)/100:L(1);
Mo=qvloc(3,1,1)*L(1)/2*sloc-qvloc(3,1,1)*sloc^2/2;   %Momento aste orizzontali
Mon=subs(Mo,sloc,slocn);
plot(slocn,Mon,[0 L(1)],[0 0])
title('Grafico Momento Aste Orizzontali')
xlabel('Ascissa s(m)')
ylabel('Momento (KNm)')
axis ij
axis([0 3 0 0.3])
subplot(2,2,2)
slocni=0:L(30)/100:L(30);
Mi=qvloc(3,1,30)*L(1)/2*sloc-qvloc(3,1,30)*sloc^2/2;  %Momento aste inclinate
Min=subs(Mi,sloc,slocni);
plot(slocni,Min);
plot(slocni,Min,[0 L(1)],[0 0])
title('Grafico Momento Aste Inclinate')
xlabel('Ascissa s(m)')
ylabel('Momento (KNm)')
axis ij
axis([0 3 0 0.3])

subplot(2,2,3)
slocn=0:L(1)/100:L(1);
Mo=qvloc(3,1,1)*L(1)/2*sloc-qvloc(3,1,1)*sloc^2/2;   
To=diff(Mo,sloc,1);                                  %Taglio aste orizzontali
Ton=subs(To,sloc,slocn);
plot(slocn,Ton,[0 L(1)],[0 0])
title('Grafico Taglio Aste Orizzontali')
xlabel('Ascissa s(m)')
ylabel('Taglio (KN)')
axis ij
axis([0 3 -0.4 0.4])

subplot(2,2,4)
slocni=0:L(30)/100:L(30);
Mi=qvloc(3,1,30)*L(30)/2*sloc-qvloc(3,1,30)*sloc^2/2;   
Ti=diff(Mi,sloc,1);                                     %Taglio aste inclinate
Tin=subs(Ti,sloc,slocni); 
plot(slocni,Tin,[0 L(30)],[0 0])
title('Grafico Taglio Aste Inclinate')
xlabel('Ascissa s(m)')
ylabel('Taglio (KN)')
axis ij
axis([0 3 -0.4 0.4])

%% Grafico Numerazione Nodi

figure;
hold on
axis equal
axis([-1 10 -1 10]) 
P=[x;y]';
scatter (P(:,1),P(:,2),'k');
for e=1:1:m
    I = conn(e,1);
    J = conn(e,2);
    line([P(I,1) P(J,1)],[P(I,2) P(J,2)],'Color','b');
end
for i=1:n
str(i)={i};
text(x(i)+0.3,y(i)+0.2,str(i),'color','k')
end
quiver(0,0,1,0,'MaxHeadSize',2,'color','r')
quiver(0,0,0,1,'MaxHeadSize',2,'color','r')
text(0.6,-0.2,'X','color','r')
text(-0.2,0.6,'Y','color','r')
text(-0.3,-0.5,'Z uscente','color','r')