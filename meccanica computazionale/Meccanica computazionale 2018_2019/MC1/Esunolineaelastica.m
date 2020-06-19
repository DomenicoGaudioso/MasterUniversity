clear all
close all

syms A B C D E F G H q a s k b L1 L2 Ma Mb Mc Na Nb Nc Ey J real                                                          

%Calcolo derivate funzione spostamento                                     

v1=q/(k*b)+exp(a*s)*(A*sin(a*s)+B*cos(a*s))+exp(-a*s)*(C*sin(a*s)+D*cos(a*s));        %Spostamento tratto s1

v11=diff(v1,s);               %Derivata prima di v1
v12=diff(v1,s,2);             %Derivata seconda di v1
v13=diff(v1,s,3);             %Derivata terza di v1

v2=q/(k*b)+exp(a*s)*(E*sin(a*s)+F*cos(a*s))+exp(-a*s)*(G*sin(a*s)+H*cos(a*s));        %Spostamento tratto s2

v21=diff(v2,s);               %Derivata prima di v2
v22=diff(v2,s,2);             %Derivata seconda di v2
v23=diff(v2,s,3);             %Derivata terza di v2

% Calcolo condizioni al bordo

Ey*J*subs(v13,s,0)==Na;
Ey*J*subs(v12,s,0)==Ma;
Ey*J*subs(v23,s,L2)==-Nc;
Ey*J*subs(v22,s,L2)==Mc;
Ey*J*(subs(v23,s,0)-subs(v13,s,L1))==Nb;
Ey*J*(subs(v22,s,0)-subs(v12,s,L1))==Mb;
subs(v1,s,L1)-subs(v2,s,0)==0;
subs(v11,s,L1)-subs(v21,s,0)==0;

% Calcolo coefficienti

eqns=[Ey*J*subs(v13,s,0)==Na,
      Ey*J*subs(v12,s,0)==Ma,
      Ey*J*subs(v23,s,L2)==-Nc,
      Ey*J*subs(v22,s,L2)==Mc,
      Ey*J*(subs(v23,s,0)-subs(v13,s,L1))==Nb,
      Ey*J*(subs(v22,s,0)-subs(v12,s,L1))==Mb,
      subs(v1,s,L1)-subs(v2,s,0)==0,
      subs(v11,s,L1)-subs(v21,s,0)==0];
  
  vars=[A B C D E F G H];
  
  [M,V]=equationsToMatrix(eqns,vars);
  
  K=inv(M)*V;

%Valori numerici

qn=19.5;                    %KN/m
bn=1.3;                     %m
kn=25000;                   %KN/m^3
En=31500000;                %KN/m^3
Jn=0.0634;                  %m^4
Nan=240;                    %KN
Nbn=845;                    %KN
Ncn=415;                    %KN
Man=65;                     %KNm
Mbn=92;                     %KNm
Mcn=120;                    %KNm  
L1n=4;                      %m
L2n=6;                      %m
an=(kn*bn/(4*En*Jn))^(1/4); %Coefficiente alpha
n1=20;                      %Discretizzazione per grafico
n2=20;                      %Discretizzazione per grafico
sn1=0:L1n/n1:L1n;           %Ascissa tratto s1 discretizzata
sn2=0:L2n/n2:L2n;           %Ascissa tratto s2 discretizzata
sntot=[sn1 (sn2+4)];        %Ascissa totale 

v1risolta=q/(k*b)+exp(a*s)*(K(1)*sin(a*s)+K(2)*cos(a*s))+exp(-a*s)*(K(3)*sin(a*s)+K(4)*cos(a*s)); %Spostamento risolto
 
v2risolta=q/(k*b)+exp(a*s)*(K(5)*sin(a*s)+K(6)*cos(a*s))+exp(-a*s)*(K(7)*sin(a*s)+K(8)*cos(a*s)); %Spostamento risolto

%% Grafico deformata

v1n=subs(v1risolta,[Ey,J,Na,Nb,Nc,Ma,Mb,Mc,a,L1,L2,b,k,q],[En,Jn,Nan,Nbn,Ncn,Man,Mbn,Mcn,an,L1n,L2n,bn,kn,qn]); %Spostamento risolto numerico funzione di s1

v1sn=double(subs(v1n,s,sn1));

v2n=subs(v2risolta,[Ey,J,Na,Nb,Nc,Ma,Mb,Mc,a,L1,L2,b,k,q],[En,Jn,Nan,Nbn,Ncn,Man,Mbn,Mcn,an,L1n,L2n,bn,kn,qn]); %Spostamento risolto numerico funzione di s2

v2sn=double(subs(v2n,s,sn2));

vntot=[v1sn v2sn];

hold on
subplot(3,1,1)
plot(sntot,vntot,[0 L1n+L2n],[0 0])
title('Grafico spostamenti')
xlabel('Ascissa(m)')
ylabel('Spostamento linea asse(m)')
lgd=legend('Deformata','Conf.iniziale','location','eastoutside');
axis ij
  
format short e
vA=double(subs(v1n,s,0))               %Spostamento in A
vBs1=double(subs(v1n,s,L1n))           %Spostamento in B (s1)
vBs2=double(subs(v2n,s,0))             %Spostamnto in B (s2)
vC=double(subs(v2n,s,L2n))             %Spostamento in C

%% Grafico Taglio

T1n=subs(-Ey*J*diff(v1risolta,s,3),[Ey,J,Na,Nb,Nc,Ma,Mb,Mc,a,L1,L2,b,k,q],[En,Jn,Nan,Nbn,Ncn,Man,Mbn,Mcn,an,L1n,L2n,bn,kn,qn]); %Spostamento risolto numerico funzione di s1

T1sn=double(subs(T1n,s,sn1));

T2n=subs(-Ey*J*diff(v2risolta,s,3),[Ey,J,Na,Nb,Nc,Ma,Mb,Mc,a,L1,L2,b,k,q],[En,Jn,Nan,Nbn,Ncn,Man,Mbn,Mcn,an,L1n,L2n,bn,kn,qn]); %Spostamento risolto numerico funzione di s2

T2sn=double(subs(T2n,s,sn2));

Tntot=[T1sn T2sn];

hold on
subplot(3,1,2)
plot(sntot,Tntot,[0 L1n+L2n],[0 0])
title('Grafico Taglio')
xlabel('Ascissa(m)')
ylabel('Taglio (KN)')
lgd=legend('Taglio','Linea asse','location','eastoutside');
axis ij

%% Grafico Momento

M1n=subs(-Ey*J*diff(v1risolta,s,2),[Ey,J,Na,Nb,Nc,Ma,Mb,Mc,a,L1,L2,b,k,q],[En,Jn,Nan,Nbn,Ncn,Man,Mbn,Mcn,an,L1n,L2n,bn,kn,qn]); %Spostamento risolto numerico funzione di s1

M1sn=double(subs(M1n,s,sn1));

M2n=subs(-Ey*J*diff(v2risolta,s,2),[Ey,J,Na,Nb,Nc,Ma,Mb,Mc,a,L1,L2,b,k,q],[En,Jn,Nan,Nbn,Ncn,Man,Mbn,Mcn,an,L1n,L2n,bn,kn,qn]); %Spostamento risolto numerico funzione di s2

M2sn=double(subs(M2n,s,sn2));

Mntot=[M1sn M2sn];

hold on
subplot(3,1,3)
plot(sntot,Mntot,[0 L1n+L2n],[0 0])
title('Grafico Momento')
xlabel('Ascissa (m)')
ylabel('Momento (KNm)')
lgd=legend('Momento','Linea asse','location','eastoutside');
axis ij
