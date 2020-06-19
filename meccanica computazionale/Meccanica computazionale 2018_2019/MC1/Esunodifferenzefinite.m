clear all
close all

% Dati numerici
                     
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
an=(kn*bn/(4*En*Jn))^(1/4);
qstar=4*an^4*qn/(kn*bn);    %KN/m

n=25;                       

n1=n*4;                     %Numero punti tratto 1
n2=n*6;                     %Numero punti tratto 2
Ds1=L1n/(n1-1);             %Lunghezza intervallo tratto 1
Ds2=L2n/(n2-1);             %Lunghezza intervallo tratto 2

K=zeros(n1+n2+8);           %Matrice coefficienti nulla
B=zeros(n1+n2+8,1);         %Vettore terini noti nullo

for i=1:n1
    K(i,i)=1;
    K(i,i+1)=-4;
    K(i,i+2)=6+4*an^4*Ds1^4;
    K(i,i+3)=-4;
    K(i,i+4)=1;
    B(i,1)=qstar*Ds1^4;
end

for j=i+1:i+n2
    K(j,j+4)=1;
    K(j,j+5)=-4;
    K(j,j+6)=6+4*an^4*Ds2^4;
    K(j,j+7)=-4;
    K(j,j+8)=1;
    B(j,1)=qstar*Ds2^4;
end

%Condizioni al bordo

%Condizione 1
K(j+1,1)=-1;
K(j+1,2)=2;
K(j+1,4)=-2;
K(j+1,5)=1;
B(j+1,1)=2*Nan*Ds1^3/(En*Jn);
%Condizione 2
K(j+2,2)=1;
K(j+2,3)=-2;
K(j+2,4)=1;
B(j+2,1)=Man*Ds1^2/(En*Jn);
%Condizione 3
K(j+3,n1+n2+8)=1;
K(j+3,n1+n2+7)=-2;
K(j+3,n1+n2+5)=2;
K(j+3,n1+n2+4)=-1;
B(j+3,1)=-2*Ncn*Ds2^3/(En*Jn);
%Condizione 4
K(j+4,n1+n2+7)=1;
K(j+4,n1+n2+6)=-2;
K(j+4,n1+n2+5)=1;
B(j+4,1)=Mcn*Ds2^2/(En*Jn);
%Condizione 5
K(j+5,n1+9)=1/Ds2^3;
K(j+5,n1+8)=-2/Ds2^3;
K(j+5,n1+6)=2/Ds2^3;
K(j+5,n1+5)=-1/Ds2^3;
K(j+5,n1+4)=-1/Ds1^3;
K(j+5,n1+3)=2/Ds1^3;
K(j+5,n1+1)=-2/Ds1^3;
K(j+5,n1)=1/Ds1^3;
B(j+5,1)=2*Nbn/(En*Jn);
%Condizione 6
K(j+6,n1+8)=1/Ds1^2;
K(j+6,n1+7)=-2/Ds1^2;
K(j+6,n1+6)=1/Ds1^2;
K(j+6,n1+3)=-1/Ds1^2;
K(j+6,n1+2)=2/Ds1^2;
K(j+6,n1+1)=-1/Ds1^2;
B(j+6,1)=Mbn/(En*Jn);
%Condizione 7
K(j+7,n1+2)=1;
K(j+7,n1+7)=-1;
B(j+7,1)=0;
%Condizione 8
K(j+8,n1+3)=1;
K(j+8,n1+6)=1;
K(j+8,n1+1)=-1;
K(j+8,n1+8)=-1;
B(j+8,1)=0;

W=inv(K)*B;


%% Grafico deformata

Wplot1=W(3:n1+2);
sn1=0:L1n/(n1-1):L1n;
Wplot2=W(n1+7:n1+n2+6);
sn2=0:L2n/(n2-1):L2n;
Wplottot=[Wplot1' Wplot2'];
sntot=[sn1 sn2+L1n];
ind=[0 0];

hold on
subplot(3,1,1)
plot(sntot,Wplottot,[0 L1n+L2n],ind)
title('Grafico Deformata')
xlabel('Ascissa (m)')
ylabel('Spostamento linea asse (m)')
lgd=legend('Deformata','Conf. iniziale','location','eastoutside');
axis ij

%% Grafico momento

W1=W(1:n1+4);
for k=1:n1
    W12(k)=(W1(k+3)-2*W1(k+2)+W1(k+1))/Ds1^2;   %Derivata seconda spostamento tratto 1
end
M1=-En*Jn*W12;                                  %Momento tratto 1

W2=W(n1+5:n1+n2+8);
for m=1:n2
    W22(m)=(W2(m+3)-2*W2(m+2)+W2(m+1))/Ds2^2;   %Derivata seconda spostamento tratto 2
end
M2=-En*Jn*W22;                                  %Momento tratto 2

Mtot=[M1 M2];
stot=[sn1 sn2+L1n];

hold on
subplot(3,1,2)
plot(stot,Mtot,[0 L1n+L2n],ind)
title('Grafico Momento')
xlabel('Ascissa (m)')
ylabel('Momento (KNm)')
lgd=legend('Momento','Linea asse','location','eastoutside');
axis ij

%% Grafico taglio

for t=1:n1
    W13(t)=(W1(t+4)-2*W1(t+3)+2*W1(t+1)-W1(t))/(2*Ds1^3);   %Derivata terza spostamento tratto 1
end
T1=-En*Jn*W13;                                  %Taglio tratto 1

for r=1:n2
    W23(r)=(W2(r+4)-2*W2(r+3)+2*W2(r+1)-W2(r))/(2*Ds2^3);   %Derivata terza spostamento tratto 2
end
T2=-En*Jn*W23;                                  %Taglio tratto 2

Mtot=[T1 T2];
stot=[sn1 sn2+L1n];

hold on
subplot(3,1,3)
plot(stot,Mtot,[0 L1n+L2n],ind)
title('Grafico Taglio')
xlabel('Ascissa (m)')
ylabel('Taglio (KN)')
lgd=legend('Taglio','Linea asse','location','eastoutside');
axis ij