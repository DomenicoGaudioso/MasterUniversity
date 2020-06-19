clear all 
close all
   
syms phi phi1 phi2 phi3 s E J q L L1 k Na Nb Nc Ma Mb Mc Nc B Q real

n=6;                         %Numero funzioni(inserire numero pari)

phi(1)=1;
phi1(1)=0;
phi2(1)=0;
phi3(1)=0;

phi(2)=s-L1;
phi1(2)=1;
phi2(2)=0;
phi3(2)=0;

for i=3:2:n
    phi(i)=heaviside(L1-s)*(s-L1)^((i+1)/2);
    phi1(i)=heaviside(L1-s)*((i+1)/2)*(s-L1)^((i+1)/2-1);
    phi2(i)=heaviside(L1-s)*((i+1)/2)*((i+1)/2-1)*(s-L1)^((i+1)/2-2);
    phi3(i)=heaviside(L1-s)*((i+1)/2)*((i+1)/2-1)*((i+1)/2-2)*(s-L1)^((i+1)/2-3);
    
    phi(i+1)=heaviside(s-L1)*(s-L1)^((i+1)/2);
    phi1(i+1)=heaviside(s-L1)*((i+1)/2)*(s-L1)^((i+1)/2-1);
    phi2(i+1)=heaviside(s-L1)*((i+1)/2)*((i+1)/2-1)*(s-L1)^((i+1)/2-2);
    phi3(i+1)=heaviside(s-L1)*((i+1)/2)*((i+1)/2-1)*((i+1)/2-2)*(s-L1)^((i+1)/2-3);
end
%%
if (-1)^n>0
for i = 1:n                            
    for j = 1:n
        B(i, j) = int(E*J*phi2(i)*phi2(j), s, 0, L)+int(k*phi(i)*phi(j), s, 0, L);
        Q(i) = int(q*phi(i), s, 0, L)+Na*subs(phi(i),s,0)+Nb*subs(phi(i),s,L1)+Nc*subs(phi(i),s,L)+Ma*subs(phi1(i),s,0)+Mb*subs(phi1(i),s,L1)+Mc*subs(phi1(i),s,L); 
    end
end
else if (-1)^n<0
        for i = 1:n+1                      
    for j = 1:n+1
        B(i, j) = int(E*J*phi2(i)*phi2(j), s, 0, L)+int(k*phi(i)*phi(j), s, 0, L);
        Q(i) = int(q*phi(i), s, 0, L)+Na*subs(phi(i),s,0)+Nb*subs(phi(i),s,L1)+Nc*subs(phi(i),s,L)+Ma*subs(phi1(i),s,0)+Mb*subs(phi1(i),s,L1)+Mc*subs(phi1(i),s,L); 
    end
        end
    end
end
   
%% Valori numerici
 
qn=19.5;                     %KN/m
kn=32500;                    %Cosante molla (KN/m^3)
En=31500000 ;                %Modulo young (KN/m^2)
Jn=0.0634;                   %Momento inerzia (m^4)
Ln=10;                       %Lunghezza (m)
L1n=4;                       %Lunghezza (m)
Nan=240;                     %KN
Nbn=845;                     %KN
Ncn=415;                     %KN
Man=-65;                     %KNm
Mbn=-92;                     %KNm
Mcn=120;                     %KNm
sn=0:Ln/1000:Ln;             %Ascissa

Bn=subs(B,[E,J,L,L1,k,Na,Nb,Nc,Ma,Mb,Mc,q],[En,Jn,Ln,L1n,kn,Nan,Nbn,Ncn,Man,Mbn,Mcn,qn]);
Qn=subs(Q,[E,J,L,L1,k,Na,Nb,Nc,Ma,Mb,Mc,q],[En,Jn,Ln,L1n,kn,Nan,Nbn,Ncn,Man,Mbn,Mcn,qn]);

C=inv(Bn)*Qn';

%% Grafici

v=double(subs(subs(phi*C, L1, L1n), s, sn));
M=double(subs(subs(-En*Jn*phi2*C, L1, L1n), s, sn));
T=double(subs(subs(-En*Jn*phi3*C, L1, L1n), s, sn));

hold on
subplot(3,1,1)
plot(sn,v,[0 Ln],[0 0]);
title('Grafico spostamenti')
xlabel('Ascissa s(m)')
ylabel('Spostamento linea asse (m)')
axis ij

subplot(3,1,2)
plot(sn,T,[0 Ln],[0 0])
title('Grafico Taglio')
xlabel('Ascissa s(m)')
ylabel('Taglio (KN)')
axis ij

subplot(3,1,3)
plot(sn,M,[0 Ln],[0 0])
title('Grafico Momento')
xlabel('Ascissa s(m)')
ylabel('Momento (KNm)')
axis ij