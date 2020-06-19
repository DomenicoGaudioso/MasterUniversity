%%  Computational Mechanics - PRIMA ESERCITAZIONE
% Triangular loaded on beam stuck to the ends  
%% Clear workspace and close any open windows
clear all
close all
%% dati profilo IPE
% nel nostro caso IPE 400
h=0.4      % (m) altezza profilato
b=0.18     % (m) lunghezza piattabanda
a=0.0086   % (m) spessore anima
e=0.0133   % (m) spessore piattabanda

A=2*b*e+(h-2*e)*a                                         % (m^2) area IPE trascurando gli smussi
Ix=2*([(1/12)*b*e^3]+b*e*[(h-e)/2]^2)+(1/12)*a*(h-2*e)^3  % (m^4) moment of inertia IPE secondo l'asse x
Iy=2*[(1/12)*e*b^3]+[(1/12)*(h-2*e)*a^3]                  % (m^4) moment of inertia IPE secondo l'asse y
I=Ix                                                       % poichè in questo caso useremo soltanto Ix
%% dati problema
E = 210e9/1000;                  % (KN/m^2) Young's modulus steel
Matricola=506682
L = Matricola/(40*1000)          % (m) length       
syms s                           % (m) abscissa (for plots)
q=Matricola/10000                % (KN/m) max forze on triangle load     
%% METODO DI RAYLEGH-RITZ 
% con funzione di base trigonometrica
digits(5); %  variable precision used

for j=1:1;
phi(j)=1-cos(2*(j)*pi*(s)/(L.^(j)))  
dphi(j)= diff(phi(j));             % derivata prima Funzioni di base (polinomiale)
d2phi(j)= diff(phi(j),2);  % derivata seconda Funzioni di base (polinomiale)
end

C=sym(['C'], [max(j) 1])

vRR=sum(phi*C);              % spostamento approssimato (v tilde)
u=(sum(d2phi*C)).^2  

U=(1/2)*E*I*int(u,s,0,L)   %Energia di deformazioneelastica approssimata (U tilde)

w=q*(s/L)*(vRR)
W=int(w,s,0,L);             % Lavoro dei carichi approssimato (W tilde)

V=(U-W);                    % Energia potenziale totale (V tilde)


for j=1:max(j);  
derivate_parziali(j)= diff(V,C(j))==0; % derivate parziali rispetto a c1,c2 e cj di V  (V tilde)
end

[MatA, b] =equationsToMatrix(derivate_parziali,[C]);
Costanti_1=vpa(linsolve(MatA,b));

vRRt=vpa(subs(vRR,C,(Costanti_1)));
%% Caratteristiche della sollecitazione
M_RRt(s)= vpa(-E*I*diff(vRRt,2)) % Bending moment
T_RRt(s)=diff(M_RRt(s))             % Shear force
%% Grafici 
% deformata
s_plot = linspace(0,L,50);
vRRt_plot=subs(vRRt,s,s_plot);
subplot(2,2,[1 2])
plot(s_plot,vRRt_plot,'LineWidth',2,'DisplayName','RRt')

xlabel('s (m)')
ylabel('w (m)')
axis ij
hold on
%CdS
%Bending moment
subplot(2,2,4)
plot(s_plot,M_RRt(s_plot),'Color', [0 0.25 0.25],'LineWidth',1.5,'DisplayName','MRRt')
alpha(.1)
axis([0 max(s_plot) -500 500])
title('Bending Moment')
xlabel('s  [m]')
ylabel('M(s)    [KN*m]' )
axis ij
hold on
legend
%Shear force
subplot(2,2,3)
plot(s_plot,T_RRt(s_plot),'Color', 'r','LineWidth',1.5,'DisplayName','TRRt')
alpha(.1)
axis([0 max(s_plot) -300 300])
title('Shear Force')
xlabel('s  [m]')
ylabel('T(s)    [KN]')
hold on
legend
%% CONFRONTO CON LA LINEA ELASTICA
Metodo_linea_elastica;
s_plot = linspace(0,L,50);
vex_plot=subs(vex,s,s_plot);

subplot(2,2,[1 2])
plot(s_plot,vex_plot,'g','LineWidth',2,'DisplayName','L.E')
title('Deflection with RAYLEGH-RITZ')
axis([0 max(s_plot) 0 0.05])
xlabel('s (m)')
ylabel('v (m)')
axis ij
hold on
legend
%CdS
%Bending moment
subplot(2,2,4)
area(s_plot,M(s_plot),'FaceColor', [0 0.25 0.25],'LineWidth',1.5,'DisplayName','Mesatto')
alpha(.1)
axis([0 max(s_plot) -500 500])
title('Bending Moment')
xlabel('s  [m]')
ylabel('M(s)    [KN*m]' )
axis ij
hold on
legend
%Shear force
subplot(2,2,3)
area(s_plot,T(s_plot),'FaceColor', 'r','LineWidth',1.5,'DisplayName','Tesatto')
alpha(.1)
axis([0 max(s_plot) -300 300])
title('Shear Force')
xlabel('s  [m]')
ylabel('T(s)    [KN]')
hold on
legend
