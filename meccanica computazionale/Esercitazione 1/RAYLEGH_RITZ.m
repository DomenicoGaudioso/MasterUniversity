%%  Computational Mechanics - PRIMA ESERCITAZIONE
% Triangular loaded on beam stuck to the ends  
%% Clear workspace and close any open windows
clear all
close all
%% Data of the IPE section
% In our case, IPE 400 section
h=0.4;      % (m) altezza 
b=0.18;     % (m) lunghezza piattabanda
a=0.0086;   % (m) spessore anima
e=0.0133;   % (m) spessore piattabanda

A=2*b*e+(h-2*e)*a;                                         % (m^2) area IPE (trascurando gli smussi)
Ix=2*([(1/12)*b*e^3]+b*e*[(h-e)/2]^2)+(1/12)*a*(h-2*e)^3;  % (m^4) moment of inertia IPE secondo l'asse x
Iy=(2*[(1/12)*e*b^3])+[(1/12)*(h-2*e)*a^3];                % (m^4) moment of inertia IPE secondo l'asse y
I=Ix                                                       % poichè in questo caso è chiamato in causa soltanto Ix
%% data of the problem
E = 21E7;                         % (KN/m^2) Young's modulus steel
Matricola=506682;
L = Matricola/(40*1000);          % (m) length 
q=Matricola/10000;                % (KN/m) max forze on triangle load     
%% METODO DI RAYLEGH-RITZ (con una funzione polinomiale)
digits(5); %  variable precision used
syms s  %  abscissa 

j=1:1 % pedice che determina il grado della funzione
phi=2*s.^(3+(j))-4*L*s.^(2+(j))+2*L^2*s.^(1+(j));   % Funzione di base che rispetta le condizioni al contorno
dphi(j)= diff(phi(j));     % derivata prima della funzioni di base (polinomiale)
d2phi(j)= diff(phi(j),2);  % derivata seconda funzioni di base (polinomiale)

C=sym(['C'], [max(j) 1]);   % vettore simbolico delle costanti c1....cj

v=sum(phi*C);               % spostamento approssimato (v tilde)

u=(sum(d2phi*C)).^2;
U=(1/2)*E*I*int(u,s,0,L);   % Energia di deformazioneelastica approssimata (U tilde)

w=q*(s/L)*(v);
W=int(w,s,0,L);             % Lavoro dei carichi approssimato (W tilde)

V=vpa(U-W) % Energia potenziale totale (V tilde)

for j=1:max(j);  
derivate_parziali(j)= diff(V,C(j))==0;  % derivate parziali rispetto a c1,c2 e cj di V  (V tilde)
end

[MatA, b] =equationsToMatrix(derivate_parziali,[C])  % Convert set of linear equations to matrix form
Costanti_1=linsolve(MatA,b)                          % Solve linear system of equations
vRR1(s)=vpa(subs(v,C,(Costanti_1)))
%% Maximum displacement point of vRR1
d=solve(diff(vRR1));
vRR1_MAX=max(subs(vRR1,s,d(3)));  % maximum displacement
%% Caratteristiche della sollecitazione
M_RR1(s)= vpa(-E*I*diff(vRR1,2)); % Bending moment 
T_RR1(s)=diff(M_RR1(s));          % Shear force
%% Graphics
% Deflection
s_plot = linspace(0,L,50)
vRR1_plot=subs(vRR1,s,s_plot)
%subplot(3,3,1)
subplot(2,2,[1 2])
plot(s_plot,vRR1_plot,'m','LineWidth',2,'DisplayName','Raylegh Ritz 1')
axis([0 max(s_plot) 0 0.05])
title('Deflection with RAYLEGH-RITZ ')
xlabel('s  [m]')
ylabel('v  [m]')
axis ij
hold on

legend
%% METODO DI RAYLEGH-RITZ (con 2 funzione polinomiale)
digits(5);
j=1:2;
procedimento_RR; %file dove ho descritto i procedimenti essenziali
%per il metodo di Raylegh-Ritz
vRR2(s)=vpa(subs(v,C,(Costanti_1)));
%% %% Maximum displacement point of vRR2
d=solve(diff(vRR2));
vRR2_MAX=max(subs(vRR2,s,d(3)));  % maximum displacement
%% Caratteristiche della sollecitazione
M_RR2(s)= vpa(-E*I*diff(vRR2,2))  % Bending moment
T_RR2(s)=diff(M_RR2(s))           % Shear force
%% Graphics
% Deflection
subplot(2,2,[1 2])
s_plot = linspace(0,L,50);
vRR2_plot=subs(vRR2,s,s_plot);
%subplot(3,3,4)

plot(s_plot,vRR2_plot,'k','LineWidth',2,'DisplayName','Raylegh Ritz 2')
axis([0 max(s_plot) 0 0.05])
title('Deflection with RAYLEGH-RITZ')
xlabel('s  [m]')
ylabel('v  [m]')
axis ij
hold on
legend

%% CdS RR1
%Bending moment
d_Mnullo=solve(M_RR1,s,'MaxDegree',3);     % per trovare la distanza dove il momento è nullo
d_Mnullo1=(d_Mnullo(1));
d_Mnullo2=(d_Mnullo(2));
s1_M=linspace(0,d_Mnullo1,10);
s2_M=linspace(d_Mnullo1,d_Mnullo2,10);
s3_M=linspace(d_Mnullo2,L,10);
%subplot(3,3,2);
subplot(2,2,3)
area(s1_M,M_RR1(s1_M),'FaceColor', 'r','LineWidth',1.5)
alpha(.5)
axis ij
hold on
area(s2_M,M_RR1(s2_M),'FaceColor', [0 0.25 0.25],'LineWidth',1.5)
alpha(.4)
axis ij
hold on
area(s3_M,M_RR1(s3_M),'FaceColor', 'r','LineWidth',1.5)
axis([0 max(s_plot) -500 500])
title('Bending Moment')
alpha(.4)
xlabel('s  [m]')
ylabel('M(s)    [KN*m]')
axis ij
hold on
%Shear force
d=solve(T_RR1);                                    
d_Tnullo=(d(1));
M_RR1_max=subs(M_RR1(s),s,d_Tnullo);
s1_taglio=linspace(0,d_Tnullo,d_Tnullo*2);
s2_taglio=linspace(d_Tnullo,L,d_Tnullo*2);
%subplot(2,2,3)
subplot(2,2,4);

area(s1_taglio,T_RR1(s1_taglio),'FaceColor', '[0 0.6 0.8]','LineWidth',2)
alpha(.5)
hold on
area(s2_taglio,T_RR1(s2_taglio),'FaceColor', '[1 0.9 0]','LineWidth',2)
axis([0 max(s_plot) -250 250])
alpha(.5)
title('Shear Force')
xlabel('s  [m]')
ylabel('T(s)    [KN]')
hold on

%% CdS 2
%Bending moment
d_Mnullo=solve(M_RR2,s,'MaxDegree',3);     % per trovare la distanza dove il momento è nullo
d_Mnullo1=(d_Mnullo(2));
d_Mnullo2=(d_Mnullo(3));
s1_M=linspace(0,d_Mnullo1,10);
s2_M=linspace(d_Mnullo1,d_Mnullo2,10);
s3_M=linspace(d_Mnullo2,L,10);
%subplot(3,3,5)
subplot(2,2,3)
area(s1_M,M_RR2(s1_M),'FaceColor', 'r','LineWidth',2)
alpha(.4)
axis ij
hold on
area(s2_M,M_RR2(s2_M),'FaceColor', [0 0.25 0.25],'LineWidth',2)
alpha(.1)
axis ij
hold on
area(s3_M,M_RR2(s3_M),'FaceColor', 'r','LineWidth',2)
alpha(.1)
axis ij
hold on
%Shear force
d=solve(T_RR2);                                    
d_Tnullo=(d(2));
M_RR_max=subs(M_RR2(s),s,d_Tnullo);
s1_taglio=linspace(0,d_Tnullo,d_Tnullo*2);
s2_taglio=linspace(d_Tnullo,L,d_Tnullo*2);
%subplot(3,3,6)
subplot(2,2,4)
area(s1_taglio,T_RR2(s1_taglio),'FaceColor', '[0 0.6 0.8]','LineWidth',2)
alpha(.1)
hold on
area(s2_taglio,T_RR2(s2_taglio),'FaceColor', '[1 0.9 0]','LineWidth',2)
alpha(.1)
hold on
%% confronto con la linea elastica
% Linea elastica 
Metodo_linea_elastica;
subplot(2,2,[1 2])
s_plot = linspace(0,L,50)
vex_plot=subs(vex,s,s_plot)
plot(s_plot,vex_plot,'--g','LineWidth',1.5,'DisplayName','v_e_x_a_c_t')
hold on
legend