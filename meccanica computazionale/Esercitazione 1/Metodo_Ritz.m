%% METODO DI RAYLEGH-RITZ (con una funzione polinomiale)
digits(5); %  variable precision used
syms s  %  abscissa 

j=1:1 % pedice che determina il grado della funzione
phi=2*s.^(3+(j))-4*L*s.^(2+(j))+2*L^2*s.^(1+(j));   % Funzione di base che rispetta le condizioni al contorno
dphi(j)= diff(phi(j));     % derivata prima della funzioni di base (polinomiale)
d2phi(j)= diff(phi(j),2);  % derivata seconda funzioni di base (polinomiale)

C=sym(['C'], [max(j) 1]);

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
%% Caratteristiche della sollecitazione
M_RR1(s)= vpa(-E*I*diff(vRR1,2)); % Bending moment 
T_RR1(s)=diff(M_RR1(s));          % Shear force
%% Graphics
% Deflection
s_plot = linspace(0,L,50)
vRR1_plot=subs(vRR1,s,s_plot)
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
%% Caratteristiche della sollecitazione
M_RR2(s)= vpa(-E*I*diff(vRR2,2))  % Bending moment
T_RR2(s)=diff(M_RR2(s))           % Shear force
%% Graphics
% Deflection
s_plot = linspace(0,L,50);
vRR2_plot=subs(vRR2,s,s_plot);
plot(s_plot,vRR2_plot,'k','LineWidth',2,'DisplayName','Raylegh Ritz 2')
axis([0 max(s_plot) 0 0.06])
title('Deflection with RAYLEGH-RITZ')
xlabel('s  [m]')
ylabel('v  [m]')
axis ij
hold on
%% CdS 1
%Bending moment
%subplot(2,2,3)
plot(s_plot,M_RR1(s_plot),'LineWidth',1.5,'DisplayName','Raylegh Ritz 1')
title('Bending Moment')
xlabel('s  [m]')
ylabel('M(s)    [KN*m]')
axis ij
%Shear force
subplot(2,2,4);
plot(s_plot,T_RR1(s_plot),'LineWidth',1.5,'DisplayName','Raylegh Ritz 1')
title('Shear Force')
xlabel('s  [m]')
ylabel('T(s)    [KN]')
%% CdS 2
%Bending moment
plot(s_plot,M_RR2(s_plot),'LineWidth',1.5,'DisplayName','Raylegh Ritz 2')
axis([0 max(s_plot) -500 500])
title('Bending Moment')
xlabel('s  [m]')
ylabel('M(s)    [KN*m]')
axis ij
hold on
legend
%Shear force
plot(s_plot,T_RR2(s_plot),'FaceColor', '[1 0.9 0]','LineWidth',2)
axis([0 max(s_plot) -250 250])
title('Shear Force')
xlabel('s  [m]')
ylabel('T(s)    [KN]')
hold on
legend