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
I=Ix;                                                      % poichè in questo caso è chiamato in causa soltanto Ix
%% data of the problem
E = 21E7;                         % (KN/m^2) Young's modulus steel
Matricola=506682;
L = Matricola/(40*1000);          % (m) length 
q=Matricola/10000;                % (KN/m) max forze on triangle load     
%% METODO DELLA LINEA ELASTICA (exact solution)
digits(5); %  variable precision used
syms s;    %  abscissa

c = sym('c', [1 4]);  % definisco un vettore di costanti
V = sym('V', [5 1]);  % definisco un vettore per le derivate, dove:
                      %V(1)=derivata quarta dello spostamento;
                      %V(2)=derivata terza dello spostamento;
                      %V(3)=derivata seconda dello spostamento;
                      %V(4)=derivata prima dello spostamento;
                      %V(5)=spostamento.

V(1)=((q*(s/L))/(E*I));
for t =1:4   
V(t+1)=(int(V(t))+c(t))
end

% Siccome non mi piace che V(5) sia lo spostamento e V(1) la sua derivata quarta, 
% definisco una matrice "i" che mi inverta questo vettore. 
i=[0 0 0 0 1; 0 0 0 1 0; 0 0 1 0 0; 0 1 0 0 0; 1 0 0 0 0]; 
v=i*V;

% condizioni al bordo
cc1=subs(v(1),s,0)==0;
cc2=subs(v(1),s,L)==0;
cc3=subs(v(2),s,0)==0;
cc4=subs(v(2),s,L)==0;

[MatA, b] = equationsToMatrix([cc1 cc2 cc3 cc4],c);  % Convert set of linear equations to matrix form
Costanti=linsolve(MatA,b)                            % Solve linear system of equations

vex(s)=vpa(subs(v(1),c,transpose(Costanti)))
%% Maximum displacement point
d=solve(diff(vex(s)));
vMAX=subs(vex,s,d(3));  % maximum displacement 
%% Caratteristiche della sollecitazione
M(s)= vpa(-E*I*diff(vex,2)) % Bending moment
T(s)=diff(M(s))             % Shear force
%% Graphics
% Deflection
s_plot = linspace(0,L,50)
vex_plot=subs(vex,s,s_plot)

subplot(2,2,[1,2])
plot(s_plot,vex_plot,'g','LineWidth',2)
axis([0 max(s_plot) 0 0.05])
title('Deflection')
xlabel('s [m]')
ylabel('v [m]')
axis ij
hold on

%CdS
%Bending moment
d_Mnullo=solve(M,s,'MaxDegree',3)     % per trovare la distanza dove il momento è nullo
d_Mnullo1=(d_Mnullo(2))
d_Mnullo2=(d_Mnullo(3))

s1_M=linspace(0,d_Mnullo1,10)
s2_M=linspace(d_Mnullo1,d_Mnullo2,10)
s3_M=linspace(d_Mnullo2,L,10)

subplot(2,2,4)
area(s1_M,M(s1_M),'FaceColor', 'r','LineWidth',1.5)
axis ij
alpha(.5)
hold on
area(s2_M,M(s2_M),'FaceColor', [0 0.25 0.25],'LineWidth',1.5)
axis ij
alpha(.5)
hold on
area(s3_M,M(s3_M),'FaceColor', 'r','LineWidth',1.5)
axis([0 max(s_plot) -500 500])
alpha(.5)
title('Bending Moment')
xlabel('s  [m]')
ylabel('M(s)    [KN*m]' )
axis ij
hold on

%Shear force
d=solve(T)                                    
d_Tnullo=(d(2))
Mmax=subs(M(s),s,d_Tnullo)
s1_taglio=linspace(0,d_Tnullo,d_Tnullo*2)
s2_taglio=linspace(d_Tnullo,L,d_Tnullo*2)

subplot(2,2,3)
area(s1_taglio,T(s1_taglio),'FaceColor', '[0 0.6 0.8]','LineWidth',1.5)
alpha(.5)
hold on
area(s2_taglio,T(s2_taglio),'FaceColor', '[1 0.9 0]','LineWidth',1.5)
axis([0 max(s_plot) -300 300])
alpha(.5)
title('Shear Force')
xlabel('s  [m]')
ylabel('T(s)    [KN]')
hold on
