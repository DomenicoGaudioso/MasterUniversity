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
%% METODO DELLE DIFFERENZE FINITE
digits(5); %  variable precision used
i=1:20;    % Dominio
delta=L/(max(i)-5);  % max(i)-5 poichè i tratti in cui è divisa la trave non sono quanto i punti del dominio   
v= sym('v', [1 (max(i))]);
Equazioni= sym('v_', [1 (max(i))]);
%% Equazione di campo
j=0:(max(i)-5);
s_diff=j*delta;  % abscissa 

% Procedimento per parametrizzare l'equazioni di campo
x=1:(max(i)-4);
cc=[(q*(s_diff(x)))/(E*I*L)]*(delta^4);
for n=3:(max(i)-2);
    E_c(n)=vpa([v(n+2)-4*v(n+1)+6*v(n)-4*v(n-1)+v(n-2)]);  
end 
u=3:(max(i)-2);
E_diff=E_c(u);
E_diff_finite=E_diff==cc;

n=3:(max(i)-2);
j=1:(max(i)-4);
Equazioni(n)= E_diff_finite(j);
%% condizioni al contorno
Equazioni(1)=v(3)==0 ;                                  % condizione al contorno 1
n=3;
Equazioni(2)=[v(n+1)-v(n-1)]/(2*delta)==0;              % condizione al contorno 2
n=(max(i)-2);
Equazioni(max(i)-1)=v(n)==0;                            % condizione al contorno 3
Equazioni(max(i))=[v(n+1)-v(n-1)]/(2*delta)==0;         % condizione al contorno 4
%% Risoluzione del sistema
[MatA, b] = equationsToMatrix([Equazioni],v);          % Convert set of linear equations to matrix form
vfinito=double(linsolve(MatA,b));                      % Solve linear system of equations     
n=3:(max(i)-2);
v_finito=double(vfinito(n));                           % Deflection
%% Graphics  Deflection
subplot(2,2,[1 2]);
plot(s_diff,v_finito,'LineWidth',2,'DisplayName','diff. finite con dominio di 12 tratti');
axis([0 max(s_diff) 0 0.05]);
title('Deflection');
xlabel('s  [m]');
ylabel('v  [m]');
axis ij
hold on
legend
%% Bending moment
n=3:(max(i)-2)
M_finito=[vfinito(n+1)+vfinito(n-1)-2*vfinito(n)]/(delta.^2)
M_diff=double(M_finito*(E*I))
% Graphics
%plot(s_diff,M_diff,'LineWidth',1,'DisplayName','diff. finite con 12 tratti')
subplot(2,2,3)
area(s_diff,M_diff,'FaceColor', [0 0.25 0.25],'LineWidth',2,'DisplayName','diff. finite con dominio di 12 tratti')
axis([0 max(s_diff) -500 500])
alpha(.2)
title('Bending Moment')
xlabel('s  [m]')
ylabel('M    [KN*m]')
hold on
legend
%% Shear force
T_finito=[vfinito(n+2)-2*vfinito(n+1)+2*vfinito(n-1)-vfinito(n-2)]/[2*(delta.^3)]
T_diff=double(T_finito*(E*I))
% Graphics
%plot(s_diff,T_diff,'LineWidth',1,'DisplayName','diff. finite con 12 tratti')
subplot(2,2,4)
area(s_diff,T_diff,'FaceColor', 'r','LineWidth',2,'DisplayName','diff. finite con dominio di 12 tratti')
axis([0 max(s_diff) -200 300])
alpha(.1)
axis ij
title('Shear force')
xlabel('s  [m]')
ylabel('T    [KN]')
hold on

legend