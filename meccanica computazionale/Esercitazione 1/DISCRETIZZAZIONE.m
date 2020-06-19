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
%% Rigid-segment discretisation
n=1:6;              % numero di tratti rigidi (discretizzazione)
a=L/max(n);         % lunghezza tratto discretizzato
s=a*n;              % (m) abscissa 


Q=[q*(s-a)/L]*a+0.5*a*[(q*s/L)-(q*(s-a)/L)];   % Risultante del carico per ogni singolo tratto

for j=1:(max(n));   
X(j)=([(q*a*(s(j)-a))/L]*a*0.5+[((q*s(j))/L)-(q*(s(j)-a))/L]*((a.^2)/3))/Q(j);   % Braccio del carico al variare di s
end 

Ko=(E*I)/a ; %rigidezza delle molle

Theta= sym('Theta', [1 max(n)]);  % vettore simbolico delle rotazioni
Theta(max(n))=-[sum(Theta)-Theta(max(n))]; % l'ultimo Theta è la somma di tutti i Theta precedenti

% Procedimento per parametrizzare l'energia di deformazione elastica
i=1:(max(n)-1);    
Uiniziale=vpa(0.5*((2*Ko)*(Theta(1).^2)));
Ufinale=0.5*(2*Ko)*Theta(max(n)).^2;
U(i)=0.5*Ko*((Theta(i+1)-Theta(i)).^2);

U=sum(U);
Utot=Ufinale+Uiniziale+U;

% Procedimento per parametrizzare il lavoro W
sum_theta= sym('sum_angle', [(max(n)-1) 1]);

for j=1:max((n)-2);
sum_theta(1)=Theta(1);
sum_theta(j+1)=sum(sum_theta(j)+Theta(j+1));
end
for j=1:max(n);
W(j)=vpa(Q(j)*Theta(j)*X(j));
end
for j=2:max((n)); 
W_intermedio(j)=vpa(a*sum_theta((j)-1)*Q(j)); 
end 
%W_intermedio è una componente in più, ma essendo zero a me non cambia nulla,
% sommare uno zero non comporta incrementi
Wtot=sum(W)+sum(W_intermedio);
V=Utot-Wtot;

for j=1:(max(n)-1);  
derivate_parziali(j)= diff(V,Theta(j))==0;  % derivate parziali rispetto a c1,c2 e cj di V  (V tilde)
end

Theta= sym('Theta', [1 (max(n)-1)]); % vettore simbolico dei risultati
[MatA,b] =equationsToMatrix(derivate_parziali,[Theta]);     % Convert set of linear equations to matrix form
Theta=double(linsolve(MatA,b))                                       % Solve linear system of equations

%% Deflection
vrig=zeros(1,max(n)+1);
vrig(1)=0;

for j=1:max(n)-1;
vrig(j+1)=vrig(j)+Theta(j)*a;
end

s_plot = linspace(0,L,max(n)+1);
%% Plot deformed shapes
subplot(2,2,[1 2])
plot(s_plot,vrig,'--','LineWidth',2,'color',[1 0.5 0],'DisplayName','rigid segments')
axis([0 max(s_plot) 0 0.06])
title('deformed shapes ')
xlabel('s  [m]')
ylabel('v    [m]')
axis ij
hold on
%% Bending moment
M=sym('M_discr',[1 (max(n)+1)])
M(1)=double(2*Ko*Theta(1))
x=2:(max(n)-1)
M(x)= double(Ko*(Theta(x)-Theta(x-1)))
M(max(n))=double(Ko*((-sum(Theta)-Theta(max(n)-1))))
M(max(n)+1)=double(2*Ko*(sum(Theta)))
%plot(s_plot,M)
subplot(2,2,3)
area(s_plot,M,'FaceColor', [0 0.25 0.25],'LineWidth',2)
alpha(.2)
axis([0 max(s_plot) -300 500])
title('Bending moment')
xlabel('s  [m]')
ylabel('M    [KN*m]')
hold on
%% Shear 
% i=1:max(n)
% dist1=L-X(i)
% i=0:(max(n)-1)
% dist2=i*a
%braccio=dist1-dist2  % braccio delle forze per l'equilibrio al momento nell'incastro finale
% for i=1:max(n)
% R(i)=(1/L)*[braccio(i)*Q(i)]
% end
%R1=sum(R)     % Reazione vincolare all'incastro iniziale 
%R2=sum(Q)-R1  % Reazione vincolare all'incastro finale

%Reazioni ottenuti per equilibrio della struttura
R2=q*L/3-M(1)/L+M(max(n)+1)/L  
R1=q*L/2-R2

T_discr=R1-0.5*((q*((s_plot).^2))/L)
%plot(s_plot,T_discr)
subplot(2,2,4)
area(s_plot,T_discr,'FaceColor', 'r','LineWidth',2)
alpha(.2)
axis([0 max(s_plot) -300 200])
title('Shear')
xlabel('s  [m]')
ylabel('T    [KN]')
hold on

legend off