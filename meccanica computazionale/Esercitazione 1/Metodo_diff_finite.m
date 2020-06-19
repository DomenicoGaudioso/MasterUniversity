digits(5); %  variable precision used
delta=L/(max(i)-5);  % max(i)-5 poichè i tratti in cui è divisa la trave non sono quanto i punti del dominio   
v= sym('v', [1 (max(i))]);
Equazioni= sym('v_', [1 (max(i))]);
%% Equazione di campo
j=0:(max(i)-5)
s_diff=j*delta  % abscissa 

% Procedimento per parametrizzare l'equazioni di campo
x=1:(max(i)-4)
cc=[(q*(s_diff(x)))/(E*I*L)]*(delta^4)
for n=3:(max(i)-2)
    E_c(n)=vpa([v(n+2)-4*v(n+1)+6*v(n)-4*v(n-1)+v(n-2)])  
end 
u=3:(max(i)-2)
E_diff=E_c(u)
E_diff_finite=E_diff==cc

n=3:(max(i)-2)
j=1:(max(i)-4)
Equazioni(n)= E_diff_finite(j)
%% condizioni al contorno
Equazioni(1)=v(3)==0                                   % condizione al contorno 1
n=3
Equazioni(2)=[v(n+1)-v(n-1)]/(2*delta)==0              % condizione al contorno 2
n=(max(i)-2)
Equazioni(max(i)-1)=v(n)==0                            % condizione al contorno 3
Equazioni(max(i))=[v(n+1)-v(n-1)]/(2*delta)==0         % condizione al contorno 4
%% Risoluzione del sistema
[MatA, b] = equationsToMatrix([Equazioni],v)          % Convert set of linear equations to matrix form
vfinito=double(linsolve(MatA,b))                      % Solve linear system of equations     
n=3:(max(i)-2)
v_finito=double(vfinito(n))                           % Deflection
%% Bending moment
n=3:(max(i)-2)
M_finito=[vfinito(n+1)+vfinito(n-1)-2*vfinito(n)]/(delta.^2)
M_diff=double(M_finito*(E*I))
%% Shear force
T_finito=[vfinito(n+2)-2*vfinito(n+1)+2*vfinito(n-1)-vfinito(n-2)]/[2*(delta.^3)]
T_diff=double(T_finito*(E*I))