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
Theta=double(linsolve(MatA,b));                              % Solve linear system of equations

%% Deflection
vrig=zeros(1,max(n)+1);
vrig(1)=0;

for j=1:max(n)-1;
vrig(j+1)=vrig(j)+Theta(j)*a;
end

s_plot = linspace(0,L,max(n)+1);

%% Bending moment
M=sym('M_discr',[1 (max(n)+1)])
M(1)=double(2*Ko*Theta(1))
x=2:(max(n)-1)
M(x)= double(Ko*(Theta(x)-Theta(x-1)))
M(max(n))=double(Ko*((-sum(Theta)-Theta(max(n)-1))))
M(max(n)+1)=double(2*Ko*(sum(Theta)))
%% Shear 
i=1:max(n)
dist1=L-X(i)
i=0:(max(n)-1)
dist2=i*a
braccio=dist1-dist2  % braccio delle forze per l'equilibrio al momento nell'incastro finale
for i=1:max(n)
R(i)=(1/L)*[braccio(i)*Q(i)]
end
R1=sum(R)     % Reazione vincolare all'incastro iniziale 
R2=sum(Q)-R1  % Reazione vincolare all'incastro finale

T_discr=R1-0.5*((q*((s_plot).^2))/L)