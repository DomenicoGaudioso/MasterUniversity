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