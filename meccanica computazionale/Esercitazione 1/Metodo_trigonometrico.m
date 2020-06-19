for j=1:max(j);
phi(j)=1-cos(2*(j)*pi*(s)/(L));  % Funzioni di base (polinomiale)
dphi(j)= diff(phi(j));             % derivata prima Funzioni di base (polinomiale)
d2phi(j)= diff(phi(j),2);  % derivata seconda Funzioni di base (polinomiale)
end

C=sym(['C'], [max(j) 1]);

vRR=sum(phi*C);              % spostamento approssimato (v tilde)
u=(sum(d2phi*C)).^2;  

U=(1/2)*E*I*int(u,s,0,L);   %Energia di deformazioneelastica approssimata (U tilde)

w=q*(s/L)*(vRR);
W=int(w,s,0,L);             % Lavoro dei carichi approssimato (W tilde)

V=(U-W);                    % Energia potenziale totale (V tilde)


for j=1:max(j);  
derivate_parziali(j)= diff(V,C(j))==0; % derivate parziali rispetto a c1,c2 e cj di V  (V tilde)
end

[MatA, b] =equationsToMatrix(derivate_parziali,[C]);
Costanti_1=linsolve(MatA,b);

vRRt=subs(vRR,C,(Costanti_1));

%% Grafici 
% deformata
s_plot = linspace(0,L,50);
vRRt_plot=subs(vRRt,s,s_plot);