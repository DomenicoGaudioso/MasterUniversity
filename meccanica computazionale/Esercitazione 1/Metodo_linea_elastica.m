syms s
c = sym('c', [1 4]);
V = sym('V', [5 1]);
V(1)=(q*(s/L))/(E*I);
for t =1:4;
V(t+1)=int(V(t))+c(t);
end
% per far tornare i pedici di v come voglio io, cioè v=v(1) e cosi via
i=[0 0 0 0 1; 0 0 0 1 0; 0 0 1 0 0; 0 1 0 0 0; 1 0 0 0 0]; 
v=i*V;
% condizioni al bordo
cc1=subs(v(1),s,0)==0;
cc2=subs(v(1),s,L)==0;
cc3=subs(v(2),s,0)==0;
cc4=subs(v(2),s,L)==0;
[MatA, b] = equationsToMatrix([cc1 cc2 cc3 cc4],c);
Costanti=linsolve(MatA,b);
vex(s)=(subs(v(1),c,transpose(Costanti)));
%% punto di massimo di vex
d=solve(diff(vex(s)));
vMAX=subs(vex,s,d(3));  % spostamento massimo 
%% Caratteristiche della sollecitazione
M(s)= vpa(-E*I*diff(vex,2));
T(s)=diff(M(s));