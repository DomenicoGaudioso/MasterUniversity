syms z p real
L=10;
f=L/2;
Y(z)=(4*f*(L-z)*z)/L^2;
N(z)=simplify(-(L^2*p*sqrt(1+(16*f^2*(L-2*z)^2)/L^4))/(8*f));
a=0:0.1:L;       
ang(z)=atan(diff(Y,z))
angle=ang(a)
NE=double(N(a)/p)
NE2=double(N(L/2)/p)
%% Plot
figure(1)
c=1/5;
pos1 = [-0.06 0.06 0.7 0.8];
subplot('Position',pos1);
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
for n=1:length(a)-1 
 B=double([a(n) a(n+1)  (a(n+1)+c*NE(n+1)*sin(angle(n+1)))  (a(n)+c*NE(n)*sin(angle(n+1)))]);
 C=double([Y(a(n)) Y(a(n+1))  (Y(a(n+1))-c*NE(n+1)*cos(angle(n+1)))  (Y(a(n))-c*NE(n)*cos(angle(n)))]);
hpatch=patch(B,C,'r','Facealpha',0.6);
axis equal
axis off
hold on  
end
%Struttura
plot(a,Y(a),'-k','LineWidth',2)
title(['\fontname{Courier}\fontsize{15}Forza Normale N_(_z_) con f=L/2'],'color','K');
axis equal
axis off
hold on
% plot per i valori
pos2 = [0.7 0.4 0.2 0.2];
subplot('Position',pos2);
area(a,NE,'FaceColor', 'r','LineWidth',1.5)
alpha(.6)
title('diagramma della forza Normale')
xlabel('z  [m]')
ylabel('N(z)/p   [KN]')
hold on
%% F=L
f=L;
Y(z)=(4*f*(L-z)*z)/L^2;
N(z)=simplify(-(L^2*p*sqrt(1+(16*f^2*(L-2*z)^2)/L^4))/(8*f));
a=0:0.05:L;       
ang(z)=atan(diff(Y,z))
angle=ang(a)
NE=double(N(a)/p)/5
NE2=double(N(L/2)/p)
%% Plot f=L
figure(2)
c=1;
pos1 = [-0.06 0.06 0.7 0.8];
subplot('Position',pos1);
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
for n=1:length(a)-1 
 B=double([a(n) a(n+1)  (a(n+1)+c*NE(n+1)*sin(angle(n+1)))  (a(n)+c*NE(n)*sin(angle(n+1)))]);
 C=double([Y(a(n)) Y(a(n+1))  (Y(a(n+1))-c*NE(n+1)*cos(angle(n+1)))  (Y(a(n))-c*NE(n)*cos(angle(n)))]);
hpatch=patch(B,C,'r','Facealpha',0.6);
axis equal
axis off
hold on  
end
%Struttura
plot(a,Y(a),'-k','LineWidth',2)
title(['\fontname{Courier}\fontsize{15}Forza Normale N_(_z_) con f=L'],'color','K');
axis equal
axis off
hold on
% plot per i valori
pos2 = [0.7 0.4 0.2 0.2];
subplot('Position',pos2);
area(a,NE,'FaceColor', 'r','LineWidth',1.5)
alpha(.6)
title('diagramma della forza Normale')
xlabel('z  [m]')
ylabel('N(z)/p   [KN]')
hold on