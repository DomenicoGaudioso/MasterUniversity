%% Esercizio 6- Teoria delle strutture
clear all
close all
%% parametri simbolici
syms  z A B C  real 
%%
L=10; % larghezza arco
f=L/2; % altezza massima dekk'arco
p=1; %carico verticale
q=p; % carico orizzontale
H=L/20; % altezza della sezione
E=3E7; % Young's modulus concrete [KN/m^2]
b=0.25 %m 
h=L/20 
digits(5)
J=vpa((b*h^3)/12)
Y(z)=(4*f*(L-z)*z)/L^2;
N(z)=simplify((6*q* (L - 2*z) - L*p* (1 + (2 - (4*z)/L)^2))/(4*sqrt(1+(2 -(4 *z)/L)^2)));
M(z)=simplify((q*z *(L^2 - 3* L* z + 2*z^2))/L);
T(z)=simplify((q *(L^2 - 8*L*z + 8*z^2))/(2*L*sqrt(1 + (2 - (4*z)/L)^2)));
a=0:0.1:L;       
ang(z)=atan(diff(Y,z));
angle=ang(a);
NE=double(N(a));
ME=double(M(a));
TE=double(T(a));
NE2=double(N(L/2));
%% Equazione dell'intradosso
y(z)=A*z^2+B*z+C;
delta=B^2-4*A*C;

E1=-B/(2*A)==L/2;                            % poiché il punto Vx del vertice L/2
E2=-delta/(4*A)==f-(L/40);                 % poiché il punto Vy del vertice è a f
E3=y((H/2)*sin(angle(1)))==-(H/2)*cos(angle(1));         % poichè la parabola passa per questo punto

S=solve([E1 E2 E3],[A B C]);     %risolvo il sistema( è un sistema non lineare)
A1=S.A(1); B1=S.B(1); C1=S.C(1);  % quindi scelgo i punti che davvero soddisfano l'equazione e non sono soluzioni banali
Y1(z)=simplify(subs(y,[A B C],[A1 B1 C1]));
a1=double((H/2)*sin(angle(1))):0.1:double(L-1*(H/2)*sin(angle(1)))
%% Equazione dell'estradosso
E4=-B/(2*A)==L/2;                                                    % poiché il punto Vx del vertice L/2
E5=-delta/(4*A)==f+(L/40);                                       % poiché il punto Vy del vertice è a f
E6=y(-(H/2)*sin(angle(1)))==(H/2)*cos(angle(1));         % poichè la parabola passa per questo punto

S2=solve([E4 E5 E6],[A B C]);     %risolvo il sistema( è un sistema non lineare)
A2=S2.A(1); B2=S2.B(1); C2=S2.C(1);  % quindi scelgo i punti che davvero soddisfano l'equazione e non sono soluzioni banali
Y2(z)=simplify(subs(y,[A B C],[A2 B2 C2]));
a2=-(H/2)*sin(angle(1)):0.1:L+1*(H/2)*sin(angle(1))
%% Equazione dell'intradosso del terzo medio
y(z)=A*z^2+B*z+C;
delta=B^2-4*A*C;

E7=-B/(2*A)==L/2;                            % poiché il punto Vx del vertice L/2
E8=-delta/(4*A)==f-(H/6);                 % poiché il punto Vy del vertice è a f
E9=y((H/6)*sin(angle(1)))==-(H/6)*cos(angle(1));         % poichè la parabola passa per questo punto

S3=solve([E7 E8 E9],[A B C]);     %risolvo il sistema( è un sistema non lineare)
A3=S3.A(1); B3=S3.B(1); C3=S3.C(1);  % quindi scelgo i punti che davvero soddisfano l'equazione e non sono soluzioni banali
Y3(z)=simplify(subs(y,[A B C],[A3 B3 C3]));
a3=double((H/6)*sin(angle(1))):0.1:double(L-1*(H/6)*sin(angle(1)))
%% Equazione dell'estradosso del terzo medio
E10=-B/(2*A)==L/2;                            % poiché il punto Vx del vertice L/2
E11=-delta/(4*A)==f+(H/6);                 % poiché il punto Vy del vertice è a f
E12=y(-(H/6)*sin(angle(1)))==(H/6)*cos(angle(1));         % poichè la parabola passa per questo punto

S4=solve([E10 E11 E12],[A B C]);     %risolvo il sistema( è un sistema non lineare)
A4=S4.A(1); B4=S4.B(1); C4=S4.C(1);  % quindi scelgo i punti che davvero soddisfano l'equazione e non sono soluzioni banali
Y4(z)=simplify(subs(y,[A B C],[A4 B4 C4]));
a4=-(H/6)*sin(angle(1)):0.1:L+1*(H/6)*sin((angle(1)))
%% Plot N
c=0.2
figure(1)
pos1 = [-0.06 0.06 0.7 0.8];
subplot('Position',pos1);
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
for n=1:length(a)-1 
 NN=NE(n);
 B=double([a(n) a(n+1)  (a(n+1)+c*NE(n+1)*sin(angle(n+1)))  (a(n)+c*NE(n)*sin(angle(n+1)))]);
 C=double([Y(a(n)) Y(a(n+1))  (Y(a(n+1))-c*NE(n+1)*cos(angle(n+1)))  (Y(a(n))-c*NE(n)*cos(angle(n)))]);
 if NN<0
hpatch=patch(B,C,'r','Facealpha',0.6);
 else
hpatch=patch(B,C,'b','Facealpha',0.8);
 end
axis equal
axis off
hold on  
end
%Struttura
plot(a,Y(a),'-k','LineWidth',3)
title(['\fontname{Courier}\fontsize{15}Forza Normale N_(_z_) con f=L/2'],'color','K');
axis off
hold on
% plot per i valori
pos2 = [0.7 0.4 0.2 0.2];
subplot('Position',pos2);
area(a,N(a),'FaceColor', 'r','LineWidth',1.5)
alpha(.6)
title({'diagramma della forza Normale'; 'p=1 & q=p'})
xlabel('z  [m]')
ylabel('N(z)   [KN]')
hold on
%% Plot M
figure(2)
pos1 = [-0.06 0.06 0.7 0.8];
subplot('Position',pos1);
c=0.2
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
for n=1:length(a)-1 
    MM=ME(n);
 B=double([a(n) a(n+1)  (a(n+1)+c*ME(n+1)*sin(angle(n+1)))  (a(n)+c*ME(n)*sin(angle(n+1)))]);
 C=double([Y(a(n)) Y(a(n+1))  (Y(a(n+1))-c*ME(n+1)*cos(angle(n+1)))  (Y(a(n))-c*ME(n)*cos(angle(n)))]);
 if MM<0
hpatch=patch(B,C,'r','Facealpha',0.8);
 else
hpatch=patch(B,C,'b','Facealpha',0.4);
 end
axis equal
axis off
hold on  
end
%Struttura
plot(a,Y(a),'-k','LineWidth',2)
title(['\fontname{Courier}\fontsize{15}Momento Flettente M_(_z_) con f=L/2'],'color','K');
axis equal
axis off
hold on
% plot per i valori
pos2 = [0.7 0.4 0.2 0.2];
subplot('Position',pos2);
area(a,M(a),'FaceColor', 'b','LineWidth',1.5)
alpha(.4)
title({'diagramma del Momento'; 'p=1 & q=p'})
xlabel('z  [m]')
ylabel('M(z)   [KN*m]')
%% Plot T
figure(3)
c=0.2
pos1 = [-0.06 0.06 0.7 0.8];
subplot('Position',pos1);
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
for n=1:length(a)-1 
    TT=TE(n);
 B=double([a(n) a(n+1)  (a(n+1)+c*TE(n+1)*sin(angle(n+1)))  (a(n)+c*TE(n)*sin(angle(n+1)))]);
 C=double([Y(a(n)) Y(a(n+1))  (Y(a(n+1))-TE(n+1)*c*cos(angle(n+1)))  (Y(a(n))-c*TE(n)*cos(angle(n)))]);
 if TT<0
hpatch=patch(B,C,'y','Facealpha',0.6);
 else
hpatch=patch(B,C,'g','Facealpha',0.5);
 end
axis equal
axis off
hold on  
end
%Struttura
plot(a,Y(a),'-k','LineWidth',2)
title(['\fontname{Courier}\fontsize{15}Forza di Taglio T_(_z_) con f=L/2'],'color','K');
axis equal
axis off
hold on
% plot per i valori
pos2 = [0.7 0.4 0.2 0.2];
subplot('Position',pos2);
area(a,T(a),'FaceColor', 'y','LineWidth',1.5)
alpha(.4)
title({'diagramma del Taglio'; 'p=1 & q=p'})
xlabel('z  [m]')
ylabel('T(z)   [KN]')
%% Plot  curva delle pressioni
e(z)=M(z)/N(z)
eE=double((e(a)));
q=p/50
N1(z)=simplify((6*q* (L - 2*z) - L*p* (1 + (2 - (4*z)/L)^2))/(4*sqrt(1+(2 -(4 *z)/L)^2)));
M1(z)=simplify((q*z *(L^2 - 3* L* z + 2*z^2))/L);
e2(z)=M1(z)/N1(z)
e2E=double((e2(a)));
figure(4)
%subplot(1,2,1)
fig = gcf; % current figure handle
fig.Color = [1 1 1];
fig.ToolBar = 'none';
%Struttura
plot(a,Y(a),'--k','LineWidth',1,'DisplayName','line d asse')
title(['\fontname{Courier}\fontsize{15}Curva delle Pressioni e_(_z_) con f=L/2'],'color','K');
hold on
plot(a1,Y1(a1),'-k','LineWidth',1,'DisplayName','intradosso') %intradosso
hold on
plot(a2,Y2(a2),'-k','LineWidth',1,'DisplayName','estradosso') %estradosso
% hold on
plot(a3,Y3(a3),'-b','LineWidth',1,'DisplayName','intradosso terzo medio') %intradosso del terzo medio
hold on
plot(a4,Y4(a4),'-b','LineWidth',1,'DisplayName','estradosso  terzo medio') %estradosso del terzo medio
hold on

for n=1:length(a) 
 ex(n)=double(a(n)-eE(n)*sin(angle(n)));
 ey(n)=double(Y(a(n))+eE(n)*cos(angle(n)));
 
 e2x(n)=double(a(n)-e2E(n)*sin(angle(n)));
 e2y(n)=double(Y(a(n))+e2E(n)*cos(angle(n)));
end
plot(ex,ey,'-g','LineWidth',2,'DisplayName','e(z) con p=1 q=p/14')
hold on
plot(e2x,e2y,'-r','LineWidth',2,'DisplayName','e(z) con p=1 q=p/50')
title({'\fontname{Courier}\fontsize{15}Curva delle Pressioni e_(_z_)'},'color','K');
axis equal
axis off
hold on
legend

figure(5)
q=p/50
c=8000 %amplifica lo spostamento
vx(1)=0;
vx(2)=L/4+c*(L^4* q*(-830*sqrt(2)+ 332*sqrt(5)+27*(-2+sqrt(10))*asinh(2)))/(32768*sqrt(5)*E*J)*sin(ang(L/4)); %a L/4
vx(3)=L/2-c*((L^4*q *(166*sqrt(5)-27*asinh(2)))/(24576*E*J))*sin(ang(L/2)); % alla chiave
vx(4)=(3/4)*L+c*(L^4*q*(-830*sqrt(2)+ 332*sqrt(5)+27*(-2+sqrt(10))*asinh(2)))/(32768*sqrt(5)*E*J)*sin(ang(3*L/4));
vx(5)=L;
vy(1)=0;
vy(2)=Y(L/4)+c*(L^4* q*(-830*sqrt(2)+ 332*sqrt(5)+27*(-2+sqrt(10))*asinh(2)))/(32768*sqrt(5)*E*J)*cos(ang(L/4)); %a L/4
vy(3)=Y(L/2)-c*((L^4*q *(166*sqrt(5)-27*asinh(2)))/(24576*E*J))*cos(ang(L/2)); % alla chiave %0.000028838 m
vy(4)=Y((3/4)*L)-c*(L^4* q*(-830*sqrt(2)+ 332*sqrt(5)+27*(-2+sqrt(10))*asinh(2)))/(32768*sqrt(5)*E*J)*cos(ang(3*L/4));
vy(5)=0;
vxq=vx(1):0.1:vx(5)
vyq=vy(1):0.1:vy(5)
vyIN= interp1(vx,vy,vxq,'cubic')
plot(vx,vy,'*g','LineWidth',2)
hold on
plot(vxq,vyIN,'.-g','LineWidth',2)
hold on
%Struttura
plot(a,Y(a),'-k','LineWidth',2)
title(['\fontname{Courier}\fontsize{15}spostamento v_(_z_) con f=L/2 e q=p/50'],'color','K');
axis equal
axis off
hold on