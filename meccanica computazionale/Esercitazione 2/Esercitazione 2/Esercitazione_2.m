%%  Computational Mechanics - SECONDA ESERCITAZIONE
% TRUSS STRUCTURE
%% Clear workspace and close any open windows
clear all
close all
%% Set Command Window output display format
format bank;
%% ALTER THE NEXT LINES TO CHOOSE AN OUTPUT FILE FOR THE RESULTS
% Open file for output of results
fid = fopen('Ese2_results.txt','w'); 
%% ALTER THE NEXT LINE TO CHOOSE AN INPUT FILE
disp('Results printed in file : truss_1_results.txt');
%% data of the problem                      
Matricola=506682.;                 % my university number
a= Matricola/(200*1000);           % [m] length 
t=Matricola/20000;                 % thermal variation [°C]
alpha=12E-6;                       %thermal expansion coefficient
E=21E7;                            % Young's modulus steel [KN/m^2]
%% data structure
nnd = 12;              % Number of nodes:
nel = 20;               % Number of elements:
nne = 2 ;              % Number of nodes for element:
nodof =2 ;             % Number of degrees of freedom for node
eldof = nne*nodof;     % Number of degrees of freedom for element
%% Nodes cFoordinates X and      Y
x(1)=0.0;  y(1)=0.0;  % X and Y coord. node 1
x(2)=a;    y(2)=0;    % X and Y coord. node 2
x(3)=0.0;  y(3)=2*a;  % X and Y coord. node 3
x(4)=a;    y(4)=2*a;  % X and Y coord. node 4
x(5)=-a;   y(5)=4*a;  % X and Y coord. node 5
x(6)=0.;   y(6)=4*a;  % X and Y coord. node 6
x(7)=a;    y(7)=4*a;  % X and Y coord. node 7
x(8)=3*a;  y(8)=4*a;  % X and Y coord. node 8
x(9)=5*a;  y(9)=4*a;  % X and Y coord. node 9
x(10)=0.0; y(10)=5*a; % X and Y coord. node 10
x(11)=a;   y(11)=5*a; % X and Y coord. node 11
x(12)=3*a; y(12)=5*a; % X and Y coord. node 12 
%% Element connectivity
conn(1,1)=1;       conn(1,2)=3;    % 1st and 2nd node of element 1
conn(2,1)=2;       conn(2,2)=4;    % 1st and 2nd node of element 2
conn(3,1)=2;       conn(3,2)=3;    % 1st and 2nd node of element 3
conn(4,1)=3;       conn(4,2)=4;    % 1st and 2nd node of element 4
conn(5,1)=4;       conn(5,2)=7;    % 1st and 2nd node of element 5
conn(6,1)=3;       conn(6,2)=6;    % 1st and 2nd node of element 6
conn(7,1)=3;       conn(7,2)=7;    % 1st and 2nd node of element 7
conn(8,1)=6;       conn(8,2)=7;    % 1st and 2nd node of element 8
conn(9,1)=5;       conn(9,2)=6;    % 1st and 2nd node of element 9
conn(10,1)=5;      conn(10,2)=10;  % 1st and 2nd node of element 10
conn(11,1)=6;      conn(11,2)=10;  % 1st and 2nd node of element 11
conn(12,1)=6;      conn(12,2)=11;  % 1st and 2nd node of element 12
conn(13,1)=7;      conn(13,2)=11;  % 1st and 2nd node of element 13
conn(14,1)=10;     conn(14,2)=11;  % 1st and 2nd node of element 14
conn(15,1)=7;      conn(15,2)=8;   % 1st and 2nd node of element 15
conn(16,1)=8;      conn(16,2)=11;  % 1st and 2nd node of element 16
conn(17,1)=11;     conn(17,2)=12;  % 1st and 2nd node of element 17
conn(18,1)=8;      conn(18,2)=12;  % 1st and 2nd node of element 18
conn(19,1)=8;      conn(19,2)=9;   % 1st and 2nd node of element 19
conn(20,1)=9;      conn(20,2)=12;  % 1st and 2nd node of element 20
%% Area section [m^2]
A_m=0.005383;    % Area montanti verticali_HE200A [m^2]
A_c=0.005381;    % Area correnti orizzontali_IPE300 [m^2]
A_d=0.005508;    % Area diagonali_2L_120x120x12 [m^2]

A=zeros(nel,1);
A=[A_m; ... %  A of element 1
   A_m; ... %  A of element 2
   A_d; ... %  A of element 3
   A_c; ... %  A of element 4
   A_m; ... %  A of element 5
   A_m; ... %  A of element 6
   A_d; ... %  A of element 7
   A_c; ... %  A of element 8
   A_c; ... %  A of element 9
   A_d; ... %  A of element 10
   A_m; ... %  A of element 11
   A_d; ... %  A of element 12
   A_m; ... %  A of element 13
   A_c; ... %  A of element 14
   A_c; ... %  A of element 15
   A_d; ... %  A of element 16
   A_c; ... %  A of element 17
   A_m; ... %  A of element 18
   A_c; ... %  A of element 19
   A_d];    %  A of element 20 

%% Element properties
for e = 1:nel;
    I = conn(e,1); % node start
    J = conn(e,2); % node end
    % Element length
    L(e) = sqrt((x(J)-x(I))^2 + (y(J)-y(I))^2);
    % angle formed by the element with the horizontal
    theta(e) = atan2(y(J)-y(I), x(J)-x(I));
    % stiffness matrix in local
    KLoc(:,:,e) =((E*A(e))/L(e))*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0] ;
    % rotation matrix
    Q = [ cos(theta(e)) -sin(theta(e)); sin(theta(e)) cos(theta(e)) ];
    % transfer matrix
    T = [ Q zeros(2,2); zeros(2,2) Q ];
    % stiffness matrix in global
    KGlob(:,:,e) = T*KLoc(:,:,e)*T';
end

%% Mass of the structure
% specific weight ordinary steel
gamma=78.5; % [kN/m3]
for i=1:nel;
M(i)=gamma*A(i).*L(i); % Mass element [KN]
end
Mtot=sum(M); % Total Mass [KN]

%% Global stiffness matrix
K = zeros(nne*nnd,nne*nnd); % dimensions of the global stiffness matrix
for e = 1:nel;
    I = conn(e,1); % node start
    J = conn(e,2); % node end
    
 % transfer of the submatrices to their position in the stiffness matrix
 % global
    K(2*I-1:2*I,2*I-1:2*I) = K(2*I-1:2*I,2*I-1:2*I) + KGlob(1:2,1:2,e);
    K(2*I-1:2*I,2*J-1:2*J) = K(2*I-1:2*I,2*J-1:2*J) + KGlob(1:2,3:4,e);
    K(2*J-1:2*J,2*I-1:2*I) = K(2*J-1:2*J,2*I-1:2*I) + KGlob(3:4,1:2,e);
    K(2*J-1:2*J,2*J-1:2*J) = K(2*J-1:2*J,2*J-1:2*J) + KGlob(3:4,3:4,e);
end
%% Load vector in [KN]
p= zeros(nnd,nne);       % concentrated forces
p(3,:)=[20.0 0.0];       %forces in X and Y directions at node 3
p(5,:)=[15.0 -20.0];     %forces in X and Y directions at node 5
p(10,:)=[5.0 0.0];       %forces in X and Y directions at node 10
p(9,:)=[0.0 -80.0];      %forces in X and Y directions at node 10

pt=zeros(nnd,nne); % thermal forces
pt(5,:)=[-E*A(9)*alpha*t  0.0];                   % thermal forces in X and Y directions at node 5
pt(6,:)=[+E*A(9)*alpha*t-E*A(8)*alpha*t  0.0];    % thermal forces in X and Y directions at node 6
pt(7,:)=[+E*A(8)*alpha*t-E*A(15)*alpha*t  0.0];   % thermal forces in X and Y directions at node 7
pt(8,:)=[+E*A(15)*alpha*t-E*A(19)*alpha*t  0.0];  % thermal forces in X and Y directions at node 8
pt(9,:)=[+E*A(19)*alpha*t  0.0];                  % thermal forces in X and Y directions at node 9

pq=zeros(nnd,nne); % distributed forces
q=50; % distributed load
pq(10,:)=[0.0 -q*L(14)/2];          % distributed load in X and Y directions at node 10
pq(11,:)=[0.0 -q*(L(14)+L(17))/2];  % distributed load in X and Y directions at node 11
pq(12,:)=[0.0 -q*(L(17))/2];        % distributed load in X and Y directions at node 12

P=(p+pq)-pt; % Load vector

% This cycle trasforms the matrix [nndx2;2] into a single vector
nf=ones(nnd,nne);
c=0; for i=1:nnd;
for j=1:2;
c=c+1;
nf(i,j)=c;
end
end

for i=1:nnd;
for j=1:2;
F(nf(i,j))=P(i,j);
end
end
%% Restraints
Kstar = K;

Kstar(1:4,:) = zeros(4,2*nnd);
Kstar(:,1:4) = zeros(2*nnd,4);

Kstar(1,1) = 1;
Kstar(2,2) = 1;
Kstar(3,3) = 1;
Kstar(4,4) = 1;

pstar =transpose(F);
%% Solution
v= inv(Kstar)*pstar;  % nodal diplacements 
f= K*v;               % reaction forces in the nodes
r =(f-pstar);         % reaction force in external constraints

%% Exact reaction force in external constraints

Ry2=20*2-20+15*4+5*5+80*5+50*(a/2)+50*4*a;
Ry1=-Ry2+20+80+50*3*a;

%% Normal forces

i=1:2:nne*nnd;
vx=v(i);  % nodal diplacements in X directions in global 
fx=f(i);  % reaction forces in the nodes in X directions in global
j=2:2:nne*nnd;
vy=v(j);  % nodal diplacements in Y directions in global
fy=f(j);  % reaction forces in the nodes in Y directions in global

for e=1:nel;
i = conn(e,1); % node start
j= conn(e,2);  % node end
vglobal(:,:,e)=[vx(i) vy(i);vx(j) vy(j)]; % matrix of the displacements of the element in global
theta(e) = atan2(y(j)-y(i), x(j)-x(i));   % angle formed by the element with the horizontal
Q=[cos(theta(e)) -sin(theta(e)); sin(theta(e)) cos(theta(e))]; % rotation matrix
vlocal(:,:,e)=vglobal(:,:,e)*Q; % matrix of the displacements of the element in local
L(e)=sqrt((x(j)-x(i))^2 + (y(j)-y(i))^2); % Element length
epslon(e)=(vlocal(2,1,e)-vlocal(1,1,e))/L(e); % epslon total 
end

% conditions for thermal variations
epslon_elastic=epslon; % epslon elastic

epslon_elastic(8)=epslon(8)-(-alpha*t);   % epslon elastic element 8
epslon_elastic(9)=epslon(9)-(-alpha*t);   % epslon elastic element 9
epslon_elastic(15)=epslon(15)-(-alpha*t); % epslon elastic element 15
epslon_elastic(19)=epslon(19)-(-alpha*t); % epslon elastic element 19

for e=1:nel
N(e)=E*A(e)*(epslon_elastic(e)); % Normal force in Local
end

%% Bending Moment and Share for elements loaded with distributed load
s=0:L(14)/20:L(14);
M14=-fy(10)*s-(q*(s.^2))/2;
T14=-fy(10)-q*s;

s=0:L(17)/20:L(17);
M17=-(fy(11)-fy(10))*s-(q*(s.^2))/2;
T17=-(fy(11)-fy(10))-q*s;
%% Print model data
print_dati_Ese2 
%%
print_results_Ese2;
%
fclose(fid);

%% Post-Processing
%% Plot structure 
figure(1)
pos0 = [0.06 0.03 0.94 0.94];
subplot('Position',pos0)
for i=1:nel
node_1=conn(i,1); % node start
node_2=conn(i,2); % node end
xx=[x(node_1),x(node_2)];
yy=[y(node_1),y(node_2)];
S(i)=plot(xx,yy,'K-o','LineWidth',1.5,'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5],'DisplayName','structure');
axis equal
axis([-5 14 -1 15])
axis off
hold on

%% rectangle
width = 0.4; % whatever
height = 0.3; % whatever...
xCenter =(x(node_1)+x(node_2))/2; % Wherever...
yCenter =(y(node_1)+y(node_2))/2; % Wherever...
xLeft = xCenter - width/2;
yBottom = yCenter - height/2;
rectangle('Position', [xLeft, yBottom, width, height], 'EdgeColor', 'b', 'FaceColor', 'w', 'LineWidth', 1);

text((x(node_1)+x(node_2))/2,(y(node_1)+y(node_2))/2,{(i)},'HorizontalAlignment','center','Color','m','FontSize',8);
axis off
for i=1:nnd
text(x(i)+0.6,y(i)+0.25,{(i)},...
    'HorizontalAlignment','right',...
    'FontSize',6)
axis off
hold on
end
end
%% Plot deformed structure
for i=1:nel;
node_1=conn(i,1); % node start
node_2=conn(i,2); % node end
xx=[x(node_1)+10*vx(node_1),x(node_2)+10*vx(node_2)];
yy=[y(node_1)+10*vy(node_1),y(node_2)+10*vy(node_2)];
D(i)=plot(xx,yy,'--o','Color',[0 1 0],'LineWidth',1,'MarkerSize',5,'MarkerEdgeColor',[0 0.3 0],'MarkerFaceColor',[0 0.8 0],'DisplayName','deformed structure');
axis([-5 14 -0.06  15]);
axis off
title(['\fontname{Courier}\fontsize{13}Plot Structure & Deformed '],'color','k');
end
%% print forces
% distributed load
rectangle('Position', [0,5.25*a, 3*a, (5.6*a-5.3*a)+0.15], 'EdgeColor', 'r', 'FaceColor', 'w', 'LineWidth', 1);
text(-1,5.4*a,'50 KN/m','HorizontalAlignment','center','Color','r','FontSize',8);
axis off
hold on
for i=0:(L(14)+L(17))/15:(L(14)+L(17))
plot(i,5.3*a,'vr','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
p0 = [i 5.6*a]; 
p1 = [i 5.3*a]; 
vectarrow(p0,p1)
axis off
hold on
end
% external constraints
plot(0.0,-0.06,'^k','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','w')
plot(a,-0.06,'^k','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','w')





% Concentrated loads in X direction
text(-a-1.5,4*a+0.25,'15 KN','HorizontalAlignment','center','Color','r','FontSize',8);
plot(-a-0.3,4*a,'>r','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
p0 = [-a-1.5 4*a]; 
p1 = [-a-0.3 4*a]; 
vectarrow(p0,p1)
axis off


text(-1.5,5*a+0.25,'5 KN','HorizontalAlignment','center','Color','r','FontSize',8);
plot(-0.3,5*a,'>r','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
p0 = [-1 5*a]; 
p1 = [-0.3 5*a]; 
vectarrow(p0,p1)
axis off

text(-1.5,2*a+0.25,'20 KN','HorizontalAlignment','center','Color','r','FontSize',8);
plot(-0.3,2*a,'>r','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
p0 = [-1.5 2*a]; 
p1 = [-0.3 2*a]; 
vectarrow(p0,p1)
axis off

% Concentrated loads in Y direction
text(-a+1,4*a-1,'20 KN','HorizontalAlignment','center','Color','r','FontSize',8);
plot(-a,3.4*a,'vr','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
p0 = [-a 3.9*a]; 
p1 = [-a 3.4*a]; 
vectarrow(p0,p1)
axis off

text(5*a-1,4*a-1,'80 KN','HorizontalAlignment','center','Color','r','FontSize',8);
plot(5*a,3.1*a,'vr','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
p0 = [5*a 3.9*a]; 
p1 = [5*a 3.1*a]; 
vectarrow(p0,p1)
axis off


% Thermal variations
width = 0.2; % whatever
height = 0.5; % whatever...
xCenter =(-a*0.7); % Wherever...
yCenter =(4*a); % Wherever...
xLeft = xCenter - width/2;
yBottom = yCenter - height/2;
rectangle('Position', [xLeft, yBottom, width, height], 'EdgeColor', 'k', 'FaceColor', 'g', 'LineWidth', 1);
text(-a*0.7,3.8*a,'-t','HorizontalAlignment','center','Color','k','FontSize',7);
axis off

width = 0.2; % whatever
height = 0.5; % whatever...
xCenter =(a*0.7); % Wherever...
yCenter =(4*a); % Wherever...
xLeft = xCenter - width/2;
yBottom = yCenter - height/2;
rectangle('Position', [xLeft, yBottom, width, height], 'EdgeColor', 'k', 'FaceColor', 'g', 'LineWidth', 1);
text(a*0.7,3.8*a,'-t','HorizontalAlignment','center','Color','k','FontSize',7);
axis off

width = 0.2; % whatever
height = 0.5; % whatever...
xCenter =(a+2*a*0.7); % Wherever...
yCenter =(4*a); % Wherever...
xLeft = xCenter - width/2;
yBottom = yCenter - height/2;
rectangle('Position', [xLeft, yBottom, width, height], 'EdgeColor', 'k', 'FaceColor', 'g', 'LineWidth', 1);
text(a+2*a*0.7,3.8*a,'-t','HorizontalAlignment','center','Color','k','FontSize',7);
axis off

width = 0.2; % whatever
height = 0.5; % whatever...
xCenter =(3*a+2*a*0.7); % Wherever...
yCenter =(4*a); % Wherever...
xLeft = xCenter - width/2;
yBottom = yCenter - height/2;
rectangle('Position', [xLeft, yBottom, width, height], 'EdgeColor', 'k', 'FaceColor', 'g', 'LineWidth', 1);
text(3*a+2*a*0.7,3.8*a,'-t','HorizontalAlignment','center','Color','k','FontSize',7);
axis off
hold on

%legend
legend([S(nel),D(nel)],'Structure','Deformed')
legend('Location',[0.58 0.53 0.03 0.03])
legend('boxoff')
%% Value Table vx,vy 
i=1:nnd;
T3=table(i',vx,vy,'VariableNames',{'node','vx','vy'});
% Get the table in string form
TString=evalc('disp(T3)');
%Use TeX Markup for bold formatting and underscores.
TString=strrep(TString,'<strong>','\bf');
TString=strrep(TString,'</strong>','\rm');
TString=strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth=get(0,'FixedWidthFontName');
% Output the table using the annotation command
annotation(gcf,'TextBox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0.53 0.06 0.18 0.42]);
hold on;
axis off
%% Plot Tension/Compression forces
figure(2)
pos2 = [0.05 0.02 0.94 0.94];
subplot('Position',pos2)
for i=1:nel
node_1=conn(i,1); % initial node of the element i
node_2=conn(i,2); % final node of element i
xx=[x(node_1),x(node_2)];
yy=[y(node_1),y(node_2)];
if N(i)>0;
plot(xx,yy,'r-o','LineWidth',2,'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5],'DisplayName','Tension');
axis equal
axis([-7 14 -0.06 13]);
hold on
axis off
else
plot(xx,yy,'b-o','LineWidth',2,'MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5],'DisplayName','Compression');
axis equal;
axis([-7 14 -0.06 13]);
title(['\fontname{Courier}\fontsize{13} Tension/Compression Forces'],'color','K');
hold on
axis off
end
%% rectangle 
width = 1; % whatever
height = 0.3; % whatever...
xCenter =(x(node_1)+x(node_2))/2; % Wherever...
yCenter =(y(node_1)+y(node_2))/2; % Wherever...
xLeft = xCenter - width/2;
yBottom = yCenter - height/2;
rectangle('Position', [xLeft, yBottom, width, height], 'EdgeColor', 'k', 'FaceColor', 'w', 'LineWidth', 1);
text((x(node_1)+x(node_2))/2,(y(node_1)+y(node_2))/2,{N(i)},'HorizontalAlignment','center','Color','k','FontSize',6);
axis off
end
% external constraints
plot(0.0,-0.06,'^k','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','w')
plot(a,-0.06,'^k','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','w')
%legend
legend('Tension','Compression')
legend('Location',[0.63 0.4 0.03 0.03])
legend('boxoff')
%% Plot Diagram Normal Forces
figure(3)
pos2 = [0.00 0.03 0.94 0.94];
subplot('Position',pos2);
for e=1:nel  
node_1=conn(e,1); % node start
node_2=conn(e,2); % node end
theta(e) = atan2(y(node_2)-y(node_1), x(node_2)-x(node_1));
Q=[cos(theta(e)), -sin(theta(e)), sin(theta(e)), cos(theta(e))];
if rad2deg(theta(e))==90 && N(e)>0;
B=[x(node_1) x(node_2) (x(node_2)+0.0015*abs(N(e))) (x(node_1)+0.0015*abs(N(e)))];
C=[y(node_1) y(node_2)  y(node_2) y(node_1)];
hpatch=patch(B,C,'r','Facealpha',0.2);
axis off;
hold on;
elseif rad2deg(theta(e))==0 && N(e)>0;
B=[x(node_1) x(node_2) x(node_2) x(node_1)];
C=[y(node_1) y(node_2)  (y(node_2)-0.003*abs(N(e))) (y(node_1)-0.003*abs(N(e)))];
patch(B,C,'r','Facealpha',0.2);
axis off;
hold on;
elseif rad2deg(theta(e))==0 && N(e)<0;
B=[x(node_1) x(node_2) x(node_2) x(node_1)];
C=[y(node_1) y(node_2)  (y(node_2)+0.003*abs(N(e))) (y(node_1)+0.003*abs(N(e)))];
patch(B,C,'b','Facealpha',0.2);
axis off;
hold on;
elseif rad2deg(theta(e))==90 && N(e)<0;
B=[x(node_1) x(node_2) (x(node_2)-0.0015*abs(N(e))) (x(node_1)-0.0015*abs(N(e)))];
C=[y(node_1) y(node_2)  y(node_2) y(node_1)];
hpatch=patch(B,C,'b','Facealpha',0.2);
rotate(hpatch,[1 0 0],theta(e));
axis equal
hold on
elseif rad2deg(theta(e))<65 ;
B=[x(node_1) x(node_2) x(node_2)+0.0015*(abs(N(e))*sin(theta(e)))  x(node_1)+0.0015*(abs(N(e))*sin(theta(e)))];
C=[y(node_1) y(node_2)  y(node_2)-0.0015*(abs(N(e))*cos(theta(e))) y(node_1)-0.0015*(abs(N(e))*cos(theta(e)))];
hpatch=patch(B,C,'r','Facealpha',0.2);
rotate(hpatch,[1 0 0],theta(e));
axis equal;
hold on;
elseif rad2deg(theta(e))>90 ; 
B=[x(node_1) x(node_2) x(node_2)-0.0015*(abs(N(e))*cos(theta(e)-pi/2))  x(node_1)-0.0015*(abs(N(e))*cos(theta(e)-pi/2))];
C=[y(node_1) y(node_2)  y(node_2)-0.0015*(abs(N(e))*sin(theta(e)-pi/2)) y(node_1)-0.0015*(abs(N(e))*sin(theta(e)-pi/2))];
hpatch=patch(B,C,'r','Facealpha',0.2);
rotate(hpatch,[1 0 0],theta(e));
axis equal;
axis([-7 14 -0.1 13.4]);
title(['\fontname{Courier}\fontsize{13}Plot Diagram Normal Forces'],'color','K');
hold on;
axis off;
end
text((x(node_1)+x(node_2))/2,(y(node_1)+y(node_2))/2,{e},'HorizontalAlignment','center','Color','k','FontSize',6);
end
%% Value Table N
i=1:nel;
T1=table(i',N','VariableNames',{'element','N'});
% Get the table in string form
TString=evalc('disp(T1)');
%Use TeX Markup for bold formatting and underscores.
TString=strrep(TString,'<strong>','\bf');
TString=strrep(TString,'</strong>','\rm');
TString=strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth=get(0,'FixedWidthFontName');
% Output the table using the annotation command
annotation(gcf,'TextBox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0.55 0.05 0.18 0.61]);
hold on;
axis off

%% Plot Bending Moment and Share
figure(4);
% Bending moment beam 14
pos3 = [0.04 0.56 0.2 0.3];
subplot('Position',pos3);
s_14_plot=0:L(14)/20:L(14);
area(s_14_plot,M14,'FaceColor','b','LineWidth',2,'FaceAlpha',0.4);
axis([0 max(s_14_plot) 0 max(M17)+10]);
title(['\fontname{Courier}\fontsize{8}Bending moment beam 14'],'color','b');
xlabel('s [m]');
ylabel('M [KN*m]');
axis off;
axis ij;
% Bending moment beam 17
pos4 = [0.24 0.56 0.4 0.3];
subplot('Position',pos4);
s_17_plot=0:L(17)/20:L(17);
area(s_17_plot,M17,'FaceColor','b','LineWidth',2,'FaceAlpha',0.4);
axis([0 max(s_17_plot) 0 max(M17)+10]);
title(['\fontname{Courier}\fontsize{8}Bending moment beam 17'],'color','b');
xlabel('s [m]');
ylabel('M [KN*m]');
axis off;
axis ij;
% Shear force beam 14
pos5 = [0.04 0.15 0.2 0.3];
subplot('Position',pos5);
s_14_plot=0:L(14)/20:L(14);
area(s_14_plot,T14,'FaceColor','r','LineWidth',2,'FaceAlpha',0.4);
axis([0 max(s_14_plot) -max(T17)-10 max(T17)+10]);
title(['\fontname{Courier}\fontsize{8}Shear force beam 14'],'color','r');
xlabel('s [m]');
ylabel('T [KN]');
axis off;
hold on;
% Shear force beam 17
pos6 = [0.24 0.15 0.4 0.3];
subplot('Position',pos6);
s_17_plot=0:L(17)/20:L(17);
area(s_17_plot,T17,'FaceColor','r','LineWidth',2,'FaceAlpha',0.4);
axis([0 max(s_17_plot) -max(T17)-10  max(T17)+10]);
title(['\fontname{Courier}\fontsize{8}Shear force beam 17'],'color','r');
xlabel('s [m]');
ylabel('T[KN]');
axis off;
hold on;
%% Value Table M & T 
T2=table(s_14_plot',M14',T14',s_17_plot',M17',T17','VariableNames',{'s14','M14','T14','s17','M17','T17'});
% Get the table in string form
TString=evalc('disp(T2)');
%Use TeX Markup for bold formatting and underscores.
TString=strrep(TString,'<strong>','\bf');
TString=strrep(TString,'</strong>','\rm');
TString=strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth=get(0,'FixedWidthFontName');
% Output the table using the annotation command
annotation(gcf,'TextBox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0.65 0.2 0.34 0.63]);
axis off