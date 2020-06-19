%% Linea elastica 
Metodo_linea_elastica;
subplot(2,2,[1 2])
s_plot = linspace(0,L,50)
vex_plot=subs(vex,s,s_plot)
plot(s_plot,vex_plot,'g','LineWidth',1.5,'DisplayName','v_e_x_a_c_t')
hold on
legend
%% CDS
% per definire il grafico di una retta
y=s_plot*0
subplot(2,2,4)
plot(s_plot,y,'K','DisplayName','asse_t_r_a_v_e')
hold on

subplot(2,2,3)
plot(s_plot,y,'K','DisplayName','asse_t_r_a_v_e')
hold on

%Bending moment
subplot(2,2,4)
plot(s_plot,-M(s_plot),'LineWidth',2,'color',[0 0.25 0.25],'DisplayName','M_e_x_a_c_t')
axis([0 max(s_plot) -500 500])
hold on

%Shear force
subplot(2,2,3)
plot(s_plot,-T(s_plot),'LineWidth',2,'color','r','DisplayName','T_e_x_a_c_t')
axis([0 max(s_plot) -300 200])
hold on
%% Rigid-segment discretisation 1
n=1:4
Metodo_discretizzazione
%% Plot deformed shapes
subplot(2,2,[1 2])
plot(s_plot,vrig,'--','LineWidth',1,'DisplayName','4 rigid segments')
axis([0 max(s_plot) 0 0.06])
title('deformed shapes ')
xlabel('s  [m]')
ylabel('v    [m]')
axis ij


subplot(2,2,4)
plot(s_plot,M,'DisplayName','4 rigid segments')


subplot(2,2,3)
plot(s_plot,-T_discr,'DisplayName','4 rigid segments')
%% Rigid-segment discretisation 2
n=1:8
Metodo_discretizzazione
%% Plot deformed shapes
subplot(2,2,[1 2])
plot(s_plot,vrig,'--','LineWidth',1,'DisplayName','8 rigid segments')
axis([0 max(s_plot) 0 0.06])
title('deformed shapes')
xlabel('s  [m]')
ylabel('v    [m]')
axis ij
hold on
legend

subplot(2,2,4)
plot(s_plot,M,'DisplayName','8 rigid segments')

subplot(2,2,3)
plot(s_plot,-T_discr,'DisplayName','8 rigid segments')
axis ij 
%% Rigid-segment discretisation 3
n=1:16
Metodo_discretizzazione
%% Plot deformed shapes
subplot(2,2,[1 2])
plot(s_plot,vrig,'--','LineWidth',1,'DisplayName','16 rigid segments')
axis([0 max(s_plot) 0 0.06])
title('deformed shapes ')
xlabel('s  [m]')
ylabel('v    [m]')
axis ij
hold on
legend

subplot(2,2,4)
plot(s_plot,M,'DisplayName','16 rigid segments')
hold on
subplot(2,2,3)
plot(s_plot,-T_discr,'DisplayName','16 rigid segments')
alpha(.1)
hold on
%% Rigid-segment discretisation 4
n=1:64
Metodo_discretizzazione
%% Plot deformed shapes
subplot(2,2,[1 2])
plot(s_plot,vrig,'--','LineWidth',1,'color','k','DisplayName','64 rigid segments')
axis([0 max(s_plot) 0 0.06])
title('deformed shapes ')
xlabel('s  [m]')
ylabel('v    [m]')
axis ij
hold on
legend

subplot(2,2,4)
plot(s_plot,M,'DisplayName','64 rigid segments')
axis([0 max(s_plot) -300 450])
title('Bending moment')
xlabel('s  [m]')
ylabel('M    [KN*m]')
hold on
legend

subplot(2,2,3)
plot(s_plot,-T_discr,'LineWidth',2,'color','k','DisplayName','64 rigid segments')
axis([0 max(s_plot) -200 300])
title('Shear')
xlabel('s  [m]')
ylabel('T    [KN]')
hold on
legend
%% Differenze finite 1
i=1:12;    % Dominio
Metodo_diff_finite
% Graphics  Deflection
subplot(2,2,[1 2])
plot(s_diff,v_finito,'-*','LineWidth',1,'DisplayName','diff. finite con dominio di 12 tratti')
hold on

% Momento
subplot(2,2,4)
plot(s_diff,M_diff,'-*','LineWidth',1,'DisplayName','diff. finite con 12 tratti')
hold on

% Taglio
subplot(2,2,3)
plot(s_diff,T_diff,'-*','LineWidth',1,'DisplayName','diff. finite con 12 tratti')
hold on
%% Differenze finite 2
i=1:24;    % Dominio
Metodo_diff_finite
% Graphics  Deflection
subplot(2,2,[1 2])
plot(s_diff,v_finito,'-*','LineWidth',1,'DisplayName','diff. finite con dominio di 24 tratti')
hold on

% Momento
subplot(2,2,4)
plot(s_diff,M_diff,'-*','LineWidth',1,'DisplayName','diff. finite con 24 tratti')
hold on

% Taglio
subplot(2,2,3)
plot(s_diff,T_diff,'-*','LineWidth',1,'DisplayName','diff. finite con 24 tratti')
hold on
%% Differenze finite 3
i=1:48;    % Dominio
Metodo_diff_finite
% Graphics  Deflection
subplot(2,2,[1 2])
plot(s_diff,v_finito,'-*','LineWidth',1,'DisplayName','diff. finite con dominio di 48 tratti')
axis([0 max(s_diff) 0 0.06])
title('Deflection')
xlabel('s  [m]')
ylabel('v  [m]')
axis ij
hold on
legend

% Momento
subplot(2,2,4)
plot(s_diff,M_diff,'-*','LineWidth',1,'DisplayName','diff. finite con 48 tratti')
axis([0 max(s_diff) -500 500])
title('Bending Moment')
xlabel('s  [m]')
ylabel('M    [KN*m]')
hold on
legend

% Taglio
subplot(2,2,3)
plot(s_diff,T_diff,'-*','LineWidth',1,'DisplayName','diff. finite con 48 tratti')
axis([0 max(s_diff) -200 300])
title('Shear force')
xlabel('s  [m]')
ylabel('T    [KN]')
hold on

legend