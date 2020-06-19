%%  Computational Mechanics - TERZA ESERCITAZIONE
% risoluzione di elementi quadratici a 9 nodi [Funzioni di forma
% 1uadratiche]
%% Clear workspace and close any open windows
clear all;
close all;
%% Shape functions
[Xi, Etha] = meshgrid (-1: 0.05: 1);

l0_2Xi=(Xi-1).*Xi
l0_2Etha=(Etha-1).*Etha
l1_2Xi=(Xi+1).*(1-Xi)
l1_2Etha=(Etha+1).*(1-Etha)
l2_2Xi=(1/2).*(Xi+1).*Xi
l2_2Etha=(1/2).*(Etha+1).*Etha

figure(1)
N1=l0_2Xi.*l0_2Etha;
s=surfc(Xi,Etha,N1,'FaceAlpha',0.8)
colormap(jet)
% s.EdgeColor = 'interp'
axis vis3d
zlabel('N1(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N1(\zeta,\eta)]=l^{2}_0(\zeta)*l^{2}_0(\eta)','LineWidth',20)

figure(2)
N2=l2_2Xi.*l0_2Etha;
s=surfc(Xi,Etha,N2,'FaceAlpha',0.8)
colormap(jet)
axis vis3d
zlabel('N2(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N2(\zeta,\eta)]=l^{2}_2(\zeta)*l^{2}_0(\eta)','LineWidth',20)

figure(3)
N3=l2_2Xi.*l2_2Etha;
s=surfc(Xi,Etha,N3,'FaceAlpha',0.8)
colormap(jet)
axis vis3d
zlabel('N3(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N3(\zeta,\eta)]=l^{2}_2(\zeta)*l^{2}_2(\eta)','LineWidth',20)

figure(4)
N4=l0_2Xi.*l2_2Etha;
s=surfc(Xi,Etha,N4,'FaceAlpha',0.8)
colormap(jet)
axis vis3d
zlabel('N3(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N4(\zeta,\eta)]=l^{2}_0(\zeta)*l^{2}_2(\eta)','LineWidth',20)

figure(5)
N5=l1_2Xi.*l0_2Etha;
s=surfc(Xi,Etha,N5,'FaceAlpha',0.8)
colormap(jet)
axis vis3d
zlabel('N4(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N5(\zeta,\eta)]=l^{2}_1(\zeta)*l^{2}_0(\eta)','LineWidth',20)

figure(6)
N6=l2_2Xi.*l1_2Etha;
s=surfc(Xi,Etha,N6,'FaceAlpha',0.8)
colormap(jet)
axis vis3d
zlabel('N4(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N6(\zeta,\eta)]=l^{2}_2(\zeta)*l^{2}_1(\eta)','LineWidth',20)

figure(7)
N7=l1_2Xi.*l2_2Etha;
s=surfc(Xi,Etha,N7,'FaceAlpha',0.8)
colormap(jet)
axis vis3d
zlabel('N4(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N7(\zeta,\eta)]=l^{2}_1(\zeta)*l^{2}_2(\eta)','LineWidth',20)

figure(8)
N8=l0_2Xi.*l1_2Etha;
s=surfc(Xi,Etha,N8,'FaceAlpha',0.8)
colormap(jet)
axis vis3d
zlabel('N4(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N8(\zeta,\eta)]=l^{2}_0(\zeta)*l^{2}_1(\eta)','LineWidth',20)

figure(9)
N9=l1_2Xi.*l1_2Etha;
s=surfc(Xi,Etha,N9,'FaceAlpha',0.8)
colormap(jet)
axis vis3d
zlabel('N4(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N8(\zeta,\eta)]=l^{2}_1(\zeta)*l^{2}_1(\eta)','LineWidth',20)
