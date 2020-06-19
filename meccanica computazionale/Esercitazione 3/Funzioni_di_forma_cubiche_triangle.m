%%  Computational Mechanics - TERZA ESERCITAZIONE
% risoluzione di elementi triangolari a 10 nodi [Funzioni di forma
% cubiche]
%% Clear workspace and close any open windows
clear all;
close all;
%% Shape functions
figure(1)
[Xi, Etha] = meshgrid (0: 0.05: 1);
lamda=1-Xi-Etha;
N1=(1/2).*lamda.*(3.*lamda-1).*(3.*lamda-2);
N1 (Xi + Etha> 1) = nan;
tri = delaunay(Xi,Etha);
trisurf(tri,Xi,Etha,N1,'FaceAlpha',0.8)
axis vis3d
zlabel('N1(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N1(\zeta,\eta)]=0.5*\lambda*(3*\lambda-1)(3*\lambda-2)','LineWidth',20)

figure(2)
N2=(1/2).*Xi.*(3.*Xi-1).*(3.*Xi-2); %mettere ilpunto prima della moltiplicazione sennò non fa la superficie per bene
N2(Xi + Etha> 1) = nan;
trisurf(tri,Xi,Etha,N2,'FaceAlpha',0.8)
axis vis3d
zlabel('N2(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N2(\zeta,\eta)]=0.5*\zeta*(3*\zeta-1)(3*\zeta-2)','LineWidth',20)

figure(3)
N3=(1/2).*Etha.*(3.*Etha-1).*(3.*Etha-2);
N3(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,N3,'FaceAlpha',0.8)
axis vis3d
zlabel('N3(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N3(\zeta,\eta)]=0.5*\eta*(3*\eta-1)(3*\eta-2)','LineWidth',20)
% imagesc(NI)
% contourf(NI)
% colormap(jet)
figure(4)
N4=(9/2).*lamda.*Xi.*(3*lamda-1);
N4(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,N4,'FaceAlpha',0.8)
axis vis3d
zlabel('N4(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N4(\zeta,\eta)]=(9/2)*\lambda*\zeta*(3*\lambda-1)','LineWidth',20)

figure(5)
N5=(9/2).*Xi.*Etha.*(3*Xi-1);
N5(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,N5,'FaceAlpha',0.8)
axis vis3d
zlabel('N5(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N5(\zeta,\eta)]=(9/2)*\zeta*\eta*(3*\zeta-1)','LineWidth',20)

figure(6)
N6=(9/2).*Etha.*lamda.*(3*Etha-1);
N6(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,N6,'FaceAlpha',0.8)
axis vis3d
zlabel('N6(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N6(\zeta,\eta)]=(9/2)*\zeta*\lambda*(3*\zeta-1)','LineWidth',20)

figure(7)
N7=(9/2).*lamda.*Xi.*(3*Xi-1);
N7(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,N7,'FaceAlpha',0.8)
axis vis3d
zlabel('N7(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N7(\zeta,\eta)]=(9/2)*\lambda*\zeta*(3*\zeta-1)','LineWidth',20)

figure(8)
N8=(9/2).*Xi.*Etha.*(3*Etha-1);
N8(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,N8,'FaceAlpha',0.8)
axis vis3d
zlabel('N8(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N8(\zeta,\eta)]=(9/2)*\zeta*\eta*(3*\eta-1)','LineWidth',20)

figure(9)
N9=(9/2).*Etha.*lamda.*(3*lamda-1);
N9(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,N9,'FaceAlpha',0.8)
axis vis3d
zlabel('N9(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N9(\zeta,\eta)]=(9/2)*\zeta*\lambda*(3*\lambda-1)','LineWidth',20)

figure(10)
N10=27.*Etha.*lamda.*Xi;
N10(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,N10,'FaceAlpha',0.8)
axis vis3d
zlabel('N10(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [N10(\zeta,\eta)]=27*\zeta*\lambda*\eta','LineWidth',20)