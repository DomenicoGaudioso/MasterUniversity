%%  Computational Mechanics - TERZA ESERCITAZIONE
% risoluzione di elementi triangolari a 6 nodi [Funzioni di forma
% quadratiche]
%% Clear workspace and close any open windows
clear all;
close all;
%% Shape functions
figure(1)
[Xi, Etha] = meshgrid (0: 0.05: 1);
lamda=1-Xi-Etha;
NI=lamda.*(2.*lamda-1);
NI (Xi + Etha> 1) = nan;
tri = delaunay(Xi,Etha);
trisurf(tri,Xi,Etha,NI,'FaceAlpha',0.8)
axis vis3d
zlabel('NI(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [NI(\zeta,\eta)]=\lambda*(2*\lambda-1)','LineWidth',20)

figure(2)
NJ=Xi.*(2.*Xi-1); %mettere ilpunto prima della moltiplicazione sennò non fa la superficie per bene
NJ(Xi + Etha> 1) = nan;
trisurf(tri,Xi,Etha,NJ,'FaceAlpha',0.8)
axis vis3d
zlabel('NJ(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [NJ(\zeta,\eta)]=\zeta*(2*\zeta-1)','LineWidth',20)

figure(3)
NK=Etha.*(2.*Etha-1);
NK(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,NK,'FaceAlpha',0.8)
axis vis3d
zlabel('NK(\zeta,\eta)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [NK(\zeta,\eta)]=\eta*(2*\eta-1)','LineWidth',20)
% imagesc(NI)
% contourf(NI)
% colormap(jet)

figure(4)
NL=4.*Etha.*lamda;
NL(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,NL,'FaceAlpha',0.8)
axis vis3d
zlabel('NK(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [NL(\zeta,\eta)]=4*\eta*\lambda','LineWidth',20)

figure(5)
NM=4.*Xi.*lamda;
NM(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,NM,'FaceAlpha',0.8)
axis vis3d
zlabel('NM(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [NM(\zeta,\eta)]=4*\zeta*\lambda','LineWidth',20)

figure(6)
NN=4.*Xi.*Etha;
NN(Xi+Etha> 1) = nan;
trisurf(tri,Xi,Etha,NM,'FaceAlpha',0.8)
axis vis3d
zlabel('NM(Xi,Etha)','LineWidth',2)
xlabel('\zeta','LineWidth',2)
ylabel('\eta','LineWidth',2)
title(' Shape functions [NN(\zeta,\eta)]=4*\zeta*\eta','LineWidth',20)