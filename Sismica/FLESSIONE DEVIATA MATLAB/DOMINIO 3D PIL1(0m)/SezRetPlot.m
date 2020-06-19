% [] = SezRetPlot(b,h,XYAArm)

function [] = SezRetPlot(b,h,XYAArm,sezname)

% vettori coordinate perimetro
x = [0 b b 0 0];
y = [0 0 h h 0];

figure(1)
hold on
xlabel({'Base  [mm]';'Area barre  [mm^2]'})
ylabel('Altezza  [mm]')
axis equal
% plottaggio perimetro
axis ([-10 b+10 -10 h+10])
title({['SEZIONE "' sezname '"'];'';'M_x ^+ , M_y ^+ : antiorari'})
% plottaggio armature
plot(x,y,'k')
plot(XYAArm(:,1),XYAArm(:,2),'ob')
% plottaggio area armature
[r,c] = size(XYAArm);
for i=1:r
area = num2str(XYAArm(i,3));
text(XYAArm(i,1)+3.5,XYAArm(i,2)+3.5,area)
end
% salvataggio
filename = ['./' sezname '/SEZIONE_' sezname];
saveas(gcf,filename,'fig');
saveas(gcf,filename,'jpg');