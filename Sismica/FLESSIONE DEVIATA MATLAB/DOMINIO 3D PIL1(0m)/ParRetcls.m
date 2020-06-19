% [ecls,scls] = ParRetcls(fcd,e2,eu)
%
% Questa funzione calcola i punti della curva tensione-deformazione
% (parabola-rettangolo) del cls.
%
% INPUT:
% fcd   valore tensione resistenza a compressione [MPa]
% e2    valore deformazione tratto parabolico
% eu    valore deformazione ultima
%
% OUTPUT:
% ecls  vettore deformazioni
% scls  vettore tensioni


function [ecls,scls] = ParRetcls(fcd,e2,eu,sezname)

e = [0:0.01:e2];

for i=1:length(e)

scls(i) = fcd * (2*e(i)/e2-(e(i)/e2)^2);

end

ecls = [e eu];
scls = [scls fcd];

figure(100)
hold on
title('Diagramma TENSIONE-DEFORMAZIONE CLS')
xlabel('Deformazione \epsilon [%]')
ylabel('Tensione \sigma  [MPa]')
plot(ecls,scls,'r')

fcd = num2str(fcd);
e2 = num2str(e2);
eu = num2str(eu);
text(0.2,5,{['f_{cd} = ' fcd];['e_2 = ' e2];['e_u = ' eu]})



% salvataggio
filename = ['./' sezname '/cls'];
saveas(gcf,filename,'fig');
saveas(gcf,filename,'tif');