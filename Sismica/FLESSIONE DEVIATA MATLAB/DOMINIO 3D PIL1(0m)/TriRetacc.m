% [eacc,sacc] = TriRetacc(fyd,ey,esu)
%
% Questa funzione calcola i punti della curva tensione-deformazione
% (triangolo-rettangolo) dell'acciaio.
%
% INPUT:
% fyd   valore tensione snervamento [MPa]
% ey    valore deformazione limite elastico
% esu   valore deformazione ultima
%
% OUTPUT:
% eacc  vettore deformazioni
% sacc  vettore tensioni


function [eacc,sacc] = TriRetacc(fyd,ey,esu,sezname)

eacc = [-esu -ey 0 ey esu];
sacc = [-fyd -fyd 0 fyd fyd];

figure(200)
hold on
title('Diagramma TENSIONE-DEFORMAZIONE ACCIAIO')
xlabel('Deformazione \epsilon [%]')
ylabel('Tensione \sigma  [MPa]')
plot(eacc(3:5),sacc(3:5),'r')

fyd = num2str(fyd);
ey = num2str(ey);
esu = num2str(esu);
text(5,100,{['f_{yd} = ' fyd];['e_y = ' ey];['e_{su} = ' esu]})

% salvataggio
filename = ['./' sezname '/acciaio'];
saveas(gcf,filename,'fig');
saveas(gcf,filename,'jpg');