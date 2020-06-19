% FUNCTION: SezRetDef
%
% [DefSez] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec,es)
%
% La funzione calcola la deformazione delle aree elementari (quadrate)
% e delle barre d'armatura,  partendo dai  valori  di  deformazione di
% riferimento.
%
%
% INPUT:
% YSezRot  matrice delle coordinate y
% ec       deformazione estremo superiore
% es       deformazione armatura di riferimento
%
% OUTPUT:
% DefSez  Matrice delle deformazioni delle aree elementari
% DefArm  vettore delle deformazioni delle barre d'armatura


function [DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec,es)


% massimo valore coordinata y
yc = max(max(YSezrifRot));

% valore coordinata barra acciaio di riferimento
ys = min(XYAArmRot(:,2));

y = [ys yc];
e = [es ec];

% calcolo delle deformazioni
DefSez = interp1(y,e,YSezRot,'linear','extrap');

DefArm = interp1(y,e,XYAArmRot(:,2),'linear','extrap');