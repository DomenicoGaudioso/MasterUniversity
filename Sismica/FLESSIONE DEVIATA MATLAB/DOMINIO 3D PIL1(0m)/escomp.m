% [esc] = escomp (YSezrifRot,XYAArmRot,eu);
%
% Questa funzione calcola la deformazione di compressione della
% barra d'armatura di riferimento (barra inferiore).
%
% INPUT:
% YSezrifRot  matrice delle coordinate y sezione di riferimento
% XYAArmRot   matrice delle coordinate e delle aree delle barre di armatura ruotate
% eu          valore di deformazione ultima cls 
%
% OUTPUT:
% esc         valore deformazione di compressione barra armatura di riferimento

function [esc] = escomp (YSezrifRot,XYAArmRot,eu);


ymax = max(max(YSezrifRot));

ys = min(XYAArmRot(:,2));

esc = eu*ys/ymax;