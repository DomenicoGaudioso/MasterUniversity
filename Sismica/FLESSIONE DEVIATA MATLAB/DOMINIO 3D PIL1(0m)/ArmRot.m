% FUNCTION: ArmRot
%
% [XYAArmRot,LArm,ANGArmRot] = ArmRot(XYAArm,angle)
%
% La funzione calcola la posizione delle barre d'armatura rispetto all'origine
% assegnata la rotazione 'angle'.
%
% INPUT:
% XYAArm     matrice delle coordinate e dalle aree delle barre di armatura [X1 Y1 A1 ; .. .. .. ; Xn Yn An] [mm][mm][mm^2]
% angle      angolo di rotazione in gradi sessa-decimali
%
% OUTPUT:
% XYAArmRot  matrice delle coordinate e delle aree delle barre di armatura ruotate
% LArm       vettore delle distanze delle barre di armatura dall'origine
% ANGArmRot  vettore degli angoli che formano le distanze LArmRot(i,j) dei vertici dall'origine rispetto all'asse delle x


function [XYAArmRot,LArm,ANGArmRot] = ArmRot(XYAArm,angle)

% Conversione angolo in radianti
angle = angle*pi/180;

[r c] = size(XYAArm);

% calcolo distanze LArm dall'origine
for i=1:r
LArm(i) = (XYAArm(i,1)^2+XYAArm(i,2)^2)^0.5;
end

% calcolo angoli di riferimento
for i=1:r
ANGArmRot(i) = angle+atan(XYAArm(i,2)/XYAArm(i,1));
end

% calcolo delle coordinate ruotate
for i=1:r
XYAArmRot(i,1) = LArm(i)*cos(ANGArmRot(i));
XYAArmRot(i,2) = LArm(i)*sin(ANGArmRot(i));
XYAArmRot(i,3) = XYAArm(i,3);
end