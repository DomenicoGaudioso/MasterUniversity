% [Nn,Mnx,Mny] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc)
%
% Questa funzione calcola le coppie N-Mrd.
%
% INPUT:
% Xg          coordinata x del baricentro
% Yg          coordinata y del baricentro
% XSez        matrice coordinate x aree elementari configurazione iniziale
% YSez        matrice coordinate y aree elementari configurazione iniziale
% XYAArm      matrice delle coordinate x-y/aree barre di armatura
% DefSez      matrice delle deformazioni aree elementari
% DefArm      matrice delle deformazioni barre di armatura
% ecls        vettore delle deformazioni digramma tensione-deformazione cls
% scls        vettore delle tensioni digramma tensione-deformazione cls
% eacc        vettore delle deformazioni digramma tensione-deformazione acciaio
% sacc        vettore delle tensioni digramma tensione-deformazione acciaio
%
% OUTPUT:
% Nn          Sforzo normale agente nella sezione [kN]
% Mnx         Momento resistente x della sezione [kNm]
% Mny         Momento resistente y della sezione [kNm]

function [Nn,Mnx,Mny] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);


% Calcolo tensioni aree sezioni elementari cls
TenSez = interp1(ecls,scls,DefSez,'linear',0);

% Calcolo vettore tensioni barre d'armatura
for i=1:length(DefArm)
TenArm(i) = XYAArm(i,3)*interp1(eacc,sacc,DefArm(i),'linear',0);
end

% Calcolo sforzo normale [kN]
Nn = (sum(TenArm) + sum(sum(TenSez)))/1000;

% Calcolo momento resistente ultimo Mnx [kNm]
[r c] = size(YSez);

for i=1:c
	for j=1:r
	Mnxcls(j,i) = TenSez(j,i) * (YSez(j,i)-Yg);
	end
end


for i=1:length(TenArm)
Mnxacc(i) = TenArm(i) * (XYAArm(i,2)-Yg);
end

Mnx = (sum(sum(Mnxcls)) + sum(Mnxacc))/10^6;

% Calcolo momento resistente ultimo Mny [kNm]
[r c] = size(XSez);

for i=1:c
	for j=1:r
	Mnycls(j,i) = TenSez(j,i) * (XSez(j,i)-Xg);
	end
end


for i=1:length(TenArm)
Mnyacc(i) = TenArm(i) * (XYAArm(i,1)-Xg);
end

Mny = (sum(sum(Mnycls)) + sum(Mnyacc))/10^6;