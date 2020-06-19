% FUNCTION: SezRetRot
%
% [XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSezRot,YSezRot,L,ANGRot] = SezRetRot(b,h,angle)
%
% La funzione suddivide la sezione rettangolare in aree elementari (quadrate)
% e ne calcola le coordinate rispetto ad un sistema di  riferimento  centrato
% sul vertice inferiore sisnistro.
% La funzione calcola le coordinate delle aree elementari (quadrate) per  una
% data rotazione rispetto ad un sistema di  riferimento  centrato sul vertice
% inferiore sisnistro.
%
% INPUT:
% b          base sezione [mm]
% h          altezza sezione [mm]
% angle      angolo di rotazione in gradi sessa-decimali
%
% OUTPUT:
% XSezrifRot  matrice delle coordinate x
% YSezrifRot  matrice delle coordinate y
% Lrif        matrice delle distanze dei vertici dall'origine
% ANGrifRot   matrice degli angoli che formano le distanze L(i,j) dei vertici dall'origine rispetto all'asse delle x
% XSezRot     matrice delle coordinate x
% YSezRot     matrice delle coordinate y
% L           matrice delle distanze dall'origine
% ANGRot      matrice degli angoli che formano le distanze L(i,j) dall'origine rispetto all'asse delle x

function [XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSezRot,YSezRot,L,ANGRot] = SezRetRot(b,h,angle)


% Conversione angolo in radianti
angle = angle*pi/180;

% determinazione coordinate baricentriche aree elementari
for i=1:b
XSez(i)=0.5+(i-1);
end

for i=1:h
YSez(i)=0.5+(i-1);
end

% Calcolo matrice distanze
for i=1:b

	for j=1:h
	
	L(j,i)=(XSez(i)^2+YSez(j)^2)^0.5;
	
	end
	
end

% Calcolo matrice angoli
for i=1:b

	for j=1:h
	
	ANGRot(j,i) = angle+atan(YSez(j)/XSez(i));
	
	end
	
end

% Calcolo matrici coordinate ruotate
for i=1:b

	for j=1:h
	
	XSezRot(j,i) = L(j,i)*cos(ANGRot(j,i));
	YSezRot(j,i) = L(j,i)*sin(ANGRot(j,i));
	
	end
	
end


% calcolo di Lrif
Lrif = [ 0 b ; h (b^2+h^2)^0.5];

% calcolo di ANGrif
ANGrifRot = [angle angle ; angle+pi/2 angle+atan(h/b)];

% calcolo matrice coordinate ruotate sezione di riferimento
for i=1:2

	for j=1:2
	
	XSezrifRot(j,i) = Lrif(j,i)*cos(ANGrifRot(j,i));
	YSezrifRot(j,i) = Lrif(j,i)*sin(ANGrifRot(j,i));
	
	end
	
end