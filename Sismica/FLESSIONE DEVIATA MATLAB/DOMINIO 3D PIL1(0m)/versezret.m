% [] = versezret (Nval,Mxval,Myval,NMxMy)
%
% Questa funzione esegue la verifica a  'presso-flessione'  deviata della
% sezione, ricavando il dominio di resistenza Mx-My, riferito allo sforzo
% normale agente, per interpolazione lineare.
%
% INPUT:
% Nval   vettore degli sforzi normali per i quali è stato calcolato il dominio 3D
% Mxval  matrice dei momenti resistenti Mx riferiti ai diversi valori Nval
% Myval  matrice dei momenti resistenti My riferiti ai diversi valori Nval
% NMxMy  matrice delle triplette N-Mx-My sollecitanti
%
% OUTPUT:
% VER    matrice risultati verifiche

function [VER] = versezret (folder,Nval,Mxval,Myval,NMxMy)

% eliminazione colonne dati doppie (limite tra i vari cicli)
[rr cc] = size(Mxval);
k = cc/4;
Mxval = [Mxval(:,1:k-1) Mxval(:,k+1:2*k-1) Mxval(:,2*k+1:3*k-1) Mxval(:,3*k+1:4*k-1)];
Myval = [Myval(:,1:k-1) Myval(:,k+1:2*k-1) Myval(:,2*k+1:3*k-1) Myval(:,3*k+1:4*k-1)];

% calcolo angolo momento totale sollecitante
[r c] = size(NMxMy);
for i=1:r
	if NMxMy(i,2)>0 & NMxMy(i,3)==0
	ANGMxMy(i) = 0;
	elseif NMxMy(i,2)==0 & NMxMy(i,3)>0
	ANGMxMy(i) = pi/2;
	elseif NMxMy(i,2)<0 & NMxMy(i,3)==0
	ANGMxMy(i) = pi;
	elseif NMxMy(i,2)==0 & NMxMy(i,3)<0
	ANGMxMy(i) = (3/2)*pi;
	elseif NMxMy(i,2)>0 & NMxMy(i,3)>0
	ANGMxMy(i) = atan(NMxMy(i,3)/NMxMy(i,2));
	elseif NMxMy(i,2)<0 & NMxMy(i,3)>0
	ANGMxMy(i) = pi+atan(NMxMy(i,3)/NMxMy(i,2));
	elseif NMxMy(i,2)<0 & NMxMy(i,3)<0
	ANGMxMy(i) = pi+atan(NMxMy(i,3)/NMxMy(i,2));
	elseif NMxMy(i,2)>0 & NMxMy(i,3)<0
	ANGMxMy(i) = 2*pi+atan(NMxMy(i,3)/NMxMy(i,2));
	end
end

% interpolazione dominio
[rr cc] = size(Mxval);
for i=1:r
% determinazione dominio riferito allo sforzo sollecitante NMxMy(i,1)
	for j=1:cc
		Mxint(j) = interp1(Nval,Mxval(:,j),NMxMy(i,1),'linear');
		Myint(j) = interp1(Nval,Myval(:,j),NMxMy(i,1),'linear');
% calcolo angolo momento totale resistente
		if Mxint(j)>0 & Myint(j)==0
		ANGMtotint(j) = 0;
		elseif Mxint(j)==0 & Myint(j)>0
		ANGMtotint(j) = pi/2;
		elseif Mxint(j)<0 & Myint(j)==0
		ANGMtotint(j) = pi;
		elseif Mxint(j)==0 & Myint(j)<0
		ANGMtotint(j) = (3/2)*pi;
		elseif Mxint(j)>0 & Myint(j)>0
		ANGMtotint(j) = atan(Myint(j)/Mxint(j));
		elseif Mxint(j)<0 & Myint(j)>0
		ANGMtotint(j) = pi+atan(Myint(j)/Mxint(j));
		elseif Mxint(j)<0 & Myint(j)<0
		ANGMtotint(j) = pi+atan(Myint(j)/Mxint(j));
		elseif Mxint(j)>0 & Myint(j)<0
		ANGMtotint(j) = 2*pi+atan(Myint(j)/Mxint(j));
		end
	end

% Ordinamento vettore angoli dominio di resistenza
	ANGMtotintrot(1) = 0;
	if Myint(1)>=0
	alpha = ANGMtotint(1);
		for j=2:length(ANGMtotint)
		ANGMtotintrot(j) = ANGMtotint(j)-alpha;
		end
	ANGMxMyrot(i) = ANGMxMy(i)-alpha;
	else
	alpha = 2*pi-ANGMtotint(1);
		for j=2:length(ANGMtotint)
		ANGMtotintrot(j) = ANGMtotint(j)+alpha;
		end
	ANGMxMyrot = ANGMxMy(i)+alpha;
	end

% calcolo intervallo intersezione dominio
% ciclo per gli n-1 intervalli
	for j=1:length(ANGMtotintrot)-1
		if ANGMxMyrot>=ANGMtotintrot(j) & ANGMxMyrot<=ANGMtotintrot(j+1)
		Mx1 = Mxint(j);
		My1 = Myint(j);
		Mx2 = Mxint(j+1);
		My2 = Myint(j+1);
		end
	end
% condizione per l'ultimo intervallo	
	if ANGMxMyrot>=ANGMtotintrot(j+1)
	Mx1 = Mxint(1);
	My1 = Myint(1);
	Mx2 = Mxint(j+1);
	My2 = Myint(j+1);
	end
	
% calcolo punto intersezione dominio
% calcolo parametri retta dominio
	q = interp1([Mx1 Mx2],[My1 My2],0,'linear','extrap');
	m = (My1-My2)/(Mx1-Mx2);
% calcolo parametri retta momento sollecitante
	ms = (NMxMy(i,3))/(NMxMy(i,2));
% calcolo momento resistente MxRd
	MxRd(i) = q/(ms-m);
% calcolo momento resistente MyRd
	MyRd(i) = ms*MxRd(i);
% calcolo coeff. di sicurezza
	r(i) = (MxRd(i)^2+MyRd(i)^2)^0.5/(NMxMy(i,2)^2+NMxMy(i,3)^2)^0.5;
% controllo verifica
	if r(i)>=1
	esito = 'OK';
	else
	esito = 'NO';
	end

	
% plottaggio figura	
	Ns = num2str(NMxMy(i,1));
	Mxs = num2str(NMxMy(i,2));
	Mys = num2str(NMxMy(i,3));
	Mxres = num2str(MxRd(i));
	Myres = num2str(MyRd(i));
	rs = num2str(r(i));
	sezione = num2str(i);
	figure(i+10000)
	hold on
	grid on
	title({['Dominio di resistenza  N = ' Ns ' [kN]'] ['Mx = ' Mxs ' [kNm]   ' 'My = ' Mys ' [kNm]'] ['MxRd = ' Mxres ' [kNm]  ' 'MyRd = ' Myres ' [kNm]'] ['\rho= ' rs]})
	xlabel('Momento flettente Mx  [kNm]')
	ylabel('Momento flettente My  [kNm]')
% plottaggio perimetro dominio
	plot([Mxint Mxint(1)],[Myint Myint(1)],'-b')
% plottaggio momenti resistenti
	plot([0 MxRd(i)],[0 MyRd(i)],'-k')
	plot(MxRd(i),MyRd(i),'ko')
% plottaggio momenti sollecitanti
    sz = 40;
	scatter(NMxMy(i,2),NMxMy(i,3),sz,'MarkerEdgeColor','k',...
              'MarkerFaceColor','r',...
              'LineWidth',1.5)
    
% salvataggio figura
	filename = [folder '/SEZ_' sezione ' ' esito];
	saveas(gcf,filename,'fig');
	saveas(gcf,filename,'jpg');
end

% matrice risultati
VER = [NMxMy MxRd' MyRd' r'];