function [Nclean] = cleancurve (N)

% pulizia vettore N da valori duplicati
Nclean(1) = N(1);
% Controllo valore di trazione pura
for i=2:length(N)
	if N(1) == N(i)
	n = i;
	end
end

% Correzzione vettore degli sforzi normali
for i=2:n
Nclean(i) = N(1)+0.0000000001*i;
end

for i=n+1:length(N)
	if N(i-1) < N(i)
	Nclean(i) = N(i);
	elseif N(i-1) >= N(i)
	Nclean(i) = N(i-1)+0.0000000001;
	end
end