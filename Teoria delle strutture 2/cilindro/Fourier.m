%% Approssimazione Serie di Fourier

syms t
w=1;
for k=0:25
    X(k+1)=(8/pi^2)*(((-1)^k)*(sin((2*k+1)*w*t)/(2*k+1)^2));
end
for j=1:(length(X)-1)
    A(1) = X(1);
    A(j+1)=X(j+1)+A(j);
    t_plot = linspace(0,2*pi,100);
    AA =subs(A(j+1),t,t_plot);
    %plot(t_plot,AA)
    %hold on
end
As(t) = A(length(X));
%axis([0 2*pi -1 1]);