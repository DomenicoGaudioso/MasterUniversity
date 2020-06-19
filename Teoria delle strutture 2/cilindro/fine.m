%% Plot
clear all
close all

syms alpha r

%% Numerico
mu = 1;
H = 10;
R = 1;

%% Fase Elastica
Tauel(alpha,r) = (alpha*mu*r)/H;
Fourier;
alpha1 = vpa(As(pi/4));  %Fine Fase Elastica
c=40;
i = linspace(0 , alpha1, c);
rn=linspace(0,R,10);

for j = 1:c;
    figure()
    alphai(j) = vpa(As(i(j)));
%     t=-pi/2:0.01:pi/2; 
%     x=cos(t); 
%     y=sin(t); 
%     plot(x,y,'k-','LineWidth',2) 
%     hold on
    b = num2str(double(alphai(j)),3);
    area(rn,Tauel(i(j),rn),'FaceColor','b','DisplayName',b);
    %title('Andamento Tensioni \tau - Fase di Carico Elastico');
    xlabel('r');
    ylabel('\tau');
    
    axis([0 1 -0.05 0.05]);
    
    g=legend('show','Location','northwest')
    title(g,'\alpha :')
    n = num2str(j);
    print(n,'-djpeg')

end

%% Fase Plastica
syms  gammaps
alpha2 = vpa(As(pi/2));  %Fine Fase Plastica
Re(alpha) = (alpha1/alpha)*R;
gammap(alpha, r) = (-R*alpha1 + alpha*r)/(2*H);
Taup(alpha, r) = 2*mu*gammap
Tauel(alpha, r) = (alpha*mu*r)/H;
i = linspace(1, 2, c); 
for j = 1:c;
    alphai(j) = vpa(As((pi/4)*i(j)));
    Ri = Re(alphai(j));
    A = linspace(Ri, R, 10);
    %figure()
    %plot(A, gammap(alphai(j), A))
    B = linspace(0, Ri, 10);
    figure()
    b = num2str(double(alphai(j)),3);
    area(B, Tauel(alphai(j), B),'Facecolor', 'b','DisplayName',b)
    hold on
    area(A, Tauel(alphai(j), A)-Taup(alphai(j), A),'Facecolor', 'r','DisplayName',b)
    hold on
    %title('Andamento Tensioni \tau - Fase di Carico Plastico');
    xlabel('r');
    ylabel('\tau');
    axis([0 1 -0.05 0.05]);
    g=legend('show','Location','northwest')
    title(g,'\alpha :')
    n = num2str(c+j);
    print(n,'-djpeg')
end
%% Scarico Elastico
R2 = Re(alpha2);
i = linspace(2, 4, (40)); 
for j = 1:(c);  
        alphai(j) = vpa(As((pi/4)*i(j)));
        gammap2(r) = gammap(alpha2, r);
        Taup2(r) = 2*mu*gammap2;
        %Ri = Re(alphai(j));
        A = linspace(R2, R, 10);
        B = linspace(0, R2, 10);
        b = num2str(double(alphai(j)),3);
        figure()
        area(B, Tauel(alphai(j), B),'Facecolor', 'b', 'DisplayName',b)
        hold on
        area(A, (Tauel(alphai(j), A) - Taup2(A)),'Facecolor', 'r','DisplayName',b)
        hold on  
        %title('Andamento Tensioni \tau - Fase di Scarico Elastico');
        xlabel('r');
        ylabel('\tau');
        axis([0 1 -0.05 0.05])
    g=legend('show','Location','northwest')
    title(g,'\alpha :')
        n = num2str(2*c+j);
        print(n,'-djpeg')
end

%% Scarico Elastico pigreco-3mezzipigreco
%
% %% Scarico Plastico
% syms r2 gammap
% alpha4 = vpa(As(5*pi/4));  %Fine Fase Scarico Plastico
% gammaps(alpha, r2) = (R*alpha1 + alpha*r2)/(2*H);
% %Re(alpha) = (alpha1/alpha)*R;
% Rp(alpha,r1) = (2*H*gammap2/alpha)+((alpha1/alpha)*R);
% Taups(alpha, r2) = 2*mu*gammaps;
% i = linspace(5, 6, 11); 
% for j = 1:10;
%     alphai(j) = vpa(As((pi/4)*i(j+1)));
%     Ri = Re(alphai(j));
%     Rp2i = Rp(alphai(j),R);
%     A = linspace(0, Ri, 10);
%     B = linspace(Ri, Rp2i, 10);
%     C = linspace(Rp2i, R, 10);
%     %figure()
%     %plot(A, gammap(alphai(j), A))
%     figure()
%     b = num2str(double(alphai(j)),3);
%     area(A, Tauel(alphai(j), A),'Facecolor', 'b','DisplayName',b)
%     hold on
%     area(B, Tauel(alphai(j),B) - Taup2(B),'Facecolor', 'r','DisplayName',b)
%     hold on
%     area(C, Tauel(alphai(j),C) - Taup2(C) - Taups(alphai(j),C),'Facecolor', 'g','DisplayName',b)
%     hold on
%     title('Andamento Tensioni \tau - Fase di Carico Plastico');
%     xlabel('r');
%     ylabel('\tau');
%     axis([0 1 -0.05 0.05]);
%     legend('show')
% end
    