close all
clear all
A=235;
E=2000;
sigma1=-A:A;
epsilon1=sigma1/E;
plot(epsilon1,sigma1, 'b','LineWidth',1)
grid on
hold on
% title({'deformata'})
xlabel('\epsilon  [-]')
ylabel('\sigma   [N/mm^2]')
epsilon_star=10*max(epsilon1);
epsilon2=max(epsilon1):epsilon_star
sigma2=E*A*(epsilon2+1)/(A+E)
plot(epsilon2,sigma2, 'r','LineWidth',1)
hold on
epsilon3=epsilon_star-2*A/E
epsilon3_1=epsilon_star:-0.1:epsilon3
x1=epsilon_star;
x2=epsilon3;
y1=E*A*(epsilon_star+1)/(A+E)
y2=y1-2*A
sigma3=(y2-y1)*(epsilon3_1-x1)/(x2-x1)+y1
plot(epsilon3_1,sigma3, 'k','LineWidth',1)
hold on



