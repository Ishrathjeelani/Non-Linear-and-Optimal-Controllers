clc; clear all; close all;
%% System Parameters
g = -9.81;
m = 0.18;
I = 0.00025;
L = 0.086;
%% Control Design
%% Solving closed loop system dynamics
% Write your code here
Xd=0.75; Yd=0.1;
% Define the closed loop dynamic equations
kp1=7; kd1=1.42; kp2i=2.5; kd2i=0.56; kp2=0.04; kd2=0.008;
%u1=m*g+kp1*(Yd-x(2))+kd1*(-x(5)); u2=kp2*((kp2i*(Xd-x(1))-kd2i*x(4))-x(3))+kd2*((-kp2i*x(4)+kd2i*sin(x(3))/m*(m*g+kp1*(Yd-x(2))-kd1*x(5)))-x(6));
func = @(t,x)[x(4);x(5);x(6);
    -(m*g+kp1*(Yd-x(2))+kd1*(-x(5)))*sin(x(3))/m;
    -g+(m*g+kp1*(Yd-x(2))+kd1*(-x(5)))*cos(x(3))/m;
    (kp2*((kp2i*(Xd-x(1))-kd2i*x(4))-x(3))+kd2*((-kp2i*x(4)+kd2i*sin(x(3))*(m*g+kp1*(Yd-x(2))-kd1*x(5))/m)-x(6)))/I];

%Lyapunov based
k=1; 
%u1=(m*g-k*x(5))/cos(x(3)); u2=(x(4)*(g-k*x(5))*x(3)*I)/(x(6)*m);
func_lya = @(t,x)[x(4);x(5);x(6);
    -((m*g-k*x(5)+k*x(4)*x(3))*sin(x(3)))/m;
    ((-m*g+(m*g-k*x(5)+k*x(4)*x(3)))*cos(x(3)))/m;
    ((k*x(4)*(m*g)*x(3)-k*x(4)*x(5)*x(3)*(1+1/m)))];
%u1=(m*g-x(5)-x(2)+x(1)+x(4))/cos(x(3)); u2=I*(-x(6));
% func_lya = @(t,x)[x(4);x(5);x(6);
%     -((m*g-x(5)-x(2)+x(1)*m/x(3)+x(4)*x(3))/cos(x(3))*sin(x(3)))/m;
%     (-m*g+(m*g-x(5)-x(2)+x(1)*m/x(3)+x(4)*x(3))/cos(x(3))*cos(x(3)))/m;
%     (-2*x(4)*x(5)*x(3)/m+x(4)*x(3)*g-x(4)*x(2)*x(3)/m-x(5)*x(1)*m/x(3))*x(6)];
tspan=[0,10];
% Initial Conditions
x0 = [0,0,0.01,0.1,0.1,0.01];
% Solve the non linear function
[ts,xs] = ode45(func,tspan,x0);
p1=plot(xs(:,4)); 
t = ['$\dot{x} and $ ','$\dot{y} $ ',' vs Time'];
title(t,'interpreter','latex')
hold on;
p2=plot(xs(:,5)); p3=plot(xs(:,6));
ylabel(['$\dot{x} ,$ ','$\dot{y} $ '],'interpreter','latex'); xlabel('t(s)');
legend('$\dot{x} ,$ ','$\dot{y} $ ','interpreter','latex');
hold off;
% Animate the Quadcopter
animateCopter(xs,ts,L);

function animateCopter(xs,ts,L)
figure
for ii=1:length(ts)
    plot(xs(1:ii,1),xs(1:ii,2),'r.'), hold on
    plot([xs(ii,1)-L*cos(xs(ii,3)) xs(ii,1)+L*cos(xs(ii,3))],...
        [xs(ii,2)-L*sin(xs(ii,3)) xs(ii,2)+L*sin(xs(ii,3))],...
        'ko-', 'linewidth',1,'markerfacecolor','k'),
    plot(xs(ii,1),xs(ii,2),'bs',...
        'markersize',10,'markerfacecolor','b')
    hold off
    axis([-1 1 -1 1]);
    box on, grid on
    drawnow
end
end