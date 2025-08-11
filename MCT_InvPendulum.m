clc; clear all; close all

%% System Parameters
m = 2;
L = 1;
g = 9.81;

%% Task 1: Open Loop Response
% Non Liner state space representation of open loop system dynamics of
% inverted pendulum

% Define the dynamic equations
func = @(t,x)[x(2);(g/L)*sin(x(1))]; % @(t,x)[f1;f2]
tspan=[0,10];
% Initial Conditions
x0 = [pi 1];
% Solve the non linear function
[ts,xs] = ode45(func,tspan,x0);

% Animate the inverted pendulum
%plotPendulum(xs,ts);

%% Task 2: Non Linear Control
% Write your code here
% Define the closed loop dynamic equations
func_lya = @(t,x)[x(2);(g/L)*sin(x(1))-(g/L)*(sin(x(1))+x(2)+x(1))]; % @(t,x)[f1;f2]
% Initial Conditions
x1 = [pi/2 0];
% Solve the non linear function
[tc,xc] = ode45(func_lya,tspan,x1);
% Animate the inverted pendulum
%plotPendulum(xc,tc);

%% Task 3: LQR
% Write your code here
%Linearized System
A=[0 1; g/L 0]; B=[0; 1]; Q=[1 0;0 1]; R=0.25;
syms p11 p12 p21 p22
P=[p11 p12; p21 p22];
eqns=P*(B*B')*P/R-Q-P*A-A'*P;
soln=solve([eqns(1,1)==0,eqns(1,2)==0,eqns(2,1)==0,eqns(2,2)==0],[p11 p12 p21 p22]);
P_val=[soln.p11(6) soln.p12(6); soln.p21(6) soln.p22(6)];
det(P_val)
U_lqr=-B*B'*P_val/R
A_new = A+U_lqr;
%eig(A_new)
% Define the closed loop dynamic equations
func_lqr = @(t,x)[x(2);(g/L)*sin(x(1))-19.822*x(1)-6.606*x(2)];
% Initial Conditions
x2 = [pi/2 0];
% Solve the non linear function
[tl,xl] = ode45(func_lqr,tspan,x2);
p1=plot(xl(:,1)); 
t = ['$\theta $ and ','$\dot{\theta} $',' vs Time'];
title(t,'interpreter','latex')
hold on;
p2=plot(xl(:,2)); ylabel(['$\theta ,$ ','$\dot{\theta} $'],'interpreter','latex'); xlabel('t(s)');
legend('$\theta $ ','$\dot{\theta} $','interpreter','latex');
hold off;
% Animate the inverted pendulum
plotPendulum(xl,tl);

function plotPendulum(xs,ts)
figure 
for ii=1:length(ts)
    plot([0 cos(0.5*pi-xs(ii))],[0 sin(0.5*pi-xs(ii))],...
        'b-', 'linewidth',1.2), hold on
    plot(cos(0.5*pi-xs(1:ii)),sin(0.5*pi-xs(1:ii)),'k')
    plot(cos(0.5*pi-xs(ii)),sin(0.5*pi-xs(ii)),'bo','markersize',7,...
        'markerfacecolor','b')
    plot(0,0,'ko','markersize',7,'markerfacecolor','k')
    hold off
    axis([-1.2 1.2 -1.2 1.2]);
    box on, grid on
    drawnow
end
end