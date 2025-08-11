clear;clc; close all;
% Parameters
m_c = 1; % Mass of the cart
m_p = 1; % Mass of the pole
l = 0.5; % Half the length of the pole
g = -9.81; % Gravitational acceleration
I = m_p*l^2/12;
% Initial conditions
x0 = [0; 0; pi; 0];
% Time span
tspan = 0:0.1:20;
% Solve ODE
u = 0; % Define input signal here
[t,x] = ode45(@(t,x) cart_pole_ode(t,x,m_c,m_p,l,g,u,I), tspan, x0);

%Animate system
fig = figure;
sgtitle("Time response","Interpreter","latex");
%Angle vs Time
subplot(2,1,1); 
plot(t,rad2deg(x(:,3)),"LineWidth",2);
grid on;
xlabel("time$(s)$","Interpreter","latex");
ylabel("Angle $\theta$","Interpreter","latex");
title("Angle $\theta$ vs Time $(s)$ ","Interpreter","latex");
%Position vs Time
subplot(2,1,2);
plot(t,x(:,1),'LineWidth',2);
grid on;
xlabel("time $(s)$","Interpreter","latex");
ylabel("Distance $x$","Interpreter","latex");
title("Position $x$ vs Time $(s)$ ","Interpreter","latex");
%Phase plot ANgular Velocity vs Angle
fig2 = figure();
plot(x(:,4),x(:,3),'LineWidth',2);
grid on;
xlabel("Angle $\theta$","Interpreter","latex");
ylabel("Angular Velocity $\dot{\theta}$","Interpreter","latex");
title("$\dot{\theta}$ vs $\theta$ phase plot","Interpreter","latex");
saveas(fig,'Passive','png');
saveas(fig2,'PassivePh','png');
%Start of Video creation
myVideo = VideoWriter('Passive');
myVideo.FrameRate = 10;
open(myVideo)
  loops = length(t);
M(loops) = struct('cdata',[],'colormap',[]);
figanim = figure;
for i = 1:length(t)
    animate_cart_pole(t(i),x(i,:)',l);
    M(i) = getframe(gcf);
     writeVideo(myVideo, M(i));
end
close(myVideo);
%End of Video creation
function dxdt = cart_pole_ode(t,x,m_c,m_p,l,g,u,I)
    dxdt = zeros(4,1);
    %% Sliding Mode control Law
%     a1=0.1;b0=2.1;a2=2;
%     u =(-((a1*(x(1)+x(2))+m_p*l*x(4)^2*sin(x(3))-m_p*g*l*sin(x(3)))/m_p*l+b0)*sign(a1*(x(4)+x(3))+x(1)+x(2)));

    %% Sontag
%     k=m_p*(x(4));
%     lf = k*(m_p*g*l*sin(x(3))-m_p*l*x(4)^2*sin(x(3))*m_p*l*cos(x(3))/(m_p+m_c))/((I+m_p*l^2)*(m_p+m_c)-m_p^2*l^2*cos(x(3))^2);
%     lg = -k*m_p*l/((I+m_p*l^2)*(m_p+m_c)-m_p^2*l^2*cos(x(3))^2) +(m_p+m_c)*x(2)/((I+m_p*l^2)*(m_p+m_c)-m_p^2*l^2*cos(x(3))^2);
%         if lg==0
%             u=0;
%         else
%             u=-(lf+(lf^2+lg^4)^0.5)/lg;
%         end
%     
    %% Passivity
    u = m_c*x(2)*sin(x(3))*(l*x(4)^2+g*cos(x(3)))/((m_c+m_p*sin(x(3))^2*x(3))*x(3)-m_c*x(2)) + m_p*l*x(3)*(cos(x(3)));
    
    %% Control Lyapunov Function
%     kp=5;kd=0.1;
%     xdesire=0;xdotdesire=0;
%     qdot =[x(2); x(4)];
%     M= [m_c+m_p m_c*l*cos(x(3));m_c*l*cos(x(3)) I+m_c*l^2];
%     C = [0 -m_p*l*cos(x(3));0 0];
%     G = [0;-m_p*g*l*sin(x(3))];
%     B= [1;0];
%     u = (pinv(M\B))*(-M\C*qdot-M\G) +kp*(xdesire-x(3))+kd*(xdotdesire-x(4));

    %% State-space Equations
    dxdt(1) = x(2);
    dxdt(2) = (u +m_p*l*x(4)^2*sin(x(3))-m_p*l*m_p*g*l*sin(x(3))*cos(x(3))/(I+m_p*l^2))/((m_p+m_c)*(I+m_p*l^2)-m_p^2*l^2*cos(x(3))^2)*(I+m_p*l^2);
    dxdt(3) = x(4);
    dxdt(4) = (m_p*g*l*sin(x(3))-m_p*l*(u+m_p*l*x(4)^2*sin(x(3)))*cos(x(3))/(m_p+m_c))*(m_p+m_c)/((I+m_p*l^2)*(m_p+m_c)-m_p^2*l^2*cos(x(3))^2);
end
function fig = animate_cart_pole(t,x,l)
    % Unpack state vector
    x_cart = x(1);
    theta = x(3);
    % Compute position of cart and pole
    x_pole = x_cart + l*sin(theta);
    y_pole = l*cos(theta);
    % Draw cart and pole
    fig = gcf;
    clf;
    hold on;
    grid on;
     title("Animation for cart-pole system(Passive-based Control)","Interpreter","latex");
    rectangle('Position',[x_cart-0.2 -0.1 0.4 0.2],'FaceColor',[0.2 0.5 0.8],'EdgeColor','k');
    plot([x_cart x_pole],[0 y_pole],'k','LineWidth',4);
    wheel(x_cart-0.15,-0.15,0.05);%Left wheel
    wheel(x_cart-0.15,-0.15,0.015);
    wheel(x_cart+0.15,-0.15,0.05);%right wheel
    wheel(x_cart+0.15,-0.15,0.015);
    plot([x_cart-10*l x_cart+10*l],[-0.2 -0.2],'k','LineWidth',2);%floor
    grid on;
            xl=-5*l;
            xh=5*l;
        if(x_cart>4*l)
            xl=x_cart-9*l;
            xh=x_cart+l;
        elseif(x_cart<-4*l)
            xl=x_cart-l;
            xh=x_cart+9*l;
           
        else
             xl=-5*l;
            xh=5*l;
        end
        text((xl+xh)/2,1.8*l,sprintf("Time $(s)$  $= %.2f s$, $\\theta = %.2f$ deg. ",t,rad2deg(theta)),"Interpreter","latex","HorizontalAlignment","center","Interpreter","latex");
        axis equal;
        xlim manual
        ylim manual
        xlim([xl xh])
        ylim([-1 1])
         drawnow;
  end
function h = wheel(x,y,r)
d = r*2;
px = x-r;
py = y-r;
h = rectangle('Position',[px py d d],'Curvature',[1,1]);
daspect([1,1,1])
end