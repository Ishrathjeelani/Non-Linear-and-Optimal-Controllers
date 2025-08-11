clc;
clear all;
close all;

m_c = 1; % Mass of the cart
m_p = 0.1; % Mass of the pole
l = 0.5; % Half the length of the pole
g = 9.81; % Gravitational acceleration
I=m_c*l^2/12;
kp=5;
kd=0.1;
% Initial conditions
x0 = [0; pi;0; 0];
%x = [1;1;1;1];
tspan = 0:0.1:20;
xdesire = 0;
xdotdesire =0;


Al = [0 0 1 0;
      0 0 0 1;
      0 g*m_p/m_c 0 0;
      0 g*(m_c+m_p)/(l*m_c) 0 0];
Bl = [0;0; 1/m_c; 1/(l*m_c)];
K = lqr(Al,Bl,diag([1 1 1 100]),1);

        
opts = odeset('MaxStep', 0.1,'RelTol',1e-6,'AbsTol',1e-6);

[t,x] = ode45(@(t,x) cart_pole_ode(t,x,l,K,kd,kp,I,m_c,m_p,g), tspan, x0,opts);



%Animate system
myVideo = VideoWriter('PFL');
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
fig = figure;
sgtitle("Time response for partial feedback linearization","Interpreter","latex");
subplot(2,1,1); 
plot(t,rad2deg(x(:,2)),"LineWidth",2);
grid on;
xlabel("time$(s)$","Interpreter","latex");
ylabel("angle $\theta$","Interpreter","latex");
title("Angle $\theta$ vs Time $(s)$ ","Interpreter","latex");

subplot(2,1,2);
plot(t,x(:,1),'LineWidth',2);
grid on;
xlabel("time $(s)$","Interpreter","latex");
ylabel("Distance $x$","Interpreter","latex");
title("Position $x$ vs Time $(s)$ ","Interpreter","latex");
saveas(fig,'timeresponsePFL','png');
fig2 = figure;
plot(x(:,4),x(:,2),'LineWidth',2);
grid on;
xlabel("Angle $\theta$","Interpreter","latex");
ylabel("Angular Velocity $\dot{\theta}$","Interpreter","latex");
title("$\dot{\theta}$ vs $\theta$ phase plot for partial feedback linearization","Interpreter","latex");
saveas(fig2,'PhaseplotPFL','png');


function dxdt = cart_pole_ode(t,x,l,K,kp,kd,I,m_c,m_p,g)

    dxdt = zeros(4,1);
    alpha =0.1;
Etilde = .5*x(4)^2 - g*cos(x(2)) - g;
    delta_theta  = x(4)^2 + (x(2)-pi)^2;
    
    if abs(Etilde) < 1 && delta_theta < 1
           u = - K*(x-[0;pi;0;0]);
    else
         u = (-kp*x(1)-kd*x(2)-alpha*(m_p*g*l-1)*x(4)+m_p*x(4)^2 +...
(m_p*g*sin(2*x(2)))/2 )/(m_c+m_p*sin(x(2))^2);
    end
       
    dxdt(1) = x(2);
dxdt(2)=((m_p*l^2 + I)*(m_p*l*sin(x(3))*x(4)^2 + u))/(m_p^2*l^2 + I*m_c + I*m_p - m_p^2*l^2*cos(x(3))^2 + m_c*m_p*l^2) - (m_p^2*g*l^2*cos(x(3))*sin(x(3)))/(m_p^2*l^2 + I*m_c + I*m_p - m_p^2*l^2*cos(x(3))^2 + m_c*m_p*l^2);
dxdt(3) = x(4);
dxdt(4)=(m_p*g*l*sin(x(3))*(m_c + m_p))/(m_p^2*l^2 + I*m_c + I*m_p - m_p^2*l^2*cos(x(3))^2 + m_c*m_p*l^2) - (m_p*l*cos(x(3))*(m_p*l*sin(x(3))*x(4)^2 + u))/(m_p^2*l^2 + I*m_c + I*m_p - m_p^2*l^2*cos(x(3))^2 + m_c*m_p*l^2);
 end

function fig = animate_cart_pole(t,x,l)
    % Unpack state vector
    x_cart = x(1);
    theta = x(2);
    % Compute position of cart and pole
    x_pole = x_cart + l*sin(theta);
    y_pole = l*cos(theta);
    % Draw cart and pole
    fig = gcf;
    clf;
    hold on;
     title("Animation for cart-pole system with partial feedback linearization","Interpreter","latex");
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
