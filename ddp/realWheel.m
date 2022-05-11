clear all;
close all;

global xd_individual ud_individual u u1d index;
u = [];
[xd, ud] = ddp_car_obst();

% boundary conditions in state space
x0 = [0.1; 65; 0; 0.1; 0];
xf = [0;0;-pi/2;0;0];
T = 20/1000;
for i = 1:length(ud)
   if i == 1
       u1d(i) = (ud(1,i+1) - ud(1,i))/T;
   elseif i == length(ud)
       u1d(i) = (- ud(1,i-1) + ud(1,i))/T;
   else
       u1d(i) = (ud(1,i+1) - ud(1,i-1))/2*T;
   end
end
% T = 2;

%%%%%%%%% TRAJECTORY GENERATION %%%%%%%%%%%%%

% norm of initial and final velocity along desired path
% it determines how much the curves will bend
% and can be freely chosen
S.u1 = 1;
S.l = 1;

% % boundary conditions in flat output space 
% % y0 = car_h(x0);
% yf = car_h(xf);
% dy0 = [x0(4) * cos(x0(3)); x0(4) * sin(x0(3))]; % desired starting velocity
% dyf = [xf(4) * cos(xf(3)); xf(4) * sin(xf(3))]; % desired end velocity
% 
% % compute path coefficients
% A = poly3_coeff(y0, dy0, yf, dyf, T);
% 
% % plot desired path
% X = A*poly3([0:.01:T]);
% plot(x(1,:), X(2,:), '-r')
hold on

%%%%%%%%% TRAJECTORY TRACKING %%%%%%%%%%%%%
% S.A = A;

% gains
S.k1 = 2000 * [1 0; 0 1];
S.k2 = 300 * [1 0; 0 1];
S.k3 = 10 * [1 0; 0 1];

% perturb initial condition
% x = x0 + [.25;.25;.1;.1]
x = x0;

% augmented state with dynamic compensator, i.e xi=u1
xa = [x;S.u1];
Tf=T;
T0=0;
xtrack=[];
% simulate system
for i = 1 : length(xd)-1
    xd_individual = xd(:,i+1);
    ud_individual = ud(:,i);
    index = i;
    [ts, xas] = ode45(@car_ode, [T0 Tf], xa, [], S);
    xtrack=[xtrack xas'];
    xa=xas(end,:)';
    err(:,i) = xas(end,1:5)' - xd_individual;
    T0=Tf;
    Tf=Tf+T
    
end

% visualize
figure(1)
subplot(1,2,1)
plot(xtrack(1,:), xtrack(2,:), '-b', 'LineWidth', 1.5);
legend('Goal', 'Obstacle','','','','','Desired','Executed','FontSize',18)
subplot(1,2,2)
hold on
plot(linspace(0,20.02,1000),u(1,:))
plot(linspace(0,20.02,1000),u(2,:))
legend('$u_{1d}$','$u_{2d}$','$u_1$', '$u_2$','Interpreter','latex','FontSize',18)
% legend('u_1', 'u_2')
figure(2)
plot(linspace(0,20.02,1000),err(1,:),'DisplayName','$e_x(t)$')
hold on
plot(linspace(0,20.02,1000),err(2,:),'DisplayName','$e_y(t)$')
plot(linspace(0,20.02,1000),err(3,:),'DisplayName','$e_\theta (t)$')
legend('Interpreter','latex','FontSize',18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = car_h(x)

y = x(1:2);

end

% function f = poly3(t)
% f = [t.^3; t.^2; t; ones(size(t))];
% end
% 
% function f = dpoly3(t)
% f = [3*t.^2; 2*t; ones(size(t)); zeros(size(t))];
% end
% 
% function f = d2poly3(t)
% f = [6*t; 2; zeros(size(t)); zeros(size(t))];
% end
% 
% function A = poly3_coeff(y0, dy0, yf, dyf, T)
% % computes cubic curve connecting (y0,dy0) and (yf, dyf) at time T
% 
% Y = [y0, dy0, yf, dyf];
% L = [poly3(0), dpoly3(0), poly3(T), dpoly3(T)]
% A = Y*inv(L);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ua = car_ctrl(t, xa, S)
global xd_individual ud_individual u1d index;
% tracking control law
xs = xa(1);
ys = xa(2);
thetas = xa(3);
vs = xa(4);
phis = xa(5);

% get desired outputs:
% yd = S.A*poly3(t);
% dyd = S.A*dpoly3(t);
% d2yd = S.A*d2poly3(t);
yd = xd_individual(1:2);
dyd = [xd_individual(4)*cos(xd_individual(3)); ...
    xd_individual(4)*sin(xd_individual(3))];
d2yd = [ud_individual(1)*cos(xd_individual(3))-xd_individual(4)^2*sin(xd_individual(3))*tan(xd_individual(5));
    ud_individual(1)*sin(xd_individual(3))+xd_individual(4)^2*cos(xd_individual(3))*tan(xd_individual(5))];
d3yd = [-3*ud_individual(1)*xd_individual(4)*sin(xd_individual(3))*tan(xd_individual(5))-xd_individual(4)^3*cos(xd_individual(3))*(tan(xd_individual(5)))^2;
       3*ud_individual(1)*xd_individual(4)*cos(xd_individual(3))*tan(xd_individual(5))-xd_individual(4)^3*sin(xd_individual(3))*(tan(xd_individual(5)))^2]...
       + [cos(xd_individual(3)), -(xd_individual(4)^2*sin(xd_individual(3)))/cos(xd_individual(5))^2; sin(xd_individual(3)),  (xd_individual(4)^2*cos(xd_individual(3)))/cos(xd_individual(5))^2] * ...
       [u1d(index);ud_individual(2)];
% get current output
y = car_h(xa);

% compensator, i.e.  xi=u1
xi = xa(end);

% current velocity
dy = [vs*cos(thetas)*cos(phis); vs*sin(thetas)*cos(phis)];
d2y = [xi*cos(thetas) - vs^2*sin(thetas)*tan(phis);
        xi*sin(thetas) + vs^2*cos(thetas)*tan(phis)];

% error state
z1 = y - yd;
z2 = dy - dyd;
z3 = d2y - d2yd;
% d2yd = [0;0;0];

% virtual inputs
v = d3yd - S.k1*z1 -S.k2*z2 - S.k3*z3;
size(v);

% augmented inputs ua=(dxi, u2)
ua = [cos(thetas) sin(thetas); sin(thetas)*((sin(phis))^2-1)/vs^2 cos(thetas)*(cos(phis))^2/vs^2] *...
        (v - [-3*xi*vs*sin(thetas)*tan(phis)-vs^3*cos(thetas)*(tan(phis))^2;
               3*xi*vs*cos(thetas)*tan(phis)-vs^3*sin(thetas)*(tan(phis))^2]);

ua(1) = ua(1) * 1.05;
ua(2) = ua(2) - 3.3 * pi / 180;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxa = car_ode(t, xa, S)
global u index;
% unicycle ODE
ua = car_ctrl(t, xa, S);

% xi = xa(end);
% dxi = ua(1);
xi = xa(end);
dxi = ua(1);
u2 = ua(2);

u(:,index) = [xa(6); u2];

dxa = [xa(4)*cos(xa(3))*cos(xa(5));
       xa(4)*sin(xa(3))*cos(xa(5));
       xa(4)*sin(xa(5)) / S.l;
       xi;
       u2;
       dxi];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
