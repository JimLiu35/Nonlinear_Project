clear all;
close all;

global xd_individual;
xd = ddp_car_obst();

% boundary conditions in state space
x0 = [0.1; 65; 0; 0; 0];
xf = [0;0;-pi/2;0;0];
T = 20/1000;
% T = 2;

%%%%%%%%% TRAJECTORY GENERATION %%%%%%%%%%%%%

% norm of initial and final velocity along desired path
% it determines how much the curves will bend
% and can be freely chosen
S.u2 = 1;
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
S.k1 = 14 * [1 0 0; 0 1 0; 0 0 1];
S.k2 = 0.6 * [1 0 0; 0 1 0; 0 0 1];

% perturb initial condition
% x = x0 + [.25;.25;.1;.1]
x = x0;

% augmented state with dynamic compensator, i.e xi=u1
xa = x;
Tf=T;
T0=0;
xtrack=[];
% simulate system
for i = 1 : length(xd)-1
    xd_individual = xd(:,i+1);
    [ts, xas] = ode45(@car_ode, [T0 Tf], xa, [], S);
    xtrack=[xtrack xas'];
    xa=xas(end,:)';
    T0=Tf;
    Tf=Tf+T
   
end

% visualize
figure(1)
subplot(1,2,1)
plot(xtrack(1,:), xtrack(2,:), '-b');

legend('desired', 'executed')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = car_h(x)

y = x(1:3);

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
global xd_individual;
% tracking control law

% get desired outputs:
% yd = S.A*poly3(t);
% dyd = S.A*dpoly3(t);
% d2yd = S.A*d2poly3(t);
yd = xd_individual(1:3);
dyd = [xd_individual(5)*cos(xd_individual(3))*cos(xd_individual(4)); ...
    xd_individual(5)*sin(xd_individual(3))*cos(xd_individual(4)); xd_individual(5)*sin(xd_individual(4)) / S.l];

% get current output
y = car_h(xa);

% compensator, i.e.  xi=u1
% xi = xa(end);

% current velocity
dy = [xa(5)*cos(xa(3))*cos(xa(4)); xa(5)*sin(xa(3))*cos(xa(4)); xa(5)*sin(xa(4)) / S.l];

% error state
z1 = y - yd;
z2 = dy - dyd;
d2yd = [0;0;0];

% virtual inputs
v = d2yd - S.k1*z1 -S.k2*z2;
size(v);

% augmented inputs ua=(dxi, u2)
ua = pinv([-cos(xa(3))*sin(xa(4))*xa(5) cos(xa(3))*cos(xa(4));...
           -sin(xa(3))*cos(xa(4))*xa(5) sin(xa(3))*cos(xa(4));...
           cos(xa(4))*xa(5)  sin(xa(4))]) * (v - ...
           [-sin(xa(3))*cos(xa(4))*sin(xa(4))*xa(5)^2; ...
            cos(xa(3))*cos(xa(4))*sin(xa(4))*xa(5)^2;...
            0]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxa = car_ode(t, xa, S)
% unicycle ODE
ua = car_ctrl(t, xa, S);

% xi = xa(end);
% dxi = ua(1);
u1 = ua(1);
u2 = ua(2);


dxa = [xa(5)*cos(xa(3))*cos(xa(4));
       xa(5)*sin(xa(3))*cos(xa(4));
       xa(5)*sin(xa(4)) / S.l;
       u1;
       u2];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
