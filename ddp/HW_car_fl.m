% boundary conditions in state space
x0 = [-5; -3; 0.5; 1];
xf = [0; 0; 1; 1];
T = 10;

%%%%%%%%% TRAJECTORY GENERATION %%%%%%%%%%%%%

% norm of initial and final velocity along desired path
% it determines how much the curves will bend
% and can be freely chosen
S.u2 = 1;
S.l = 1;

% boundary conditions in flat output space 
y0 = car_h(x0);
yf = car_h(xf);
dy0 = [x0(4) * cos(x0(3)); x0(4) * sin(x0(3))]; % desired starting velocity
dyf = [xf(4) * cos(xf(3)); xf(4) * sin(xf(3))]; % desired end velocity

% compute path coefficients
A = poly3_coeff(y0, dy0, yf, dyf, T);

% plot desired path
X = A*poly3([0:.01:T]);
plot(X(1,:), X(2,:), '-r')
hold on

%%%%%%%%% TRAJECTORY TRACKING %%%%%%%%%%%%%
S.A = A;

% gains
S.k1 = 5 * eye(2);
S.k2 = 2 * eye(2);

% perturb initial condition
x = x0 + [.25;.25;.1;.1]

% augmented state with dynamic compensator, i.e xi=u1
xa = x;

% simulate system
[ts, xas] = ode45(@car_ode, [0 T], xa, [], S);

% visualize
plot(xas(:,1), xas(:,2), '-b');

legend('desired', 'executed')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = car_h(x)

y = x(1:2);

end

function f = poly3(t)
f = [t.^3; t.^2; t; ones(size(t))];
end

function f = dpoly3(t)
f = [3*t.^2; 2*t; ones(size(t)); zeros(size(t))];
end

function f = d2poly3(t)
f = [6*t; 2; zeros(size(t)); zeros(size(t))];
end

function A = poly3_coeff(y0, dy0, yf, dyf, T)
% computes cubic curve connecting (y0,dy0) and (yf, dyf) at time T

Y = [y0, dy0, yf, dyf];
L = [poly3(0), dpoly3(0), poly3(T), dpoly3(T)]
A = Y*inv(L);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ua = car_ctrl(t, xa, S)
% tracking control law

% get desired outputs:
yd = S.A*poly3(t);
dyd = S.A*dpoly3(t);
d2yd = S.A*d2poly3(t);

% get current output
y = car_h(xa);

% compensator, i.e.  xi=u1
% xi = xa(end);

% current velocity
dy = [cos(xa(3)); sin(xa(3))]*xa(4);

% error state
z1 = y - yd;
z2 = dy - dyd;

% virtual inputs
v = d2yd - S.k1*z1 -S.k2*z2;
size(v)

% augmented inputs ua=(dxi, u2)
ua = [v(2)*cos(xa(3)) * S.l/xa(4)^2 - v(1)*sin(xa(3)*S.l/xa(4)^2);
     (v(2)*cos(xa(3)) + v(1)*sin(xa(3)))];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxa = car_ode(t, xa, S)
% unicycle ODE
ua = car_ctrl(t, xa, S);

% xi = xa(end);
% dxi = ua(1);
u1 = ua(1);
u2 = ua(2);


dxa = [cos(xa(3))*xa(4);
       sin(xa(3))*xa(4);
       u1 * xa(4) / S.l;
       u2];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
