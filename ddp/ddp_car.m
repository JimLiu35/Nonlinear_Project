function f = ddp_car()

% time horizon and segments
tf = 20;
S.N = 32;
S.h = tf/S.N;

% cost function parameters
S.Q = .0*diag([5, 5, 1, 1, 1]);
S.R = diag([1, 5]);
S.Qf = diag([5, 5, 1, 1, 1]);

S.f = @car_f;
S.L = @car_L;
S.Lf = @car_Lf;
S.mu = 1;

% initial state
x0 = [-5; -2; -1.2; 0; 0];

% initial control sequence
us = [repmat([.1;0.1], 1, S.N/2), repmat(-[.1;0.1], 1, S.N/2)]/5;
%us = [repmat([.1;0], 1, N/2), repmat(-[.1;0], 1, N/2)]/5;

xs = ddp_traj(x0, us, S);

Jmin = ddp_cost(xs, us,  S)

subplot(1,2,1)

plot(xs(1,:), xs(2,:), '-b')
hold on

S.a = 1;

for i=1:20
  [dus, V, Vn, dV, a] = ddp(x0, us, S);

  % update controls
  us = us + dus;
  
  S.a = a;   % reuse step-size for efficiency
  
  % update trajectory
  xs = ddp_traj(x0, us, S);

  plot(xs(1,:), xs(2,:), '-b');
end

plot(xs(1,:), xs(2,:), '-g', 'LineWidth',3);

J = ddp_cost(xs, us, S)

xlabel('x')
ylabel('y')
axis equal

subplot(1,2,2)

plot(0:S.h:tf-S.h, us(1,:),0:S.h:tf-S.h, us(2,:));
xlabel('sec.')
legend('u_1','u_2')



function [L, Lx, Lxx, Lu, Luu] = car_L(k, x, u, S)
% car cost (just standard quadratic cost)

if (k == S.N+1)
  L = x'*S.Qf*x/2;
  Lx = S.Qf*x;
  Lxx = S.Qf;
  Lu = [];
  Luu = [];
else
  L = S.h/2*(x'*S.Q*x + u'*S.R*u);
  Lx = S.h*S.Q*x;
  Lxx = S.h*S.Q;
  Lu = S.h*S.R*u;
  Luu = S.h*S.R;
end


function [x, A, B] = car_f(k, x, u, S)
% car dynamics and jacobians

h = S.h;
c = cos(x(3));
s = sin(x(3));
v = x(4);
w = x(5);

A = [1 0 -h*s*v h*c 0;
     0 1 h*c*v h*s 0;
     0 0 1 0 h;
     0 0 0 1 0;
     0 0 0 0 1];

B = [0 0;
     0 0;
     0 0;
     h 0;
     0 h];

x = [x(1) + h*c*v;
     x(2) + h*s*v;
     x(3) + h*w;
     v + h*u(1);
     w + h*u(2)];