function f = ddp_arm()

% model parameters
S.m1 = 1;
S.m2 = 1;
S.l1 = .5;
S.l2 = .5;
S.lc1 = .25;
S.lc2 = .25;
S.I1 = S.m1*S.l1/12;
S.I2 = S.m2*S.l2/12;
S.g = 9.8;

% time horizon and number of segments
tf = 2;
S.N = 128;
S.h = tf/S.N;

% Cost function parameters
S.Q = diag([1, 1, .1, .1]);
S.R = diag([.05, .05]);
S.Pf = diag([5, 5, 1, 1]);

S.f = @arm_f;
S.L = @arm_L;
S.mu = 0;

% initial state 
x0 = [0; 0; 0; 0];
% final desired state
S.xf = [pi/2; pi/2; 0; 0];

% initial controls 
us = zeros(2, S.N);

%us = [repmat([.1;0.1], 1, S.N/2), repmat(-[.1;0.1], 1, S.N/2)]/5;
%us = [repmat([.1;0], 1, N/2), repmat(-[.1;0], 1, N/2)]/5;

xs = ddp_traj(x0, us, S);

Jmin = ddp_cost(xs, us,  S)

subplot(1,2,1)

plot(xs(1,:), xs(2,:), '-b')
hold on
plot(S.xf(1),S.xf(2),'*g')
axis equal

S.a = 1;

for i=1:20
  [dus, V, Vn, dV, a] = ddp(x0, us, S);

  % update controls
  us = us + dus;

  S.a = a;   %reuse old step-size for efficiency
  
  % update trajectory using new controls
  xs = ddp_traj(x0, us, S);

  plot(xs(1,:), xs(2,:), '-b');
  
%  keyboard
end
plot(xs(1,:), xs(2,:), '-g','LineWidth',3);

J = ddp_cost(xs, us, S)

xlabel('q_1')
ylabel('q_2')

subplot(1,2,2)

plot(0:S.h:tf-S.h, us(1,:),0:S.h:tf-S.h, us(2,:));
xlabel('sec.')
ylabel('Nm')
legend('u_1','u_2')


function [L, Lx, Lxx, Lu, Luu] = arm_L(k, x, u, S)
% arm cost function

dx = x - S.xf;

if (k==S.N+1)
  L = dx'*S.Pf*dx/2;
  Lx = S.Pf*dx;
  Lxx = S.Pf;
  Lu = [];
  Luu = [];
else
  L = S.h/2*(dx'*S.Q*dx + u'*S.R*u);
  Lx = S.h*S.Q*dx;
  Lxx = S.h*S.Q;
  Lu = S.h*S.R*u;
  Luu = S.h*S.R;
end


function [x, A, B] = arm_f(k, x, u, S)
% arm discrete dynamics
% set jacobians A, B to [] if unavailable
% the following parameters should be set:
% S.m1  : mass of first body
% S.m2  : mass of second body
% S.l1  : length of first body
% S.l2  : length of second body
% S.lc1 : distance to COM
% S.lc2 : distance to COM
% S.I1  : inertia 1
% S.I2  : inertia 2
% S.g   : gravity
%
% S.h : time-step

q = x(1:2);
v = x(3:4);

c1 = cos(q(1));
c2 = cos(q(2));
s2 = sin(q(2));
c12 = cos(q(1) + q(2));

% coriolis matrix
C = -S.m2*S.l1*S.lc2*s2*[v(2), v(1) + v(2);
                    -v(1), 0] + diag([.2;.2]);

% mass elements
m11 = S.m1*S.lc1^2 + S.m2*(S.l1^2 + S.lc2^2 + 2*S.l1*S.lc2*c2) + ...
      S.I1 + S.I2;

m12 = S.m2*(S.lc2^2 + S.l1*S.lc2*c2) + S.I2;

m22 = S.m2*S.lc2^2 + S.I2;

% mass matrix
M = [m11, m12;
     m12, m22];

% gravity vector
fg = [(S.m1*S.lc1 + S.m2*S.l1)*S.g*c1 + S.m2*S.lc2*S.g*c12;
      S.m2*S.lc2*S.g*c12];

% acceleration
a = inv(M)*(u - C*v - fg);
v = v + S.h*a;

x = [q + S.h*v;
     v];

% leave empty to use finite difference approximation
A= [];
B= [];
