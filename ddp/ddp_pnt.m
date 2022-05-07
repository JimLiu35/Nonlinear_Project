function f = ddp_pnt()

% time-step and # of segments
S.h = 1;
S.N = 10;

% system mass
S.m = 2;

% cost function specification
S.Q = 0.01*diag([1,1,.5,.5]);
S.R = diag([.1, .1]);
S.Pf = diag([1,1,5,5]);

S.f = @pnt_f;
S.L = @pnt_L;
S.Lf = @pnt_Lf;

S.mu = 0;

% initial state
x0 = [-1; -1; .1; 0];

% initial control 
us = [repmat([.1;.05], 1, S.N/2), repmat(-[.1;.05], 1, S.N/2)]/2;

xs = ddp_traj(x0, us, S);
V = ddp_cost(xs, us, S);

plot(xs(1,:), xs(2,:), '-b')
hold on

% two iteration should be enough for a linear-quadratic problem
for i=1:2
  [dus, V, Vn, dV, a] = ddp(x0, us, S);
  % update control
  us = us + dus;
  S.a = a;
  xs = ddp_traj(x0, us, S);
  plot(xs(1,:), xs(2,:), '-g');
end
plot(xs(1,:), xs(2,:), '-m');



function [L, Lx, Lxx, Lu, Luu] = pnt_L(k, x, u, S)

if k < S.N+1
  L = S.h*(x'*S.Q*x/2 + u'*S.R*u/2);
  Lx = S.h*S.Q*x;
  Lxx = S.h*S.Q;
  Lu = S.h*S.R*u;
  Luu = S.h*S.R;
else
  L = x'*S.Pf*x/2;
  Lx = S.Pf*x;
  Lxx = S.Pf;
  Lu = zeros(S.m,1);
  Luu = zeros(S.m,S.m);
end


function [x, fx, fu] = pnt_f(k, x, u, S)

A = [eye(2), eye(2)*S.h;
     zeros(2), eye(2)];

B = [zeros(2); 
     S.h*eye(2)];

x = A*x + B*u;

fx = A;

fu = B;