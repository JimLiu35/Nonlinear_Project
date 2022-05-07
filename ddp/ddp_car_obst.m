function output = ddp_car_obst()

% time horizon and segments
tf = 20;
S.N = 1000;
S.h = tf/S.N;

% cost function parameters
S.Q = 0.0*diag([5, 5, 1, 1, 1]);
S.R = diag([1, 5]);
S.Qf = diag([5, 5, 1, 1, 1]);

S.f = @car_f;
S.L = @car_L;
S.Lf = @car_Lf;
S.mu = 0;

% initial state
x0 = [0.1; 65; 0; 0; 0];
S.xf = [0;0;-pi/2;0;0];

S.os(1).p = [-10; 50];
S.os(1).r = 15;
S.os(2).p = [10; 20];
S.os(2).r = 15;
S.os(3).p = [-10; 0];
S.os(3).r = 3;
S.os(4).p = [10; 0];
S.os(4).r = 3;
S.os(5).p = [0;-10];
S.os(5).r = 3;
% S.os(6).p = [-5.5;1.5];
% S.os(6).r = 1;

S.ko = 10;

% initial control sequence
% us = [repmat([.1;0.1], 1, S.N/2), repmat(-[.1;0.1], 1, S.N/2)]/5;
%us = [repmat([.1;0], 1, N/2), repmat(-[.1;0], 1, N/2)]/5;
us = zeros(2,S.N);

xs = ddp_traj(x0, us, S);

Jmin = ddp_cost(xs, us, S)
figure(1)
subplot(1,2,1)

% plot(xs(1,:), xs(2,:), '-b')
hold on
plot(S.xf(1),S.xf(2),'*g')
axis equal

if isfield(S, 'os')
  da = .1;
  a = -da:da:2*pi;
  for i=1:length(S.os)
    % draw obstacle
    plot(S.os(i).p(1) + cos(a)*S.os(i).r,  S.os(i).p(2) + sin(a)*S.os(i).r, ...
         '-r','LineWidth',2);
  end
  axis equal
end



S.a = 1;

for i=1:50
  [dus, V, Vn, dV, a] = ddp(x0, us, S);

  % update controls
  us = us + dus;
  
  S.a = a;   % reuse step-size for efficiency
  
  % update trajectory
  xs = ddp_traj(x0, us, S);

%   plot(xs(1,:), xs(2,:), '-b');
%   pause(1)
end

plot(xs(1,:), xs(2,:), '-g');
fprintf("-------------------------------\n")
disp(xs(1,end))
disp(xs(2,end))
fprintf("-------------------------------\n")
J = ddp_cost(xs, us, S)
output = xs;
xlabel('x')
ylabel('y')

subplot(1,2,2)

plot(0:S.h:tf-S.h, us(1,:),0:S.h:tf-S.h, us(2,:));
xlabel('sec.')
legend('u_1','u_2')



function [L, Lx, Lxx, Lu, Luu] = car_L(k, x, u, S)
% car cost (just standard quadratic cost)
dx = x - S.xf;
if (k == S.N+1)
  L = dx'*S.Qf*dx/2;
  Lx = S.Qf*dx;
  Lxx = S.Qf;
  Lu = [];
  Luu = [];
else
  L = S.h/2*(dx'*S.Q*dx + u'*S.R*u);
  Lx = S.h*S.Q*dx;
  Lxx = S.h*S.Q;
  Lu = S.h*S.R*u;
  Luu = S.h*S.R;
end

% quadratic penalty term
if isfield(S, 'os')
  for i=1:length(S.os)
    g = x(1:2) - S.os(i).p;
    c = S.os(i).r - norm(g);
    if c < 0
      continue
    end
    
    L = L + S.ko/2*c^2;
    v = g/norm(g);
    Lx(1:2) = Lx(1:2) - S.ko*c*v;
    Lxx(1:2,1:2) = Lxx(1:2,1:2) + S.ko*v*v';  % Gauss-Newton appox
  end
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