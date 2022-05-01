function f = cem_test

% the most basic test of optimizing a quadratic over R^2

opts.v = .9;

opts.iter = 1;

opts.N = 500;

opts.sigma = 0;

d = 3;

x0 = .5*randn(d,1);
C0 = 5*eye(d);

iters = 40;
cs=zeros(iters,1);
x=x0;
opts.C = C0;

opts.H = diag(randn(d));% - .2*randn(d));

global ws
ws = ones(opts.N,1)/opts.N;

for i=1:iters  
  % optimize a simple quadratic
  [x, c, mu, C] = cem(@objfun, x, opts, opts);
  x = mu;
  cs(i) = c;
  opts.C = C;
end
plot(cs);
hold on

opts.tilt = 1;
x=x0;
opts.C = C0;

for i=1:iters
  % optimize a simple quadratic
  [x, c, mu, C] = cem(@objfun, x, opts, opts);
  cs(i) = c;
  opts.C = C;
end
plot(cs,'g');
legend('CE', 'CEt')

function f = objfun(x, opts)

c = .1*ones(length(x),1);

f = x'*opts.H*x + (x-1)'*opts.H*(x-2)/2;
