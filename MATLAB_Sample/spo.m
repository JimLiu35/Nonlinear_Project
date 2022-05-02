function [pxs, pk, ks, xss, Js, S] = spo(px0, pk, pg, S)
% Stochastic Policy Optimization using reward-weighted-regression (RWR)
% or the cross-entropy (CE) method
%
% @param px0 distribution over the initial state
% @param pk initial distribution over the policy parameters k
% @param pg distribution over the context variables g
% @param S: required settings:
%    S.N : number of trajectory segments
%    S.c : control dimension
%    S.M : number of samples per iteration
%    S.f : discrete-time dynamics function x = f(x, u, S)
%    S.L : trajectory cost function Lk = L(x, u, S)
%    S.Lf : terminal cost function Lf = Lf(x, S)
%    S.phi: control policy basis function, Phi = phi(x, g, S), used to
%           set the control according to  u = Phi(x, g)*k for state x
%           and context g
%    S.valid : state constraint function b = valid(x, S) return
%              true/false, (optional)
%    S.usece: whether to use CE or RWR-type method (optional,
%             default is 0), 
%    S.v: use low-pass smoothing param when changing distribution (default
%         is 0.9)
%
% @return pxs resulting estimated distribution over trajectories
% @return pk updated policy distribution
% @return ks sampled policy parameters
% @param xss the trajectory samples
% @param Js corresponding costs
% @param S updated settings
%
% Author: Marin Kobilarov, marin(at)jhu.edu


% set the same random seed
if ~isfield(S,'rs')  
  S.rs = rng;
else
  rng(S.rs);
end

% set the same random seed
if ~isfield(S,'usece')  
  S.usece = 0;
end

% set smoothing param v
if ~isfield(S,'v')  
  S.v = 0.9;
end


N = S.N;
n = length(px0.m);
c = S.c;
ng = length(pg.m);

nk = length(pk.m);

M = S.M;

S.n = n;
S.nk = nk;

xss = zeros(n, N+1, M); % states
uss = zeros(c, N, M);   % controls
gs = zeros(ng, M);      % context (goal, obstacles, etc...)
ks = zeros(nk, M);      % gains

Js = zeros(1,M);        % costs


for j=1:M,
  % sample trajectory and contex

  while(1)
    
    [xs, us, g, k] = sample_traj(px0, pk, pg, S.f, S);
    
    xss(:,:,j) = xs;
    uss(:,:,j) = us;
    gs(:,j) = g;
    ks(:,j) = k;
    
    S.g = g;
    Js(j) = 0;
    for i=1:N,
      Js(j) = Js(j) +  S.L(xss(:,i,j), uss(:,i,j), S);
    end
    Js(j) = Js(j) +  S.Lf(xss(:,N+1,j), S);
    
    if Js(j) < S.Jmax  % sample until lower than Jmax
      break
    else
      Js(j) = S.Jmax;
      break;
    end
    
  end
  
  %  if (S.n>5)
  %    plot3(xs(1,:),xs(2,:),xs(3,:),'-b','LineWidth',1);
  %  else 
  plot(xs(1,:),xs(2,:),'-b','LineWidth',1);
  %  end
  hold on
end

for i=2:N,
  xis = reshape(xss(:,i,:), n, M);
  pxs(i).m = mean(xis')';
  pxs(i).S = cov(xis');
%  plotcov2(pxs(i).m(1:2), pxs(i).S(1:2,1:2), 'plot-opts',{'g'}); 
end

% update params k

cs = Js;

[cmin, imin] = min(cs);

usemin = 0

if usemin
  cmin = min(cs);
  cmax = max(cs);
  
  cvar = cov(cs);

  %  b = 1/mean(cs);
  ws = exp(-(cs - cmin)/(cmin*cvar));
  
else

  if (S.usece)
    Jse = sort(Js);
    qnt = Jse(floor(.1*S.M));
    ws = (Js<qnt);
  else
    b = max(1/min(cs), 1e-3);
    b = 10/max(cs);
    ws = exp(-b*cs);
  end  
  
end
ws = ws/sum(ws);
%disp('ws'), ws'

Sk = wcov(ks', ws);

disp(['det(Sk)=' num2str(det(Sk))]);
pk.m = (1 - S.v).*pk.m + S.v.*(ks*ws');

%  pus(i).S = (1 - v).*pus(i).S + v.*wcov(dus', ws);% + 1e-6*eye(c);    
pk.S = (1 - S.v).*pk.S + S.v.*Sk;
  


function [xs, us, g, k] = sample_traj(px0, pk, pg, f, S)

n = length(px0.m);
c = S.c;
N = S.N;

xs = zeros(n, N+1);
us = zeros(c, N);
while(1)
  
  xs(:,1) = mvnrnd(px0.m, px0.S);
  %xs(:,1) = pxs(1).m;
  
  if (norm(diag(pg.S)) > 1e-16)
    g = mvnrnd(pg.m, pg.S)';
  else
    g = pg.m;
  end
  
  while(1)
    k = mvnrnd(pk.m, pk.S)';
    if ~length(find(k<0))
      break
    end
  end
  
  for i=1:N
    us(:,i) = S.phi(xs(:,i), g, S)'*k;
    xs(:,i+1) = S.f(xs(:,i), us(:,i), S);
  end
  
  if isfield(S, 'valid') & S.valid(xs, [], g, S)
    break
  end
end

