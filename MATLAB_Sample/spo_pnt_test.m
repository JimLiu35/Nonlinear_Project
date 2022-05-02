function f = spo_pnt_test

% whether in 3d or 2d
S.full3d = 0;

sepfig = 1;
saveFrames = 0;

S.M = 100;
S.n = 6;
S.c = 3;

S.usece = 2;

S.N = 100;

S.f = @pnt_f;
S.L = @pnt_L;
S.Lf = @pnt_Lf;
S.phi = @pnt_phi;
S.valid = @pnt_valid;

S.Q = 0.*diag([5, 5, 5, 1, 1, 1]);
S.R = .1*diag([1, 1, 1]);
S.Qf = diag([5, 5, 5, 1, 1, 1]);

% obstacles

if S.full3d
  S.os(1).p = [-1.5; -4.5; 0];
  S.os(1).r = 1.5;
  
  S.os(2).p = [-5; -5; 0];
  S.os(2).r = 2;

  S.os(3).p = [-5; 0; 0];
  S.os(3).r = 1;  
else
  S.os(1).p = [-1.5; -4.5; 0];
  S.os(1).r = 1.5;
  
  S.os(2).p = [-5; -5; 0];
  S.os(2).r = 2;

  S.os(3).p = [-5; 0; 0];
  S.os(3).r = 1;
end

S.ko = 1;  % obstacle gain
S.dr = 10; % detection radius
%S.os = [];

% initial state distribution
if S.full3d
  px0.m = [-10; -10; 0; 0; 0; 0];;
  px0.S = 4*diag([1; 1; 1; .0001; .0001; .0001])
else
  px0.m = [-10; -10; 0; 0; 0; 0];;
  px0.S = 4*diag([1; 1; .0001; .0001; .0001; .0001])
end

% goal distribution
pg.m = [0; 0; 0; 0; 0; 0];
pg.S = .0*diag([1; 1; 1; 0.0001; 0.0001; 0.0001]);

% initial policy parameters (gains) distribution
pk.m = [1; 1];
pk.S = 16*diag([1; 1]);

% add obstacle avoidance gains if obstacles are present
if isfield(S, 'os') & ~isempty(S.os)
  pk.m = [pk.m; S.ko];
  pk.S = diag([diag(pk.S); 20]);
end

N = S.N;
tf = 16;
S.dt = tf/N;

Jms = [];
Jubs = [];
Jebs = [];
Jrobs = [];
S.Jmax = 100;

Pcm = [];  % Pc metrics [Pce, Pch, PcB, Pcr]
kms = [];
kss = [];

f1 = figure
hold on

f2 = figure
f3 = figure
f4 = figure

for it=1:30
  figure(f1)
  
  kms = [kms, pk.m];
  kss = [kss, sqrt(diag(pk.S))];
  
  pk0 = pk;
  [pxs, pk, ks, xss, Js, S] = spo(px0, pk, pg, S);

  % collect expected costs
  Jm = mean(Js);
  Jms = [Jms; Jm];

  % draw obstacles
  if isfield(S, 'os')
    for j=1:length(S.os)
      
      if S.full3d
        draw_obst([], S.os(j), 'r')
      else
        draw_circle(S.os(j), '-r');
      end
    %  draw_circle(S.os(j), '-r');
    
    end
  end
  axis equal
  if S.full3d
    axis([-15 5 -15 5 -5 5])
    view(3)
  else    
    axis([-15 5 -15 5])
  end
    
  if saveFrames
    saveas(gca,['sopt/frames/path' num2str(it) '.eps'],'epsc');
  end
    hold off
    
% Robust bound
  ps = mvnpdf(ks', pk.m', pk.S)';
  p0s = mvnpdf(ks', pk0.m', pk0.S)';
  ws = ps./p0s;
  ws = S.M*ws/sum(ws);
  
  bcs = ones(1,S.M); %flag whether colliding
  for j=1:S.M
    c = pnt_obst(xss(:,N,j), S);
    bcs(j) = (length(find(c > 0)) > 0);
  end
  Mc = sum(bcs);
  Pce = Mc/S.M;  %empirical prob of collision
  
  % optionally we can limit sampling to trajectory with costs less
  % than Jmax: here we set Jmax to the maximum cost of non-colliding
  % trajectories from current iteration
  Ifs = find(~bcs);
  Jfs = Js(Ifs);
  S.Jmax = max(Jfs);

  
  % the following is only for analysis and does not affect the algorithm
  wfs = ws(Ifs);
  smax = rda(pk, pk0, 2)  
  [Jh, Jha, Jrob, a] = bnd_robust(1, pk, pk0, wfs, Jfs, smax, S.Jmax, .9);  
  Jrobs = [Jrobs, [Jh; Jha; Jrob]];  
  [Pch, Pcha, Pcrob, a] = bnd_robust(1, pk, pk0, ws, bcs, smax, 1, .9);
  Pcm = [ Pcm, [Pch; Pcha; Pcrob]];
  
  
  % display
  figure(f2)
  if ~sepfig
    subplot(1,3,1)  
  end
  plot(1:it,kms(1,:), 'r-.', 1:it,kms(2,:), 'g:', 1:it,kms(3,:),'b','LineWidth',3)
  hold on
  if (it >1)
    shadedErrorBar(1:it,kms(1,:),kss(1,:),{'r-.','LineWidth',3},1);
    shadedErrorBar(1:it,kms(2,:),kss(2,:),{'g:','LineWidth',3},1);
    shadedErrorBar(1:it,kms(3,:),kss(3,:),{'b','LineWidth',3},1);
  end
  legend('k_p', 'k_d', 'k_o')
  title('Gains')
  if saveFrames
    saveas(gca,['frames/plot_pk' num2str(it) '.eps'],'epsc');
  end
  hold off

  if sepfig
    figure(f3)
  else
    subplot(1,3,2)
  end
  h=plot(1:it,Jrobs(1,:), 1:it,Jrobs(2,:),1:it,Jrobs(3,:),'LineWidth',3)
  legend('empirical J', 'truncated J', 'robust J^+')
  title('Cost')
  set(h(1),'LineStyle', '-.')
  set(h(2),'LineStyle', ':')
  if saveFrames
    saveas(gca,['frames/plot_J' num2str(it) '.eps'],'epsc');
  end
  
  if sepfig
    figure(f4)
  else
    subplot(1,3,3)
  end
  h=plot(1:it, Pcm(1,:), 1:it, Pcm(2,:),1:it, Pcm(3,:),'LineWidth',3) 
  title('Probability of Collision')
  set(h(1),'LineStyle', '-.')
  set(h(2),'LineStyle', ':')

  legend('empirical P_c', 'truncated P_c', 'robust P^+_c')
  drawnow
  if saveFrames
    saveas(gca,['frames/plot_Pc' num2str(it) '.eps'],'epsc');
  end
  %  input('Press Enter to continue')
  
end




function f = pnt_phi(x, g, S)

d = S.n/2;

f = [ (-x(1:d) + g(1:d))';
      (-x(d+1:end) + g(d+1:end))'];


if ~isempty(S.os)
  G = zeros(3,1);

  p = x(1:d);
  v = x(d+1:end);  
  
  for i=1:length(S.os)
    o = S.os(i);
    N = o.p - p;
    d = norm(N);  
    
    od = abs(d - o.r);
    N = N/d;
    
    if (od < S.dr),
      Sv = cross(N, v/norm(v));
      
      b = sign(asin(norm(Sv)))*acos(dot(N,v/norm(v)));      
      if (norm(b) <= pi/2),
        % disp('obstacle')
        Sv = 1/od*Sv/norm(Sv);
        %ps =  p + R*Sv;
        %plot3([p(1),ps(1)],[p(2),ps(2)],[p(3),ps(3)],'-b');
      else
        Sv = [0;0;0];
      end
      G = G + Sv;
    end
  end  
  
  f = [f;
       (so3_hat(G)*v)'];  
end



function L = pnt_L(x, u, S)
% car cost (just standard quadratic cost)
dx = x - S.g;
L = S.dt/2*(dx'*S.Q*dx + u'*S.R*u);

if 0
if isfield(S, 'os')
  d = S.n/2;
  p = x(1:d);
  for i=1:length(S.os)
    g = p - S.os(i).p;
    c = S.os(i).r - norm(g);
    if c < 0
      continue
    end    
    L = L + 1000/c^2;
  end
end
end


function L = pnt_Lf(x, S)
dx = x - S.g;
L = dx'*S.Qf*dx/2;

function f = pnt_obst(x, S)

f = [];
% returns obstacle function c(x)<0
if isfield(S, 'os')
  d = S.n/2;
  p = x(1:d);
  no = length(S.os);
  if no
    f = zeros(no,1);
  else
    f = [];
  end
  
  for i=1:no
    f(i) = S.os(i).r - norm(p - S.os(i).p);
  end
end



function x = pnt_f(x, u, S)
% car dynamics and jacobians

% if colliding then stop
c = pnt_obst(x, S);

if length(find(c > 0)) > 0
  return
end

d = S.n/2;

% control bounds
umax = 5;
u = min(u, umax);

% add disturbance
u = u + randn(d,1)*.25;  

if ~S.full3d
  u(3) = 0;
end

x = x + S.dt*[x(d+1:end);
              u];


function f = pnt_valid(x, u, g, S)
% box state bounds
f = ~length(find(abs(x)>20));
f = f & (pnt_obst(g, S) <= 0);


function G = draw_circle(o,ls)
da = .1;
a = -da:da:2*pi;
plot(o.p(1) + cos(a)*o.r,  o.p(2) + sin(a)*o.r, ls,'LineWidth',2);
  


%%%%%%%%%%%%%%%%%% below is just for analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = normal_KL(p, q)
k = length(p.m);
f = 1/2*(trace(inv(q.S)*p.S) + (q.m - p.m)'*inv(q.S)*(q.m - p.m) - ...
         k - log(det(p.S)/det(q.S)));




function [Jh, Jha, Jrob, a] = bnd_robust(a, pk, pk0, ws, Js, smax, ...
                                         Jmax, delta)

M = length(Js);

a = sqrt(smax*log(1/delta)/M/smax);

Jh = sum(Js.*ws)/M;

Jha = 1/(M*a)*sum(logtrunc(a*Js.*ws));

Jrob = Jha + Jmax*(a/2*smax + 1/a/M*log(1/delta));



function [Le, Ve, Lb1, Lb2] = emp_bounds(Js, M, d, Jmax)

Le = mean(Js)
Ve = var(Js)

%Ve=0;
Lb1 = Le + Jmax*sqrt(log(1/d)/2/M);
Lb2 = Le + sqrt(2*Ve*log(2/d)/M) + Jmax*7*log(2/d)/3/(M-1);


function y = logtrunc(x)
y = log(1 + x + x.^2/2);


function f = Da(p0, p1, a)
% Renyi divergence

%dS0=det(p0.S)
%det(p1.S)

Sa = (1-a)*p0.S + a*p1.S;
%det(Sa)

dm = p1.m - p0.m;
f = a/2*dm'*inv(Sa)*dm + log(det(Sa)/(det(p0.S)^(1-a)*det(p1.S)^a))/(1-a)/2;


function f = rda(p0, p1, a)
f = exp((a-1)*Da(p0,p1,a));


function G = draw_obst(G, o, c)
  
if isempty(c),
  c = 'b';
end

if (isempty(G))
  
  G.hs = [];
  
  [Xs,Ys,Zs]=ellipsoid(0, 0, 0, o.r, o.r, o.r, 20);
  
  G.Xs{1}=Xs;
  G.Ys{1}=Ys;
  G.Zs{1}=Zs;
  
  G.hs(1) = surf(Xs, Ys, Zs);
  set(G.hs(1),'EdgeColor','none', ...
              'FaceColor', c, ...
              'FaceAlpha', .3, ...
              'FaceLighting','phong', ...
              'AmbientStrength',0.3, ...
              'DiffuseStrength',0.8, ...
              'SpecularStrength',0.9, ...
              'SpecularExponent',25, ...
              'BackFaceLighting','lit');
  
  camlight left;
  hidden off
end

R = eye(3);
x = o.p;

e = so3_log(R,1);

for i=1:length(G.hs),
  set(G.hs(i),'XData', G.Xs{i} + x(1));
  set(G.hs(i),'YData', G.Ys{i} + x(2));
  set(G.hs(i),'ZData', G.Zs{i} + x(3));

  if (norm(e) > .0001)
    rotate(G.hs(i), e/norm(e), norm(e)*180/pi, x);
  end
end

%set(h,'Visible','off')
axis square
axis equal


function f = so3_log(R, varargin)

if (nargin>1)
  if (norm(R-eye(3),'fro') < 2*eps)
    f = zeros(3,1);
    return
  end
end

phi = acos(1/2*(trace(R)-1));

if (nargin>1)
  if (abs(phi) < 1e-10)
    f = zeros(3,1);
    return
  end
end

f = so3_hatinv(phi/(2*sin(phi))*(R-R'));


function f = so3_hat(w)
f=[0, -w(3), w(2);
   w(3), 0, -w(1);
   -w(2), w(1), 0];


function v = so3_hatinv(xi)
v = [xi(3,2); xi(1,3); xi(2,1)];