function f = cem_planning_test()
% EN530.678: demonstration of the cross-entropy method (CEM)
% for simple path planning with obstacles
%
% Marin Kobilarov, marin(at)jhu.edu


sn = 6;  % total intermediate knots
ng = 1;  % number of Gaussian mixture models
iter = 30; % total iterations

S = si_init(sn, ng);

N = 50;  % total samples
p =.1;    % rho-quantile
nf = round(p*N);

S.ss = [];

s = feval(S.f_sinit, S);

crs = inf*ones(1,iter);

video =0;  %set to 1 to record a video
mov = 0;
if video,
  mov = avifile('ce.avi');
  mov.Quality = 100;
  mov.fps = 30;
end

uf = 0;   % set to 1 to generate and save separate plots for each iteration

S.uf = uf;

f1 = [];
f2 = [];

if uf
  f1 = figure
else
  figure
  set(gcf,'Position', [0 0 1280 480])
  subplot(1,3,1)
  subplot('Position',[0.02 0.05 .35 .9])
end

draw_env(S, 0);
view(2)
axis([0 20 0 20])
axis equal
drawnow

if uf
  f3 = figure
else
  subplot(1,3,3)
  subplot('Position',[0.8 0.2 .18 .48])
end
plot(1:iter, crs)

drawnow


if uf
  f2 = figure
else
  subplot(1,3,2)
  subplot('Position',[0.4 0.05 .35 .9])
end

feval(S.f_draws, s, S);

drawnow

pss = zeros(N, S.n*S.sn);

for k=1:iter,

  S.k = k;  
  if uf
    figure(f1);
  else    
    subplot(1,3,1)
    subplot('Position',[0.02 0.05 .35 .9])
  end

  hold off

  draw_env(S, 0);
  view(2)
  axis([0 20 0 20])
  axis equal
  
  disp('sampling')
  tic
  
  j0 = 2;
  if k==1
    j0=1;
  end
    
  for j=j0:N
    pss(j,:) = feval(S.f_sample, s, 1, S);

    xs = feval(S.f_traj, pss(j,:), S);

    S.ss(j).xs = xs;
    S.ss(j).c = feval(S.f_cost, pss(j,:), S);
        
    draw_path(xs, S, 'm', 1, 1);

    draw_bc(S, 0)

    axis([0 20 0 20])
    view(2)    
    drawnow
    if video,
      %    mov = addframe(mov,getframe(gca));
      mov = addframe(mov,gcf);
    end
    
  end
  toc
  
  if k==1
    Js1 = [S.ss.c]  % Question: what is S.ss.c?
    save('ce_Js.txt', 'Js1', '-ascii');
  end

  
  disp('fitting')
  tic  
  [pss, cs] = ce_sort(pss, S);

  a = .9;
  sn = ce_fit(pss(1:nf,:), s, S);

  s.mu = a*sn.mu + (1-a)*s.mu;
  s.Sigma = a*sn.Sigma + (1-a)*s.Sigma;

  toc 
 
  % draw optimal
  xs = feval(S.f_traj, pss(1,:), S);
  draw_path(xs, S, 'k-', 3, 10)

%  for j=1:nf
%    draw_path(xss(:,:,j), S, 'b', 3);
%  end
  drawnow

  if uf
    saveas(gca,[num2str(k) '_1.eps'],'psc2');
  end
  
  crs(k) = cs(1);
  
  if uf
    figure(f3)
  else
    subplot(1,3,3)
    subplot('Position',[0.80 0.2 .18 .48])
  end

  if 1
    hold off   
    if (S.sys == 'ui')
      semilogy(1:k, crs(1:k), '-o', 'LineWidth',3)
    else
      plot(1:k, crs(1:k), '-o', 'LineWidth',3)
    end
      %  set(gca,'XTick',1:k)
    xlabel('iterations', 'FontSize', 15)
    ylabel('cost', 'FontSize', 15)
    
    set(gca, 'FontSize',15)
  end

  if uf
    saveas(gca,[num2str(k) '_3.eps'],'psc2');
  end
  
  if uf
    figure(f2);
  else
    subplot(1,3,2)
    subplot('Position',[0.4 0.05 .35 .9])
  end
  feval(S.f_draws, s, S);
  drawnow
  
  if uf
    saveas(gca,[num2str(k) '_2.eps'],'psc2');
  else
    %uncomment to save a new figure at each iteration
    % saveas(gca,[num2str(k) '.eps'],'psc2');
  end
    
  if video,
    %    mov = addframe(mov,getframe(gca));
    mov = addframe(mov,gcf);
  end
  
end

if video,
  mov = close(mov);
  saveas(gca,['all.eps'],'psc2');
end


function [pss, cs] = ce_sort(pss, S)
N = size(pss, 1); % number of samples
cs = zeros(N, 1);
for i=1:N,
  cs(i) = feval(S.f_cost, pss(i,:), S);
end
[cs,is] = sort(cs, 1, 'ascend');
pss = pss(is,:);


function s = ce_fit(pss, s, S)

% fit independently each point
if (S.ng == 1) 
  for i=1:S.sn,
    % extract all samples at i-th point
    ind = (i-1)*S.n + (1:S.n);
    s.mu(ind) = mean(pss(:,ind));
    s.Sigma(ind,ind) = cov(pss(:,ind)) + diag(S.rf.*rand(S.n,1));
  end
else
  if S.k > inf
    s = gmdistribution.fit(pss, S.ng, 'Start', s);
  else
    s = gmdistribution.fit(pss, S.ng, 'Regularize',.01);
  end
end


%%%%%%% problem specific %%%%%%%%%%

function S = si_init(sn, ng)

S.sys = 'si';

% problem specific params

% general parameters
S.f_sinit = @si_sinit;
S.f_traj = @si_traj;
S.f_draws = @si_draws;
S.f_sample = @si_sample;
S.f_valid = @si_valid;
S.f_cost = @si_cost;

S.n = 2;            % Question: What is S.n? system dimension?
S.sn = sn;
S.ng = ng;

circ = 0            % Circular obstacle

% Here defines obstacle
if ~circ

S.Ps(1).xv = [10, 11, 13, 14, 12, 10];
S.Ps(1).yv = [15, 13, 13, 15, 17, 15];

S.Ps(2).xv = [5.5, 10, 9.2, 5, 7.5, 5.5];
S.Ps(2).yv = [4, 5.5, 10, 10, 7, 4];

S.Ps(3).xv = [14, 15, 16, 17, 14, 14];
S.Ps(3).yv = [5, 4.5, 5, 8, 10, 5];

S.xi = [2.5; 6.5];
%S.xi = [2.5; 1];
S.xf = [18; 15];

else

  S.Os(1).q = [10,10];
  S.Os(1).r =4;

  S.xi = [2; 2];
  %S.xi = [2.5; 1];
  S.xf = [18; 18];


end

%S.xi = [5.5; 1];

S.snf = 6*S.sn;
S.ph = 0;

S.rf = [.1; .1];


function xs = stline(xi, xf, sn)        % Question: stline = straight line?
n = length(xi);
xs = zeros(n, sn + 2);
for i=1:n
  xs(i,:) = linspace(xi(i), xf(i), sn+2);
end


function s = si_sinit(S)
s.mu = zeros(S.ng, S.n*S.sn);
s.Sigma = zeros(S.n*S.sn, S.n*S.sn, S.ng);

xs = stline(S.xi,S.xf, S.sn); 

for j=1:S.ng,
  s.mu(j,:) = reshape(xs(:,2:end-1), S.n*S.sn, 1);
  s.Sigma(:,:,j) = 100*eye(S.n*S.sn);
end

s.PComponents = 1/S.sn*ones(S.ng,1);



function xs = si_traj(ps, S)
xs = [S.xi, reshape(ps, S.n, S.sn), S.xf];
xs = interp1(linspace(0, 1, size(xs,2)), xs', ...
             linspace(0, 1, S.snf), 'cubic')';



function ps = si_sample(s, c, S)

ps = zeros(1, S.n*S.sn);
while(1)
  if (S.ng > 1)
    obj = gmdistribution(s.mu, s.Sigma, s.PComponents);
    ps = random(obj,1);
  else
    for i=1:S.sn,
      ind = (i-1)*S.n + (1:S.n);
      ps(ind) = mvnrnd(s.mu(ind), s.Sigma(ind,ind), 1);
    end
  end
  if ~c || feval(S.f_valid, ps, S)
    break
  end
end



function f = si_cost(ps, S)
ps = reshape(ps, S.n, S.sn);
xsa = [S.xi, ps];
xsb = [ps, S.xf];
dxs = xsb - xsa;
f = sum(sqrt(sum(dxs.*dxs, 1)));


function f = si_draws(s, S)
hold off
axis([0 20 0 20])
draw_env(S, 0);

Z = [];

xr = 0:.5:20;
yr = 0:.5:20;

[X,Y] = meshgrid(xr, yr);
Xr = reshape(X, size(X,1)*size(X,2), 1);
Yr = reshape(Y, size(Y,1)*size(Y,2), 1);

for i=1:S.sn,  
  ind = (i-1)*S.n + (1:2);
  obj = gmdistribution(s.mu(:,ind), s.Sigma(ind,ind,:), ...
                       s.PComponents);

  p = pdf(obj, [Xr, Yr]);
  if i==1
    Z = p;
  else
    Z = Z + p;
  end
  
  
%  ezcontour(@(x,y)pdf(obj,[x y]),[0 20],[0 20])
end


ps = 0

if ps
%  alpha(0);

  h = surf(xr, yr, reshape(Z, length(xr), length(yr)));

  set(h,'FaceAlpha',0);
  view(-20,48)


else
  view(2)
  hold on
  [C,h] = contour(xr, yr, reshape(Z, length(xr), length(yr)));
  set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
  text_handle = clabel(C,h);
  set(text_handle,'FontSize',15)
end
set(gca,'FontSize',15);
drawnow

hold on

for i=1:S.ng,  
  p = reshape(s.mu(i,:), S.n, S.sn);
  %  z = interp2(X,Y,reshape(Z, length(xr), length(yr)), p(1,:),
  %  p(2,:));
  z = 0*ones(size(p,2));
  plot3(p(1,:), p(2,:), z, 'r--', 'LineWidth', 2)
  plot3(p(1,:), p(2,:), z, 'r*', 'MarkerSize', 10, 'LineWidth', 2)

end



function f = draw_path(xs, S, c, lw, ms)

plot3(xs(1,:), xs(2,:), S.ph*ones(size(xs,2),1), [c '-'], ...
      'LineWidth', lw, 'MarkerSize', ms);



function f = draw_env(S, z)

plot3([0,0,20,20], [0,20,20,0], [z z z z], '-k')
hold on

if isfield(S, 'Os')
  for i=1:length(S.Os),
    q = S.Os(i).q;
    r = S.Os(i).r;
    a = 0:.1:2*pi;
    plot3(q(1) + cos(a)*r, q(2) + sin(a)*r, z*ones(size(a)), '-k', 'LineWidth',3)
  end
end

if isfield(S, 'Ps')
  for i=1:length(S.Ps),    
    length(S.Ps(i).xv)
    plot3(S.Ps(i).xv, S.Ps(i).yv, z*ones(size(S.Ps(i).xv)), '-k', 'LineWidth',3);
  end
end

draw_bc(S,z)
set(gca, 'FontSize',15)


function f = draw_bc(S,z)
plot3(S.xi(1), S.xi(2), z, 'og','MarkerSize',10, 'LineWidth',3) 
h = text(S.xi(1)-1.2, S.xi(2)-.8, 'S');
set(h, 'FontSize',20);
plot3(S.xf(1), S.xf(2), z, 'or','MarkerSize',10, 'LineWidth',3) 
h = text(S.xf(1)+.8, S.xf(2)+.4, 'G');
set(h, 'FontSize',20);

set(gca, 'FontSize',15)


function f = si_valid(ps, S)

if si_inside(reshape(ps,2,length(ps)/2), S)
  f = 0
  return
end
xs = si_traj(ps, S);
f = ~si_inside(xs, S);


function f = si_inside(xs, S)

dxs = zeros(size(xs));

if isfield(S, 'Os')
  for o=1:length(S.Os),
    dxs(1,:) = xs(1,:) - S.Os(o).q(1);
    dxs(2,:) = xs(2,:) - S.Os(o).q(2);
    ds = sqrt(sum(dxs.*dxs,1));
    r = S.Os(o).r;
    %  if (ds < r)
    if ds < r
      f = 1;
      return
    end
    
    for i=1:size(xs,2)-1
      a = xs(:,i+1) - xs(:,i);
      a = a/norm(a);
      b = S.Os(o).q(1:2)' - xs(:,i);
      d = b'*a;
      if (d > 0 && norm(d*a - b) < S.Os(o).r)
        f = 1;
        return
      end
    end
      
  end
end

if isfield(S, 'Ps')
  for i=1:length(S.Ps),
    if inpolygon(xs(1,:), xs(2,:), S.Ps(i).xv, S.Ps(i).yv)
      f = 1
      return
    end
    
    [xi, yi] = polyxpoly(xs(1,:), xs(2,:), S.Ps(i).xv, S.Ps(i).yv);
    if ~isempty(xi)
      f = 1;
      return
    end
  end
end


f = 0;