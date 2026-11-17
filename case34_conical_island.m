% CASE 4 -- Solitary wave runup around a conical island (Briggs et al. 1995)
%
% Choose the incident wave with VARIANT below:
%   'A' : a0 = 0.045   'B' : a0 = 0.091   'C' : a0 = 0.181
% (a0 = normalized wave height, lengths scaled by the still-water depth).
% An initial hump at rest splits into two solitary waves; the left-going
% one is discarded by the domain, the right-going one hits the island.
% Elevation is recorded at the 8 gauges (WG1,2,3,4,6,9,16,22) of the
% experiment. Output: results_conical_island_<VARIANT>.mat
%
% Requires: SWE2D_engine_v2.m, Pnonhydros1Nwd_new.m, closest.m, inv_label.m

clear; close all;

VARIANT = 'A';                       % <-- 'A', 'B' or 'C'

switch VARIANT
    case 'A', a0 = 0.045;
    case 'B', a0 = 0.091;
    case 'C', a0 = 0.181;
    otherwise, error('VARIANT must be ''A'', ''B'' or ''C''.');
end
CASE_NAME = ['conical_island_' VARIANT];

%% Parameters (normalized by the still-water depth)
Lx = 120;
Ly = 70;
h  = 0.32;                           % physical depth [m], used only for scaling
Dz = 0.32/h;                         % normalized still-water depth (=1)
g  = 9.81;
dx = 0.25; dy = 0.25;
dt = 0.005;
thres = 5e-3;                        % wet/dry threshold

options.g     = g;
options.dx    = dx;
options.dy    = dy;
options.dt    = dt;
options.thres = thres;

x = 0:dx:Lx;  y = 0:dy:Ly;
[xx,yy] = meshgrid(x,y);
Nx = length(x); Ny = length(y); Nxy = Nx*Ny;
dz = Dz*ones(Nx,Ny);
T  = 25; Nt = floor(T/dt);

% Solitary wave
gamm = sqrt(0.75*a0);
LL   = (2*h/gamm)*acosh(sqrt(20))/h;   % effective wave length
X0   = 40;                             % position of the initial hump

%% Bottom: conical island centred at (xp,yp)
r1 = 3.6/h;  r2 = 1.1/h;               % base and top radius
xp = X0+LL+r1+5;
yp = 35;
d  = zeros(Nx,Ny);                     % island height above the flat bottom
for k = 1:Ny
    for j = 1:Nx
        r = sqrt((xx(k,j)-xp)^2+(yy(k,j)-yp)^2);
        if r <= r1 && r >= r2
            d(j,k) = 0.25*(r1-r);
        elseif r < r2
            d(j,k) = 0.625/h;
        end
    end
end
d = dz-d;                              % water depth

figure
surf(xx,yy,d','edgecolor','none')
title('Depth')

figure
surf(xx,yy,-d','edgecolor','none')
title('Bottom Topography')

%% Preallocation
E0 = zeros(Nx,Ny);   E1 = zeros(Nx,Ny);   hh = zeros(Nx,Ny);
U0 = zeros(Nx+1,Ny); U1 = zeros(Nx+1,Ny);
V0 = zeros(Nx,Ny+1); V1 = zeros(Nx,Ny+1);
qu = zeros(Nx+1,Ny); qv = zeros(Nx,Ny+1);
W0 = zeros(Nx,Ny);   W1 = zeros(Nx,Ny);    % vertical velocity at the surface
Wb0 = zeros(Nx,Ny);  Wb1 = zeros(Nx,Ny);   % vertical velocity at the bed
pd = zeros(Nxy,1);   P = zeros(Nx,Ny);

%% Initial condition: solitary hump at rest
for k = 1:Ny
    for j = 1:Nx
        E0(j,k) = max(2*a0*(sech(gamm*(j*dx-X0)))^2,-d(j,k)+thres);
        hh(j,k) = E0(j,k)+d(j,k);
    end
end

%% Gauges (Briggs et al. numbering)
gauge_id = [1 2 3 4 6 9 16 22];
gauge_xy = [xp-r1-LL/2, yp+7;
            xp-r1-LL/2, yp+2.34;
            xp-r1-LL/2, yp-2.34;
            xp-r1-LL/2, yp-7;
            xp-r1,      yp;
            xp-8.125,   yp;
            xp,         yp-8.125;
            xp+8.125,   yp];
Ngauge   = numel(gauge_id);
gauge_ix = closest(x,gauge_xy(:,1));
gauge_iy = closest(y,gauge_xy(:,2));

figure
surf(xx,yy,E0','edgecolor','none')
title('Initial eta \eta_0')
view(3,25)
xlabel('x'); ylabel('y'); zlabel('z')
figure
surf(xx,yy,E0','edgecolor','none')
title('Initial eta \eta_0')
view(0,1)
xlabel('x'); zlabel('z')
figure
surf(xx,yy,E0','edgecolor','none')
hold on
for ig = 1:Ngauge
    plot3(x(gauge_ix(ig)),y(gauge_iy(ig)),E0(gauge_ix(ig),gauge_iy(ig)),'.r','markersize',10)
end
view(0,90)
xlabel('x'); ylabel('y')

%% Diagnostics (saved in the results file)
Volume_hist   = zeros(Nt+1,1);
minH_raw_hist = zeros(Nt+1,1);
minH_hist     = zeros(Nt+1,1);
nClipped_hist = zeros(Nt+1,1);
t_fact_hist   = zeros(Nt,1);
t_solve_hist  = zeros(Nt,1);
eta_gauge     = zeros(Nt+1,Ngauge);     % row n+1 = time n*dt; row 1 (t=0) is 0

Volume_hist(1)   = sum(hh(:))*dx*dy;
minH_raw_hist(1) = min(hh(:));
minH_hist(1)     = min(hh(:));

t_total_start = tic;

%% Time loop: finite volume with flooding & drying
colormap('winter'); caxis([-2, 2]);
fig = figure;
set(fig,'Position',[600 200 1200 400],'color',[1 1 1])
for n = 1:(Nt-1)

    % Face depths and fluxes
    for j = 2:Nx
        for k = 1:Ny
            if U0(j,k) > 0
                hh(j,k) = d(j-1,k)+E0(j-1,k);       % upwind: take depth from the left
            elseif U0(j,k) < 0
                hh(j,k) = d(j,k)+E0(j,k);           % upwind: take depth from the right
            else
                hh(j,k) = min(d(j-1,k),d(j,k))+max(E0(j-1,k),E0(j,k));
            end
            qu(j,k) = U0(j,k)*hh(j,k);
        end
    end
    qu(1,1:Ny)    = hh(1,1:Ny).*U0(1,1:Ny);
    qu(Nx+1,1:Ny) = hh(Nx,1:Ny).*U0(Nx+1,1:Ny);
    qv(1:Nx,1)    = hh(1:Nx,1).*V0(1:Nx,1);
    qv(1:Nx,Ny+1) = hh(1:Nx,Ny).*V0(1:Nx,Ny+1);
    for k = 2:Ny
        for j = 1:Nx
            if V0(j,k) > 0
                hh(j,k) = d(j,k-1)+E0(j,k-1);
            elseif V0(j,k) < 0
                hh(j,k) = d(j,k)+E0(j,k);
            else
                hh(j,k) = min(d(j,k-1),d(j,k))+max(E0(j,k-1),E0(j,k));
            end
            qv(j,k) = V0(j,k)*hh(j,k);
        end
    end

    % Surface elevation (mass equation)
    for k = 1:Ny
        for j = 1:Nx
            E1(j,k) = E0(j,k)-dt/dx*(qu(j+1,k)-qu(j,k))-dt/dy*(qv(j,k+1)-qv(j,k));
            hh(j,k) = d(j,k)+E1(j,k);
        end
    end

    % Hydrostatic momentum step and vertical velocity at the bed
    [V1,U1] = SWE2D_engine_v2(E1,U0,V0,hh,qu,qv,options);
    for k = 2:Ny-1
        for j = 2:Nx-1
            if hh(j,k) > thres || hh(j-1,k) > thres
                Wb1(j,k) = -max(0,U1(j,k))*(d(j,k)-d(j-1,k))/dx-min(0,U1(j+1,k))*(d(j+1,k)-d(j,k))/dx ...
                    -max(0,V1(j,k))*(d(j,k)-d(j,k-1))/dy-min(0,V1(j,k+1))*(d(j,k+1)-d(j,k))/dy;
                W1(j,k) = W0(j,k)+(Wb1(j,k)-Wb0(j,k));   % W_t = -Wb_t
            else
                Wb1(j,k) = 0; W1(j,k) = 0;
            end
        end
    end

    % Clip negative depths
    h_before_clip = hh;
    nClipped      = nnz(h_before_clip < 0);
    hh            = max(hh,0);

    % Non-hydrostatic pressure
    [pd,t_fact_hist(n),t_solve_hist(n)] = Pnonhydros1Nwd_new(hh,U1,V1,W1,Wb1,Wb0,options);

    Volume_hist(n+1)   = sum(hh(:))*dx*dy;
    minH_raw_hist(n+1) = min(h_before_clip(:));
    minH_hist(n+1)     = min(hh(:));
    nClipped_hist(n+1) = nClipped;

    for m = 1:Nxy
        [i,j] = inv_label(m,Nx);
        if hh(i,j) > thres
            P(i,j) = pd(m);
        else
            P(i,j) = 0;
        end
    end

    % Pressure correction of U, V, W; update of W
    for k = 1:Ny
        for j = 2:Nx
            U1(j,k) = U1(j,k)-dt/dx/2*(P(j,k)-P(j-1,k));
        end
    end
    for k = 2:Ny
        for j = 1:Nx
            V1(j,k) = V1(j,k)-dt/dy/2*(P(j,k)-P(j,k-1));
        end
    end
    for k = 1:Ny
        for j = 1:Nx
            if hh(j,k) > thres
                W1(j,k) = W1(j,k)+2*dt/hh(j,k)*P(j,k);
            else
                W1(j,k) = 0;
            end
        end
    end
    for k = 2:Ny-1
        for j = 2:Nx-1
            W1(j,k) = W1(j,k)-dt*(U1(j,k)+U1(j+1,k))/2*(W1(j,k)-W1(j-1,k))/dx ...
                -dt*(V1(j,k)+V1(j,k+1))/2*(W1(j,k)-W1(j,k-1))/dy;   % advection -u*w_x - v*w_y
        end
    end
    Wb1(1,1:Ny) = Wb1(2,1:Ny); Wb1(Nx,1:Ny) = Wb1(Nx-1,1:Ny);

    % Plotting
    if mod(n,100) == 1
        figure(6)
        colormap('winter')
        F1 = surf(xx,yy,-d'); hold on;     % bottom
        F2 = surf(xx,yy,E1'); hold off;    % surface elevation
        set(F1,'FaceColor','None','EdgeColor','k');
        set(F2,'FaceColor','interp');
        view(10,15)
        zlim([-0.05,2*a0]);
        title(sprintf('surface profile at t= %1.4s',num2str(n*dt)))
        pause(0.01);
    end

    % Record gauges
    for ig = 1:Ngauge
        eta_gauge(n+1,ig) = E1(gauge_ix(ig),gauge_iy(ig));
    end

    % Advance in time
    E0 = E1;
    U0 = U1; V0 = V1;
    W0 = W1; Wb0 = Wb1;
end

%% Save results and print a summary
T_total = toc(t_total_start);
t_hist  = (0:Nt)*dt;

results = struct();
results.case_name     = CASE_NAME;
results.Nx = Nx; results.Ny = Ny; results.Nt = Nt;
results.dx = dx; results.dy = dy; results.dt = dt;
results.t_hist        = t_hist;
results.eta_gauge     = eta_gauge;       % raw units, columns = gauges 1,2,3,4,6,9,16,22
results.Volume_hist   = Volume_hist;     % mass conservation
results.minH_raw_hist = minH_raw_hist;   % positivity
results.minH_hist     = minH_hist;
results.nClipped_hist = nClipped_hist;
results.T_total       = T_total;         % computational cost
results.t_fact_hist   = t_fact_hist;
results.t_solve_hist  = t_solve_hist;
save(sprintf('results_%s.mat',CASE_NAME),'-struct','results');

Vh = Volume_hist(1:Nt);
fprintf('--- %s ---\n',CASE_NAME);
fprintf('  Total time: %.2f s (%d steps, %dx%d grid, %.2f s/step avg)\n', ...
    T_total,Nt,Nx,Ny,T_total/Nt);
fprintf('  Max relative volume drift: %.3e %%\n', ...
    max(abs((Vh-Vh(1))/Vh(1)*100)));
fprintf('  Clipping triggered on %.3f %% of (cell,step) pairs\n', ...
    100*sum(nClipped_hist(1:Nt))/(Nx*Ny*Nt));
fprintf('  Saved to results_%s.mat\n\n',CASE_NAME);
