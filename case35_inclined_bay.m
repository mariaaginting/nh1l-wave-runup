% CASE 5 -- Solitary wave runup in an inclined bay (Ozkan-Haller & Kirby 1997;
%           Hubbard & Dodd 2002)
%
% A solitary wave is imposed at the left boundary (x=-20) and runs up a bay
% whose shoreline is x = -(L/pi)*cos(pi*y/L). The shoreline elevation is
% tracked at five alongshore sections y/L = 0, 0.25, 0.5, 0.75, 1.
% Output: results_inclined_bay.mat
%
% Set use_nonhydrostatic = false for a plain shallow-water run.
%
% Requires: SWE2D_engine_v2.m, Pnonhydros1Nwd_new.m, closest.m, inv_label.m

clear; close all;

CASE_NAME = 'inclined_bay';

%% Parameters
Lx = 30;
Ly = 10;
h  = 1.273;                    % offshore still-water depth
g  = 9.81;
dx = 0.2; dy = 0.2;
dt = 0.005;
hthr = 5e-3;                   % wet/dry threshold

use_nonhydrostatic = true;     % false = shallow-water run (no pressure correction)

options.g     = g;
options.dx    = dx;
options.dy    = dy;
options.dt    = dt;
options.thres = hthr;
options.wet_wet_coupling = true;   % drop wet-dry couplings in the pressure system

x = -20:dx:Lx;  y = -Ly:dy:Ly;
[xx,yy] = meshgrid(x,y);
Nx = length(x); Ny = length(y); Nxy = Nx*Ny;
T  = 60; Nt = floor(T/dt);
t  = 0:dt:T;

%% Bottom: inclined bay
L  = Ly;
fc = 3*L/pi;
Nc = closest(x,fc);
d  = zeros(Nx,Ny);
for k = 1:Ny
    for j = 1:Nx
        if j > Nc
            d(j,k) = h*(1+(x(j)-fc)/(-L/pi*cos(pi/L*y(k))-fc));
        else
            d(j,k) = h;
        end
    end
end

figure
surf(xx,yy,d','edgecolor','none')
title('Depth')

figure
surf(xx,yy,-d','edgecolor','none')
title('Bottom Topography')

%% Preallocation
E0 = zeros(Nx,Ny);   E1 = zeros(Nx,Ny);
U0 = zeros(Nx+1,Ny); U1 = zeros(Nx+1,Ny);
V0 = zeros(Nx,Ny+1); V1 = zeros(Nx,Ny+1);
qu = zeros(Nx+1,Ny); qv = zeros(Nx,Ny+1);
W0 = zeros(Nx,Ny);   W1 = zeros(Nx,Ny);    % vertical velocity at the surface
Wb0 = zeros(Nx,Ny);  Wb1 = zeros(Nx,Ny);   % vertical velocity at the bed
pd = zeros(Nxy,1);   P = zeros(Nx,Ny);
ss = zeros(Nt,Ny);                         % normalized shoreline elevation

%% Incident solitary wave (time series imposed at x = x(1))
alp = 0.02; Amp = alp*h;
c0  = sqrt(g*h);
gm  = c0/L*(0.75*alp*(1+alp))^0.5*pi/0.4;
X0  = 3.5;

E_in = Amp*sech(gm*(1:Nt)'*dt-X0).^2;

figure
plot(t(1:Nt),E_in,'-r')
xlabel('time'); ylabel('\eta (m)')

%% Initial condition: still water
E0 = max(E0,-d);
hh = E0+d;

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

%% Run-up sections: y/L = 0, 0.25, 0.5, 0.75, 1
%  (the shoreline is assumed on the increasing-x side of the domain)
yoL      = [0 0.25 0.5 0.75 1];
gauge_y  = yoL*Ly;
Ngauge   = numel(gauge_y);
gauge_iy = closest(y,gauge_y);          % y-index of each section

%% Diagnostics (saved in the results file)
Volume_hist   = zeros(Nt+1,1);
minH_raw_hist = zeros(Nt+1,1);
minH_hist     = zeros(Nt+1,1);
nClipped_hist = zeros(Nt+1,1);
t_fact_hist   = zeros(Nt,1);
t_solve_hist  = zeros(Nt,1);
eta_gauge     = zeros(Nt+1,Ngauge);     % run-up (bed elevation at the shoreline)

Volume_hist(1)   = sum(hh(:))*dx*dy;
minH_raw_hist(1) = min(hh(:));
minH_hist(1)     = min(hh(:));
for ig = 1:Ngauge
    eta_gauge(1,ig) = runup_at(hh,d,gauge_iy(ig),hthr);
end

t_total_start = tic;

%% Time loop: finite volume with flooding & drying
colormap('winter'); caxis([-2, 2]);
fig = figure;
set(fig,'Position',[600 200 1200 400],'color',[1 1 1])
for n = 1:(Nt-1)

    % Incident wave at the left boundary
    E0(1,:) = E_in(n);
    U0(1,:) = sqrt(g/h)*E0(1,:);

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

    % Shoreline elevation at every y, normalized by the incident amplitude
    for j = 1:Ny
        xm = find(hh(:,j) > hthr,1,'last');
        if isempty(xm)
            ss(n,j) = 0;               % wave has not reached this row yet
        else
            ss(n,j) = E1(xm,j)/Amp;
        end
    end

    % Clip negative depths and pin eta to the bed on dry land. Negative hh
    % and dry-face accelerations feed Inf/NaN into the pressure system.
    nClipped = nnz(hh < 0);
    minH_raw = min(hh(:));
    E1 = max(E1, -d);
    hh = E1 + d;

    % Hydrostatic momentum step
    [V1,U1] = SWE2D_engine_v2(E1,U0,V0,hh,qu,qv,options);
    U1(1,:) = U0(1,:);                 % incident-wave face (SWE leaves U1(1,:)=0)

    if use_nonhydrostatic
        [U1,V1] = zero_dry_faces(U1,V1,hh,hthr);

        % Vertical velocity at the bed and surface (wet cells only)
        Wb1(:) = 0; W1(:) = 0;
        for k = 2:Ny-1
            for j = 2:Nx-1
                if hh(j,k) > hthr
                    Wb1(j,k) = -max(0,U1(j,k))*(d(j,k)-d(j-1,k))/dx-min(0,U1(j+1,k))*(d(j+1,k)-d(j,k))/dx ...
                        -max(0,V1(j,k))*(d(j,k)-d(j,k-1))/dy-min(0,V1(j,k+1))*(d(j,k+1)-d(j,k))/dy;
                    W1(j,k) = W0(j,k)+(Wb1(j,k)-Wb0(j,k));   % W_t = -Wb_t
                end
            end
        end

        % Non-hydrostatic pressure (wet-wet couplings only)
        [pd,t_fact_hist(n),t_solve_hist(n)] = Pnonhydros1Nwd_new(hh,U1,V1,W1,Wb1,Wb0,options);
        pd(~isfinite(pd)) = 0;
        for m = 1:Nxy
            [i,j] = inv_label(m,Nx);
            if hh(i,j) > hthr
                P(i,j) = pd(m);
            else
                P(i,j) = 0;
            end
        end
        P(~isfinite(P)) = 0;

        % Pressure correction of U, V, W; update of Wb and W
        for k = 1:Ny
            for j = 2:Nx
                if hh(j,k) > hthr || hh(j-1,k) > hthr
                    U1(j,k) = U1(j,k)-dt/dx/2*(P(j,k)-P(j-1,k));
                end
            end
        end
        for k = 2:Ny
            for j = 1:Nx
                if hh(j,k) > hthr || hh(j,k-1) > hthr
                    V1(j,k) = V1(j,k)-dt/dy/2*(P(j,k)-P(j,k-1));
                end
            end
        end
        for k = 1:Ny
            for j = 1:Nx
                if hh(j,k) > hthr
                    W1(j,k) = W1(j,k)+2*dt/hh(j,k)*P(j,k);
                else
                    W1(j,k) = 0;
                end
            end
        end
        [U1,V1] = zero_dry_faces(U1,V1,hh,hthr);
        U1(1,:) = U0(1,:);
        for k = 2:Ny-1
            for j = 2:Nx-1
                if hh(j,k) > hthr
                    Wb1(j,k) = -max(0,U1(j,k))*(d(j,k)-d(j-1,k))/dx-min(0,U1(j+1,k))*(d(j+1,k)-d(j,k))/dx ...
                        -max(0,V1(j,k))*(d(j,k)-d(j,k-1))/dy-min(0,V1(j,k+1))*(d(j,k+1)-d(j,k))/dy;
                    W1(j,k) = W1(j,k)-dt*(U1(j,k)+U1(j+1,k))/2*(W1(j,k)-W1(j-1,k))/dx ...
                        -dt*(V1(j,k)+V1(j,k+1))/2*(W1(j,k)-W1(j,k-1))/dy;   % advection -u*w_x - v*w_y
                else
                    Wb1(j,k) = 0; W1(j,k) = 0;
                end
            end
        end
        Wb1(1,1:Ny) = Wb1(2,1:Ny); Wb1(Nx,1:Ny) = Wb1(Nx-1,1:Ny);
        Wb1(1:Nx,1) = Wb1(1:Nx,2); Wb1(1:Nx,Ny) = Wb1(1:Nx,Ny-1);
        W1(hh <= hthr) = 0;
        Wb1(hh <= hthr) = 0;
    end

    % Diagnostics
    Volume_hist(n+1)   = sum(hh(:))*dx*dy;
    minH_raw_hist(n+1) = minH_raw;
    minH_hist(n+1)     = min(hh(:));
    nClipped_hist(n+1) = nClipped;
    for ig = 1:Ngauge
        eta_gauge(n+1,ig) = runup_at(hh,d,gauge_iy(ig),hthr);
    end

    % Plotting
    if mod(n,100) == 1
        figure(5)
        colormap('winter')
        F1 = surf(xx,yy,-d'); hold on;     % bottom
        F2 = surf(xx,yy,E1'); hold off;    % surface elevation
        set(F1,'FaceColor','None','EdgeColor','k');
        set(F2,'FaceColor','interp','EdgeColor','None');
        view(-39,31)
        zlim([-2*Amp,2*Amp]);
        title(sprintf('surface profile at t= %1.4s',num2str(n*dt)))
        pause(0.01);
        
        figure(6)
        plot(x,E1(:,floor(Ny/2)+1))
        grid on
        xlabel('x[m]'); ylabel('\eta[m]')
        title([' \eta [m] at t = ', num2str(t(n)) ,' s'])
        axis([-20 Lx-5 -0.01 2*Amp])
        pause(0.01)
    end

    figure(6)
    plot(x,E1(:,floor(Ny/2)+1))
    grid on
    xlabel('x[m]'); ylabel('\eta[m]')
    title([' \eta [m] at t = ', num2str(t(n)) ,' s'])
    axis([-20 Lx-5 -0.01 2*Amp])
    pause(0.01)

    % Advance in time
    E0 = E1;
    U0 = U1; V0 = V1;
    W0 = W1; Wb0 = Wb1;
end

%% Shoreline time series at the five sections
s1 = ss(:,floor(Ny/2)+1);       % y = 0
s2 = ss(:,floor(5*Ny/8)+1);     % y = 2.5
s3 = ss(:,floor(6*Ny/8)+1);     % y = 5
s4 = ss(:,floor(7*Ny/8)+1);     % y = 7.5
s5 = ss(:,floor(Ny));           % y = 10

max_runup = max(ss,[],1)';
min_runup = min(ss,[],1)';

figure
plot(y/L,max_runup,'r',y/L,min_runup,'r');
xlabel('y/L'); ylabel('\eta/A');

%% Save results and print a summary
T_total = toc(t_total_start);
t_hist  = (0:Nt)*dt;

results = struct();
results.case_name     = CASE_NAME;
results.Nx = Nx; results.Ny = Ny; results.Nt = Nt;
results.dx = dx; results.dy = dy; results.dt = dt;
results.t_hist        = t_hist;
results.eta_gauge     = eta_gauge;       % raw units; normalize by the incident amplitude when comparing
results.shoreline_ss  = ss;              % shoreline elevation / Amp at every y, size Nt x Ny
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

%% Local functions
function r = runup_at(hh,d,jcol,hthr)
% Bed elevation at the last wet cell of column jcol (shoreline position).
wet = find(hh(:,jcol) > hthr);
if isempty(wet), i_shore = 1; else, i_shore = max(wet); end
r = -d(i_shore,jcol);
end

function [U1,V1] = zero_dry_faces(U1,V1,hh,hthr)
% Zero face velocities whose upwind cell is dry (no flux, no Wb source).
[Nx,Ny] = size(hh);
Uc = U1(2:Nx,:);
Uc((Uc < 0 & hh(2:Nx,:) <= hthr) | (Uc > 0 & hh(1:Nx-1,:) <= hthr)) = 0;
U1(2:Nx,:) = Uc;
Vc = V1(:,2:Ny);
Vc((Vc < 0 & hh(:,2:Ny) <= hthr) | (Vc > 0 & hh(:,1:Ny-1) <= hthr)) = 0;
V1(:,2:Ny) = Vc;
end
